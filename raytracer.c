// CS 430 Project 4 RayTracer
// Charles Beck
// 11/15/16
//========= OVERVIEW =========
//	This program was designed to read in a JSON file, parse it, store the data for later use.
//	It will then build a scene based of off the JSON values. The scene can be complex or simple,
//	as in it the program can handle spahes, color, shadows, illumination, reflection, and refraction.
//	
//	The program has a usage after compiling that goes as:
// 		
//			gcc raytracer.c -o raytrace
//
//			raytrace width height output.ppm input.ppm
//	
//	The program will use the data stored to build each pixel and then use ray shooting to detect 
//	instances of interaction and the placement of the pixels. The rays can be used to create shadows, color
//	, reflection, refraction, and the shape of the objects in pixels. After the program has laid out the shadows
//	, illumination, and colors of the objects it will then add the reflection values and refraction values to the 
//	objects in the image. The program outputs a file in P3 ppm format, GIMP is recomended to view image.
//==========================================================================================================
//============ IMPORTS =================
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "parser.c"
#include "ppmwrite.c"
#include "vector_math.h"

#define MAX_DEPTH 7 // max depth for recursion



Object objects[128];
// === STATIC METHODS ==============
static inline double sqr(double v)
{
  return v*v;
}
static inline double v3_len(double* x)
{
    return sqrt(sqr(x[0]) + sqr(x[1]) + sqr(x[2]));
}

static inline void normalize(double* v)
{
  double len = sqrt(sqr(v[0]) + sqr(v[1]) + sqr(v[2]));
  v[0] /= len;
  v[1] /= len;
  v[2] /= len;
}

void v3_reflect(double* x, double* y, double* z)
{
    double scalar = 2.0 * v3_dot(x, y);

    double tmp[3];

    //printf("not working");
    v3_scale(y, scalar, tmp);

    v3_subtract(x, tmp, z);

}

// ====== CLAMP FOR ILL =================
double clamp(double colorVal){
    double max = 1;
    double min = 0;

	if(colorVal > max) return max;

    else if (colorVal < min) return min;

    else   return colorVal;
}

static inline double fres_calc(double a, double b, double m){
    return b * m + a * (1 - m);
}


// ================ PLANE INTERSECT =====================================
double plane_intersect(double* p, double* n, double* Rd, double* Ro){
    double alph,delta;

    normalize(n);

    alph = v3_dot(n, Rd);

    if (fabs(alph) <0.0001) return -1;

    double incident_vector[3];
    v3_subtract(p, Ro, incident_vector);
    delta = v3_dot(incident_vector, n);

    double t = delta/alph; 

    if (t<0.0) return -1;

    return t; 

}
// ======================== SPHERE INTERSECT =============================
double sphere_intersect(double* p, double r, double* Rd, double* Ro)
{

     double a, b;
    double vector_diff[3];
    v3_subtract(Ro, p, vector_diff);

    a = 2 * (Rd[0]*vector_diff[0] + Rd[1]*vector_diff[1] + Rd[2]*vector_diff[2]);
    b = sqr(vector_diff[0]) + sqr(vector_diff[1]) + sqr(vector_diff[2]) - sqr(r);


    double disc = sqr(a) - 4*b;
    double t;
    if (disc < 0) {
		return -1;
	}
    disc = sqrt(disc);
    t = (-a - disc) / 2.0;
    if (t < 0.0) {
		t = (-a + disc) / 2.0;
	}
    if (t < 0.0) {
		return -1;
	}
    return t;
}

// ================ SHADOWS SHOWN =================================================
int shadows(Object objects[], double* newRd, double* newRo, int items, int closestObject, double maxDistance){
    int k;
    int newBest_o = -1;
    normalize(newRd);
    double newBestT = INFINITY;

    for(k = 0; k < items; k++){   
        double newT = 0;
        if (k == closestObject){ 
            continue;
        }
        else if(strcmp(objects[k].type, "sphere") == 0){
            newT = sphere_intersect(objects[k].structures.sphere.position, objects[k].structures.sphere.radius, newRd, newRo);

        } 
		else if(strcmp(objects[k].type, "plane") == 0){
            newT = plane_intersect(objects[k].structures.plane.position, objects[k].structures.plane.normal, newRd, newRo);
        }

        if (maxDistance != INFINITY && newT > maxDistance){
			continue;}

        if (newT > 0 && newT < newBestT){
            newBest_o = k;
        }

    }
    return newBest_o;

}

// ================== DIFFUSE HANDLE =======================================
void diffuseHandle(double *objNormal, double *light, double *illumColor, double *objDiffuse, double *outColor) {

    double normDotLight = v3_dot(objNormal, light);

    if (normDotLight > 0)
    {
        double diffuseProduct[3];

        diffuseProduct[0] = objDiffuse[0] * illumColor[0];
        diffuseProduct[1] = objDiffuse[1] * illumColor[1];
        diffuseProduct[2] = objDiffuse[2] * illumColor[2];

        v3_scale(diffuseProduct, normDotLight, outColor);
    }
    else{
        outColor[0] = 0;
        outColor[1] = 0;
        outColor[2] = 0;
    }
}

// ================= SPECULAR HANDLE =============================================
void specularHandle(double ns, double *light, double *lightRef, double *objNormal, double *V, double *objSpecular, double *illumColor, double *outColor) {

    double rayDotLight = v3_dot(V, lightRef);
    double normDotLight = v3_dot(objNormal, light);

    if (rayDotLight > 0 && normDotLight > 0){
        double rayPowerOfNS = pow(rayDotLight, ns);
        double specularProduct[3];
        specularProduct[0] = objSpecular[0] * illumColor[0];
        specularProduct[1] = objSpecular[1] * illumColor[1];
        specularProduct[2] = objSpecular[2] * illumColor[2];

        v3_scale(specularProduct, rayPowerOfNS, outColor);
    }

    else{
        outColor[0] = 0;
        outColor[1] = 0;
        outColor[2] = 0;
    }
}

// ==================== ANGULAR ATTENUATION ===================================================
double angular_attenuation(Object objects[], double intersection[3], int items, int currLight){
    double thetaRad = objects[currLight].structures.light.theta * (M_PI / 180);
    double cosTheta = cos(thetaRad);
    double directZero = objects[currLight].structures.light.direction[0];
    double directOne = objects[currLight].structures.light.direction[1];
    double directTwo = objects[currLight].structures.light.direction[2];
   if (objects[currLight].structures.light.theta == 0){
	   return 1.0;}
   if (directOne == 0 && directTwo == 0 && directZero == 0){
	   return 1.0;}
   
    double cosAlph= v3_dot(objects[currLight].structures.light.direction, intersection);
    if (cosAlph < cosTheta) return 0.0;
    return pow(cosAlph, objects[currLight].structures.light.angular_a0);
}

// =================== RADIAL ATTENUATION ===================================================
double radial_attenuation(double aOne, double aTwo, double aZero, double distance){
    if(distance == INFINITY){
        return 1.0;
    }

    return 1/(aTwo*pow(distance,2) + aOne*distance + aZero);
}

// ===================== RAYTRACE ================================================================
// This method is the main recursive ray tracing method.
// This is how we acheive reflection and refraction.
// This is done witht the information retreived from the JSON file 
// and the information is used to trace rays through out the image finding interactions
// with objects of a certain type.

void ray_trace(Object objects[], double* Ro, double* Rd, int items, int currDepth, int currObject, double* newColor){
    int i;
    double best_t = INFINITY;
    int best_i = 0;
    for (i = 0; i < items; i++){
        double t = 0;
        if(currObject == i) continue;
        else{
            if(strcmp(objects[i].type, "sphere") == 0){
                t = sphere_intersect(objects[i].structures.sphere.position, objects[i].structures.sphere.radius,Rd, Ro);
           } 
		else if(strcmp(objects[i].type, "plane") == 0){
                t = plane_intersect(objects[i].structures.plane.position, objects[i].structures.plane.normal, Rd, Ro);
            }
        }
        if (t > 0 && t < best_t)
        {
            best_t = t;
            best_i = i;
        }
        int l,k;

        if(best_t > 0 && best_t != INFINITY && best_t != -1){
			
            if(strcmp(objects[best_i].type, "sphere") == 0){
                double reflect = objects[best_i].structures.sphere.reflectivity;
                double refract = objects[best_i].structures.sphere.refractivity;
                double indexOfRefraction = objects[best_i].structures.sphere.ior;
                if((reflect >= 0 || refract >= 0 || indexOfRefraction >=0) && currDepth < 7){
                    double temp[3];
                    double Ron[3];
                    double Rdn[3];
                    v3_scale(Rd, best_t, temp);
                    v3_add(temp, Ro, Ron);
                    v3_subtract(Ron, objects[best_i].structures.sphere.position, Rdn);

                    normalize(Rdn);
                    int inside = 0;
                    double bias = 0.0001;
                    if (v3_dot(Rd,Rdn) > 0)
                    {
                        v3_scale(Rdn, -1, Rdn);
                        inside = 1;
                    }

                    double rayDir[3] = {-Rd[0], -Rd[1], -Rd[2]}; 
                    double facingratio = v3_dot(rayDir,Rdn);
                    double fresneleffect = fres_calc(pow(1 - facingratio, 3), 1, 0.1);
                    double refldir[3];
                    v3_reflect(Rd, Rdn, refldir);
                    normalize(refldir);
                    double reflection[3];
                    double nextRayO[3];
                    v3_scale(Rdn, bias, nextRayO);
                    v3_add(Ron, nextRayO, nextRayO);
                    ray_trace(objects,nextRayO, refldir, items, currDepth+1, best_i, reflection);
                    double refraction[3] = {0,0,0};
                    if(refract > 0){
                        double eta = 0;
                        if(inside) eta = indexOfRefraction;
                        else eta = 1/indexOfRefraction;
                        double inRdn[3] = {-Rdn[0], -Rdn[1], -Rdn[2]};
                        double cosi = v3_dot(inRdn, Rd);
                        double k = 1 - pow(eta, eta) * (1 - cosi * cosi);
                        double temp1 = eta * cosi - sqrt(k);
                        double refractD[3];
                        double temp2[3];
                        v3_scale(Rdn, temp1, temp2);
                        v3_scale(Rd, eta, refractD);
                        v3_add(refractD, temp2, refractD);
                        double refractO[3];
                        v3_scale(Rdn, bias, refractO);
                        v3_subtract(Ron, refractO, refractO);
                        normalize(refractD);
                        ray_trace(objects, refractO, refractD, items, currDepth+1, best_i, refraction);
                    }
                        newColor[0] += reflection[0] * fresneleffect + refraction[0] * (1 - fresneleffect) * refract;
                        newColor[1] += reflection[1] * fresneleffect + refraction[1] * (1 - fresneleffect) * refract;
                        newColor[2] += reflection[2] * fresneleffect + refraction[2] * (1 - fresneleffect) * refract;
				}
                else{
                    for (l = 0; l < items; l++){
                        if(strcmp(objects[l].type, "light") == 0){   
                            double temp[3];
                            double Ron[3];
                            double Rdn[3];
                            v3_scale(Rd, best_t, temp);
                            v3_add(temp, Ro, Ron);
                            v3_subtract(objects[l].structures.light.position, Ron, Rdn);
                            double distanceTLight = v3_len(Rdn);
                            normalize(Rdn);
                            double shadow_intersect = shadows(objects, Rdn, Ron, items, best_i, distanceTLight);
                            if(shadow_intersect != -1){  
                                continue;
                            }
                            else{
                               double sphere_position[3] = {objects[best_i].structures.sphere.position[0],objects[best_i].structures.sphere.position[1],objects[best_i].structures.sphere.position[2]};
                                double n[3] = {Ron[0] - sphere_position[0], Ron[1]-sphere_position[1], Ron[2]-sphere_position[2]}; 
                                normalize(n);
                                double vector_L[3] = {Rdn[0], Rdn[1], Rdn[2]}; 
                                normalize(vector_L);

                                double reflection_L[3];
                                normalize(reflection_L);
                                double V[3] = {Rd[0], Rd[1], Rd[2]}; 
                                double diffuseColor[3];
                                double specularColor[3];
                                double diffuseSpec[3];
                                double object_light_range[3];
                                v3_reflect(vector_L, n, reflection_L);
                                diffuseHandle(n, vector_L, objects[l].structures.light.color, objects[best_i].structures.sphere.diffuseColor, diffuseColor);
                                specularHandle(20, vector_L, reflection_L, n, V, objects[best_i].structures.sphere.specularColor, objects[l].structures.light.color, specularColor);
                                v3_add(diffuseColor, specularColor, diffuseSpec);
                                v3_scale(Rdn, -1, object_light_range);
                                double fang = angular_attenuation(objects, Ron, items, l);
                                double frad = radial_attenuation(objects[l].structures.light.radial_a1, objects[l].structures.light.radial_a2, objects[l].structures.light.radial_a0, distanceTLight);
								newColor[0] += frad * fang * diffuseSpec[0];
                                newColor[1] += frad * fang * diffuseSpec[1];
                                newColor[2] += frad * fang * diffuseSpec[2];
                            }
                        }
                    }
                }
            }
            else if(strcmp(objects[best_i].type, "plane") == 0){

                    for (l = 0; l < items; l++){
                        if(strcmp(objects[l].type, "light") == 0){ 
                            double temp[3];
                            double Ron[3];
                            double Rdn[3];
                            v3_scale(Rd, best_t, temp);
                            v3_add(temp, Ro, Ron);
                            v3_subtract(objects[l].structures.light.position, Ron, Rdn);
                            double distanceTLight = v3_len(Rdn);
                            normalize(Rdn);
                            double shadow_intersect = shadows(objects, Rdn, Ron, items, best_i, distanceTLight);
                            if(shadow_intersect != -1){
                                continue;
                            }
                            else{
                                double n[3] = {objects[best_i].structures.plane.normal[0],objects[best_i].structures.plane.normal[1],objects[best_i].structures.plane.normal[2]}; //objects[best_i]->normal;
                                normalize(n);
                                double vector_L[3] = {Rdn[0], Rdn[1], Rdn[2]}; 
                                double reflection_L[3];
                                double V[3] = {-Rd[0], -Rd[1], -Rd[2]};
                                v3_reflect(vector_L, n, reflection_L); 
                                normalize(reflection_L);
                                double diffuseColor[3];
                                double specularColor[3];
                                double diffuseSpec[3];
                                double object_light_range[3];
                                diffuseHandle(n, vector_L, objects[l].structures.light.color, objects[best_i].structures.plane.diffuseColor, diffuseColor);
                                specularHandle(20, vector_L, V, n, reflection_L, objects[best_i].structures.plane.specularColor, objects[l].structures.light.color, specularColor);
                                v3_add(diffuseColor, specularColor, diffuseSpec);
                                v3_scale(Rdn, -1, object_light_range);
                                double fang = angular_attenuation(objects, Ron, items, l);
                                double frad = radial_attenuation(objects[l].structures.light.radial_a1, objects[l].structures.light.radial_a2, objects[l].structures.light.radial_a0, distanceTLight);
                                newColor[0] += frad * fang * diffuseSpec[0];
                                newColor[1] += frad * fang * diffuseSpec[1];
                                newColor[2] += frad * fang * diffuseSpec[2];
                            }
                        }
                    }
            }
        }
        else{
            newColor[0] = 0*255 ;
            newColor[1] = 0*255;
            newColor[2] = 0*255;

        }
    }
}

// ================ RAYCAST ====================================================
//
// This method handles the heavy work of the image generation
// It will shoot rays to determine the position of lights and objects.
// It will also handle putting all the calculations together that allow us
// to recreate illumination, shadows, and color to create three dimensional objects.
// This method will use the data that was paresd and use the pixeldata to colot objects,
// as well as the position of objects to determine shadows and color changes
//

int ray_cast(Object objects[], Pixmap * buffer, double width, double height, int items){
	double cx, cy, h, w, pixelHeight, pixelWidth;
	int i, x, y;
    double Ro[3] = {0, 0, 0};
	double Rd[3] = {0, 0, 0};
	double point[3] = {0,0, 1};
	double view[2] = {0,0};

	cx = 0;
	cy = 0;

	buffer->width = width;
	buffer->height = height;
	buffer->color = 255;

	for(i = 0; i < items; i++){
		if(strcmp(objects[i].type, "camera") == 0){
			h = (double)objects[i].structures.camera.height;
			w = (double)objects[i].structures.camera.width;
		}
	}
	pixelHeight = h / height;
    pixelWidth = w / width;

	for (y = 0; y < width; y++){
		for (x = 0; x < height; x++){
			point[1] = -(view[1] - h/2.0 + pixelHeight*(y + 0.5));
            point[0] = view[0] - w/2.0 + pixelWidth*(x + 0.5);
            normalize(point);
			Rd[0] = point[0];
            Rd[1] = point[1];
            Rd[2] = point[2];
            double best_t = INFINITY;

			int best_i = 0;
			for (i = 0; i < items; i++)
            {
				double t = 0;
				if(strcmp(objects[i].type, "sphere") == 0){

					t = sphere_intersect(objects[i].structures.sphere.position, objects[i].structures.sphere.radius,Rd, Ro);
				} 
				else if(strcmp(objects[i].type, "plane") == 0){
					t = plane_intersect(objects[i].structures.plane.position, objects[i].structures.plane.normal, Rd, Ro);
				}
				if (t > 0 && t < best_t){
					best_t = t;
                    best_i = i;
				}
				int l,k;
                double* color = malloc(sizeof(double)*3);

				if(best_t > 0 && best_t != INFINITY && best_t != -1){
					color[0] = 0;
                    color[1] = 0;
                    color[2] = 0;

                    if(strcmp(objects[best_i].type, "sphere") == 0){
                        double reflect = objects[best_i].structures.sphere.reflectivity;
                        double refract = objects[best_i].structures.sphere.refractivity;
                        double indexOfRefraction = objects[best_i].structures.sphere.ior;
                        if(reflect >= 0 || refract >= 0 || indexOfRefraction >=0) {
                            double temp[3];
                            double Ron[3];//phit
                            double Rdn[3];//nhit
                            v3_scale(Rd, best_t, temp);
                            v3_add(temp, Ro, Ron);
                            v3_subtract(Ron, objects[best_i].structures.sphere.position, Rdn);

                            normalize(Rdn);
                            int inside = 0;
                            double bias = 0.0001;
                            if (v3_dot(Rd,Rdn) > 0)
                            {
                                v3_scale(Rdn, -1, Rdn);
                                inside = 1;
                            }

                            double rayDir[3] = {-Rd[0], -Rd[1], -Rd[2]}; 
                            double facingratio = v3_dot(rayDir,Rdn);
                            double fresneleffect = fres_calc(pow(1 - facingratio, 3), 1, 0.1);
                            double refldir[3];
                            v3_reflect(Rd, Rdn, refldir);
                            normalize(refldir);
                            double reflection[3];

                            double nextRayO[3];
                            v3_scale(Rdn, bias, nextRayO);
                            v3_add(Ron, nextRayO, nextRayO);

                            ray_trace(objects, nextRayO, refldir, items, 1, best_i, reflection);

                            double refraction[3] = {0,0,0};
                            if(refract > 0)
                            {
                                double eta = 0;
                                if(inside) eta = indexOfRefraction;
                                else eta = 1/indexOfRefraction;
                                double inRdn[3] = {-Rdn[0], -Rdn[1], -Rdn[2]};
                                double cosi = v3_dot(inRdn, Rd);
                                double k = 1 - pow(eta, eta) * (1 - cosi * cosi);
                                double temp1 = eta * cosi - sqrt(k);
                                double refractD[3];
                                double temp2[3];
                                v3_scale(Rdn, temp1, temp2);
                                v3_scale(Rd, eta, refractD);
                                v3_add(refractD, temp2, refractD);
                                double refractO[3];
                                v3_scale(Rdn, bias, refractO);
                                v3_subtract(Ron, refractO, refractO);
                                normalize(refractD);
                                ray_trace(objects, refractO, refractD, items, 1, best_i, refraction);
                            }
                                color[0] += reflection[0] * fresneleffect + refraction[0] * (1 - fresneleffect) * refract;
                                color[1] += reflection[1] * fresneleffect + refraction[1] * (1 - fresneleffect) * refract;
                                color[2] += reflection[2] * fresneleffect + refraction[2] * (1 - fresneleffect) * refract;
                        }
                       for (l = 0; l < items; l++){
                            if(strcmp(objects[l].type, "light") == 0){   
                                double temp[3];
                                double Ron[3];
                                double Rdn[3];
                                v3_scale(Rd, best_t, temp);
                                v3_add(temp, Ro, Ron);
                                v3_subtract(objects[l].structures.light.position, Ron, Rdn);
                                double distanceTLight = v3_len(Rdn);
                                normalize(Rdn);
                                double shadow_intersect = shadows(objects, Rdn, Ron, items, best_i, distanceTLight);
                                if(shadow_intersect != -1)
                                {   

                                    continue;
                                }
                                
                                else{
                                    double sphere_position[3] = {objects[best_i].structures.sphere.position[0],objects[best_i].structures.sphere.position[1],objects[best_i].structures.sphere.position[2]};
                                    double n[3] = {Ron[0] - sphere_position[0], Ron[1]-sphere_position[1], Ron[2]-sphere_position[2]};
                                    normalize(n);
                                    double vector_L[3] = {Rdn[0], Rdn[1], Rdn[2]}; 
                                    normalize(vector_L);
                                    double reflection_L[3];
                                    normalize(reflection_L);
                                    double V[3] = {Rd[0], Rd[1], Rd[2]};
                                    double diffuseColor[3];
                                    double specularColor[3];
                                    double diffuseSpec[3];
                                    double object_light_range[3];
                                    v3_reflect(vector_L, n, reflection_L);
                                    diffuseHandle(n, vector_L, objects[l].structures.light.color, objects[best_i].structures.sphere.diffuseColor, diffuseColor);
                                    specularHandle(20, vector_L, reflection_L, n, V, objects[best_i].structures.sphere.specularColor, objects[l].structures.light.color, specularColor);
                                    v3_add(diffuseColor, specularColor, diffuseSpec);
                                    v3_scale(Rdn, -1, object_light_range);
                                    double fang = angular_attenuation(objects, Ron, items, l);
                                    double frad = radial_attenuation(objects[l].structures.light.radial_a1, objects[l].structures.light.radial_a2, objects[l].structures.light.radial_a0, distanceTLight);
                                    color[0] += frad * fang * diffuseSpec[0];
                                    color[1] += frad * fang * diffuseSpec[1];
                                    color[2] += frad * fang * diffuseSpec[2];
                               }

                            }

                        }

                        buffer->image[y*3 * buffer->width + x*3].r = clamp(color[0]) *255;
                        buffer->image[y*3 * buffer->width + x*3+1].g = clamp(color[1]) *255;
                        buffer->image[y*3 * buffer->width + x*3+2].b = clamp(color[2]) *255;

                    }
                    else if(strcmp(objects[best_i].type, "plane") == 0){ 
                        double reflect = objects[best_i].structures.plane.reflectivity;
                        double refract = objects[best_i].structures.plane.refractivity;
                        double indexOfRefraction = objects[best_i].structures.plane.ior;
                        if(reflect >= 0 || refract >= 0 || indexOfRefraction >=0){
                            double temp[3];
                            double Ron[3];
                            double Rdn[3];
                            v3_scale(Rd, best_t, temp);
                            v3_add(temp, Ro, Ron);
                            v3_subtract(Ron, objects[best_i].structures.plane.position, Rdn);
                            normalize(Rdn);
                            int inside = 0;
                            double bias = 0.0001;
                            if (v3_dot(Rd,Rdn) > 0){
                                v3_scale(Rdn, -1, Rdn);
                                inside = 1;
                            }
                            double rayDir[3] = {-Rd[0], -Rd[1], -Rd[2]}; 
                            double facingratio = v3_dot(rayDir,Rdn);
                            double fresneleffect = fres_calc(pow(1 - facingratio, 3), 1, 0.1);
                            double refldir[3];
                            v3_reflect(Rd, Rdn, refldir);
                            normalize(refldir);
                            double reflection[3];
                            double nextRayO[3];
                            v3_scale(Rdn, bias, nextRayO);
                            v3_add(Ron, nextRayO, nextRayO);
                            ray_trace(objects, nextRayO, refldir, items, 1, best_i, reflection);
                            double refraction[3] = {0,0,0};
                            if(refract > 0){
                                double eta = 0;
                                if(inside) eta = indexOfRefraction;
                                else eta = 1/indexOfRefraction;
                                double inRdn[3] = {-Rdn[0], -Rdn[1], -Rdn[2]};
                                double cosi = v3_dot(inRdn, Rd);
                                double k = 1 - pow(eta, eta) * (1 - cosi * cosi);
                                double temp1 = eta * cosi - sqrt(k);
                                double refractD[3];
                                double temp2[3];
                                v3_scale(Rdn, temp1, temp2);
                                v3_scale(Rd, eta, refractD);
                                v3_add(refractD, temp2, refractD);
                                double refractO[3];
                                v3_scale(Rdn, bias, refractO);
                                v3_subtract(Ron, refractO, refractO);
                                normalize(refractD);
                                ray_trace(objects, refractO, refractD, items, 1, best_i, refraction);
                            }
                                color[0] += reflection[0] * fresneleffect + refraction[0] * (1 - fresneleffect) * refract;
                                color[1] += reflection[1] * fresneleffect + refraction[1] * (1 - fresneleffect) * refract;
                                color[2] += reflection[2] * fresneleffect + refraction[2] * (1 - fresneleffect) * refract;
                        }
      
                        for (l = 0; l < items; l++){
                            if(strcmp(objects[l].type, "light") == 0){   
                                double temp[3];
                                double Ron[3];
                                double Rdn[3];
                                v3_scale(Rd, best_t, temp);
                                v3_add(temp, Ro, Ron);
                                v3_subtract(objects[l].structures.light.position, Ron, Rdn);
                                double distanceTLight = v3_len(Rdn);
                                normalize(Rdn);
                                double shadow_intersect = shadows(objects, Rdn, Ron, items, best_i, distanceTLight);
                                if(shadow_intersect != -1){
									continue;
                                }
                                else{
                                    double n[3] = {objects[best_i].structures.plane.normal[0],objects[best_i].structures.plane.normal[1],objects[best_i].structures.plane.normal[2]}; //objects[best_i]->normal;
                                    normalize(n);
                                    double vector_L[3] = {Rdn[0], Rdn[1], Rdn[2]};
                                    normalize(vector_L);
                                    double reflection_L[3];
                                    double V[3] = {-Rd[0], -Rd[1], -Rd[2]};
                                    v3_reflect(vector_L, n, reflection_L); 
                                    normalize(reflection_L);
                                    double diffuseColor[3];
                                    double specularColor[3];
                                    double diffuseSpec[3];
                                    double object_light_range[3];
                                    diffuseHandle(n, vector_L, objects[l].structures.light.color, objects[best_i].structures.plane.diffuseColor, diffuseColor);
                                    specularHandle(20, vector_L, V, n, reflection_L, objects[best_i].structures.plane.specularColor, objects[l].structures.light.color, specularColor);
                                    v3_add(diffuseColor, specularColor, diffuseSpec);
                                    v3_scale(Rdn, -1, object_light_range);
                                    double fang = angular_attenuation(objects, Ron, items, l);
                                    double frad = radial_attenuation(objects[l].structures.light.radial_a1, objects[l].structures.light.radial_a2, objects[l].structures.light.radial_a0, distanceTLight);
                                    color[0] += frad * fang * diffuseSpec[0];
                                    color[1] += frad * fang * diffuseSpec[1];
                                    color[2] += frad * fang * diffuseSpec[2];
                                }
                            }

                        }
						// add color to image
                        buffer->image[y*3 * buffer->width + x*3].r = clamp(color[0]) *255;
                        buffer->image[y*3 * buffer->width + x*3+1].g = clamp(color[1]) *255;
                        buffer->image[y*3 * buffer->width + x*3+2].b = clamp(color[2]) *255;
                    }
				}
				else{
					// If all else fails color one color
                    buffer->image[y*3 * buffer->width + x*3].r = 0*255 ;
                    buffer->image[y*3 * buffer->width + x*3+1].g = 0*255;
                    buffer->image[y*3 * buffer->width + x*3+2].b = 0*255;
				}
				// free up color space
            free(color);
			}
		}
	}
	return 0;
}

int main(int argc, char *argv[]){
	FILE *json;
	int items, i;
	int ppmFormat = 3;
	double width, height;
	width = atof(argv[1]);
	height = atof(argv[2]);

	Pixmap picbuffer;
    picbuffer.image = (PixelColor*)malloc(sizeof(PixelColor)*width* (height*3));
	json = fopen(argv[4], "r");
	if(json == NULL){
		fprintf(stderr, "Error: could not open file.\n");
		fclose(json);
		exit(-1);
	}
	else{
		items = read_scene(json, objects);
        ray_cast(objects, &picbuffer, width, height, items);

        int size = height * width;

        ppmWriter(&picbuffer, argv[3], size , ppmFormat);

	}
    fclose(json);
    free(picbuffer.image);

	return 0;
}
