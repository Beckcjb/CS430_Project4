// Parser
// Charles Beck
// CS 430
// Project 4 parser
//
//===== IMPORTS ==========
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include "objectsStructs.h"

int lineNumber = 1;

// ================== NEXT Character =================
int next_c(FILE* json) {
    int c = fgetc(json);
    // If new line add one to the line number counter
    if (c == '\n')
    {
        lineNumber += 1;
    }

    if (c == EOF)
    {
        fprintf(stderr, "Error: Unexpected end of file on line number %d.\n", lineNumber);
        fclose(json);
        exit(1);
    }

    return c;
}

// =============== EXPECT C NEXT ================
void expect_c(FILE* json, int d) {
    int c = next_c(json);
    if (c == d){
		return;
	}
    fprintf(stderr, "Error: Expected '%c' on line %d.\n", d, lineNumber);
    fclose(json);
    exit(1);
}

//================== SKIP WHITE SPACE =================
void skip_ws(FILE* json) {
    int c = next_c(json);

    while (isspace(c))
    {
      c = next_c(json);
    }
    ungetc(c, json);
}

// =============== NEXT STRING ================

char* next_str(FILE* json)
{
    char buffer[129];
    int c = next_c(json), i = 0;

    if (c != '"')
    {
        fprintf(stderr, "Error: Expected string on line %d.\n", lineNumber);
        fclose(json);
        exit(-1);
    }
    c = next_c(json);

    while (c != '"')
    {
        if (i >= 128)
        {
          fprintf(stderr, "Error: Strings longer than 128 characters in length are not supported.\n");
          fclose(json);
          exit(-1);
        }
        if (c == '\\')
        {
          fprintf(stderr, "Error: Strings with escape codes are not supported.\n");
          fclose(json);
          exit(-1);
        }
        if (c < 32 || c > 126)
        {
          fprintf(stderr, "Error: Strings may contain only ascii characters.\n");
          fclose(json);
          exit(-1);
        }
        buffer[i] = c;
        i += 1;
        c = next_c(json);
    }
    buffer[i] = 0;
    return strdup(buffer);
}

// =================== NEXT NUMBER ======================
double next_num(FILE* json)
{
    double value;
    if(fscanf(json, "%lf", &value) == 0)
    {
        fprintf(stderr, "Error, line number %d; expected numeric value.\n", lineNumber);
        fclose(json);
        exit(1);
    }
    return value;
}

// ============ NEXT VECTOR================
double* next_vec(FILE* json){
    double* vector = malloc(3*sizeof(double));

    expect_c(json, '[');
    skip_ws(json);
    vector[0] = next_num(json);

    skip_ws(json);
    expect_c(json, ',');
    skip_ws(json);
    vector[1] = next_num(json);

    skip_ws(json);
    expect_c(json, ',');
    skip_ws(json);
    vector[2] = next_num(json);

    skip_ws(json);
    expect_c(json, ']');

    return vector;
}
// ============ READ SCENE ===================
/* In the Read scene method the program will parse the json file
*  and store the data based on the name of the values. The method
*  will save any data that is valid for later use.
*/
int read_scene(FILE *json, Object objects[]){
   
    int c, items;
	double *vector;
	char *name, *value;
	items = 0;

	skip_ws(json);

	c = next_c(json);

	if(c != '[')
    {
		fprintf(stderr, "Error, line number %d; invalid scene definition '%c'\n", lineNumber, c);
		fclose(json);
		exit(-1);
	}

	skip_ws(json);
	c = next_c(json);

	if(c != ']') ungetc(c, json);

	while(c != ']')
    {
		skip_ws(json);
		c = next_c(json);

		if(c != '{')
        {
			fprintf(stderr, "Error, line number %d; invalid object definition '%c'\n", lineNumber, c);
			fclose(json);
			exit(-1);
		}

		skip_ws(json);
		c = next_c(json);

		while(c != '}')
        {
			if(c == '"') ungetc(c, json);

			name = next_str(json);

			if(strcmp(name, "type") == 0)
            {
				skip_ws(json);
				c = next_c(json);

				if(c != ':')
                {
					fprintf(stderr, "Error, line number %d; invalid separator '%c'.\n", lineNumber, c);
					fclose(json);
					exit(-1);

				}
				else
                {
					skip_ws(json);
					value = next_str(json);
					objects[items].type = value;
					if(strcmp(value, "light") == 0)
                    {
                        objects[items].structures.light.angular_a0 = 0;
                        objects[items].structures.light.radial_a1 = 0;
                        objects[items].structures.light.radial_a2 = 0;
                        objects[items].structures.light.radial_a0 = 0;

                        objects[items].structures.light.direction[0] = 0;
                        objects[items].structures.light.direction[1] = 0;
                        objects[items].structures.light.direction[2] = 0;

                    }
				}

			} else if(strcmp(name, "width") == 0) {

				skip_ws(json);
				c = next_c(json);

				if(c != ':')
                {
					fprintf(stderr, "Error, line number %d; invalid separator '%c'.\n", lineNumber, c);
					fclose(json);
					exit(-3);
				}
				else
                {
					skip_ws(json);
					objects[items].structures.camera.width = next_num(json);
				}

			} else if(strcmp(name, "height") == 0) {

				skip_ws(json);
				c = next_c(json);

				if(c != ':')
                {
					fprintf(stderr, "Error, line number %d; invalid separator '%c'.\n", lineNumber, c);
					fclose(json);
					exit(-1);
				}
				else
                {
					skip_ws(json);
					objects[items].structures.camera.height = next_num(json);
				}

			} else if(strcmp(name, "radius") == 0) {

				skip_ws(json);
				c = next_c(json);

				if(c != ':')
                {
					fprintf(stderr, "Error, line number %d; invalid separator '%c'.\n", lineNumber, c);
					fclose(json);
					exit(-1);
				}
				else
                {
					skip_ws(json);
					objects[items].structures.sphere.radius = next_num(json);
				}

			} else if(strcmp(name, "reflectivity") == 0) {

				skip_ws(json);
				c = next_c(json);

				if(c != ':')
                {
					fprintf(stderr, "Error, line number %d; invalid separator '%c'.\n", lineNumber, c);
					fclose(json);
					exit(-1);
				}
				else
                {
					skip_ws(json);
                    if(strcmp(objects[items].type, "plane") == 0)
                        objects[items].structures.plane.reflectivity  = next_num(json);
                    else if(strcmp(objects[items].type, "sphere") == 0)
                        objects[items].structures.sphere.reflectivity  = next_num(json);


				}

			} else if(strcmp(name, "refractivity") == 0) {

				skip_ws(json);
				c = next_c(json);

				if(c != ':')
                {
					fprintf(stderr, "Error, line number %d; invalid separator '%c'.\n", lineNumber, c);
					fclose(json);
					exit(-1);
				}
				else
                {
					skip_ws(json);
                    if(strcmp(objects[items].type, "plane") == 0)
                        objects[items].structures.plane.refractivity  = next_num(json);
                    else if(strcmp(objects[items].type, "sphere") == 0)
                        objects[items].structures.sphere.refractivity  = next_num(json);				}

			} else if(strcmp(name, "ior") == 0) {

				skip_ws(json);
				c = next_c(json);

				if(c != ':')
                {
					fprintf(stderr, "Error, line number %d; invalid separator '%c'.\n", lineNumber, c);
					fclose(json);
					exit(-1);
				}
				else
                {
					skip_ws(json);
                    if(strcmp(objects[items].type, "plane") == 0)
                        objects[items].structures.plane.ior  = next_num(json);
                    else if(strcmp(objects[items].type, "sphere") == 0)
                        objects[items].structures.sphere.ior  = next_num(json);				}

			} else if(strcmp(name, "radial-a2") == 0) {

				skip_ws(json);
				c = next_c(json);

				if(c != ':')
                {
					fprintf(stderr, "Error, line number %d; invalid separator '%c'.\n", lineNumber, c);
					fclose(json);
					exit(-1);
				}
				else
                {
					skip_ws(json);
					objects[items].structures.light.radial_a2 = next_num(json);
				}

			} else if(strcmp(name, "radial-a1") == 0) {

				skip_ws(json);
				c = next_c(json);

				if(c != ':')
                {
					fprintf(stderr, "Error, line number %d; invalid separator '%c'.\n", lineNumber, c);
					fclose(json);
					exit(-1);
				}
				else
                {
					skip_ws(json);
					objects[items].structures.light.radial_a1 = next_num(json);
				}

			} else if(strcmp(name, "radial-a0") == 0) {

				skip_ws(json);
				c = next_c(json);

				if(c != ':')
                {
					fprintf(stderr, "Error, line number %d; invalid separator '%c'.\n", lineNumber, c);
					fclose(json);
					exit(-1);
				}
				else
                {
					skip_ws(json);
					objects[items].structures.light.radial_a0 = next_num(json);
				}

			} else if(strcmp(name, "angular-a0") == 0) {

				skip_ws(json);
				c = next_c(json);

				if(c != ':')
                {
					fprintf(stderr, "Error, line number %d; invalid separator '%c'.\n", lineNumber, c);
					fclose(json);
					exit(-1);
				}
				else
                {
					skip_ws(json);
					objects[items].structures.light.angular_a0 = next_num(json);
				}

			} else if(strcmp(name, "theta") == 0) {

				skip_ws(json);
				c = next_c(json);

				if(c != ':')
                {
					fprintf(stderr, "Error, line number %d; invalid separator '%c'.\n", lineNumber, c);
					fclose(json);
					exit(-1);
				}
				// If it is the angular_a0 then lets place it in the camera structure in the objects array
				else
                {
					skip_ws(json);
					objects[items].structures.light.theta = next_num(json);
				}

			}else if(strcmp(name, "color") == 0) {

				skip_ws(json);
				c = next_c(json);

				if(c != ':')
                {
					fprintf(stderr, "Error, line number %d; invalid separator '%c'.\n", lineNumber, c);
					fclose(json);
					exit(-1);
				}
				else
                {
					skip_ws(json);
					vector = next_vec(json);

                        if(strcmp(objects[items].type, "light") == 0){
                        objects[items].structures.light.color[0] = vector[0];
						objects[items].structures.light.color[1] = vector[1];
						objects[items].structures.light.color[2] = vector[2];
					}

				}

			}else if(strcmp(name, "diffuse_color") == 0) {

				skip_ws(json);
				c = next_c(json);
				if(c != ':'){
					fprintf(stderr, "Error, line number %d; invalid separator '%c'.\n", lineNumber, c);
					fclose(json);
					exit(-1);
				}
				else{
					skip_ws(json);
					vector = next_vec(json);
					if(strcmp(objects[items].type, "plane") == 0) {

						objects[items].structures.plane.diffuseColor[0] = vector[0];
						objects[items].structures.plane.diffuseColor[1] = vector[1];
						objects[items].structures.plane.diffuseColor[2] = vector[2];

					}
					else if(strcmp(objects[items].type, "sphere") == 0) {

                    objects[items].structures.sphere.diffuseColor[0] = vector[0];
                    objects[items].structures.sphere.diffuseColor[1] = vector[1];
                    objects[items].structures.sphere.diffuseColor[2] = vector[2];
					}
                }

            } 
			else if(strcmp(name, "specular_color") == 0) {

				skip_ws(json);
				c = next_c(json);

				if(c != ':'){
					fprintf(stderr, "Error, line number %d; invalid separator '%c'.\n", lineNumber, c);
					fclose(json);
					exit(-1);
				}
				else{
					skip_ws(json);
					vector = next_vec(json);
					if(strcmp(objects[items].type, "plane") == 0) {

						objects[items].structures.plane.specularColor[0] = vector[0];
						objects[items].structures.plane.specularColor[1] = vector[1];
						objects[items].structures.plane.specularColor[2] = vector[2];

					}
					else if(strcmp(objects[items].type, "sphere") == 0) {

                        objects[items].structures.sphere.specularColor[0] = vector[0];
                        objects[items].structures.sphere.specularColor[1] = vector[1];
                        objects[items].structures.sphere.specularColor[2] = vector[2];
					}
                }

            }
			else if(strcmp(name, "position") == 0) {

				skip_ws(json);
				c = next_c(json);

				if(c != ':'){
					fprintf(stderr, "Error, line number %d; invalid separator '%c'.\n", lineNumber, c);
					fclose(json);
					exit(-1);
				}
				else
                {
					skip_ws(json);
					vector = next_vec(json);

					if(strcmp(objects[items].type, "sphere") == 0){
						objects[items].structures.sphere.position[0] = vector[0];
						objects[items].structures.sphere.position[1] = vector[1];
						objects[items].structures.sphere.position[2] = vector[2];

					} else if(strcmp(objects[items].type, "plane") == 0) {
						objects[items].structures.plane.position[0] = vector[0];
						objects[items].structures.plane.position[1] = vector[1];
						objects[items].structures.plane.position[2] = vector[2];

					} else if(strcmp(objects[items].type, "light") == 0){
                        objects[items].structures.light.position[0] = vector[0];
						objects[items].structures.light.position[1] = vector[1];
						objects[items].structures.light.position[2] = vector[2];
					}

				}

			} else if(strcmp(name, "normal") == 0) {

				skip_ws(json);
				c = next_c(json);

				if(c != ':'){
					fprintf(stderr, "Error, line number %d; unexpected character '%c'.\n", lineNumber, c);
					fclose(json);
					exit(-1);
				}
				else
                {
					skip_ws(json);
					vector = next_vec(json);

					objects[items].structures.plane.normal[0] = vector[0];
					objects[items].structures.plane.normal[1] = vector[1];
					objects[items].structures.plane.normal[2] = vector[2];
				}

			} else if(strcmp(name, "direction") == 0) {

				skip_ws(json);
				c = next_c(json);

				if(c != ':')
                {
					fprintf(stderr, "Error, line number %d; unexpected character '%c'.\n", lineNumber, c);
					fclose(json);
					exit(-1);
				}
				else
                {
					skip_ws(json);
					vector = next_vec(json);
					objects[items].structures.light.direction[0] = vector[0];
					objects[items].structures.light.direction[1] = vector[1];
					objects[items].structures.light.direction[2] = vector[2];
				}

			}
			else
            {
				fprintf(stderr, "Error, line number %d; invalid type '%s'.\n", name);
				fclose(json);
				exit(-1);
			}

			skip_ws(json);
			c = next_c(json);

			if(c == ',')
            {
				skip_ws(json);
                c = next_c(json);
			}
		}  

		skip_ws(json);
		c = next_c(json);

		if(c == '{') {
			ungetc(c, json);
		}

		if(c == ','){
			skip_ws(json);
			c = next_c(json);

			if(c == '{'){
				ungetc(c, json);
			}
			}
		items += 1;
	}
// end of read scene
	return items;
}