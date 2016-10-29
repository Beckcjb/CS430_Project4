#ifndef OBECTSTRUCTS_H_INCLUDED
#define OBECTSTRUCTS_H_INCLUDED

typedef struct Light
{
    double color[3];
    double position[3];
    double direction[3];
    double radial_a2, radial_a1, radial_a0;
    double angular_a0;
    double theta;

} Light;

typedef struct Camera
{
	double width;
	double height;
} Camera;


typedef struct Plane
{
	double diffuseColor[3];
	double specularColor[3];
    double position[3];
	double normal[3];

} Plane;


typedef struct Sphere
{
	double diffuseColor[3];
	double specularColor[3];
    double position[3];
	double radius;

} Sphere;

typedef struct Object
{
	char *type;

	union structures
	{
		Camera camera;
        Light light;
		Plane plane;
		Sphere sphere;

	} structures;

} Object;

#endif // OBECTSTRUCTS_H_INCLUDED