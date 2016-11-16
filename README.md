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
