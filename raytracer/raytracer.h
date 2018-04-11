/***********************************************************
	
	Starter code for Assignment 3
	
	This file contains the interface and datastructures of the raytracer.  
	Simple traversal and addition code to the datastructures are given to you.

***********************************************************/
#pragma once

#include "util.h"
#include "scene_object.h"
#include "light_source.h"

class Raytracer {
public:
	// Renders 3D scene to an image given camera and lights setup.
	void render(Camera& camera, Scene& scene, LightList& light_list, Image& image);
		
private:

	// Return the color of the ray after intersection and shading, call 
	// this function recursively for reflection and refraction.  
	Color shadeRay(Ray3D& ray, Scene& scene, LightList& light_list,
                          int depth,
                          bool env_map_loaded, unsigned long int widthBMP, long int heightBMP, unsigned char *rbuffer, unsigned char *gbuffer, unsigned char *bbuffer);

	// Traversal code for the scene, the ray is transformed into 
	// the object space of each node where intersection is performed.
	void traverseScene(Scene& scene, Ray3D& ray);

	// After intersection, calculate the color of the ray by shading it
	// with all light sources in the scene.
	void computeShading(Ray3D& ray, LightList& light_list, Scene& scene);

	// Precompute the modelToWorld and worldToModel transformations for each
    // object in the scene.
	void computeTransforms(Scene& scene);

	Vector3D refract(const Vector3D &I, const Vector3D &N, const float &ior);

	double clamp(double value, double lo, double hi);

	Vector3D convert_xyz_to_cube_uv(double x, double y, double z);

	Color getEnvMapPixelColour(Vector3D dir, unsigned long int width, long int height, unsigned char *rbuffer, unsigned char *gbuffer, unsigned char *bbuffer);
};
