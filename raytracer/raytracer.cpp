/***********************************************************
	
	Starter code for Assignment 3

	Implementations of functions in raytracer.h, 
	and the main function which specifies the scene to be rendered.	

***********************************************************/


#include "raytracer.h"
#include <cmath>
#include <iostream>
#include <cstdlib>

void Raytracer::traverseScene(Scene& scene, Ray3D& ray)  {
	for (size_t i = 0; i < scene.size(); ++i) {
		SceneNode* node = scene[i];

		if (node->obj->intersect(ray, node->worldToModel, node->modelToWorld)) {
			ray.intersection.mat = node->mat;
		}
	}
}

void Raytracer::computeTransforms(Scene& scene) {
	// right now this method might seem redundant. But if you decide to implement 
	// scene graph this is where you would propagate transformations to child nodes
		
	for (size_t i = 0; i < scene.size(); ++i) {
		SceneNode* node = scene[i];

		node->modelToWorld = node->trans;
		node->worldToModel = node->invtrans; 
	}
}

void Raytracer::computeShading(Ray3D& ray, LightList& light_list, Scene& scene) {
	for (size_t  i = 0; i < light_list.size(); ++i) {
		LightSource* light = light_list[i];
		
    //Create Shadow Ray
    Ray3D shadow_ray;
    //Origin is at intersection
    shadow_ray.origin = ray.intersection.point;
    //Direction is towards light
    Point3D light_origin(light->get_position());
    Vector3D light_dir(
      light_origin[0] - ray.intersection.point[0],
      light_origin[1] - ray.intersection.point[1],
      light_origin[2] - ray.intersection.point[2]
    );
    shadow_ray.dir = light_dir;
    shadow_ray.dir.normalize();
    shadow_ray.origin = ray.intersection.point;

    bool process_ray = true;
    
    //Find shadow intersections
    traverseScene(scene, shadow_ray);

    //Checking if ray met another object
    if(!shadow_ray.intersection.none){
      double light_distance = (light_origin - shadow_ray.intersection.point).length();
      double intsc_distance = (shadow_ray.intersection.point - ray.intersection.point).length();
      //Don't process_ray if it is being blocked by another object
      //Need to check t_value to see which direction, wrt the light, the other object is in
      if(intsc_distance < light_distance && shadow_ray.intersection.t_value >= 0){
        process_ray = false;
        light->shade(ray, 1);
      }
      light->shade(ray, 1);
    }
    
    if(process_ray){
      //default lighting
      light->shade(ray, 0);
    }


		// Each lightSource provides its own shading function.
		// Implement shadows here if needed.
		//light->shade(ray);        
	}
}

Color Raytracer::shadeRay(Ray3D& ray, Scene& scene, LightList& light_list) {
	Color col(0.0, 0.0, 0.0); 
	traverseScene(scene, ray); 

	// Don't bother shading if the ray didn't hit 
	// anything.
	if (!ray.intersection.none) {
		computeShading(ray, light_list, scene); 
		col = ray.col;  
	}

	// You'll want to call shadeRay recursively (with a different ray, 
	// of course) here to implement reflection/refraction effects.  

	return col; 
}	

void Raytracer::render(Camera& camera, Scene& scene, LightList& light_list, Image& image) {
	computeTransforms(scene);

	Matrix4x4 viewToWorld;
	double factor = (double(image.height)/2)/tan(camera.fov*M_PI/360.0);

	viewToWorld = camera.initInvViewMatrix();

	// Create more rays for each pixel here

	// Divide pixel into 3x3 grid and fire a ray at a random location in each segment


	// Construct a ray for each pixel.
	for (int i = 0; i < image.height; i++) {
		for (int j = 0; j < image.width; j++) {


			//
			// Anti Aliasing Implementation
			//

			// inside each pixel, create 3x3 grid:
		    // properties: grid boundaries are separated by a factor of 1/3 = 0.333

			// factor 0.5 puts the ray in the middle, so I want to pick randoms within the following boundaries:
			// ROW 1
			// x (0.000, 0.333) y (0.0, 0.333)
			// x (0.333, 0.666) y (0.0, 0.333)
			// x (0.666, 0.999) y (0.0, 0.333)

			// ROW 2
			// x (0.000, 0.333) y (0.333, 0.666)
			// x (0.333, 0.666) y (0.333, 0.666)
			// x (0.666, 0.999) y (0.333, 0.666)

			// ROW 3
			// x (0.000, 0.333) y (0.666, 0.999)
			// x (0.333, 0.666) y (0.666, 0.999)
			// x (0.666, 0.999) y (0.666, 0.999)

			//const double one_third = 0.333;
			Color avg_color_for_pixel;
			int num_rays_fired_per_pixel = 9;
			const double one_third = 1 / sqrt(num_rays_fired_per_pixel);
      const double sq_root = sqrt(num_rays_fired_per_pixel);

			for (double x_boundary_min = 0.000; x_boundary_min  < sq_root*one_third - 0.001; x_boundary_min += one_third){ // Do not want it to run with x_bound_min = 0.999, it will be in another pixel
				for( double y_boundary_min = 0.000; y_boundary_min  < sq_root*one_third - 0.001; y_boundary_min += one_third){
					// Sets up ray origin and direction in view space, 
					// image plane is at z = -1.
					Point3D origin(0, 0, 0);
					Point3D imagePlane;
					// Modify to be in some segment and obtain a random offset within it
					// Can obtain a random integer up to 1000 and divide it to get a 3 digit double


					//int x_rand_int = std::rand() % 1000;

					//double x_rand_double = (double)std::rand() / RAND_MAX;


					double x_rand_double =  std::abs( (double) (std::rand() / (RAND_MAX + 1.0)) );

					//printf("%d \n", x_rand_double);

					//double y_rand = 


					//imagePlane[0] = (-double(image.width)/2 + 0.5 + j)/factor;
					//imagePlane[1] = (-double(image.height)/2 + 0.5 + i)/factor;

					// Ideally would like to change the added one_third/2 to a random double d :  0 < d < one_third

					imagePlane[0] = (-double(image.width)/2 + (x_boundary_min + one_third/2) + j)/factor;
					imagePlane[1] = (-double(image.height)/2 + (y_boundary_min + one_third/2) + i)/factor;
					imagePlane[2] = -1;

			
			
			
					// TODO: Convert ray to world space  
					Vector3D direction = imagePlane - origin; // Will always just be imagePlane since origin is (0,0,0)
					Vector3D ray_direction( viewToWorld * direction);
					Point3D ray_origin(viewToWorld*origin);
					Ray3D ray(ray_origin, ray_direction);
			

				  //Color col = shadeRay(ray, scene, light_list);
					avg_color_for_pixel = avg_color_for_pixel + shadeRay(ray, scene, light_list);
					// Need to take average of colour from each ray

							
				}
			}
		    avg_color_for_pixel[0] = avg_color_for_pixel[0] / num_rays_fired_per_pixel;
			avg_color_for_pixel[1] = avg_color_for_pixel[1] / num_rays_fired_per_pixel;
			avg_color_for_pixel[2] = avg_color_for_pixel[2] / num_rays_fired_per_pixel;

			avg_color_for_pixel.clamp();
			image.setColorAtPixel(i, j, avg_color_for_pixel);
		}
	}
}

