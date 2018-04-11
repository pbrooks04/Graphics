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
  const double epsilon = 0.01; //used to get rid of noise in hard shadow
	for (size_t  i = 0; i < light_list.size(); ++i) {
   
    //size of area light. Larger val = larger area light
    const int light_size = 5;
    //light definition. Larger val = more definition/smaller space between shadows
    const double light_def = 5;
    Color avg_shadow;
    avg_shadow[0] = 0;
    avg_shadow[1] = 0;
    avg_shadow[2] = 0;

    #pragma omp parallel for
    for (int j=0; j<(2*light_size); j++){
      for (int k=0; k<(2*light_size); k++){
		
        LightSource* light = light_list[i];
          
        //Create Shadow Ray
        Ray3D shadow_ray;
        //Origin is at intersection
        shadow_ray.origin = ray.intersection.point;
        //Direction is towards light
        Point3D light_origin(light->get_position());
        Vector3D light_dir(
          (light_origin[0] - ((light_size - j)/light_def)) - ray.intersection.point[0],
          (light_origin[1] - ((light_size - k)/light_def)) - ray.intersection.point[1],
          (light_origin[2]) - ray.intersection.point[2]
        );
        shadow_ray.dir = light_dir;
        shadow_ray.dir.normalize();
        shadow_ray.origin = ray.intersection.point + (epsilon * shadow_ray.dir);

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
        avg_shadow = avg_shadow + ray.col; 
      }
    }

    ray.col[0] = avg_shadow[0]/pow((light_size*2), 2.0);
    ray.col[1] = avg_shadow[1]/pow((light_size*2), 2.0);
    ray.col[2] = avg_shadow[2]/pow((light_size*2), 2.0);
	}
}

Color Raytracer::shadeRay(Ray3D& ray, Scene& scene, LightList& light_list, int depth) {
	Color col(0.0, 0.0, 0.0); 
	traverseScene(scene, ray); 

	// Don't bother shading if the ray didn't hit 
	// anything.
	if (!ray.intersection.none) {
		if (!ray.intersection.mat->refraction_enabled){

			computeShading(ray, light_list, scene); 
			col = ray.col;  
  
			if (depth>0){
			  Ray3D reflect_ray;
			  reflect_ray.origin = ray.intersection.point;
			  reflect_ray.dir = ray.dir - 2*(ray.dir.dot(ray.intersection.normal))*ray.intersection.normal;
			  reflect_ray.origin = ray.intersection.point + (0.01 * reflect_ray.dir);
			  reflect_ray.dir.normalize();
      
			  depth = depth - 1;
			  col = col + shadeRay(reflect_ray, scene, light_list, depth);
			}

		} else { // Refract ray: no shading 
			
			// col = ray.col;

			if (depth>0){
				Vector3D refracted_ray_dir( refract(ray.dir, ray.intersection.normal, ray.intersection.mat->ior ) );
				Ray3D refracted_ray;
				refracted_ray.origin = ray.intersection.point;
				refracted_ray.dir = refracted_ray_dir;
				refracted_ray.dir.normalize();

				depth = depth - 1;
				col = col + shadeRay(refracted_ray, scene, light_list, depth);
			}
		}
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

	// Divide pixel into nxn grid and fire a ray at a random location in each segment

	

	// Construct a ray for each pixel.
  #pragma omp parallel for
	for (int i = 0; i < image.height; i++) {
		for (int j = 0; j < image.width; j++) {


			//
			// Anti Aliasing Implementation
			//

			// inside each pixel, create nxn grid:
		   
			Color avg_color_for_pixel;
			int num_rays_fired_per_pixel = 9;
			const double one_nth = 1 / sqrt(num_rays_fired_per_pixel);
			const double sq_root = sqrt(num_rays_fired_per_pixel);

			for (double x_boundary_min = 0.000; x_boundary_min  < sq_root*one_nth - 0.001; x_boundary_min += one_nth){ // Do not want it to run with x_bound_min = 0.999, it will be in another pixel
				for( double y_boundary_min = 0.000; y_boundary_min  < sq_root*one_nth - 0.001; y_boundary_min += one_nth){
					// Sets up ray origin and direction in view space, 
					// image plane is at z = -1.
					Point3D origin(0, 0, 0);
					Point3D imagePlane;
					

					double x_rand_double = (double) (std::rand() / (RAND_MAX + 1.0));
					double y_rand_double = (double) (std::rand() / (RAND_MAX + 1.0));


					//imagePlane[0] = (-double(image.width)/2 + 0.5 + j)/factor;
					//imagePlane[1] = (-double(image.height)/2 + 0.5 + i)/factor;

					imagePlane[0] = (-double(image.width)/2 + (x_boundary_min + one_nth*x_rand_double) + j)/factor;
					imagePlane[1] = (-double(image.height)/2 + (y_boundary_min + one_nth*y_rand_double) + i)/factor;
					imagePlane[2] = -1;

			
			
			
					// TODO: Convert ray to world space  
					Vector3D direction = imagePlane - origin; // Will always just be imagePlane since origin is (0,0,0)
					Vector3D ray_direction( viewToWorld * direction);
					Point3D ray_origin(viewToWorld*origin);
					Ray3D ray(ray_origin, ray_direction);
			

				  //Color col = shadeRay(ray, scene, light_list);
					avg_color_for_pixel = avg_color_for_pixel + shadeRay(ray, scene, light_list, 2); // **********************************
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

Vector3D Raytracer::refract(const Vector3D &I, const Vector3D &N, const float &ior)
{
	double cosi = clamp( I.dot(N), -1, 1 );
	float etai = 1, etat = ior;
	Vector3D n = N;
	if (cosi < 0) { 
		cosi = -cosi; 
	} else { 
		std::swap(etai, etat);
		n= -N;
	}
	float eta = etai / etat;
	float k = 1 - eta * eta * (1 - cosi * cosi);

	return k < 0 ? Vector3D(0,0,0) : eta * I + (eta * cosi - std::sqrt(k)) * n;
} 

double Raytracer::clamp(double value, double lo, double hi){
	if (value < lo)
		return lo;
	else if (value > hi)
		return hi;
	return value;
}