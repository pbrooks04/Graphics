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

Color Raytracer::shadeRay(Ray3D& ray, Scene& scene, LightList& light_list,
                          int depth,
                          bool env_map_loaded, unsigned long int widthBMP, long int heightBMP, unsigned char *rbuffer, unsigned char *gbuffer, unsigned char *bbuffer) {

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
			  col = col + shadeRay(reflect_ray, scene, light_list, depth, 
			                       env_map_loaded, widthBMP, heightBMP, rbuffer, gbuffer, bbuffer);
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
				col = col + shadeRay(refracted_ray, scene, light_list, depth,
			                         env_map_loaded, widthBMP, heightBMP, rbuffer, gbuffer, bbuffer);
			}
		}
	} else if (env_map_loaded){
		col = getEnvMapPixelColour(ray.dir, widthBMP, heightBMP, rbuffer, gbuffer, bbuffer);
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

	//
	// ***** Environment Mapping *****
	// Load image for missed rays
	//
	//
	// Image source https://raptor.developpez.com/tutorial/opengl/skybox/
	//
	char fileName[] = "cubemaps/sky.bmp"; 
	unsigned long int widthBMP;
	long int heightBMP;
	unsigned char *rbuffer;
	unsigned char *gbuffer;
	unsigned char *bbuffer;
	bool image_read_failed;

	image_read_failed = bmp_read(fileName, &widthBMP, &heightBMP, &rbuffer, &gbuffer, &bbuffer);

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
					avg_color_for_pixel = avg_color_for_pixel + shadeRay(ray, scene, light_list, 4,
																		 !image_read_failed, widthBMP, heightBMP, rbuffer, gbuffer, bbuffer); // *****Env Mapping*****
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
//
//
// Implented under the guidance of https://www.scratchapixel.com/lessons/3d-basic-rendering/introduction-to-shading/reflection-refraction-fresnel
//
//
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
//
// From recommended wiki page https://en.wikipedia.org/wiki/Cube_mapping
//
Vector3D Raytracer::convert_xyz_to_cube_uv(double x, double y, double z)
{
  float absX = std::abs(x);
  float absY = std::abs(y);
  float absZ = std::abs(z);
  
  int isXPositive = x > 0 ? 1 : 0;
  int isYPositive = y > 0 ? 1 : 0;
  int isZPositive = z > 0 ? 1 : 0;
  
  float maxAxis, uc, vc;

  double face_index, u, v;
  
  // POSITIVE X
  if (isXPositive && absX >= absY && absX >= absZ) {
    // u (0 to 1) goes from +z to -z
    // v (0 to 1) goes from -y to +y
    maxAxis = absX;
    uc = -z;
    vc = y;
    face_index = 0;
  }
  // NEGATIVE X
  if (!isXPositive && absX >= absY && absX >= absZ) {
    // u (0 to 1) goes from -z to +z
    // v (0 to 1) goes from -y to +y
    maxAxis = absX;
    uc = z;
    vc = y;
    face_index = 1;
  }
  // POSITIVE Y
  if (isYPositive && absY >= absX && absY >= absZ) {
    // u (0 to 1) goes from -x to +x
    // v (0 to 1) goes from +z to -z
    maxAxis = absY;
    uc = x;
    vc = -z;
    face_index = 2;
  }
  // NEGATIVE Y
  if (!isYPositive && absY >= absX && absY >= absZ) {
    // u (0 to 1) goes from -x to +x
    // v (0 to 1) goes from -z to +z
    maxAxis = absY;
    uc = x;
    vc = z;
    face_index = 3;
  }
  // POSITIVE Z
  if (isZPositive && absZ >= absX && absZ >= absY) {
    // u (0 to 1) goes from -x to +x
    // v (0 to 1) goes from -y to +y
    maxAxis = absZ;
    uc = x;
    vc = y;
    face_index = 4;
  }
  // NEGATIVE Z
  if (!isZPositive && absZ >= absX && absZ >= absY) {
    // u (0 to 1) goes from +x to -x
    // v (0 to 1) goes from -y to +y
    maxAxis = absZ;
    uc = -x;
    vc = y;
    face_index = 5;
  }

  // Convert range from -1 to 1 to 0 to 1
  u = 0.5f * (uc / maxAxis + 1.0f);
  v = 0.5f * (vc / maxAxis + 1.0f);

  return Vector3D(face_index, u, v);
}

Color Raytracer::getEnvMapPixelColour(Vector3D dir, unsigned long int width, long int height, unsigned char *rbuffer, unsigned char *gbuffer, unsigned char *bbuffer){

	// Find corresponding uv corrdinates and face index
	Vector3D values = convert_xyz_to_cube_uv(dir[0], dir[1], dir[2]);
	
	int face_index = (int) values[0];
	float u = values[1];
	float v = values[2];

	int face_origin_x, face_origin_y;

	// Assume image origin (0, 0) is at top left corner
	//	Set face origin to bottom left corner of face
	switch(face_index){
		case 0:
			face_origin_x = (int) (width / 2);
			face_origin_y = (int) (height*2/3);
			break;
		case 1:
			face_origin_x = 0;
			face_origin_y = (int) (height*2/3);
			break;
		case 2:
			face_origin_x = (int) (width/4);
			face_origin_y = (int) (height/3);
			break;
		case 3:
			face_origin_x = (int) (width/4);
			face_origin_y = (int) height;
			break;
		case 4:
			face_origin_x = (int) (width/4);
			face_origin_y = (int) (height*2/3);
			break;
		case 5:
			face_origin_x = (int) (width*3/4);
			face_origin_y = (int) (height*2/3);
			break;
	}

	int pixel_x = clamp(face_origin_x + (int) (u * width/4), 0 , width-1);
	int pixel_y = clamp(face_origin_y - (int) (v * height/3), 0, std::abs(height)-1); // Negative y direction for image origin at top left

	/*
      Reminder
	  // Set color of pixel (i,j) to col
		void setColorAtPixel(int i, int j, Color& col) {
			rbuffer[i*width+j] = int(col[0]*255);
			gbuffer[i*width+j] = int(col[1]*255);
			bbuffer[i*width+j] = int(col[2]*255);
		}
	 */


	 //
	 // ***** Environment Mapping *****
	 //   Partial Implementation
	 //
	 //		This is where the code fails, it will run for some time then give a segfault for the code below.
	 //		I don't understand why since x and y values are bounded to the size of the image and values are 
	 //		accessed just as in the starter code.
	 //


	/*
	Color env_pixel_colour( rbuffer[pixel_x*width+pixel_y] / 255,
							gbuffer[pixel_x*width+pixel_y] / 255,
							bbuffer[pixel_x*width+pixel_y] / 255);
	*/

	Color env_pixel_colour(0, 0, 0); 
	return env_pixel_colour;
}