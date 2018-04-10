/***********************************************************
	
	Starter code for Assignment 3

	Implements light_source.h

***********************************************************/

#include <cmath>
#include "light_source.h"

void PointLight::shade(Ray3D& ray, int mode) {
	// TODO: implement this function to fill in values for ray.col 
	// using phong shading.  Make sure your vectors are normalized, and
	// clamp colour values to 1.0.
	//
	// It is assumed at this point that the intersection information in ray 
	// is available.  So be sure that traverseScene() is called on the ray 
	// before this function.  


	//
	// Scene Signature
	//ray.col =  ray.intersection.mat->diffuse;

	//
	// Phong: diffuse & ambient
	//const int x=0;
	//const int y=1;
	//const int z=2;


	//Vector3D light_direction(ray.intersection.point[x] - this->pos[x], ray.intersection.point[y] - this->pos[y], ray.intersection.point[z] - this->pos[z]); // light source location to intersection point: intersection_point - light_pos
	//light_direction.normalize();
	//double lambert_factor = std::max(ray.intersection.normal.dot(light_direction), 0.0);
	//ray.col =  lambert_factor*ray.intersection.mat->diffuse*this->col_diffuse + ray.intersection.mat->ambient*this->col_ambient;
	//ray.col.clamp();

  const int x=0;
  const int y=1;
  const int z=2;

  //Number of lights per line of area light
  Vector3D light_direction(this->pos[x] - ray.intersection.point[x],
                           this->pos[y] - ray.intersection.point[y],
                           this->pos[z] - ray.intersection.point[z]); // light source location to intersection point: intersection_point - light_pos
  light_direction.normalize();

  // Diffuse
  double lambert_factor = std::max(ray.intersection.normal.dot(light_direction), 0.0);
  Color diffuse_component(lambert_factor * 
                          ray.intersection.mat->diffuse * 
                          this->col_diffuse); 

  // Ambient
  Color ambient_component(ray.intersection.mat->ambient * 
                          this->col_ambient);

  // Specular
  // r = -s + 2*(n [dot] s)*n
  Vector3D reflection( -light_direction + 
                       2 * ray.intersection.normal.dot(light_direction) * 
                       ray.intersection.normal );
  reflection.normalize();
  // Specular value = Intensity(specular)*max(0, r [dot] b)^shiny_exponent
  // where b is point to camera vector.... so the ray normalized and inverted?
  
  Vector3D b_vector(-ray.dir); 
                
  b_vector.normalize();
  double specular_value = pow( std::max(0.0, reflection.dot(b_vector)) , ray.intersection.mat->specular_exp ); 
  Color specular_component(specular_value*ray.intersection.mat->specular*this->col_specular);

  if(mode == 0){
    ray.col = diffuse_component + ambient_component + specular_component;
  } else if (mode == 1){
    ray.col = ambient_component;
    //ray.col[0] = 0.0;
    //ray.col[1] = 1.0;
    //ray.col[2] = 0.0;
  }
  ray.col.clamp();
}

