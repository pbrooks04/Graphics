/**********************************************************
	
	Starter code for Assignment 3
	
	Implements scene_object.h

***********************************************************/

#include <cmath>
#include "scene_object.h"

bool UnitSquare::intersect(Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld) {
	// TODO: implement intersection code for UnitSquare, which is
	// defined on the xy-plane, with vertices (0.5, 0.5, 0), 
	// (-0.5, 0.5, 0), (-0.5, -0.5, 0), (0.5, -0.5, 0), and normal
	// (0, 0, 1).
	//
	// Your goal here is to fill ray.intersection with correct values
	// should an intersection occur.  This includes intersection.point, 
	// intersection.normal, intersection.none, intersection.t_value.   
	//
	// HINT: Remember to first transform the ray into object space  
	// to simplify the intersection test.

	if (!ray.intersection.none && ray.intersection.t_value > 0){
		return false;
	}


	// So I guess right now it's in world space...

	Ray3D objectSpaceRay;
	objectSpaceRay.origin = worldToModel * ray.origin;
	objectSpaceRay.dir = worldToModel * ray.dir;
    objectSpaceRay.dir.normalize();

	// If the ray intersects with the sqaure, it must:

	// 1. intersect the x,y plane (i.e. z = 0 )
	// 2. have x in [-0.5, 0.5] and y in [-0.5, 0.5]

	const int x = 0;
	const int y = 1;
	const int z = 2;
	
	if (objectSpaceRay.dir[z] != 0){ // the ray has a z direction
		// I want to find where the ray has a z value of zero
		// So I want (x, y, 0) = ray_origin + t*ray direction
		// t is an unkown so I'll find it using:    0 = origin.z + t * dir.z

		double t_val = -(objectSpaceRay.origin[z] / objectSpaceRay.dir[z]);
		if (t_val > 0){ // If it is less than zero, the intersection was backwards, therefore the ray is moving away from the unit square
			double x_intersection = objectSpaceRay.origin[x] + t_val* objectSpaceRay.dir[x];
			double y_intersection = objectSpaceRay.origin[y] + t_val* objectSpaceRay.dir[y];

			// check that the x & y points are within range
			if (x_intersection <= 0.5 && x_intersection >= -0.5){
				if(y_intersection <= 0.5 && y_intersection >= -0.5){
					// hit

					// still object (model) space
					Point3D intersection_point(x_intersection, y_intersection, 0.0);

					
					ray.intersection.point = modelToWorld*intersection_point;
					ray.intersection.t_value = t_val;
					ray.intersection.normal = transNorm(worldToModel ,Vector3D(0.0, 0.0, 1.0)); // Vector3D: formed from inverse transpose
					ray.intersection.normal.normalize();
					ray.intersection.none = false;

					return true;
				}
			}
		}
	} else { // What if the ray is parallel to the plane?
		// ray must be on the plane, so any point on the ray must be on the plane 
		// so:  normal [dot] (p - p_not) = 0
		Vector3D normal(0.0, 0.0, 1.0);
		Vector3D point_on_ray(objectSpaceRay.origin[x], objectSpaceRay.origin[y], objectSpaceRay.origin[z]);
		if (normal.dot(point_on_ray) == 0){ // on plane -- could also just check if origin[z] == 0
			
			// ray = origin + t*dir must be a soln to 
			// n [dot] (p - p_not) = 0
			//	Point is on the plane if and only if z = 0 
			//	Since p_not = 0
			//	 => n [dot] p = 0
			//	 => 0*p_x + 0*p_y + 1*p_z = 0

			//
			//
			// TODO: finish this (probably silly) case
			//
			// return true;
		}
	}	

	return false;
}

bool UnitSphere::intersect(Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld) {
	// TODO: implement intersection code for UnitSphere, which is centred 
	// on the origin.  
	//
	// Your goal here is to fill ray.intersection with correct values
	// should an intersection occur.  This includes intersection.point, 
	// intersection.normal, intersection.none, intersection.t_value.   
	//
	// HINT: Remember to first transform the ray into object space  
	// to simplify the intersection test.

	if (!ray.intersection.none && ray.intersection.t_value > 0){
		return false;
	}

	Ray3D objectSpaceRay;
	objectSpaceRay.origin = worldToModel * ray.origin;
	objectSpaceRay.dir = worldToModel * ray.dir;
	objectSpaceRay.dir.normalize();

	const int x = 0;
	const int y = 1;
	const int z = 2;

  double t_val = -(objectSpaceRay.origin[z] / objectSpaceRay.dir[z]);
	if (t_val > 0){ // If it is less than zero, the intersection was backwards, therefore the ray is moving away from the unit square
	
  // If the ray intersects with the sphere, it must:

	// 1. be a root of the sphere x^2 + y^2 + z^2 = 1
	// 2. be the first (closer to ray origin) root

	// So to be a root, I need to find t in:
	//			(origin.x + t*dir.x)^2 + (origin.y + t*dir.y)^2 + (origin.z + t*dir.z)^2 = 1
	//  which is equivalent to:
	//			t^2 * (dir.x^2 + dir.y^2 + dir.z^2) + t*( 2*(origin.x*dir.x + origin.y*dir.y + origin.z*dir.z) ) + (origin.x^2 + origin.y^2 + origin.z^2 - 1) = 0
	// then sub into quadratic equation 

	double a = objectSpaceRay.dir[x]*objectSpaceRay.dir[x] + objectSpaceRay.dir[y]*objectSpaceRay.dir[y] + objectSpaceRay.dir[z]*objectSpaceRay.dir[z];
	double b = 2*(objectSpaceRay.origin[x]*objectSpaceRay.dir[x] + objectSpaceRay.origin[y]*objectSpaceRay.dir[y] + objectSpaceRay.origin[z]*objectSpaceRay.dir[z]);
	double c = objectSpaceRay.origin[x]*objectSpaceRay.origin[x] + objectSpaceRay.origin[y]*objectSpaceRay.origin[y] + objectSpaceRay.origin[z]*objectSpaceRay.origin[z] - 1.0;

	double b_squared_minus_4_a_c = b*b - 4*a*c; 

	if (b_squared_minus_4_a_c >= 0 && a != 0){ // check for non negative square root and non zero dividend 
		double t1 = (-b + sqrt(b_squared_minus_4_a_c)) / 2*a;
		double t2 = (-b - sqrt(b_squared_minus_4_a_c)) / 2*a;

		double t_val = std::min(t1, t2); // first hit will have a smaller t value, if they're the same it should just return that
		// Maybe I should also check that t_val > 0 --- no backward intersection

		// I have roots so this is a hit

		double x_intersection = objectSpaceRay.origin[x] + t_val*objectSpaceRay.dir[x];
		double y_intersection = objectSpaceRay.origin[y] + t_val*objectSpaceRay.dir[y];
		double z_intersection = objectSpaceRay.origin[z] + t_val*objectSpaceRay.dir[z];
 
		Point3D intersection_point(x_intersection, y_intersection, z_intersection);

		// I need the normal, which is the vector from the centre of the circle to the intersection point, which I guess is the point normalized...
		// but that should just be the point because the unit circle has radius 1 
		Vector3D intersection_normal(x_intersection, y_intersection, z_intersection);
		intersection_normal.normalize();


		
		ray.intersection.point = modelToWorld*intersection_point;
		ray.intersection.t_value = t_val;
		ray.intersection.normal = transNorm(worldToModel , intersection_normal); // Vector3D: formed from inverse transpose
		ray.intersection.normal.normalize();
		ray.intersection.none = false;

		return true;
	}
  }
	return false;
}

bool UnitCylinder::intersect(Ray3D& ray, const Matrix4x4& worldToModel,
		const Matrix4x4& modelToWorld) {
	// implement intersection code for UnitCylinder, which is centred 
	// on the origin.  
	//
	
	// HINT: Remember to first transform the ray into object space  
	// to simplify the intersection test.


	if (!ray.intersection.none && ray.intersection.t_value > 0){
		return false;
	}

	Ray3D objectSpaceRay;
	objectSpaceRay.origin = worldToModel * ray.origin;
	objectSpaceRay.dir = worldToModel * ray.dir;
	objectSpaceRay.dir.normalize();

	const int x = 0;
	const int y = 1;
	const int z = 2;

	const double z_max = 0.5;
	const double z_min = -0.5;

	// Plot inf cylinder
	//	Ray is point on x,y unit circle

	double a = objectSpaceRay.dir[x]*objectSpaceRay.dir[x] + objectSpaceRay.dir[y]*objectSpaceRay.dir[y];
	double b = 2*(objectSpaceRay.origin[x]*objectSpaceRay.dir[x] + objectSpaceRay.origin[y]*objectSpaceRay.dir[y]);
	double c = objectSpaceRay.origin[x]*objectSpaceRay.origin[x] + objectSpaceRay.origin[y]*objectSpaceRay.origin[y] - 1.0;

	double b_squared_minus_4_a_c = b*b - 4*a*c;
	
	if (b_squared_minus_4_a_c >=0 && a!= 0){
		double t1 = (-b + sqrt(b_squared_minus_4_a_c)) / (2*a);
		double t2 = (-b - sqrt(b_squared_minus_4_a_c)) / (2*a);

		double t_val = std::min(t1, t2);

		double x_intersection = objectSpaceRay.origin[x] + t_val*objectSpaceRay.dir[x];
		double y_intersection = objectSpaceRay.origin[y] + t_val*objectSpaceRay.dir[y];
		double z_intersection = objectSpaceRay.origin[z] + t_val*objectSpaceRay.dir[z];

		if (z_intersection <= z_max && z_intersection >= z_min){
			Point3D intersection_point(x_intersection, y_intersection, z_intersection);
			Vector3D intersection_normal(x_intersection, y_intersection, 0.0);
			intersection_normal.normalize();

			ray.intersection.point = modelToWorld*intersection_point;
			ray.intersection.t_value = t_val;
			ray.intersection.normal = transNorm(worldToModel , intersection_normal); // Vector3D: formed from inverse transpose
			ray.intersection.normal.normalize();
			ray.intersection.none = false;

			return true;
		}
	}
	if (objectSpaceRay.dir[z] != 0){
		double t_top = (z_max - objectSpaceRay.origin[z])/objectSpaceRay.dir[z];
		double t_bottom = (z_min - objectSpaceRay.origin[z])/objectSpaceRay.dir[z];

		double x_top = objectSpaceRay.origin[x] + t_top*objectSpaceRay.dir[x];
		double y_top = objectSpaceRay.origin[y] + t_top*objectSpaceRay.dir[y];

		double x_bottom = objectSpaceRay.origin[x] + t_bottom*objectSpaceRay.dir[x];
		double y_bottom = objectSpaceRay.origin[y] + t_bottom*objectSpaceRay.dir[y];

		bool top_intersection = (x_top*x_top + y_top*y_top) < 1;
		bool bottom_intersection = (x_bottom*x_bottom + y_bottom*y_bottom) < 1;

		if (top_intersection && bottom_intersection){
			double t_val = std::min(t_top, t_bottom);

			Vector3D intersection_normal;
			if (t_top < t_bottom){
				intersection_normal = Vector3D(0, 0, 1);
			} else {
				intersection_normal = Vector3D(0, 0, -1);
			}
			
			double x_intersection = objectSpaceRay.origin[x] + t_val*objectSpaceRay.dir[x];
			double y_intersection = objectSpaceRay.origin[y] + t_val*objectSpaceRay.dir[y];
			double z_intersection = objectSpaceRay.origin[z] + t_val*objectSpaceRay.dir[z];

			Point3D intersection_point(x_intersection, y_intersection, z_intersection);

			ray.intersection.point = modelToWorld*intersection_point;
			ray.intersection.t_value = t_val;
			ray.intersection.normal = transNorm(worldToModel , intersection_normal); 
			ray.intersection.normal.normalize();
			ray.intersection.none = false;

			return true;

		} else if (top_intersection){

			double t_val = t_top;

			double x_intersection = objectSpaceRay.origin[x] + t_val*objectSpaceRay.dir[x];
			double y_intersection = objectSpaceRay.origin[y] + t_val*objectSpaceRay.dir[y];
			double z_intersection = objectSpaceRay.origin[z] + t_val*objectSpaceRay.dir[z];

			Point3D intersection_point(x_intersection, y_intersection, z_intersection);
			Vector3D intersection_normal(0, 0, 1);

			ray.intersection.point = modelToWorld*intersection_point;
			ray.intersection.t_value = t_val;
			ray.intersection.normal = transNorm(worldToModel , intersection_normal); 
			ray.intersection.normal.normalize();
			ray.intersection.none = false;

			return true;

		} else if (bottom_intersection){

			double t_val = t_bottom;

			double x_intersection = objectSpaceRay.origin[x] + t_val*objectSpaceRay.dir[x];
			double y_intersection = objectSpaceRay.origin[y] + t_val*objectSpaceRay.dir[y];
			double z_intersection = objectSpaceRay.origin[z] + t_val*objectSpaceRay.dir[z];

			Point3D intersection_point(x_intersection, y_intersection, z_intersection);
			Vector3D intersection_normal(0, 0, -1);

			ray.intersection.point = modelToWorld*intersection_point;
			ray.intersection.t_value = t_val;
			ray.intersection.normal = transNorm(worldToModel , intersection_normal); 
			ray.intersection.normal.normalize();
			ray.intersection.none = false;

			return true;
		}
	}
	
	return false;
}

void SceneNode::rotate(char axis, double angle) {
	Matrix4x4 rotation;
	double toRadian = 2*M_PI/360.0;
	int i;
	
	for (i = 0; i < 2; i++) {
		switch(axis) {
			case 'x':
				rotation[0][0] = 1;
				rotation[1][1] = cos(angle*toRadian);
				rotation[1][2] = -sin(angle*toRadian);
				rotation[2][1] = sin(angle*toRadian);
				rotation[2][2] = cos(angle*toRadian);
				rotation[3][3] = 1;
			break;
			case 'y':
				rotation[0][0] = cos(angle*toRadian);
				rotation[0][2] = sin(angle*toRadian);
				rotation[1][1] = 1;
				rotation[2][0] = -sin(angle*toRadian);
				rotation[2][2] = cos(angle*toRadian);
				rotation[3][3] = 1;
			break;
			case 'z':
				rotation[0][0] = cos(angle*toRadian);
				rotation[0][1] = -sin(angle*toRadian);
				rotation[1][0] = sin(angle*toRadian);
				rotation[1][1] = cos(angle*toRadian);
				rotation[2][2] = 1;
				rotation[3][3] = 1;
			break;
		}
		if (i == 0) {
			this->trans = this->trans*rotation; 	
			angle = -angle;
		} 
		else {
			this->invtrans = rotation*this->invtrans; 
		}	
	}
}

void SceneNode::translate(Vector3D trans) {
	Matrix4x4 translation;
	
	translation[0][3] = trans[0];
	translation[1][3] = trans[1];
	translation[2][3] = trans[2];
	this->trans = this->trans*translation; 	
	translation[0][3] = -trans[0];
	translation[1][3] = -trans[1];
	translation[2][3] = -trans[2];
	this->invtrans = translation*this->invtrans; 
}

void SceneNode::scale(Point3D origin, double factor[3] ) {
	Matrix4x4 scale;
	
	scale[0][0] = factor[0];
	scale[0][3] = origin[0] - factor[0] * origin[0];
	scale[1][1] = factor[1];
	scale[1][3] = origin[1] - factor[1] * origin[1];
	scale[2][2] = factor[2];
	scale[2][3] = origin[2] - factor[2] * origin[2];
	this->trans = this->trans*scale; 	
	scale[0][0] = 1/factor[0];
	scale[0][3] = origin[0] - 1/factor[0] * origin[0];
	scale[1][1] = 1/factor[1];
	scale[1][3] = origin[1] - 1/factor[1] * origin[1];
	scale[2][2] = 1/factor[2];
	scale[2][3] = origin[2] - 1/factor[2] * origin[2];
	this->invtrans = scale*this->invtrans; 
}


