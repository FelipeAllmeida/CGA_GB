#pragma once

#include <math.h>
#include "Vector3.h"


class Sphere
{
public:	
	Vector3f center() const;
	float radius() const;
	float radius2() const;
	Vector3f surfaceColor() const;
	Vector3f emissionColor() const;
	float transparency() const;
	float reflection() const;
	Sphere(
		const Vector3f &p_center,
		const float &p_radius,
		const Vector3f &p_surfaceColor,
		const float &p_reflection = 0,
		const float &p_transparency = 0,
		const Vector3f &p_emissionColor = (float)0) :
		_center(p_center), _radius(p_radius), _radius2(p_radius * p_radius), _surfaceColor(p_surfaceColor), _emissionColor(p_emissionColor),
		_transparency(p_transparency), _reflection(p_reflection) {}
	//[comment]
	// Compute a ray-sphere intersection using the geometric solution
	//[/comment]
	bool intersect(const Vector3f &rayorig, const Vector3f &raydir, float &t0, float &t1) const;
private:
	Vector3f _center;                          /// position of the sphere
	float _radius, _radius2;                   /// sphere radius and radius^2
	Vector3f _surfaceColor, _emissionColor;    /// surface color and emission (light)
	float _transparency, _reflection;          /// surface transparency and reflectivity

};

