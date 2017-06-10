#include "Sphere.h"

Vector3f Sphere::center() const
{
	return _center;
}

float Sphere::radius() const
{
	return _radius;
}

float Sphere::radius2() const
{
	return _radius2;
}

Vector3f Sphere::surfaceColor() const
{
	return _surfaceColor;
}

Vector3f Sphere::emissionColor() const
{
	return _emissionColor;
}

float Sphere::transparency() const
{
	return _transparency;
}

float Sphere::reflection() const
{
	return _reflection;
}

bool Sphere::intersect(const Vector3f & rayorig, const Vector3f & raydir, float & t0, float & t1) const
{
	Vector3f l = _center - rayorig;
	float tca = l.dot(raydir);
	if (tca < 0) return false;
	float d2 = l.dot(l) - tca * tca;
	if (d2 > _radius2) return false;
	float thc = sqrt(_radius2 - d2);
	t0 = tca - thc;
	t1 = tca + thc;

	return true;
}


