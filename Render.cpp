#include "Render.h"



Render::Render()
{
}


Render::~Render()
{
}

//[comment]
// This is the main trace function. It takes a ray as argument (defined by its origin
// and direction). We test if this ray intersects any of the geometry in the scene.
// If the ray intersects an object, we compute the intersection point, the normal
// at the intersection point, and shade this point using this information.
// Shading depends on the surface property (is it transparent, reflective, diffuse).
// The function returns a color for the ray. If the ray intersects an object that
// is the color of the object at the intersection point, otherwise it returns
// the background color.
//[/comment]

float Render::Mix(const float & p_a, const float & p_b, const float & p_mix)
{
	return p_b * p_mix + p_a * (1 - p_mix);
}

Vector3f Render::Trace(const Vector3f &p_rayOrigin, const Vector3f &p_rayDirection,	const std::vector<Sphere> &p_spheres, const int &p_depth)
{
	//if (raydir.length() != 1) std::cerr << "Error " << raydir << std::endl;
	float __tNear = INFINITY;
	const Sphere* __sphere = NULL;
	// find intersection of this ray with the sphere in the scene
	for (unsigned i = 0; i < p_spheres.size(); ++i) 
	{
		float t0 = INFINITY, t1 = INFINITY;
		if (p_spheres[i].intersect(p_rayOrigin, p_rayDirection, t0, t1)) 
		{
			if (t0 < 0) t0 = t1;
			if (t0 < __tNear) 
			{
				__tNear = t0;
				__sphere = &p_spheres[i];
			}
		}
	}
	// if there's no intersection return black or background color
	if (!__sphere) return Vector3f(1);

	Vector3f __surfaceColor = 0; // color of the ray/surfaceof the object intersected by the ray
	Vector3f __phit = p_rayOrigin + p_rayDirection * __tNear; // point of intersection
	Vector3f __nhit = __phit - __sphere->center(); // normal at the intersection point
	__nhit.normalize(); // normalize normal direction
					  // If the normal and the view direction are not opposite to each other
					  // reverse the normal direction. That also means we are inside the sphere so set
					  // the inside bool to true. Finally reverse the sign of IdotN which we want
					  // positive.
	float __bias = 1e-4; // add some bias to the point from which we will be tracing
	bool __inside = false;

	if (p_rayDirection.dot(__nhit) > 0) __nhit = -__nhit, __inside = true;

	if ((__sphere->transparency() > 0 || __sphere->reflection() > 0) && p_depth < MAX_RAY_DEPTH)
	{
		float __facingratio = -p_rayDirection.dot(__nhit);
		// change the mix value to tweak the effect
		float __fresneleffect = Mix(pow(1 - __facingratio, 3), 1, 0.1);
		// compute reflection direction (not need to normalize because all vector
		// are already normalized)
		Vector3f __reflDirection = p_rayDirection - __nhit * 2 * p_rayDirection.dot(__nhit);
		__reflDirection.normalize();
		Vector3f __reflection = Trace(__phit + __nhit * __bias, __reflDirection, p_spheres, p_depth + 1);
		Vector3f __refraction = 0;

		// if the sphere is also transparent compute refraction ray (transmission)

		if (__sphere->transparency()) 
		{
			float ior = 1.1, eta = (__inside) ? ior : 1 / ior; // are we inside or outside the surface?
			float cosi = -__nhit.dot(p_rayDirection);
			float k = 1 - eta * eta * (1 - cosi * cosi);
			Vector3f refrdir = p_rayDirection * eta + __nhit * (eta *  cosi - sqrt(k));
			refrdir.normalize();
			__refraction = Trace(__phit - __nhit * __bias, refrdir, p_spheres, p_depth + 1);
		}
		// the result is a mix of reflection and refraction (if the sphere is transparent)
		__surfaceColor = (__reflection * __fresneleffect + __refraction * (1 - __fresneleffect) * __sphere->transparency()) * __sphere->surfaceColor();


	}
	else 
	{
		// it's a diffuse object, no need to raytrace any further
		for (unsigned i = 0; i < p_spheres.size(); i++) 
		{
			if (p_spheres[i].emissionColor().x > 0)
			{
				// this is a light
				Vector3f __transmission = 1;
				Vector3f __lightDirection = p_spheres[i].center() - __phit;
				__lightDirection.normalize();
				for (unsigned j = 0; j < p_spheres.size(); j++)
				{
					if (i != j)
					{
						float __t0, __t1;
						if (p_spheres[j].intersect(__phit + __nhit * __bias, __lightDirection, __t0, __t1))
						{
							__transmission = 0;
							break;
						}
					}
				}
				__surfaceColor += __sphere->surfaceColor() * __transmission * max(float(0), __nhit.dot(__lightDirection)) * p_spheres[i].emissionColor();
			}
		}
	}

	return __surfaceColor + __sphere->emissionColor();
}

//[comment]
// Main rendering function. We compute a camera ray for each pixel of the image
// trace it and return a color. If the ray hits a sphere, we return the color of the
// sphere at the intersection point, else we return the background color.
//[/comment]

void Render::RenderScene(const std::vector<Sphere> &p_spheres, Vector3f& p_ptrImage, unsigned  p_width, unsigned  p_height)
{
	Vector3f __cameraPosition = Vector3f(0.f, 0.f, 20.f);

	p_ptrImage = *new Vector3f[p_width * p_height];
	Vector3f *_pixel = &p_ptrImage;
	//p_ptrImage = *__ptrImage;
	float __invWidth = 1 / float(p_width), invHeight = 1 / float(p_height);
	float __fieldOfView = 30.0, _aspectRatio = float(p_width) / float(p_height);
	float __angle = tan(M_PI * 0.5 * __fieldOfView / 180.);

	// Trace rays
	for (unsigned y = 0; y < p_height; ++y)
	{
		for (unsigned x = 0; x < p_width; ++x, ++_pixel)
		{
			float __newRayPosX = (2.0 * ((x + 0.25) * __invWidth) - 1.0) * __angle * _aspectRatio;
			float __newRayPosY = (1.0 - 2.0 * ((y + 0.25) * invHeight)) * __angle;
			Vector3f __rayDirection(__newRayPosX, __newRayPosY, -1);
			__rayDirection.normalize();
			*_pixel = Trace(__cameraPosition, __rayDirection, p_spheres, 0);
		}
	}
}

//void RenderScenePTR(const std::vector<Sphere> &spheres, Vec3f& p_ptrImage, unsigned  p_width, unsigned  p_height)
//{
//	Vec3f __cameraPosition = Vec3f(0.f, 0.f, 20.f);
//
//	p_ptrImage = *new Vec3f[p_width * p_height];
//	Vec3f *_pixel = &p_ptrImage;
//	//p_ptrImage = *__ptrImage;
//	float __invWidth = 1 / float(p_width), invHeight = 1 / float(p_height);
//	float __fieldOfView = 30.0, _aspectRatio = float(p_width) / float(p_height);
//	float __angle = tan(M_PI * 0.5 * __fieldOfView / 180.);
//
//	// Trace rays
//	for (unsigned y = 0; y < p_height; ++y)
//	{
//		for (unsigned x = 0; x < p_width; ++x, ++_pixel)
//		{
//			float __newRayPosX = (2.0 * ((x + 0.25) * __invWidth) - 1.0) * __angle * _aspectRatio;
//			float __newRayPosY = (1.0 - 2.0 * ((y + 0.25) * invHeight)) * __angle;
//			Vec3f __rayDirection(__newRayPosX, __newRayPosY, -1);
//			__rayDirection.normalize();
//			*_pixel = trace(__cameraPosition, __rayDirection, spheres, 0);
//		}
//	}
//}

