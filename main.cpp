// Editado por Felipe Rodrigues de Almeida
// [header]
// A very basic raytracer example.
// [/header]
// [compile]
// c++ -o raytracer -O3 -Wall raytracer.cpp
// [/compile]
// [ignore]
// Copyright (C) 2012 www.scratchapixel.com
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.
// [/ignore]
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <vector>
#include <iostream>
#include <cassert>
#include <cstdlib>
#include <string>
#include <math.h>

#include <sstream>
#include <chrono>

#include <random> 
#if defined __linux__ || defined __APPLE__
// "Compiled for Linux
#else
// Windows doesn't define these values by default, Linux does
#define M_PI 3.141592653589793
#define INFINITY 1e8
#endif
template<typename T>
class Vec3
{
public:
	T x, y, z;
	Vec3() : x(T(0)), y(T(0)), z(T(0)) {}
	Vec3(T xx) : x(xx), y(xx), z(xx) {}
	Vec3(T xx, T yy, T zz) : x(xx), y(yy), z(zz) {}
	Vec3& normalize()
	{
		T nor2 = length2();
		if (nor2 > 0) {
			T invNor = 1 / sqrt(nor2);
			x *= invNor, y *= invNor, z *= invNor;
		}
		return *this;
	}
	Vec3<T> operator * (const T &f) const { return Vec3<T>(x * f, y * f, z * f); }
	Vec3<T> operator * (const Vec3<T> &v) const { return Vec3<T>(x * v.x, y * v.y, z * v.z); }
	T dot(const Vec3<T> &v) const { return x * v.x + y * v.y + z * v.z; }
	Vec3<T> operator - (const Vec3<T> &v) const { return Vec3<T>(x - v.x, y - v.y, z - v.z); }
	Vec3<T> operator + (const Vec3<T> &v) const { return Vec3<T>(x + v.x, y + v.y, z + v.z); }
	Vec3<T>& operator += (const Vec3<T> &v) { x += v.x, y += v.y, z += v.z; return *this; }
	Vec3<T>& operator *= (const Vec3<T> &v) { x *= v.x, y *= v.y, z *= v.z; return *this; }
	Vec3<T> operator - () const { return Vec3<T>(-x, -y, -z); }
	T length2() const { return x * x + y * y + z * z; }
	T length() const { return sqrt(length2()); }
	friend std::ostream & operator << (std::ostream &os, const Vec3<T> &v)
	{
		os << "[" << v.x << " " << v.y << " " << v.z << "]";
		return os;
	}
};
typedef Vec3<float> Vec3f;
class Sphere
{
public:
	Vec3f center; /// position of the sphere
	float radius, radius2; /// sphere radius and radius^2
	Vec3f surfaceColor, emissionColor; /// surface color and emission (light)
	float transparency, reflection; /// surface transparency and reflectivity
	Sphere(const Vec3f &p_center, const float &p_radius, const Vec3f &p_surfaceColor, const float &p_reflection = 0, const float &p_alpha = 0, const Vec3f &p_emissionColor = 0) :
		center(p_center), radius(p_radius), radius2(p_radius * p_radius), surfaceColor(p_surfaceColor), emissionColor(p_emissionColor), transparency(p_alpha), reflection(p_reflection)
	{ /* empty */
	}
	//[comment]
	// Compute a ray-sphere intersection using the geometric solution
	//[/comment]
	bool Intersect(const Vec3f &p_rayOrigin, const Vec3f &p_rayDirection, float &t0, float &t1) const
	{
		Vec3f l = center - p_rayOrigin;
		float tca = l.dot(p_rayDirection);
		if (tca < 0) return false;
		float d2 = l.dot(l) - tca * tca;
		if (d2 > radius2) return false;
		float thc = sqrt(radius2 - d2);
		t0 = tca - thc;
		t1 = tca + thc;
		return true;
	}
	//Pra luz
	/*
	bool intersect(const Vec3f &orig, const Vec3f &dir, float &tNear, uint32_t &triIndex, Vec2f &uv) const
	{
		float t0, t1; // solutions for t if the ray intersects
					  // analytic solution
		Vec3f L = orig - center;
		float a = dir.dotProduct(dir);
		float b = 2 * dir.dotProduct(L);
		float c = L.dotProduct(L) - radius2;
		if (!solveQuadratic(a, b, c, t0, t1)) return false;

		if (t0 > t1) std::swap(t0, t1);

		if (t0 < 0) {
			t0 = t1; // if t0 is negative, let's use t1 instead
			if (t0 < 0) return false; // both t0 and t1 are negative
		}

		tNear = t0;

		return true;
	}
	*/
};
//[comment]
// This variable controls the maximum recursion depth
//[/comment]
#define MAX_RAY_DEPTH 5
float mix(const float &a, const float &b, const float &mix)
{
	return b * mix + a * (1 - mix);
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
Vec3f Trace(const Vec3f &p_ptrRayOrigin, const Vec3f &p_ptrRayDirection, const std::vector<Sphere> &p_spheres, const int &p_depth)
{
	//if (raydir.length() != 1) std::cerr << "Error " << raydir << std::endl;
	float __near = INFINITY;
	const Sphere* __sphere = NULL;

	// find intersection of this ray with the sphere in the scene
	for (unsigned i = 0; i < p_spheres.size(); ++i)
	{
		float __firstTouch = INFINITY, __secondTouch = INFINITY;
		if (p_spheres[i].Intersect(p_ptrRayOrigin, p_ptrRayDirection, __firstTouch, __secondTouch) == true) {
			if (__firstTouch < 0) __firstTouch = __secondTouch;
			if (__firstTouch < __near)
			{
				__near = __firstTouch;
				__sphere = &p_spheres[i];
			}
		}
	}
	// if there's no intersection return black or background color
	if (!__sphere) return Vec3f(1.0);//Vec3f(0.0, 0.49, 0.74);

	Vec3f surfaceColor = 0; // color of the ray/surfaceof the object intersected by the ray
	Vec3f phit = p_ptrRayOrigin + p_ptrRayDirection * __near; // point of intersection
	Vec3f nhit = phit - __sphere->center; // normal at the intersection point
	nhit.normalize(); // normalize normal direction
					  // If the normal and the view direction are not opposite to each other
					  // reverse the normal direction. That also means we are inside the sphere so set
					  // the inside bool to true. Finally reverse the sign of IdotN which we want
					  // positive.
	float bias = 1e-4; // add some bias to the point from which we will be tracing
	bool inside = false;
	if (p_ptrRayDirection.dot(nhit) > 0) nhit = -nhit, inside = true;
	if ((__sphere->transparency > 0 || __sphere->reflection > 0) && p_depth < MAX_RAY_DEPTH)
	{
		float facingratio = -p_ptrRayDirection.dot(nhit);
		// change the mix value to tweak the effect
		float fresneleffect = mix(pow(1 - facingratio, 3), 1, 0.1);
		// compute reflection direction (not need to normalize because all vectors
		// are already normalized)
		Vec3f refldir = p_ptrRayDirection - nhit * 2 * p_ptrRayDirection.dot(nhit);
		refldir.normalize();
		Vec3f reflection = Trace(phit + nhit * bias, refldir, p_spheres, p_depth + 1);
		Vec3f refraction = 0;
		// if the sphere is also transparent compute refraction ray (transmission)
		if (__sphere->transparency)
		{
			float ior = 1.1, eta = (inside) ? ior : 1 / ior; // are we inside or outside the surface?
			float cosi = -nhit.dot(p_ptrRayDirection);
			float k = 1 - eta * eta * (1 - cosi * cosi);
			Vec3f refrdir = p_ptrRayDirection * eta + nhit * (eta * cosi - sqrt(k));
			refrdir.normalize();
			refraction = Trace(phit - nhit * bias, refrdir, p_spheres, p_depth + 1);
		}
		// the result is a mix of reflection and refraction (if the sphere is transparent)
		surfaceColor = (
			reflection * fresneleffect +
			refraction * (1 - fresneleffect) * __sphere->transparency) * __sphere->surfaceColor;
	}
	else
	{
		// it's a diffuse object, no need to raytrace any further
		for (unsigned i = 0; i < p_spheres.size(); ++i)
		{
			if (p_spheres[i].emissionColor.x > 0)
			{
				// this is a light
				Vec3f transmission = 1;
				Vec3f lightDirection = p_spheres[i].center - phit;
				lightDirection.normalize();
				for (unsigned j = 0; j < p_spheres.size(); ++j)
				{
					if (i != j) {
						float t0, t1;
						if (p_spheres[j].Intersect(phit + nhit * bias, lightDirection, t0, t1))
						{
							transmission = 0;
							break;
						}
					}
				}
				surfaceColor += __sphere->surfaceColor * transmission * std::fmaxf(float(0), nhit.dot(lightDirection)) * p_spheres[i].emissionColor;
			}
		}
	}
	Vec3f __result = surfaceColor + __sphere->emissionColor;
	return __result;
}
//[comment]
// Main rendering function. We compute a camera ray for each pixel of the image
// trace it and return a color. If the ray hits a sphere, we return the color of the
// sphere at the intersection point, else we return the background color.
//[/comment]
float ClampValue(float p_value, float p_min, float p_max)
{
	if (p_value > p_max)
	{
		return p_max;
	}

	if (p_value < p_min)
	{
		return p_min;
	}

	return p_value;
}

float DepthBuffer(const Vec3f &p_ptrRayOrigin, const Vec3f &p_ptrRayDirection, const std::vector<Sphere> &p_spheres)
{
	float __rayLenght = INFINITY;
	const Sphere* __sphere = NULL;

	for (unsigned i = 0; i < p_spheres.size(); ++i)
	{
		float __firstTouch = INFINITY, __secondTouch = INFINITY;
		if (p_spheres[i].Intersect(p_ptrRayOrigin, p_ptrRayDirection, __firstTouch, __secondTouch) == true) {
			if (__firstTouch < 0) __firstTouch = __secondTouch;
			if (__firstTouch < __rayLenght)
			{
				__rayLenght = __firstTouch;
				__sphere = &p_spheres[i];
			}
		}
	}
	// horizontal � vermelho
	// vertical � verde
	// profundidade � azul
	if (!__sphere) return 0;

	float __renderPosLimitZ = -10;
	float __depth = 0;

	if (__rayLenght < (p_ptrRayOrigin.z - __renderPosLimitZ))
	{
		__depth = 1 - ((abs(p_ptrRayDirection.z * __rayLenght) + __renderPosLimitZ) / (p_ptrRayOrigin.z - __renderPosLimitZ));
	}

	return __depth;
}

Vec3f NormalMapBuffer(const Vec3f &p_ptrRayOrigin, const Vec3f &p_ptrRayDirection, const std::vector<Sphere> &p_spheres)
{
	float __rayLenght = INFINITY;
	const Sphere* __sphere = NULL;

	for (unsigned i = 0; i < p_spheres.size(); ++i)
	{
		float __firstTouch = INFINITY, __secondTouch = INFINITY;
		if (p_spheres[i].Intersect(p_ptrRayOrigin, p_ptrRayDirection, __firstTouch, __secondTouch) == true) {
			if (__firstTouch < 0) __firstTouch = __secondTouch;
			if (__firstTouch < __rayLenght)
			{
				__rayLenght = __firstTouch;
				__sphere = &p_spheres[i];
			}
		}
	}

	if (!__sphere) return Vec3f(0);

	Vec3f surfaceColor = 0; // color of the ray/surfaceof the object intersected by the ray
	Vec3f phit = p_ptrRayOrigin + p_ptrRayDirection * __rayLenght; // point of intersection
	Vec3f nhit = phit - __sphere->center; // normal at the intersection point
	nhit.normalize(); // normalize normal direction


	float __colorR = (nhit.x + 1) / 2;
	float __colorG = (nhit.y + 1) / 2;
	float __colorB = (nhit.z + 1) * 0.5; 

	Vec3f __surfaceColor = Vec3f(__colorR, __colorG, __colorB);
	return __surfaceColor;
}
// targetBuffer 
// 0 = normal
// 1 = depth
// 2 = normal
void RenderScenePTR(const std::vector<Sphere> &spheres, Vec3f& p_ptrImage, unsigned p_width, unsigned p_height, int p_targetBuffer)
{
	Vec3f __cameraPosition = Vec3f(0.f, 0.f, 20.f);
	p_ptrImage = *new Vec3f[p_width * p_height];
	Vec3f *__pixelColor = &p_ptrImage;
	//p_ptrImage = *__ptrImage;
	float __invWidth = 1 / float(p_width), invHeight = 1 / float(p_height);
	float __fieldOfView = 30.0, _aspectRatio = float(p_width) / float(p_height);
	float __angle = tan(M_PI * 0.5 * __fieldOfView / 180.);
	// Trace rays

	float __lastDistanceHorizontal = 0;
	float __lastDistanceVertical = 0;

	float *lastPixelDistances = new float[p_width];

	for (unsigned y = 0; y < p_height; ++y)
	{
		float *currentPixelDistances = new float[p_width];

		for (unsigned x = 0; x < p_width; ++x, ++__pixelColor)
		{
			float __newRayPosX = (2.0 * ((x + 0.25) * __invWidth) - 1.0) * __angle * _aspectRatio;
			float __newRayPosY = (1.0 - 2.0 * ((y + 0.25) * invHeight)) * __angle;
			Vec3f __rayDirection(__newRayPosX, __newRayPosY, -1);
			__rayDirection.normalize();

			if (p_targetBuffer == 0)
			{
				*__pixelColor = Trace(__cameraPosition, __rayDirection, spheres, 0);				
			}
			else if (p_targetBuffer == 1)
			{
				/*float __currentDistance =*/ *__pixelColor = DepthBuffer(__cameraPosition, __rayDirection, spheres);
				/*float __maxDistanceDifference = 0.25f;
				if (abs(__currentDistance - __lastDistanceHorizontal) > __maxDistanceDifference)
				{
					*__pixelColor = Vec3f(1);
				}
				else
				{
					if (y > 0)
					{
						float lastDistanceUpper = lastPixelDistances[x];
						if (abs(__currentDistance - lastDistanceUpper) > __maxDistanceDifference)
						{
							*__pixelColor = Vec3f(1);
						}
						else
						{
							*__pixelColor = Vec3f(0);
						}
					}
					else
					{
						*__pixelColor = Vec3f(0);

					}
				}

				currentPixelDistances[x] = __currentDistance;
				__lastDistanceHorizontal = __currentDistance;	*/		
			}
			else if (p_targetBuffer == 2)
			{
				*__pixelColor = NormalMapBuffer(__cameraPosition, __rayDirection, spheres);
			}
		}
		lastPixelDistances = currentPixelDistances;
	}
}
#define filterMeanWidth 3
#define filterMeanHeight 3
double filterMean[filterMeanHeight][filterMeanWidth] =
{
	1, 1, 1,
	1, 1, 1,
	1, 1, 1
};
double factorMean = 1.0 / 9.0;
double biasMean = 0.0;
#define filterBlurWidth 9
#define filterBlurHeight 9
double filterBlur[filterBlurHeight][filterBlurWidth] =
{
	1, 0, 0, 0, 0, 0, 0, 0, 0,
	0, 1, 0, 0, 0, 0, 0, 0, 0,
	0, 0, 1, 0, 0, 0, 0, 0, 0,
	0, 0, 0, 1, 0, 0, 0, 0, 0,
	0, 0, 0, 0, 1, 0, 0, 0, 0,
	0, 0, 0, 0, 0, 1, 0, 0, 0,
	0, 0, 0, 0, 0, 0, 1, 0, 0,
	0, 0, 0, 0, 0, 0, 0, 1, 0,
	0, 0, 0, 0, 0, 0, 0, 0, 1,
};
double blurFactor = 1.0 / 9.0;
double blurBias = 0.0;
#define filterFindEdgesWidth 9
#define filterFindEdgesHeight 9
double filterFindEdges[filterFindEdgesHeight][filterFindEdgesWidth] =
{
	0, 0, 3, 2, 2, 2, 3, 0, 0,
	0, 2, 3, 5, 5, 5, 3, 2, 0,
	3, 3, 5, 3, 0, 3, 5, 3, 3,
	2, 5, 3,-12,-23, -12, 3, 5, 2,
	2, 5, 0,-23,-40,-23, 0, 5, 2,
	2, 5, 3,-12,-23,-12, 3, 5, 2,
	3, 3, 5, 3, 0, 3, 5, 3, 3,
	0, 2, 3, 5, 5, 5, 3, 2, 0,
	0, 0, 3, 2, 2, 2, 3, 0, 0
};
double factorFindEdges = 1.0 / 36.0;
double biasFindEdges = 0.0;

float clamp(float p_value, float p_a, float p_b)
{
	if (p_value < p_a)
	{
		return p_a;
	}

	if (p_value > p_b)
	{
		return p_b;
	}

	return p_value;
}

float min(float p_left, float p_right)
{
	if (p_left < p_right)
	{
		return p_left;
	}
	return p_right;
}
float max(float p_left, float p_right)
{
	if (p_left > p_right)
	{
		return p_left;
	}
	return p_right;
}
template <size_t rows, size_t cols>
void FilterImage(Vec3f* p_ptrImage, int p_width, int p_height, int p_filterWidth, int p_filterHeight, double(&array)[rows][cols], double p_factor, double p_bias)
{
	const int __filterHeight = p_filterWidth;
	const int __filterWidth = p_filterHeight;
	//std::vector<ColorRGB> image;
	//for (int i = 0; i < p_width; i++)
	//{
	// for (int j = 0; j < p_height; j++)
	// {
	// if (i == p_width-1 && j == p_height-1)
	// {
	// std::cout << "Color x: " << p_ptrImage[i * p_width + j].x << std::endl;
	// std::cout << "Color y: " << p_ptrImage[i * p_width + j].y << std::endl;
	// std::cout << "Color z: " << p_ptrImage[i * p_width + j].z << std::endl;
	// }
	// p_ptrImage[j * p_width + i].x *= 255;
	// p_ptrImage[j * p_width + i].y *= 255;
	// p_ptrImage[j * p_width + i].z *= 255;
	// }
	//}
	for (int x = 0; x < p_width; x++)
	{
		for (int y = 0; y < p_height; y++)
		{
			double red = 0.0, green = 0.0, blue = 0.0;
			//multiply every value of the filter with corresponding image pixel
			for (int filterY = 0; filterY < __filterHeight; filterY++)
			{
				for (int filterX = 0; filterX < __filterWidth; filterX++)
				{
					int imageX = (x - __filterWidth / 2 + filterX + p_width) % p_width;
					int imageY = (y - __filterHeight / 2 + filterY + p_height) % p_height;
					red += p_ptrImage[imageY * p_width + imageX].x * array[filterY][filterX];
					green += p_ptrImage[imageY * p_width + imageX].y * array[filterY][filterX];
					blue += p_ptrImage[imageY * p_width + imageX].z * array[filterY][filterX];
					//std::cout << "Color x: " << p_ptrImage[filterY * p_width + filterX].x << " | Color y: " << p_ptrImage[filterY * p_width + filterX].y << " | Color z: " << p_ptrImage[filterY * p_width + filterX].z << std::endl;
					//std::cout << "red: " << red << " |green: " << green << " | Color z: " << blue << std::endl;
				}
			}
			//truncate values smaller than zero and larger than 255
			p_ptrImage[y * p_width + x].x = min(max(float(p_factor * red + p_bias), 0), 255);
			p_ptrImage[y * p_width + x].y = min(max(float(p_factor * green + p_bias), 0), 255);
			p_ptrImage[y * p_width + x].z = min(max(float(p_factor * blue + p_bias), 0), 255);
			//take absolute value and truncate to 255
			/*p_ptrImage[y * p_width + x].x = min(std::abs(int(p_factor * red + p_bias)), 255);
			p_ptrImage[y * p_width + x].y = min(std::abs(int(p_factor * green + p_bias)), 255);
			p_ptrImage[y * p_width + x].z = min(std::abs(int(p_factor * blue + p_bias)), 255);*/
			//std::cout << "r: " << p_ptrImage[y * p_width + x].x << " | g: " << p_ptrImage[y * p_width + x].y << " | b: " << p_ptrImage[y * p_width + x].z << std::endl;
		}
	}
	/*for (int i = 0; i < p_width; i++)
	{
	for (int j = 0; j < p_height; j++)
	{
	p_ptrImage[j * p_width + i].x /= 255;
	p_ptrImage[j * p_width + i].y /= 255;
	p_ptrImage[j * p_width + i].z /= 255;
	}
	}*/
	////draw the result buffer to the screen
	//for (int y = 0; y < p_height; y++)
	// for (int x = 0; x < p_width; x++)
	// {
	// pset(x, y, p_ptrImage[y * w + x]);
	// }
	////redraw & sleep
	//redraw();
	//sleep();
}
void UseToonColor(Vec3f* p_ptrImage, int p_width, int p_height, int p_multiplier = 10)
{
	for (int i = 0; i < p_width; i++)
	{
		for (int j = 0; j < p_height; j++)
		{
			float __r = p_ptrImage[j * p_width + i].x * 255;
			float __g = p_ptrImage[j * p_width + i].y * 255;
			float __b = p_ptrImage[j * p_width + i].z * 255;
			p_ptrImage[j * p_width + i].x = (float)(__r + (p_multiplier - (int)__r % p_multiplier)) / 255;
			p_ptrImage[j * p_width + i].y = (float)(__g + (p_multiplier - (int)__g % p_multiplier)) / 255;
			p_ptrImage[j * p_width + i].z = (float)(__b + (p_multiplier - (int)__b % p_multiplier)) / 255;
		}
	}
}

void SaveImagePPM(Vec3f* p_ptrImage, float p_width, float p_height)
{
	// Save result to a PPM image (keep these flags if you compile under Windows)
	std::ofstream __offStream("./untitled.ppm", std::ios::out | std::ios::binary);
	__offStream << "P6\n" << p_width << " " << p_height << "\n255\n";
	for (unsigned i = 0; i < p_width * p_height; ++i) {
		__offStream << (unsigned char)(min(float(1), p_ptrImage[i].x) * 255) <<
			(unsigned char)(min(float(1), p_ptrImage[i].y) * 255) <<
			(unsigned char)(min(float(1), p_ptrImage[i].z) * 255);
	}
	__offStream.close();
	delete[] p_ptrImage;
}

std::vector<Sphere> InitializeSpheres()
{
	srand(13);
	std::vector<Sphere> spheres;
	// position, radius, surface color, reflectivity, transparency, emission color

	//center red sphere
	spheres.push_back(Sphere(Vec3f(0.0, 0, 0), 3.0, Vec3f(1.00, 0.32, 0.36), 1, 0.5));

	//circular spheres
	spheres.push_back(Sphere(Vec3f(5.5, 0.0, 0.0), 1.5, Vec3f(0.0, 0.32, 1.0), 1, 0.0));
	spheres.push_back(Sphere(Vec3f(-5.5, 0.0, 0.0), 1.5, Vec3f(1.0, 0.32, 0.0), 1, 0.0));
	spheres.push_back(Sphere(Vec3f(-2.5, 2.0, -5.5), 1.5, Vec3f(0.0, 1.0, 0.32), 1, 0.0));
	spheres.push_back(Sphere(Vec3f(2.5, -2.0, 5.5), 1.5, Vec3f(1.0, 0.32, 1.0), 1, 0.0));

	//background giant spheres
	//spheres.push_back(Sphere(Vec3f(0.0, -10004, -20), 10000.0, Vec3f(0.0, 0.0, 0.0), 0.0, 0.0));

	// local light
	//spheres.push_back(Sphere(Vec3f(5.5, 20, -5.5), 3, Vec3f(0.00, 0.00, 0.00), 0, 0.0, Vec3f(50)));
	spheres.push_back(Sphere(Vec3f(10, 10, 30), 3, Vec3f(0.00, 0.00, 0.00), 0, 0.0, Vec3f(20)));
	return spheres;
}

int main(int argc, char **argv)
{
	float __width = 300, __height = 300;
	std::cout << "Set width and heigth? (Default is 300x 300y), type 'y' or 'n'" << std::endl;
	std::string __response = "";
	std::cin >> __response;

	int __targetRenderingMode = 0;
	bool __useToon = false;
	bool __findEdges = false;
	bool __blur = false;
	bool __verticalFindEdges = false;

	if (__response == "y")
	{
		std::cout << "type width" << std::endl;
		std::cin >> __width;
		std::cout << "type height" << std::endl;
		std::cin >> __height;
	}
	else
	{
		std::cout << "Using default..." << std::endl;
	}

	system("cls");

	std::cout << "0 - default" << std::endl;
	std::cout << "1 - depth" << std::endl;
	std::cout << "2 - normal map" << std::endl;
	std::cout << "Choose target rendering mode: " << std::endl;
	std::cin >> __response;
	if (__response == "0")
	{
		__targetRenderingMode = 0;
	}
	else if (__response == "1")
	{
		__targetRenderingMode = 1;
	}
	else if (__response == "2")
	{
		__targetRenderingMode = 2;
	}

	system("cls");

	std::cout << "Use toon colors? , type 'y' or 'n'" << std::endl;
	std::cin >> __response;
	if (__response == "y")
		__useToon = true;

	system("cls");

	std::cout << "Use findEdges filter? , type 'y' or 'n'" << std::endl;
	std::cin >> __response;
	if (__response == "y")
		__findEdges = true;

	system("cls");

	std::cout << "Use blur filter? , type 'y' or 'n'" << std::endl;
	std::cin >> __response;
	if (__response == "y")
		__blur = true;

	system("cls");

	std::cout << "Rendering Scene..." << std::endl;
	Vec3f *__ptrImage = new Vec3f[__width * __height];
	RenderScenePTR(InitializeSpheres(), *__ptrImage, __width, __height, __targetRenderingMode);


	std::cout << "Applying Filters..." << std::endl;
	if (__useToon)
		UseToonColor(__ptrImage, __width, __height, 50);
	if (__findEdges)
		FilterImage(__ptrImage, __width, __height, filterFindEdgesWidth, filterFindEdgesHeight, filterFindEdges, factorFindEdges, biasFindEdges);
	if (__blur)
		FilterImage(__ptrImage, __width, __height, filterBlurWidth, filterBlurWidth, filterBlur, blurFactor, blurBias);

	std::cout << "Saving image.ppm..." << std::endl;
	SaveImagePPM(__ptrImage, __width, __height);

	std::cout << "Image Finished!" << std::endl;

	system("Pause");
	return 0;
}

/*
class Object
{
public:
	Object(const Matrix44f &o2w) : objectToWorld(o2w), worldToObject(o2w.inverse()) {}
	virtual ~Object() {}
	virtual bool intersect(const Vec3f &, const Vec3f &, float &, uint32_t &, Vec2f &) const = 0;
	virtual void getSurfaceProperties(const Vec3f &, const Vec3f &, const uint32_t &, const Vec2f &, Vec3f &, Vec2f &) const = 0;
	Matrix44f objectToWorld, worldToObject;
	MaterialType type = kDiffuse;
	Vec3f albedo = 0.18;
	float Kd = 0.8; // phong model diffuse weight
	float Ks = 0.2; // phong model specular weight
	float n = 10; // phong specular exponent
};
*/

bool solveQuadratic(const float &a, const float &b, const float &c, float &x0, float &x1)
{
	float discr = b * b - 4 * a * c;
	if (discr < 0) return false;
	else if (discr == 0)
	{
		x0 = x1 = -0.5 * b / a;
	}
	else
	{
		float q = (b > 0) ?	-0.5 * (b + sqrt(discr)) : -0.5 * (b - sqrt(discr));
		x0 = q / a;
		x1 = c / q;
	}

	return true;
}

bool solveQuadratic(const float &a, const float &b, const float &c, float &x0, float &x1) 
{ 
    float discr = b * b - 4 * a * c; 
    if (discr < 0) return false; 
    else if (discr == 0) { 
        x0 = x1 = - 0.5 * b / a; 
    } 
    else { 
        float q = (b > 0) ? 
            -0.5 * (b + sqrt(discr)) : 
            -0.5 * (b - sqrt(discr)); 
        x0 = q / a; 
        x1 = c / q; 
    } 
 
    return true; 
}

/*
class Light
{
public:
	Light(const Matrix44f &l2w, const Vec3f &c = 1, const float &i = 1) : lightToWorld(l2w), color(c), intensity(i) {}
	virtual ~Light() {}
	virtual void illuminate(const Vec3f &P, Vec3f &, Vec3f &, float &) const = 0;
	Vec3f color;
	float intensity;
	Matrix44f lightToWorld;
};

class DistantLight : public Light //Directional?????
{
	Vec3f dir;
public:
	DistantLight(const Matrix44f &l2w, const Vec3f &c = 1, const float &i = 1) : Light(l2w, c, i)
	{
		l2w.multDirMatrix(Vec3f(0, 0, -1), dir);
		dir.normalize(); // in case the matrix scales the light
	}
	void illuminate(const Vec3f &P, Vec3f &lightDir, Vec3f &lightIntensity, float &distance) const
	{
		lightDir = dir;
		lightIntensity = color * intensity;
		distance = kInfinity;
	}
};
*/
/*
Step 1: create a Cartesian coordinate system in which the up vector is oriented along the shaded point normal N(the shaded point normal N
and the up vector of the coordinate system are aligned).
Step2 : create a sample using the spherical to Cartesian coordinates equations.We will show in this chapter how this can be done in practice.
Step 3 : transform the sample direction from the original coordinate system to the shaded point coordinate system.
Step 4 : trace a ray in the scene in the sampled direction.
Step 5 : if the ray intersects an object, compute the color of that object at the intersection point and add this result to a temporary variable.Because the surface is diffuse, don't forget to multiply the light intensity returned along each ray by the dot product between the ray direction (the light direction) and the shaded normal N
(see below).
Step 6 : repeat step 2 to 5 N - times.
Step 7 : divide the temporary variables that holds all the results of all the sampled rays by N, the total number of samples used.The final value is an approximation of the shaded point indirect diffuse illumination.The illumination of P
by other diffuse surfaces in the scene.
*/

/*
Vec3f castRay(Vec3f &orig;, Vec3f &dir;, const uint32_t &depth;, ...)
{
if (depth > options.maxDepth) return 0;
Vec3f hitPointColor = 0;
// compute direct ligthing
...
// step1: compute shaded point coordinate system using normal N.
...
// number of samples N
uint32_t N = 16;
Vec3f indirectDiffuse = 0;
for (uint32_t i = 0; i < N; ++i) {
// step 2: create sample in world space
Vec3f sample = ...;
// step 3: transform sample from world space to shaded point local coordinate system
sampleWorld = ...;
// step 4 & 5: cast a ray in this direction
indirectDiffuse += N.dotProduct(sampleWorld) * castRay(P, sampleWorld, depth + 1, ...);
}
// step 7: divide the sum by the total number of samples N
hitPointColor += (indirectDiffuse / N) * albedo;
...
return hitPointColor;
} */