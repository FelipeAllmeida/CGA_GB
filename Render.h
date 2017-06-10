#pragma once

//[comment]
// This variable controls the maximum recursion depth
//[/comment]
#define MAX_RAY_DEPTH 5
#define min(a,b)            (((a) < (b)) ? (a) : (b)) //O visual studio/pc não reconheceu o min e max padrão do std, então criei um
#define max(a,b)            (((a) > (b)) ? (a) : (b))
#define M_PI 3.141592653589793

#include "Sphere.h" //Sphere já tem Vector3 incluido
#include <vector>

class Render
{
public:
	Render();
	~Render();
	float Mix(const float &p_a, const float &p_b, const float &p_mix);
	Vector3f Trace(const Vector3f &p_rayorig, const Vector3f &p_raydir,	const std::vector<Sphere> &p_spheres, const int &p_depth);
	void RenderScene(const std::vector<Sphere> &p_spheres, Vector3f& p_ptrImage, unsigned  p_width, unsigned  p_height);
};

