//#include "Filter.h"
#include "Input.h"
#include "Render.h"
#include "Sphere.h"
#include "Vector3.h"

void main()
{
	Render camera;
	std::fstream input("./Input.txt", std::ios_base::in);
	float posX, posY, posZ, radius, red, green, blue, reflection, transp; // Sphere

	input >> posX;
	input >> posY;
	input >> posZ;
	input >> radius;
	input >> red;
	input >> green;
	input >> blue;
	input >> reflection;
	input >> transp;

	float lPosX, lPosY, lPosZ, lRadius, lRed, lGreen, lBlue, lReflection, lTransp, lEColor; //Light

	input >> lPosX;
	input >> lPosY;
	input >> lPosZ;
	input >> lRadius;
	input >> lRed;
	input >> lGreen;
	input >> lBlue;
	input >> lReflection;
	input >> lTransp;
	input >> lEColor;

	std::vector<Sphere> spheres;
	// position, radius, surface color, reflectivity, transparency, emission color
	spheres.push_back(Sphere(Vector3f(posX, posY, posZ), radius, Vector3f(red, green, blue), reflection, transp));
	spheres.push_back(Sphere(Vector3f(0.0, 0, -20), 4, Vector3f(1.00, 0.32, 0.36), 1, 0.5));
	spheres.push_back(Sphere(Vector3f(5.0, -1, -15), 2, Vector3f(0.90, 0.76, 0.46), 1, 0.0));
	spheres.push_back(Sphere(Vector3f(5.0, 0, -25), 3, Vector3f(0.65, 0.77, 0.97), 1, 0.0));
	spheres.push_back(Sphere(Vector3f(-5.5, 0, -15), 3, Vector3f(0.90, 0.90, 0.90), 1, 0.0));
	// light
	spheres.push_back(Sphere(Vector3f(lPosX, lPosY, lPosZ), lRadius, Vector3f(lRed, lGreen, lBlue), lReflection, lTransp, Vector3f(lEColor)));
	//spheres.push_back(Sphere(Vec3f(0.0, 20, -30), 3, Vec3f(0.00, 0.00, 0.00), 0, 0.0, Vec3f(3)));
	input.close();
	camera.RenderScene(spheres, 640, 480);
	
}