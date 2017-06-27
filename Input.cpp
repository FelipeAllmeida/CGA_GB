#include "Input.h"

std::fstream input("./Input.txt", std::ios_base::in);

Input::Input()
{
}


Input::~Input()
{
}

void Input::Read()
{
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
}