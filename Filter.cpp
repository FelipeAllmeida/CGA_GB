#include "Filter.h"



Filter::Filter()
{
}


Filter::~Filter()
{
}

/*
int** filter = new int*[filterHeight];
for (int i = 0; i < filterHeight; ++i)
	filter[i] = new int[filterWidth];

unsigned width = 640, height = 480;
Vec3f *image = new Vec3f[width * height], *pixel = image;

if (whichFilter == 0)
{
	int filterExample[3][3] =
	{
		0, 0, 0,
		0, 1, 0,
		0, 0, 0
	};
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			filter[i][j] = filterExample[i][j];
		}
	}
}
if (whichFilter == 1) //MotionBlur
{
	int filterExample2[9][9] =
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
	for (int i = 0; i < 9; i++)
	{
		for (int j = 0; j < 9; j++)
		{
			filter[i][j] = filterExample2[i][j];
			//std::cout << filter[i][j] << " ";
		}
		//std::cout << std::endl;
	}
}
if (whichFilter == 2)//FindEdge
{

	int filterExample3[5][5] =
	{
		-1,  0,  0,  0,  0,
		0, -2,  0,  0,  0,
		0,  0,  6,  0,  0,
		0,  0,  0, -2,  0,
		0,  0,  0,  0, -1,
	};
	for (int i = 0; i < 5; i++)
	{
		for (int j = 0; j < 5; j++)
		{
			filter[i][j] = filterExample3[i][j];
			//std::cout << filter[i][j] << " ";
		}
		//std::cout << std::endl;
	}
}
if (whichFilter == 3)
{
	int filterExample4[3][3] =
	{
		-1, -1, 0,
		-1, 0, 1,
		0, 1, 1
	};
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			filter[i][j] = filterExample4[i][j];
			//std::cout << filter[i][j] << " ";
		}
		//std::cout << std::endl;
	}
}
if (whichFilter == 4) //MotionBlur
{
	int filterExample5[9][9] =
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
	for (int i = 0; i < 9; i++)
	{
		for (int j = 0; j < 9; j++)
		{
			filter[i][j] = filterExample5[i][j];
			//std::cout << filter[i][j] << " ";
		}
		//std::cout << std::endl;
	}
}

if (pass == false)
{
	float invWidth = 1 / float(width), invHeight = 1 / float(height);
	float fov = 30, aspectratio = width / float(height);
	float angle = tan(M_PI * 0.5 * fov / 180.);
	// Trace rays
	for (unsigned y = 0; y < height; ++y)
	{
		for (unsigned x = 0; x < width; ++x, ++pixel)
		{
			float xx = (2 * ((x + 0.5) * invWidth) - 1) * angle * aspectratio;
			float yy = (1 - 2 * ((y + 0.5) * invHeight)) * angle;
			Vec3f raydir(xx, yy, -1);
			raydir.normalize();
			*pixel = trace(Vec3f(0), raydir, spheres, 0);
		}
	}
	pass = true;
}

//apply the filter
for (int w = 0; w < width; w++)
	for (int h = 0; h < height; h++)
	{
		double red = 0.0, green = 0.0, blue = 0.0;

		//multiply every value of the filter with corresponding image pixel
		for (int filterY = 0; filterY < filterHeight; filterY++)
			for (int filterX = 0; filterX < filterWidth; filterX++)
			{
				int imageX = (w - filterWidth / 2 + filterX + width) % width;
				int imageY = (h - filterHeight / 2 + filterY + height) % height;
				red += image[imageY * width + imageX].x * filter[filterY][filterX];
				green += image[imageY * width + imageX].y * filter[filterY][filterX];
				blue += image[imageY * width + imageX].z * filter[filterY][filterX];

				if (red > 255)
					red = 255;
				if (red < 0)
					red = 0;
				if (green > 255)
					green = 255;
				if (green < 0)
					green = 0;
				if (blue > 255)
					blue = 255;
				if (blue < 0)
					blue = 0;
			}

		//truncate values smaller than zero and larger than 255
		image[h * width + w].x = factor * red + bias; //R
		image[h * width + w].y = factor * green + bias; //G
		image[h * width + w].z = factor * blue + bias; //B
	}

// Save result to a PPM image (keep these flags if you compile under Windows)
std::ofstream ofs("./untitled.ppm", std::ios::out | std::ios::binary);
ofs << "P6\n" << width << " " << height << "\n255\n";
for (unsigned i = 0; i < width * height; ++i)
{
	ofs << (unsigned char)(min(float(1), image[i].x) * 255) <<
		(unsigned char)(min(float(1), image[i].y) * 255) <<
		(unsigned char)(min(float(1), image[i].z) * 255);
}
ofs.close();
//if (whichFilter == finalFilter)
//{
delete[] image;
for (int i = 0; i < filterHeight; ++i)
{
	delete[] filter[i];
}
//}
delete[] filter;
*/