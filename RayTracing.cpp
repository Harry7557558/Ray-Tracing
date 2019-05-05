#include "World.h"
using namespace std;




//#include "RayTracing2D.cpp"
#include "Test.cpp"

#include <stdlib.h>
#include <crtdbg.h>

int main()
{
	//system("del /f /s /q D:\\Coding\\AboutPhysics\\RayTracing\\RayTracing\\Debug\\*.* /s");
	//system("del /f /s /q D:\\Coding\\AboutPhysics\\RayTracing\\RayTracing\\Release\\*.* /s");

	fout << fixed << showpoint << setprecision(5);

	_CrtSetBreakAlloc(-1);
	Render_CTest03();
	_CrtDumpMemoryLeaks();	// Visual Studio check memory leak

	return 0;

}

// https://github.com/miloyip/light2d
