//#define _4_Threads_Rendering

#include "World.h"
using namespace std;

#include "GlassMan.h"


#include "Test.cpp"
#include "CPT.cpp"

#include <stdlib.h>
#include <crtdbg.h>

#include "HumanWalkingData.h"

int main()
{
	//WalkingMan::drawTrace();
	//WalkingMan::VisualizeFourierSeries(WalkingMan::v_foot_l, 6); return 0;

	fout << fixed << showpoint << setprecision(5);

	//system("del /f /s /q C:\\Users\\harry\\Desktop\\RayTracing\\RayTracing\\Animation\\*.bmp");

	auto t0 = NTime::now();

	_CrtSetBreakAlloc(0);
	CPT_Animation_R();
	_CrtDumpMemoryLeaks();	// Visual Studio check memory leak

	auto t1 = NTime::now();
	fsec fs = t1 - t0;

	cout << "Total " << fs.count() << "s elapsed. \n\n";

	//pause;

	return 0;

}

