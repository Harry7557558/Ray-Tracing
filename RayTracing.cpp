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
	//cout << (GlassMan_std(RunningAttitude_Constructor(21.6 - 21, vec3(point(1.8, 2.3), point(6.2, 7.2)))).Pos + point(1.8, 2.3)) << endl; pause;
	//RunningAttitude::visualize(RunningAttitude::v_shank_l); return 0;
	//RunningAttitude::VisualizeFourier(RunningAttitude::v_shank_r); return 0;
	//RunningAttitude::VisualizeFourier(); return 0;
	//RunningAttitude::visualize(RunningAttitude::upperarm_r);
	//RunningAttitude::visualize(RunningAttitude::shoulder_r, RunningAttitude::elbow_r);
	//pause; return 0;


	fout << fixed << showpoint << setprecision(5);

	//system("del /f /s /q C:\\Users\\harry\\Desktop\\RayTracing\\RayTracing\\Animation\\*.bmp");

	auto t0 = NTime::now();

	_CrtSetBreakAlloc(0);
	CPT_Animation_R();
	_CrtDumpMemoryLeaks();	// Visual Studio check memory leak

	auto t1 = NTime::now();
	fsec fs = t1 - t0;

	cout << "Total " << fs.count() << "s elapsed. \n\n";

	pause;

	return 0;

}

