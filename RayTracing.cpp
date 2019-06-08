//#define _4_Threads_Rendering

#include "World.h"
using namespace std;

#include "GlassMan.h"


#include "Test.cpp"
#include "CPT.cpp"

#include <stdlib.h>
#include <crtdbg.h>


int main()
{
	cout << thread::hardware_concurrency() << endl;

	fout << fixed << showpoint << setprecision(5);

	_CrtSetBreakAlloc(0);
	CPT_Animation_T1();
	_CrtDumpMemoryLeaks();	// Visual Studio check memory leak

	//pause;

	return 0;

}

