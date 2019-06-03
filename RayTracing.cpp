#include "World.h"
using namespace std;

#include "GlassMan.h"


#include "Test.cpp"
#include "CPT.cpp"

#include <stdlib.h>
#include <crtdbg.h>

int main()
{
	
	fout << fixed << showpoint << setprecision(5);

	_CrtSetBreakAlloc(-1);
	CPT_T1();
	_CrtDumpMemoryLeaks();	// Visual Studio check memory leak

	//pause;

	return 0;

}

