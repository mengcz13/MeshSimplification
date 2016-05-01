// #include <opencv2/opencv.hpp>
#include "parser/SimpleObject.h"
#include "parser/Vec3f.h"
#include <iostream>
#include <cstdlib>
#include "SimplifierHE.h"

using namespace std;
using namespace SimpleOBJ;

int main(int argc, char** argv) {
	if (argc != 4) {
		cout << "Usage: ./<main.exe> <input.obj> <output.obj> rate" << endl;
		return -1;
	}
	SimplifierHE csobj;
	csobj.set_alpha(atof(argv[3]));
	csobj.LoadFromObj(argv[1]);
	csobj.initialize();
	csobj.simplify();
	csobj.SaveToObj(argv[2]);
	return 0;
}