// #include <opencv2/opencv.hpp>
#include "parser/SimpleObject.h"
#include "parser/Vec3f.h"
#include <iostream>
#include "SimplifierHE.h"

using namespace std;
using namespace SimpleOBJ;

int main(int argc, char** argv) {
	if (argc != 3) {
		cout << "Usage: ./<main.exe> <input.obj> <output.obj>" << endl;
		return -1;
	}
	SimplifierHE csobj;
	csobj.LoadFromObj(argv[1]);
	csobj.initialize();
	return 0;
}