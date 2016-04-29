#pragma once
#include "parser/SimpleObject.h"
#include <vector>
#include <map>
#include <utility>

using SimpleOBJ::CSimpleObject;
using SimpleOBJ::Array;
using SimpleOBJ::Vec3f;

namespace std {
	typedef size_t VR;
	typedef size_t FR;
	typedef size_t ER;

	struct HalfEdge {
		VR end;
		ER paire;
		FR face;
		ER next;
		HalfEdge() : end(-1), paire(-1), face(-1), next(-1) {}
	};

	struct Vertex {
		Vec3f pos;
		ER edge;
		Vertex() : edge(-1) {}
	};

	struct Face {
		ER edge;
		double series[4];
		Face() : edge(-1) {}
	};

	struct VPair {
		VR v[2];
		Vec3f vbar;
		double cost;
	};

	class SimplifierHE : public CSimpleObject {
	public:
		void initialize();
		void simplify();

	private:
		vector<Vertex> vlist;
		vector<Face> flist;
		vector<HalfEdge> elist;
		map<pair<int, int>, ER> edgemap;

		void set_face_orientation();
	};
}