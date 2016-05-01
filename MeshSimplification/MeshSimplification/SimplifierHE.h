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
	typedef size_t PR;

	struct HalfEdge {
		VR end;
		ER paire;
		FR face;
		ER next;
		bool left;
		HalfEdge() : end(-1), paire(-1), face(-1), next(-1), left(true) {}
	};

	struct Vertex {
		Vec3f pos;
		ER edge;
		bool left;
		int final_seq;
		Vertex() : edge(-1), left(true), final_seq(0) {}
	};

	struct Face {
		ER edge;
		double series[4];
		bool left;
		Face() : edge(-1), left(true) {}
	};

	struct VPair {
		VR v[2];
		Vec3f vbar;
		double cost;
		VPair() : cost(0) {}
	};

	struct VPairHeap {
		map<pair<int, int>, PR> pairmap;
		vector<VPair> pairlist;
		vector<int> rankinheap;
		vector<PR> pairheap;
		void up(PR pr);
		void down(PR pr);
		VPair top();
	};

	class SimplifierHE : public CSimpleObject {
	public:
		void initialize();
		void simplify();
		void set_alpha(double al) { alpha = al; }

	private:
		vector<Vertex> vlist;
		vector<Face> flist;
		vector<HalfEdge> elist;
		map<pair<int, int>, ER> edgemap;
		VPairHeap vpairheap;

		int facenum;
		double alpha;

		void set_face_orientation();
		void solve_series(FR fr);
		void solve_vpair(PR pr);
		void cutpair(VPair& vpair);
	};
}