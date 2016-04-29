#include "SimplifierHE.h"
#include <algorithm>
#include <iostream>
#include <queue>
#include <cassert>
using namespace std;

void SimplifierHE::initialize() {
	for (int i = 0; i < this->m_nVertices; ++i) {
		Vertex v;
		v.pos = this->m_pVertexList[i];
		vlist.push_back(v);
	}
	for (int i = 0; i < this->m_nTriangles; ++i) {
		Face f;
		f.edge = elist.size();
		ER fr = flist.size();
		flist.push_back(f);
		VR vl[3];
		for (int j = 0; j < 3; ++j)
			vl[j] = this->m_pTriangleList[i][j];
		for (int j = 0; j < 3; ++j) {
			for (int k = j + 1; k < 3; ++k) {
				HalfEdge he;
				he.face = fr;
				pair<int, int> link(vl[j], vl[k]);
				if (vl[j] > vl[k])
					link = pair<int, int>(vl[k], vl[j]);
				map<pair<int, int>, ER>::iterator itf = edgemap.find(link);
				if (itf == edgemap.end()) {
					edgemap.insert(pair<pair<int, int>, ER>(link, elist.size()));
				}
				else {
					he.paire = (*itf).second;
					elist[(*itf).second].paire = elist.size();
				}
				elist.push_back(he);
			}
		}
	}
	// Arbitrary setting orientation of first face
	ER e0 = flist[0].edge; // 0-1
	ER e1 = e0 + 1; // 2-0
	ER e2 = e1 + 1; // 1-2
	elist[e0].end = m_pTriangleList[0][1]; elist[e0].next = e2; vlist[elist[e0].end].edge = e2;
	elist[e1].end = m_pTriangleList[0][0]; elist[e1].next = e0; vlist[elist[e1].end].edge = e0;
	elist[e2].end = m_pTriangleList[0][2]; elist[e2].next = e1; vlist[elist[e2].end].edge = e1;
	// First face orientation set!
	set_face_orientation();

	//for (int i = 0; i < vlist.size(); ++i) {
	//	cout << i << ' ' << vlist[i].pos.x << ' ' << vlist[i].pos.y << ' ' << vlist[i].pos.z << endl;
	//	ER e = vlist[i].edge;
	//	while (1) {
	//		cout << elist[e].face << endl;
	//		cout << this->m_pTriangleList[elist[e].face][0] << ' ' << this->m_pTriangleList[elist[e].face][1] << ' ' << this->m_pTriangleList[elist[e].face][2] << ' ' << endl;
	//		e = elist[e].paire;
	//		e = elist[e].next;
	//		if (e == vlist[i].edge)
	//			break;
	//	}
	//}
}

void SimplifierHE::simplify() {

}

void SimplifierHE::set_face_orientation() {
	queue<FR> frqueue;
	frqueue.push(0);
	while (!frqueue.empty()) {
		FR fr = frqueue.front();
		frqueue.pop();
		ER e = flist[fr].edge;
		for (int i = 0; i < 3; ++i) {
			VR end = elist[e].end;
			VR start = elist[elist[elist[e].next].next].end;
			if (elist[e].paire != -1) {
				FR nf = elist[elist[e].paire].face;
				if (elist[flist[nf].edge].next == -1) {
					int sr = 0, er = 0, th = 0;
					for (int k = 0; k < 3; ++k) {
						if (this->m_pTriangleList[nf][k] == start)
							sr = k;
						else if (this->m_pTriangleList[nf][k] == end)
							er = k;
						else
							th = k;
					}
					ER ees = flist[nf].edge + (er + sr - 1);
					ER est = flist[nf].edge + (sr + th - 1);
					ER ete = flist[nf].edge + (th + er - 1);
					elist[ees].end = start; elist[ees].next = est;
					elist[est].end = this->m_pTriangleList[nf][th]; elist[est].next = ete;
					elist[ete].end = end; elist[ete].next = ees;
					vlist[this->m_pTriangleList[nf][er]].edge = ees;
					vlist[this->m_pTriangleList[nf][sr]].edge = est;
					vlist[this->m_pTriangleList[nf][th]].edge = ete;
					// Finish orientation
					frqueue.push(nf);
				}
			}
			e = elist[e].next;
		}
	}
}