#include "SimplifierHE.h"
#include <algorithm>
#include <iostream>
#include <queue>
#include <cassert>
#include <opencv2/opencv.hpp>
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

	for (int i = 0; i < flist.size(); ++i) {
		solve_series(i);
	}

	facenum = flist.size();
}

void SimplifierHE::simplify() {
	for (int i = 0; i < flist.size(); ++i) {
		for (int j = 0; j < 2; ++j) {
			VPair vpair;
			vpair.v[0] = m_pTriangleList[i][j];
			vpair.v[1] = m_pTriangleList[i][(j + 1) % 3];
			if (vpair.v[0] > vpair.v[1])
				swap(vpair.v[0], vpair.v[1]);
			map<pair<int, int>, PR>::iterator ite = vpairheap.pairmap.find(pair<int, int>(vpair.v[0], vpair.v[1]));
			if (ite == vpairheap.pairmap.end()) {
				vpairheap.pairmap.insert(pair<pair<int, int>, PR>(pair<int, int>(vpair.v[0], vpair.v[1]), vpairheap.pairlist.size()));
				vpairheap.pairlist.push_back(vpair);
				solve_vpair(vpairheap.pairlist.size() - 1);
				vpairheap.rankinheap.push_back(vpairheap.pairheap.size());
				vpairheap.pairheap.push_back(vpairheap.pairlist.size() - 1);
			}
		}
	}
	for (int i = vpairheap.pairheap.size() - 1; i >= 0; --i) {
		vpairheap.down(i);
	}

	while ((double)facenum / (double)flist.size() > alpha) {
		VPair cut = vpairheap.top();
		cutpair(cut);
		//if ((flist.size() - facenum) % 100 == 0)
		//	cout << facenum << endl;
	}

	// Reform to vertices list and triangles list
	m_nVertices = 0;
	for (int i = 0; i < vlist.size(); ++i) {
		if (vlist[i].left) {
			m_pVertexList[m_nVertices] = vlist[i].pos;
			vlist[i].final_seq = m_nVertices;
			++m_nVertices;
		}
	}
	m_nTriangles = 0;
	for (int i = 0; i < flist.size(); ++i) {
		if (flist[i].left) {
			ER e = flist[i].edge;
			bool valid = true;
			for (int j = 0; j < 3; ++j) {
				m_pTriangleList[m_nTriangles][j] = vlist[elist[e].end].final_seq;

				if (!vlist[elist[e].end].left)
					valid = false;

				e = elist[e].next;
			}
			if (valid)
				++m_nTriangles;
		}
	}
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

void SimplifierHE::solve_series(FR fr) {
	VR vl[3];
	ER e = flist[fr].edge;
	for (int j = 0; j < 3; ++j) {
		vl[j] = elist[e].end;
		e = elist[e].next;
	}
	Vec3f vec1 = vlist[vl[1]].pos - vlist[vl[0]].pos;
	Vec3f vec2 = vlist[vl[2]].pos - vlist[vl[0]].pos;
	// cross
	double x1 = vec1[0], y1 = vec1[1], z1 = vec1[2];
	double x2 = vec2[0], y2 = vec2[1], z2 = vec2[2];
	flist[fr].series[0] = y1 * z2 - y2 * z1;
	flist[fr].series[1] = x2 * z1 - x1 * z2;
	flist[fr].series[2] = x1 * y2 - x2 * y1;
	flist[fr].series[3] = -(flist[fr].series[0] * vlist[vl[0]].pos[0] + flist[fr].series[1] * vlist[vl[0]].pos[1] + flist[fr].series[2] * vlist[vl[0]].pos[2]);
	double a = flist[fr].series[0], b = flist[fr].series[1], c = flist[fr].series[2];
	double norm2 = sqrt(a * a + b * b + c * c);
	for (int i = 0; i < 4; ++i) {
		flist[fr].series[i] /= norm2;
	}
}

void SimplifierHE::solve_vpair(PR pr) {
	cv::Matx44d Q = cv::Matx44d::zeros();
	VPair& vp = vpairheap.pairlist[pr];
	for (int i = 0; i < 2; ++i) {
		ER e = vlist[vp.v[i]].edge;
		bool border = false;
		int max = 0;
		do {
			// Ugly method to get away from dead loop...
			++max;
			cv::Matx41d ser;
			for (int j = 0; j < 4; ++j)
				ser(j) = flist[elist[e].face].series[j];
			cv::Matx44d Qv = ser * ser.t();
			Q += Qv;
			e = elist[e].paire;
			if (e == -1) {
				border = true;
				break;
			}
			if (max > 10) {
				break;
			}
			e = elist[e].next;
		} while (e != vlist[vp.v[i]].edge);
		if (border) {
			e = vlist[vp.v[i]].edge;
			e = elist[elist[e].next].next;
			e = elist[e].paire;
			while (e != -1) {
				cv::Matx41d ser;
				for (int j = 0; j < 4; ++j)
					ser(j) = flist[elist[e].face].series[j];
				cv::Matx44d Qv = ser * ser.t();
				Q += Qv;
				e = elist[elist[e].next].next;
				e = elist[e].paire;
			}
		}
	}
	cv::Matx41d vm;
	cv::Matx41d ans = cv::Matx41d::zeros(); ans(3) = 1;
	cv::Matx44d dQ;
	for (int i = 0; i < 4; ++i)
		for (int j = 0; j < 4; ++j) {
			if (j >= i)
				dQ(i, j) = Q(i, j);
			else
				dQ(i, j) = Q(j, i);
		}
	dQ(3, 0) = dQ(3, 1) = dQ(3, 2) = 0; dQ(3, 3) = 1;
	int res = cv::solve(dQ, ans, vm);
	if (res == 0) { // To be fixed...
		Vec3f threep[3];
		threep[0] = vlist[vp.v[0]].pos;
		threep[1] = vlist[vp.v[1]].pos;
		threep[2] = (vlist[vp.v[0]].pos + vlist[vp.v[1]].pos);
		threep[2] /= 2;
		int choose = 0; double min = 1e10;
		cv::Matx41d vvmm;
		for (int q = 0; q < 3; ++q) {
			vvmm = cv::Matx41d(threep[q][0], threep[q][1], threep[q][2], 1);
			double costt = (vvmm.t() * Q * vvmm)(0);
			if (costt < min) {
				min = costt;
				choose = q;
			}
		}
		vp.vbar = threep[choose];
	}
	else if (res == 1) {
		for (int i = 0; i < 3; ++i)
			vp.vbar[i] = vm(i);
	}
	for (int i = 0; i < 3; ++i)
		vm(i) = vp.vbar[i];
	vm(3) = 1;
	vp.cost = (vm.t() * Q * vm)(0);
}

void SimplifierHE::cutpair(VPair& vpair) {
	if (!vlist[vpair.v[0]].left || !vlist[vpair.v[1]].left) {
		return;
	}
	if (!flist[elist[vlist[vpair.v[0]].edge].face].left || !flist[elist[vlist[vpair.v[1]].edge].face].left)
		return;
	bool border[2] = { false, false };
	/*
	inner - border nocut!
	vector<VR>[2]
	vector<FR>[2]
	v -> vbar
	u -> no left
	for f in vector<FR>: solve_series
	deleted face: change paire
	add v[0]-x(x-v[1]), change v-x pair x-y pair
	*/
	vector<VR> adjv[2];
	vector<FR> adjf[2];
	vector<FR> trif;
	for (int i = 0; i < 2; ++i) {
		ER e = vlist[vpair.v[i]].edge;
		do {
			adjv[i].push_back(elist[e].end);
			if (!vlist[elist[e].end].left) {
				border[i] = true;
				break;
			}
			adjf[i].push_back(elist[e].face);
			if (elist[e].end == vpair.v[1 - i])
				trif.push_back(elist[e].face);
			e = elist[e].paire;
			if (e == -1) {
				border[i] = true;
				break;
			}
			e = elist[e].next;
		} while (e != vlist[vpair.v[i]].edge);
	}
	//if ((border[0] && !border[1]) || (border[1] && !border[0]))
	//	return;
	if (border[0] || border[1])
		return;
	for (int i = 0; i < 2; ++i) {
		if (border[i]) {
			ER e = vlist[vpair.v[i]].edge;
			e = elist[elist[e].next].next;
			e = elist[e].paire;
			while (e != -1) {
				adjv[i].push_back(elist[e].end);
				adjf[i].push_back(elist[e].face);
				if (elist[e].end == vpair.v[1 - i])
					trif.push_back(elist[e].face);
				e = elist[elist[e].next].next;
				e = elist[e].paire;
			}
		}
	}
	vlist[vpair.v[0]].pos = vpair.vbar;
	vlist[vpair.v[1]].left = false;

	// modify halfedges
	for (int i = 0; i < trif.size(); ++i) {
		flist[trif[i]].left = false;
	}
	vector<VR> relav;
	for (int i = 0; i < trif.size(); ++i) {
		vector<ER> olde; olde.assign(3, 0);
		int xo = 0;
		ER e = flist[trif[i]].edge;
		for (int j = 0; j < 3; ++j) {
			VR end = elist[e].end;
			VR start = elist[elist[elist[e].next].next].end;
			relav.push_back(end);
			if ((start != vpair.v[0] && start != vpair.v[1]) || (end != vpair.v[0] && end != vpair.v[1])) {
				olde[xo++] = e;
			}
			else
				olde[2] = e;
			e = elist[e].next;
		}
		for (int j = 0; j < 2; ++j) {
			ER ep = elist[olde[j]].paire;
			if (ep != -1) {
				elist[ep].paire = elist[olde[1 - j]].paire;
			}
		}
	}
	for (int j = 0; j < relav.size(); ++j) {
		VR end = relav[j];
		if (flist[elist[vlist[end].edge].face].left == false) {
			ER ef = vlist[end].edge;
			bool border = false;
			do {
				if (flist[elist[ef].face].left == true) {
					vlist[end].edge = ef;
					break;
				}
				ef = elist[ef].paire;
				if (ef == -1) {
					border = true;
					break;
				}
				ef = elist[ef].next;
			} while (ef != vlist[end].edge);
			if (border) {
				ef = vlist[end].edge;
				while (ef != -1) {
					if (flist[elist[ef].face].left == true) {
						vlist[end].edge = ef;
						break;
					}
					ef = elist[elist[ef].next].next;
					ef = elist[ef].paire;
				}
			}
		}
	}

	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < adjf[i].size(); ++j) {
			if (flist[adjf[i][j]].left) {
				ER e = flist[adjf[i][j]].edge;
				for (int k = 0; k < 3; ++k) {
					if (elist[e].end == vpair.v[1])
						elist[e].end = vpair.v[0];
					e = elist[e].next;
				}
				solve_series(adjf[i][j]);
			}
		}
	}

	for (int i = 0; i < 2; ++i) {
		for (int j = 0; j < adjv[i].size(); ++j) {
			if (vlist[adjv[i][j]].left == false)
				continue;
			ER e = vlist[adjv[i][j]].edge;
			bool borderv = false;
			do {
				VR o1 = elist[e].end;
				VR o2 = adjv[i][j];
				if (vlist[o1].left && vlist[o2].left) {
					
					if (o1 > o2)
						swap(o1, o2);
					map<pair<int, int>, PR>::iterator ite = vpairheap.pairmap.find(pair<int, int>(o1, o2));
					if (ite == vpairheap.pairmap.end()) {
						VPair vpp;
						vpp.v[0] = o1, vpp.v[1] = o2;
						vpairheap.pairlist.push_back(vpp);
						vpairheap.pairmap.insert(pair<pair<int, int>, PR>(pair<int, int>(o1, o2), vpairheap.pairlist.size() - 1));
						vpairheap.pairheap.push_back(vpairheap.pairlist.size() - 1);
						vpairheap.rankinheap.push_back(vpairheap.pairheap.size() - 1);
						solve_vpair(vpairheap.pairlist.size() - 1);
						vpairheap.up(vpairheap.pairheap.size() - 1);
					}
					else {
						PR target = (*ite).second;
						if (vpairheap.rankinheap[target] != -1) {
							double origin = vpairheap.pairlist[target].cost;
							solve_vpair(target);
							if (vpairheap.pairlist[target].cost < origin) {
								vpairheap.up(vpairheap.rankinheap[target]);
							}
							else {
								vpairheap.down(vpairheap.rankinheap[target]);
							}
						}
					}
				}
				e = elist[e].paire;
				if (e == -1) {
					borderv = true;
					break;
				}
				e = elist[e].next;
			} while (e != vlist[adjv[i][j]].edge);
			if (borderv) {
				e = vlist[adjv[i][j]].edge;
				while (e != -1) {
					VR o1 = elist[e].end;
					if (vlist[o1].left) {
						VR o2 = adjv[i][j];
						if (o1 > o2)
							swap(o1, o2);
						map<pair<int, int>, PR>::iterator ite = vpairheap.pairmap.find(pair<int, int>(o1, o2));
						if (ite == vpairheap.pairmap.end()) {
							VPair vpp;
							vpp.v[0] = o1, vpp.v[1] = o2;
							vpairheap.pairlist.push_back(vpp);
							vpairheap.pairmap.insert(pair<pair<int, int>, PR>(pair<int, int>(o1, o2), vpairheap.pairlist.size() - 1));
							vpairheap.pairheap.push_back(vpairheap.pairlist.size() - 1);
							vpairheap.rankinheap.push_back(vpairheap.pairheap.size() - 1);
							solve_vpair(vpairheap.pairlist.size() - 1);
							vpairheap.up(vpairheap.pairheap.size() - 1);
						}
						else {
							PR target = (*ite).second;
							if (vpairheap.rankinheap[target] != -1) {
								double origin = vpairheap.pairlist[target].cost;
								solve_vpair(target);
								if (vpairheap.pairlist[target].cost < origin) {
									vpairheap.up(vpairheap.rankinheap[target]);
								}
								else {
									vpairheap.down(vpairheap.rankinheap[target]);
								}
							}
						}
					}
					e = elist[elist[elist[e].next].next].paire;
				}
			}
		}
	}

	facenum -= trif.size();
}

void VPairHeap::up(PR pr) {
	while (pr > 0) {
		PR parent = ((pr - 1) >> 1);
		if (pairlist[pairheap[parent]].cost > pairlist[pairheap[pr]].cost) {
			swap(pairheap[parent], pairheap[pr]);
			swap(rankinheap[pairheap[parent]], rankinheap[pairheap[pr]]);
		}
		else
			break;
		pr = parent;
	}
}

void VPairHeap::down(PR pr) {
	while (pr < pairheap.size()) {
		PR lc = (pr << 1) + 1;
		PR rc = lc + 1;
		PR choose = pr; double min = pairlist[pairheap[pr]].cost;
		if (lc < pairheap.size() && pairlist[pairheap[lc]].cost < min) {
			choose = lc;
			min = pairlist[pairheap[lc]].cost;
		}
		if (rc < pairheap.size() && pairlist[pairheap[rc]].cost < min) {
			choose = rc;
			min = pairlist[pairheap[rc]].cost;
		}
		if (pr == choose)
			break;
		else {
			swap(pairheap[pr], pairheap[choose]);
			swap(rankinheap[pairheap[pr]], rankinheap[pairheap[choose]]);
			pr = choose;
		}
	}
}

VPair VPairHeap::top() {
	VPair ret = pairlist[pairheap[0]];
	swap(pairheap[0], pairheap[pairheap.size() - 1]);
	swap(rankinheap[pairheap[0]], rankinheap[pairheap[pairheap.size() - 1]]);
	rankinheap[pairheap[pairheap.size() - 1]] = -1;
	pairheap.pop_back();
	down(0);
	return ret;
}

