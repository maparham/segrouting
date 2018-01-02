#ifndef GUROBI_FLOWUPDATE_CPP_
#define GUROBI_FLOWUPDATE_CPP_

#include <algorithm>    // std::reverse
#include <vector>
#include <map>
#include "Snap.h"

using namespace std;
using namespace TSnap;

//#define __DEBUG__ 1
#ifdef __DEBUG__
#define PRINTF printf
#else
#define PRINTF(format, args...) ((void)0)
#endif

#define ANSI_COLOR_RESET "\033[0m"
#define RED(TXT)     MyStr("\033[1;31m" TXT ANSI_COLOR_RESET)
#define GREEN(TXT)   MyStr("\033[1;32m" TXT ANSI_COLOR_RESET)
#define YELLOW(TXT)  MyStr("\033[1;33m" TXT ANSI_COLOR_RESET)
#define BLUE(TXT)    MyStr("\033[1;34m" TXT ANSI_COLOR_RESET)

#define INF numeric_limits<int>::max()/3
#define WEIGHTATTR "weight"

//// what type of graph do you want to use?
typedef PNEANet MyGraph; // undirected graph
//typedef PNGraph PGraph;  //   directed graph
//typedef PNEGraph PGraph;  //   directed multigraph
//typedef TPt<TNodeNet<TInt> > PGraph;
//typedef TPt<TNodeEdgeNet<TInt, TInt> > PGraph;
typedef vector<int> Path;
typedef const int MyNode;
typedef const TNEANet::TEdgeI EdgeI;

typedef pair<int, int> Segment;
typedef vector<Segment> SegmentStack;
typedef map<int, SegmentStack> SSMap; // key is EId
typedef map<int, SSMap> Dest_Link_Table_2D; // key is destination id
// table for double links
typedef map<int, SSMap> SSMap_2D; // key is (Eid1,EId2)
typedef map<int, SSMap_2D> Dest_Link_Table_3D; // key is destination id

struct MyStr: string {
	const string s1;
	MyStr(const char* s) :
			s1(s) {
	}
	const char* operator+(const char* s2) {
		return (s1 + s2).c_str();
	}
};

struct Result {
	const int fail = 0;
	const int success = 0;
	const int maxStackSize = 0;
	Result(int a1, int a2, int a3) :
			fail(a1), success(a2), maxStackSize(a3) {
	}
};
struct Link {
	const int tail;
	const int head;
	const int id;
	Link(int tail, int head) :
			tail(tail), head(head), id(-1) {
	}
	Link(const TNEANet::TEdgeI& EI) :
			tail(EI.GetSrcNId()), head(EI.GetDstNId()), id(EI.GetId()) {
	}
	Link(int tail, int head, const MyGraph G) :
			tail(tail), head(head), id(G->GetEId(tail, head)) {
	}
	bool inPath(const Path& p) {
		for (int i = 1; i < p.size(); ++i) {
			if (p[i - 1] == tail && p[i] == head) {
				return true;
			}
		}
		return false;
	}
};

class TILFA {
protected:
	MyGraph G;
	Dest_Link_Table_2D table2D;
	Dest_Link_Table_3D table3D;
	TIntStrH NodeLabelH;

	double edgeWeight(const int u, const int v) {
		return G->GetFltAttrDatE(G->GetEId(u, v), WEIGHTATTR);
	}
	double edgeWeight(const int EId) {
		return G->GetFltAttrDatE(EId, WEIGHTATTR);
	}

	double setEdgeWeight(const int EId, const double val) {
		return G->AddFltAttrDatE(EId, val, WEIGHTATTR);
	}

	char* nodeLabel(const int NId) {
		return G->GetStrAttrDatN(NId, "label").CStr();
	}

	double getSP(MyNode src, MyNode dest, Path& sp) {

		sp.clear();
		if (src == dest) {
			sp.push_back(src);
			return 0;
		}

		TIntH SP_tree;
		GetWeightedShortestPathTree(G, src, SP_tree, WEIGHTATTR);

		if (!SP_tree.IsKey(dest)) { // no path found
			return INF;
		}

		int N = dest;
		double spLen = 0;
		while (SP_tree.IsKey(N) && N != src) {
			sp.push_back(N); // save the shortest path
			int parent = SP_tree.GetDat(N);
			int EID = G->GetEI(parent, N).GetId();
			double w = edgeWeight(parent, N);
//			PRINTF("w(%d,%d)=%f", parent, N, w);
//			PRINTF("isEdge(%s, %s)=%d", nodeLabel(parent), nodeLabel(N),
//					G->IsEdge(parent, N, true));

			Assert(w > 0);
			spLen += edgeWeight(EID);
			N = parent;
		}
		sp.push_back(src);

		std::reverse(sp.begin(), sp.end());
		return spLen;
	}

	bool inPSpace(MyNode node, const vector<Link>& FailLinks) {
		Path sp_;
		double len = getSP(FailLinks[0].tail, node, sp_);
		for (int i = 1; i < sp_.size(); ++i) {
			for (int j = 0; j < FailLinks.size(); ++j) {
				if (sp_[i - 1] == FailLinks[j].tail && sp_[i] == FailLinks[j].head) {
					PRINTF("%d not in PSpace of %d; ", node, FailLinks[0].tail);
//					PRINTF("SP(S=%d, %d)=%f: ", FailLinks[0].tail, node, len);
					print_path(sp_);
					return false;
				}
			}
		}
		PRINTF("%d in PSpace of %d;\n", node, FailLinks[0].tail);
		return true;
	}

	bool inQSpace(MyNode node, const vector<Link>& SF, MyNode D) {
		Path sp_;
		double len = getSP(node, D, sp_);

		for (int i = sp_.size() - 2; i > -1; --i) {
			for (int j = 0; j < SF.size(); ++j) {
				if (sp_[i] == SF[j].tail && sp_[i + 1] == SF[j].head) {
					PRINTF("%d not in QSpace of %d; ", node, D);
//					PRINTF("SP(%d,D=%d)=%f: ", node, D, len);
					print_path(sp_);
					return false;
				}
			}
		}
		PRINTF("%d in QSpace of %d;\n", node, D);
		return true;
	}

	Path SP(MyNode src, MyNode dest) {
		Path p;
		PRINTF("getSP=%f\n", getSP(src, dest, p));
		return p;
	}

	double backup_path(const Link& L, const SegmentStack& stack, Path& bp,
			vector<int>& destinationMap) {
		// invariant: destination map always covers bp before its last node
		SegmentStack stack_cpy = stack;
		PRINTF("constructing BP for L=(%d,%d), stack.size()=%d\n", L.tail, L.head, stack.size());
		bp.push_back(L.tail);	// the first SP node
		double bpLen = 0;
		while (!stack_cpy.empty()) {
			Segment seg = stack_cpy.back();
			PRINTF("seg=(%d,%d)\n", seg.first, seg.second);
			stack_cpy.pop_back();
			if (seg.first == seg.second) { // an intermediate destination (repair node)
				int S_ = bp.back();
				Path sp;
				bpLen += getSP(S_, seg.first, sp);
				PRINTF("getSP(%d, %d)=%f\n", S_, seg.first, bpLen);
				bp.insert(bp.end(), sp.begin() + 1, sp.end()); // append the rest of the shortest path
				destinationMap.resize(bp.size() - 1, seg.second); // set destination for all the SP out_links

			} else { // a tunnel link
				PRINTF("seg.first=%d, bp.back()=%d\n", seg.first, bp.back());
				Assert(seg.first == bp.back());	// the link's tail must be the last node already
				bp.push_back(seg.second);	// the head node (q)
				destinationMap.push_back(seg.second);	// for the tail node (p)
				bpLen += edgeWeight(seg.first, seg.second);
				PRINTF("edgeWeight(%d,%d)=%f\n", seg.first, seg.second, edgeWeight(seg.first, seg.second));

			}
		}
		Assert(bp.size() != 1);
		return bpLen;
	}

public:
	int nofLinks;
	int nofNodes;

	void print_path(Path& sp) {
		for (int i = 0; i < sp.size() - 1; ++i) {
			PRINTF("%d,", sp[i]);
		}
		if (sp.size()) {
			PRINTF("%d\n", sp[sp.size() - 1]);
		}
	}
	void load_graph(const string path) {
		// load graph from file
		TStrHash<TInt> InStrToNIdH;
		TStr pathobj(path.c_str());
		G = LoadWeightedEdgeListStr<MyGraph>(pathobj, 0, 1, 2, WEIGHTATTR,
				InStrToNIdH);

		for (auto NI = G->BegNI(); NI < G->EndNI(); NI++) {
			NodeLabelH.AddDat(NI.GetId(), InStrToNIdH.GetKey(NI.GetId()));
		}
		nofLinks = G->GetEdges();
		nofNodes = G->GetNodes();
	}

	void drawGraph(string path) {
		TStr pathobj(path.c_str());

		PRINTF("drawing the graph in %s\n", path.c_str());
		DrawGViz(G, TGVizLayout::gvlDot,
				pathobj.GetFPath() + pathobj.GetFMid() + ".png", "Loaded Graph",
				true);
		DrawGViz(G, TGVizLayout::gvlDot,
				pathobj.GetFPath() + pathobj.GetFMid() + "_labeled.png",
				"Loaded Graph", NodeLabelH);
	}

	void printInfo() {
		PRINTF("%d nodes, %d links\n", nofNodes, nofLinks);
	}

	// computes the segment stack w.r.t D and the first link in the given list
	// knowing the other links in the list fail as well
	bool compute_SegmentStack(vector<TNEANet::TEdgeI> EIList,
			SegmentStack& stack, const int D) {

		int S1 = EIList[0].GetSrcNId();
		int F1 = EIList[0].GetDstNId();
//		int S2 = EIList[1].GetSrcNId();
//		int F2 = EIList[1].GetDstNId();
		if (S1 == D) { // skip this destination
			return false;
		}

		Path sp;
		double len = getSP(S1, D, sp);
		if (len >= INF) {
			return false;
		}

		vector<Link> failed_links;
		PRINTF("compute_SegmentStack: ");
		for (int i = 0; i < EIList.size(); ++i) {
			failed_links.push_back(Link(EIList[i]));
			PRINTF("(%d, %d), ", failed_links[i].tail, failed_links[i].head);
		}
		PRINTF("D=%d\n", D);

		PRINTF("first shortest path SP(%d,%d).weight=%f: ", S1, D, len);
		print_path(sp);
		if (sp[1] != F1) {
			PRINTF("SP not using the links\n");
			return false;
		}

		vector<double> weights(EIList.size());

		for (int i = 0; i < EIList.size(); ++i) {
			weights[i] = edgeWeight(EIList[i].GetId());
//			PRINTF("disable link (%d,%d)\n", EIList[i].GetSrcNId(), EIList[i].GetDstNId());
			setEdgeWeight(EIList[i].GetId(), INF); // disable the link=
		}

		len = getSP(S1, D, sp); // shortest path without failed links,i.e., post convergence path

		for (int i = 0; i < EIList.size(); ++i) {
			setEdgeWeight(EIList[i].GetId(), weights[i]); // enable the link
//			PRINTF("enable link (%d,%d), w=%f\n", EIList[i].GetSrcNId(), EIList[i].GetDstNId(),
//					edgeWeight(EIList[i].GetId()));
		}

		PRINTF("second shortest path weight=%f: ", len);
		print_path(sp);
		if (len >= INF) {
			PRINTF("no backup path\n");
			return false;
		}
		Assert(sp.size() > 1);
		int PSpace_sup = S1, QSpace_inf = D;
		for (int i = 0; i < sp.size() && inPSpace(sp[i], failed_links);
				++i) {
			PSpace_sup = i;
		}
		for (int i = sp.size() - 1;
				i >= PSpace_sup && inQSpace(sp[i], failed_links, D);
				--i) { // search backward for the infimum of QSpace, stop when it touches the PSpace
			QSpace_inf = i;
		}

		MyNode p = sp[PSpace_sup], q = sp[QSpace_inf];
		PRINTF("p=%d q=%d\n", p, q);

		if (q != D) {
			stack.push_back(Segment(D, D));
		}
		int segDest = D;
		// connect p to q with segments
		for (int x = QSpace_inf - 1;
				x >= PSpace_sup; --x) {
			if (!inQSpace(sp[x], failed_links, segDest)) {
				stack.push_back(Segment(sp[x], sp[x + 1]));
				if (x > 0) { // don't push S1
					stack.push_back(Segment(sp[x], sp[x]));
				}
				segDest = stack.back().first;	// destination for the next segment
			}
		}
		if (p != S1 && stack.back().first != p) {	// p must be the top of stack unless...
			stack.push_back(Segment(p, p));
		}

		PRINTF("validating stack for (%d,%d), D=%d\n", S1, F1, D);
		Assert(stack.size() > 1 || stack[0].first == S1 && stack[0].second == D);
		for (int i = stack.size() - 1; i > -1; --i) {
			Segment seg = stack[i];
			PRINTF("checking segment (%d,%d)\n", seg.first, seg.second);
			if (seg.first != seg.second) {
				Assert(
						seg.first == S1
								|| seg.first == stack[i + 1].second);
			}
		}

		return true;
	}

	void compute_SR_table() {
		for (TNEANet::TEdgeI EI = G->BegEI(); EI < G->EndEI(); EI++) {

		} // next link
	}

	// helper
	int forEachLinkonBP(const TNEANet::TEdgeI EI,
			MyNode D,
			function<void(EdgeI&, const int)> fn) {

		// segment stack for single-link failure
		SegmentStack stack = table2D[D][EI.GetId()];
		bool possible = compute_SegmentStack( { EI },
				stack, D);
		if (!possible) { // graph disconnected
			return INF;
		}

		int S = EI.GetSrcNId();
		int F = EI.GetDstNId();
		Link L1 = Link(EI);
		Path bp;
		vector<int> destMap;
		PRINTF(">>>> L1=(%d,%d), D=%d: ", S, F, stack[0].first);
		double len = backup_path(L1, stack, bp, destMap);
		print_path(bp);
		Assert(L1.inPath(bp) == false);

		// all second failures on the first backup path
		for (int i = 1; i < bp.size(); ++i) {
			EdgeI EI2 = G->GetEI(bp[i - 1], bp[i]);
			fn(EI2, destMap[i - 1]); // the link on BP and its (possibly intermediate) destination
		}
		return stack.size();
	}

	Result eval_double_failure() {
		int fail = 0, success = 0, maxSS = 0;
		for (TNEANet::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
			int D = NI.GetId();

			for (TNEANet::TEdgeI EI1 = G->BegEI(); EI1 < G->EndEI();
					EI1++) {
				Link L1 = Link(EI1);
				double len = 0;
				// all second failures on the first backup path
				int SS = forEachLinkonBP(EI1, D,
						[&](EdgeI& EI2, const int D1) {
							Link L2=Link(EI2);
							PRINTF(">> L2=(%d,%d), D=%d: ", L2.tail, L2.head,D1);

							// segment stack for double-link failure
							//SegmentStack stack=table2D[D1][L2.id];	// case1: keep the current stack
						SegmentStack stack=table2D[D][L2.id];// case2: flush the stack down to D (means possibly more rout options)

						if( compute_SegmentStack( {EI2, EI1}, stack, D1)==false) {
							return; // nonsense cases
						}
						Path bp2;
						vector<int> destinationMap2;
						len += backup_path(L2, stack, bp2,
								destinationMap2); // case1: keep the current stack
						print_path(bp2);
						PRINTF("weight=%f, INF=%d \n",len, INF);
						Assert(len<INF );
						if (L1.inPath(bp2)) {
							++fail;
							PRINTF("++fail=%d\n", fail);
						} else {
							Assert(bp2.size()>1);
							++success;
							PRINTF("++success=%d\n", success);
						}
						if( maxSS < stack.size()) {
							maxSS = stack.size();
						}
					});
				if (SS < INF && SS > maxSS) {
					maxSS = SS;
				}
			}
		}
		return Result(fail, success, maxSS);
	}
};

#endif
