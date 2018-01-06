#ifndef GUROBI_FLOWUPDATE_CPP_
#define GUROBI_FLOWUPDATE_CPP_

//#define __DEBUG__ 1

#include <algorithm>    // std::reverse
#include <vector>
#include <map>
#include "Snap.h"
#include <utils.h>
#include <ctime>

using namespace std;
using namespace TSnap;

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
typedef map<string, SegmentStack> Dest_Link_Table;
typedef map<string, bool> TraceMap;

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
	double failRatio() {
		double ratio = (double) fail / (fail + success);
		return ratio;
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
	int inPath(const Path& p) {
		for (int i = 1; i < p.size(); ++i) {
			if (p[i - 1] == tail && p[i] == head) {
				return i - 1;
			}
		}
		return -1;
	}
};

class TILFA {
protected:
	MyGraph G;
	Dest_Link_Table table;
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
		int w = getSP(src, dest, p);
		PRINTF("getSP=%f\n", w);
		return p;
	}

	double backup_path(const EdgeI& L, const SegmentStack& stack, Path& bp,
			vector<int>& destinationMap) {
		// invariant: destination map always covers bp before its last node
		SegmentStack stack_cpy = stack;
		PRINTF("constructing BP for L=(%d,%d), stack.size()=%d\n", L.GetSrcNId(), L.GetDstNId(), stack.size());
		bp.push_back(L.GetSrcNId());	// the first SP node
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
		printf("%d nodes, %d links\n", nofNodes, nofLinks);
	}

	// computes the segment stack w.r.t D and the first link in the given list
	// knowing the other links in the list fail as well
	bool compute_SegmentStack(vector<EdgeI> EIList, EdgeI& current,
			SegmentStack& stack) {

		Assert(stack.size() == 1); // the destination must be there already
		const int dest = stack[0].second;

		int S1 = current.GetSrcNId();
		int F1 = current.GetDstNId();
		if (S1 == dest) { // skip this destination
			return false;
		}

		vector<double> weights(EIList.size());
		double len;
		PRINTF("compute_SegmentStack wrt (%d,%d), D=%d\n", current.GetSrcNId(), current.GetDstNId(), dest);

		vector<Link> failed_links;
		for (int i = 0; i < EIList.size(); ++i) {
			failed_links.push_back(Link(EIList[i]));
			PRINTF("(%d, %d), ", failed_links[i].tail, failed_links[i].head);
		}
		PRINTF("D=%d\n", dest);

		for (int i = 0; i < EIList.size(); ++i) {
			weights[i] = edgeWeight(EIList[i].GetId());
			PRINTF("disable link (%d,%d)\n", EIList[i].GetSrcNId(), EIList[i].GetDstNId());
			setEdgeWeight(EIList[i].GetId(), INF); // disable the link
		}

		Path sp;
		len = getSP(S1, dest, sp); // shortest path without all failed links,i.e., post convergence path

		for (int i = 0; i < EIList.size(); ++i) {
			setEdgeWeight(EIList[i].GetId(), weights[i]); // enable the link
			PRINTF("enable link (%d,%d), w=%f\n", EIList[i].GetSrcNId(), EIList[i].GetDstNId(),
					edgeWeight(EIList[i].GetId()));
		}

		PRINTF("second shortest path weight=%f: ", len);
		print_path(sp);
		if (len >= INF) {
			PRINTF("no backup path\n");
			return false;
		}
		Assert(sp.size() > 1);

		if (sp[0] == S1 && sp[1] == dest) {	// sp is only a single link => don't need the destination segment
			stack.clear(); // remove the destination
		}

		int PSpace_sup = S1, QSpace_inf = dest;
		for (int i = 0; i < sp.size() && inPSpace(sp[i], failed_links);
				++i) {
			PSpace_sup = i;
		}
		for (int i = sp.size() - 1;
				i >= PSpace_sup && inQSpace(sp[i], failed_links, dest);
				--i) { // search backward for the infimum of QSpace, stop when it touches the PSpace
			QSpace_inf = i;
		}

		MyNode p = sp[PSpace_sup], q = sp[QSpace_inf];
		PRINTF("p=%d q=%d\n", p, q);

//		if (q != dest) {
//			stack.push_back(Segment(dest, dest));
//		}

		int segDest = dest;
		// decide segments from q to p, p < q
		for (int x = QSpace_inf - 1;
				x >= PSpace_sup; --x) {
			if (!inQSpace(sp[x], failed_links, segDest)) {
				stack.push_back(Segment(sp[x], sp[x + 1]));
				if (x > 0) { // don'4 push S1
					stack.push_back(Segment(sp[x], sp[x]));
				}
				segDest = stack.back().first;	// destination for the next segment
			}
		}

		if (p == dest) {	// in case of an inconsistent shortest path algorithm, ie, sp is a first shortest path
//			Assert(stack.size() == 1);
				//stack.push_back(Segment(dest, dest));
			stack.push_back(Segment(sp[0]/*S1*/, sp[1]));	// force the packet at S1 take this alternative SP
			PRINTF("inconsistent SP=> pushing D=%d and first link\n", dest);

		} else if (p != S1 && stack.back().first != p) {	// p must be the top of stack unless...
			stack.push_back(Segment(p, p));
		}

		Assert(stack.size() > 0);

		PRINTF("validating stack of size %d for (%d,%d), D=%d\n", stack.size(), S1, F1, dest);
		Assert(stack.size() > 1 || stack[0].first == S1 && stack[0].second == dest);
		for (int i = stack.size() - 1; i > -1; --i) {
			Segment seg = stack[i];
			PRINTF("checking segment (%d,%d)\n", seg.first, seg.second);
			if (seg.first != seg.second) {
				Assert(
						seg.first == S1
								|| seg.first == stack[i + 1].second);
			}
		}
		Assert(p == S1 || p == dest || stack.size() > 1);

		return true;
	}

	string genKey(const vector<EdgeI>& links, const int& dest) {
		stringstream jId;
		jId << dest << '_';
		for (int i = links.size() - 1; i > -1; --i) {
			jId << links[i].GetId() << "_";
		}
		return jId.str();
	}

	// compute or retrieve stack w.r.t. given known failure
	const SegmentStack& getSegmentStack(const vector<EdgeI>& links, EdgeI& current, const int& dest) {
		PRINTF("getSegmentStack wrt (%d,%d),(%d,%d), D=%d\n", current.GetSrcNId(), current.GetDstNId(),
				(links.size() > 1) ? links[1].GetSrcNId() : -1, (links.size() > 1) ? links[1].GetDstNId() : -1, dest);

#ifndef doubleTILFA
		Assert(links.size() == 1);
#endif
		Assert(links[0].GetId() == current.GetId());

		// segment stack for single-link failure
		SegmentStack& stack = table[genKey(links, dest)];	// retrieve a cached stack or instantiate a new one

		if (stack.size() > 0) {
			PRINTF("already computed, stack.size()=%d\n", stack.size());
			return stack;
		}

		stack.push_back(Segment(dest, dest));
		bool possible = compute_SegmentStack(links, current, stack);

		if (!possible) {
			PRINTF("not possible!\n");
			Assert(stack.size() == 1);
			stack.clear();
		}
		return stack;
	}

	// must be called when the current link is on the BP of the other failed link.
	// traces a possible loop starting at the "current" link, w.r.t dest.
	// this check matter only when the knowledge of other failures is not available.
	//
	bool isFail(const vector<EdgeI>& fails, EdgeI& current, int dest, int depth = 0, TraceMap& trace =
			*(new TraceMap())) {

		PRINTF("isFail? current=(%d,%d),dest=%d\n", current.GetSrcNId(), current.GetDstNId(), dest);
		bool& seen = trace[genKey( { current }, dest)];
		if (seen) {
			PRINTF("packet already has seen the link and the destination => is loop\n");
			return true;
		}
		seen = true;

#ifdef FLUSH_STACK
		return true;
#endif
#ifdef doubleTILFA
		exit(0); // should never happen
#endif

		SegmentStack stack = getSegmentStack( { current }, current, dest);
		if (stack.size() == 0) { // no backup rout => drop packet
			return true;
		}
		if (depth > 10) {
			PRINTF("recursion too deep: %d\n => fail", depth);
			exit(11);
			return true;
		}
		Path bp;
		vector<int> destinationMap;
		double len = backup_path(current, stack, bp, destinationMap);
		Assert(len <INF && bp.size()>1);

		bool isLoop = false;
		int pos;
		for (EdgeI L : fails) {
			pos = Link(L).inPath(bp);
			if (pos >= 0) {
				isLoop = isLoop || isFail(fails, L, destinationMap[pos], ++depth, trace);
			}
		}
		//if(!isLoop) exit(0);
		return isLoop;
	}

// helper; returns number of stack items
	int forEachLinkonBP(const TNEANet::TEdgeI EI,
			MyNode D,
			function<void(EdgeI&, const int)> fn) {

		int S = EI.GetSrcNId();
		int F = EI.GetDstNId();
		Link L1 = Link(EI);
		Path sp;

		double len = getSP(S, D, sp); // SP in the original graph

		if (!(len < INF)) {
			return INF;
		}
		PRINTF("first shortest path SP(%d,%d).weight=%f: ", S, D, len);
		print_path(sp);
		Assert(sp.size() > 1);
		if (sp[1] != F) {	// the the link is irrelevant to D
			PRINTF("first SP not using the link => irrelevant case\n");
			return INF;
		}

		SegmentStack stack = getSegmentStack( { EI }, EI, D);
		if (stack.size() == 0) { // no backup path
			PRINTF("empty stack returned\n");
			return INF;
		}

		Path bp;
		vector<int> destMap;
		PRINTF(">> L1=(%d,%d), D=%d: ", S, F, stack[0].first);
		len = backup_path(EI, stack, bp, destMap);
		print_path(bp);
		Assert(L1.inPath(bp) < 0);

		// all second failures on the first backup path
		for (int i = 1; i < bp.size(); ++i) {
			EdgeI EI2 = G->GetEI(bp[i - 1], bp[i]);
			fn(EI2, destMap[i - 1]); // the link on BP and its (possibly intermediate) destination
		}
		return stack.size();
	}

	Result eval_double_failure(Reporter& report) {
		int fail = 0, success = 0, maxSS = 0;
		for (TNEANet::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
			int D = NI.GetId();

			for (TNEANet::TEdgeI EI1 = G->BegEI(); EI1 < G->EndEI();
					EI1++) {

				if (EI1.GetSrcNId() == D) {
					continue;
				}

				Link L1 = Link(EI1);
				double len = 0;

				// all second failures on the first backup path
				int SS = forEachLinkonBP(EI1, D,
						[&](EdgeI& EI2, const int D1) {
							Link L2=Link(EI2);
							PRINTF(">> L2=(%d,%d), D=%d,D1=%d: ", L2.tail, L2.head,D,D1);

#ifdef FLUSH_STACK
						const int dest = D;	// case1: flush the stack down to D (means possibly more rout options)
#else
						const int dest = D1;	// case2: keep the current stack
#endif

						vector<EdgeI> fails= {EI2, EI1}; // in reverse order of failure time

						// failures known at this point
#ifdef doubleTILFA
						vector<EdgeI> knownFails= fails;	// current and previous fail
#else
						vector<EdgeI> knownFails = {fails[0]}; // only the current fail
#endif

						SegmentStack stack = getSegmentStack(knownFails,fails[0],dest);

						if(stack.size()==0) { // graph disconnected or D=S
							return;
						}
						Path bp2;
						vector<int> destinationMap2;
						len += backup_path(EI2, stack, bp2,
								destinationMap2); // case1: keep the current stack
						Assert(len<INF );

						print_path(bp2);
						PRINTF("weight=%f, INF=%d \n",len, INF);

						if (L1.inPath(bp2)>=0 && isFail(fails,EI1,dest )) { // never true for double-link TILFA
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
						report << stack.size()-1 << '\n';
					});

				if (SS < INF) {
					report << SS - 1 << '\n';
					if (SS > maxSS) {
						maxSS = SS;
					}
				}
			} // next E1
		} //next dest
		return Result(fail, success, maxSS - 1);
	}
};

#endif
