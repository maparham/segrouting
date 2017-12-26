#include <execinfo.h>
#include <signal.h>
#include <iostream>
#include "Snap.h"
#include <fstream>
//#include <unordered_map>
#include <algorithm>    // std::reverse
#include <vector>
#include <map>
#include <time.h>

using namespace std;
using namespace TSnap;

#define __DEBUG__ 1
#define INF numeric_limits<int>::max()/3
#define WEIGHTATTR "weight"

#ifdef __DEBUG__
#define PRINTF printf
#else
#define PRINTF(format, args...) ((void)0)
#endif

void handler(int sig) {
	void* array[10];
	size_t size;

	// get void*'s for all entries on the stack
	size = backtrace(array, 10);

	// print out all the frames to stderr
	fprintf(stderr, "Error: signal %d:\n", sig);
	backtrace_symbols_fd(array, size, STDERR_FILENO);
	exit(1);
}

//// what type of graph do you want to use?
typedef PNEANet MyGraph; // undirected graph
//typedef PNGraph PGraph;  //   directed graph
//typedef PNEGraph PGraph;  //   directed multigraph
//typedef TPt<TNodeNet<TInt> > PGraph;
//typedef TPt<TNodeEdgeNet<TInt, TInt> > PGraph;
typedef vector<int> Path;
typedef const int MyNode;
//// what type of graph do you want to use?
typedef PNEANet MyGraph; // undirected graph
//typedef PNGraph PGraph;  //   directed graph
//typedef PNEGraph PGraph;  //   directed multigraph
//typedef TPt<TNodeNet<TInt> > PGraph;
//typedef TPt<TNodeEdgeNet<TInt, TInt> > PGraph;

typedef pair<int, int> Segment;
typedef vector<Segment> SegmentStack;
typedef map<int, SegmentStack> BackupPaths; // key is link id
typedef map<int, BackupPaths> Dest_Link_Table; // key is destination id

struct Logger: std::ostream {
	template<typename T>
	Logger& operator <<(const T& x) {
#ifdef DEBUG
		std::cout << x;
		std::cout.flush();
#endif
		return *this;
	}
} mylog;

struct Link {
	const int tail;
	const int head;
	const int id;
	Link(int tail, int head) :
			tail(tail), head(head), id(-1) {
	}
	Link(int tail, int head, const MyGraph G) :
			tail(tail), head(head), id(G->GetEId(tail, head)) {
	}
	bool inPath(const Path& p) {
		mylog << "inPath?\n";
		for (int i = 1; i < p.size(); ++i) {
			if (p[i - 1] == tail && p[i] == head) {
				return true;
			}
		}
		return false;
	}
};

template<class Graph>
int getSP(Graph G, MyNode src, MyNode dest, Path& sp) {

	sp.clear();
	TIntH SP_tree;
	GetWeightedShortestPathTree(G, src, SP_tree, WEIGHTATTR);

//	TIntFltH NIdDistH;
//	TFltV Attr;
//	Attr.Reserve(99999, G->GetEdges());
//	Attr[0] = Attr[1] = Attr[2] = Attr[3] = Attr[4] = Attr[5] = 1;

//	int EId = G->GetEI(1, dest).GetId();
//	mylog << "EID found=" << EId;
//	Attr[EId] = 11;
//	GetWeightedShortestPath(G, src, NIdDistH, Attr);
//	mylog << "\nSP=" << NIdDistH.GetDat(dest) << "\n";

	int N = dest;
	int spLen = 0;
	while (SP_tree.IsKey(N) && N != src) {
		sp.push_back(N); // save the shortest path
		int parent = SP_tree.GetDat(N);
		int EID = G->GetEI(parent, N).GetId();
		spLen += G->GetIntAttrDatE(EID, WEIGHTATTR);
		N = parent;
	}
	sp.push_back(src);

	std::reverse(sp.begin(), sp.end());
	return spLen;
}

template<class Graph>
Path SP(Graph G, MyNode src, MyNode dest) {
	Path p;
	PRINTF("getSP=%d\n", getSP(G, src, dest, p));
	return p;
}
void print_path(Path& sp) {
	for (int i = 0; i < sp.size() - 1; ++i) {
		PRINTF("%d,", sp[i]);
	}
	if (sp.size()) {
		PRINTF("%d\n", sp[sp.size() - 1]);
	}
}

bool inPSpace(MyGraph G, MyNode node, MyNode S, MyNode F) {
	Path sp_;
	int len = getSP(G, S, node, sp_);
	for (int i = 1; i < sp_.size(); ++i) {
		if (sp_[i - 1] == S && sp_[i] == F) {
			PRINTF("%d not in PSpace; ", node);
			PRINTF("SP(S, %d)=%d: ", node, len);
			print_path(sp_);
			return false;
		}
	}
	return true;
}

bool inQSpace(MyGraph G, MyNode node, MyNode S, MyNode F, MyNode D) {
	Path sp_;
	int len = getSP(G, node, D, sp_);
	for (int i = sp_.size() - 2; i > -1; --i) {
		if (sp_[i] == S && sp_[i + 1] == F) {
			PRINTF("%d not in QSpace; ", node);
			PRINTF("SP(%d,D)=%d: ", node, len);
			print_path(sp_);
			return false;
		}
	}
	return true;
}

void load_graph(const string path, MyGraph& G) {
	// load graph from file
	TStrHash<TInt> InStrToNIdH;
	TStr pathobj(path.c_str());
	G = LoadWeightedEdgeListStr<MyGraph>(pathobj, 0, 1, 2, WEIGHTATTR,
			InStrToNIdH);

	TIntStrH NodeLabelH;
	for (auto NI = G->BegNI(); NI < G->EndNI(); NI++) {
		NodeLabelH.AddDat(NI.GetId(), InStrToNIdH.GetKey(NI.GetId()));
	}

//	PRINTF("drawing the graph in %s\n",pathobj.GetFPath().CStr());
//	DrawGViz(G, TGVizLayout::gvlDot, pathobj.GetFPath() + "graph.png",
//			"Loaded Graph", true);
//	DrawGViz(G, TGVizLayout::gvlDot, pathobj.GetFPath() + "graph_labeled.png",
//			"Loaded Graph", NodeLabelH);
}

void compute_SR_table(MyGraph G, Dest_Link_Table& table) {
	//  iterate over all destinations
	for (TNEANet::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
		int dest = NI.GetId();
		BackupPaths bp; // for a fixed destination
		table[dest] = bp;
		for (TNEANet::TEdgeI EI = G->BegEI(); EI < G->EndEI(); EI++) {
			int S = EI.GetSrcNId();
			int F = EI.GetDstNId();
			if (S == dest) { // nothing to do
				continue;
			}
			int D = dest; // iterate on all possible destinations
			PRINTF("link (%d, %d), D=%d\n", S, F, D);

			Path sp;
			int len = getSP(G, S, D, sp);
			if (len >= INF) {
				continue;
			}
			PRINTF("first shortest path weight=%d: ", len);
			print_path(sp);
			if (sp[1] != F) {
				mylog << "SP not using the link\n";
				continue;
			}

			const int w = G->GetIntAttrDatE(EI.GetId(), WEIGHTATTR);
			G->AddIntAttrDatE(EI.GetId(), INF, WEIGHTATTR); // cripple the link (S,F)
			len = getSP(G, S, D, sp); // shortest path without (S,F) => post convergence path
			G->AddIntAttrDatE(EI.GetId(), w, WEIGHTATTR);
			PRINTF("second shortest path weight=%d: ", len);
			print_path(sp);
			if (len >= INF) {
				mylog << "no backup path\n";
				continue;
			}

			int PSpace_sup = S, QSpace_inf = D;
			for (int i = 0; i < sp.size() && inPSpace(G, sp[i], S, F); ++i) {
				PRINTF("%d in PSpace\n", sp[i]);
				PSpace_sup = i;
			}
			for (int i = sp.size() - 1;
					i >= PSpace_sup && inQSpace(G, sp[i], S, F, D); --i) { // search backward for the infimum of QSpace, stop when it touches the PSpace
				PRINTF("%d in QSpace\n", sp[i]);
				QSpace_inf = i;
			}

			MyNode p = sp[PSpace_sup], q = sp[QSpace_inf];
			PRINTF("p=%d q=%d\n", p, q);
			//	Assert(QSpace_inf == PSpace_sup + 1 || PSpace_sup == QSpace_inf);

			// deduce the necessary segments by comparing p and q
			SegmentStack stack = { Segment(D, D), Segment(q, q) };
			if (QSpace_inf > PSpace_sup) { // also push a "tunnel link" and p
				stack.push_back(Segment(sp[PSpace_sup], sp[PSpace_sup + 1]));
			}
			stack.push_back(Segment(p, p));

			PRINTF("validating stack for (%d,%d), D=%d\n", S, F, D);
			for (int i = stack.size() - 1; i > -1; --i) {
				Segment seg = stack[i];
				PRINTF("checking segment (%d,%d)\n", seg.first, seg.second);
				if (seg.first != seg.second) {
					Assert(seg.first == S || seg.first == stack[i + 1].second);
				}
			}

			table[dest][EI.GetId()] = stack;
			PRINTF("stack size=%d\n", stack.size());
			Assert(Link(S, F).inPath(SP(G, sp[PSpace_sup + 1], q)) == false); // requires a rigorous proof
		} // next link
	} // next destination
}

int backup_path(MyGraph G, const Link& L, const int D, Dest_Link_Table& table,
		Path& bp, vector<int>& destinationMap) {
	PRINTF("computing BP for L=(%d,%d)\n", L.tail, L.head);
	bp.push_back(L.tail);	// the first SP node
	int bpLen = 0;
	SegmentStack stack = table[D][L.id];
	while (!stack.empty()) {
		Segment seg = stack.back();
		PRINTF("seg=(%d,%d)\n", seg.first, seg.second);
		stack.pop_back();
		if (seg.first == seg.second) { // an intermediate destination (repair node)
			int S_ = bp.back();
			Path sp;
			bpLen += getSP(G, S_, seg.second, sp);
			bp.insert(bp.end(), sp.begin() + 1, sp.end()); // append the rest of the shortest path
			destinationMap.resize(bp.size() - 1, seg.second); // set destination for all the SP out_links

		} else { // a tunnel link
			PRINTF("seg.first=%d, bp.back()=%d\n", seg.first, bp.back());
			Assert(seg.first == bp.back());	// the link's tail must be the last node already
			bp.push_back(seg.second);	// the head node
//			destinationMap.push_back(seg.second);	// for the tail node
			bpLen += G->GetIntAttrDatE(G->GetEId(seg.first, seg.second),
			WEIGHTATTR);
		}
	}
}

void TILFA_double_failure(MyGraph G, Dest_Link_Table& table) {
	int fail = 0, success = 0;
	for (TNEANet::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
		int D = NI.GetId();
		for (TNEANet::TEdgeI EI1 = G->BegEI(); EI1 < G->EndEI(); EI1++) {
			int S = EI1.GetSrcNId();
			int F = EI1.GetDstNId();
			Link L1 = Link(S, F, G);
			Path bp1;
			vector<int> destinationMap1;
			PRINTF("L1=(%d,%d), D=%d: ", S, F, D);
			int len = backup_path(G, L1, D, table, bp1, destinationMap1);
			print_path(bp1);
			Assert(L1.inPath(bp1) == false);
			// all second failures on the first backup path
			for (int i = 1; i < bp1.size(); ++i) {
				Link L2 = Link(bp1[i - 1], bp1[i], G);
				Path bp2;
				vector<int> destinationMap2;
				PRINTF("L2=(%d,%d), D=%d: ", L2.tail, L2.head,
						destinationMap1[S]);
				len += backup_path(G, L2, destinationMap1[S], table, bp2,
						destinationMap2);	// case1: keep the current stack
				//len += backup_path(G, e2, D, table, bp); // case2: flush the stack down to D
				print_path(bp2);
				if (L1.inPath(bp2)) {
					++fail;
					PRINTF("++fail=%d\n", fail);
				} else {
					++success;
					PRINTF("++success=%d\n", success);
				}
			}
		}
	}
	printf("fail=%d, success=%d\n", fail, success);
}

int main() {
	signal(SIGSEGV, handler); // install our handler

	/* initialize random seed: */
	srand(time(NULL));
//	ifstream topgen;
//	topgen.open("/Users/mahmoud/Downloads/weights-dist/1221/latencies.intra");

//	if (!topgen.is_open()) {
//		mylog("file not opened.");
//		return nodes;
//	}

//	G= TSnap::GenRndGnm<MyGraph>(10,20);

	std::ifstream filelist("files.txt");
	string path;
	while (filelist >> path) {
		PRINTF(" Loading %s\n", path.c_str());
		MyGraph G;
		load_graph(path, G);
		int nofLinks = G->GetEdges();
		int nofNodes = G->GetNodes();
		PRINTF("Loaded graph: %d nodes, %d links\n", nofNodes, nofLinks);

		Dest_Link_Table table;
		compute_SR_table(G, table);
		PRINTF("table size=%u\n", table.size());

		TILFA_double_failure(G, table);
	}

	const TStr inFile =
			"/Users/mahmoud/Downloads/weights-dist/3967/latencies.intra",
			inFile1 = "in.txt";

//	int randNode = sp[rand() % (sp.size() - 1)];
//	PRINTF("randNode=%d", randNode);
//	int failedLink = G->GetEI(sp[randNode], sp[randNode + 1]).GetId();
//	PRINTF("failedLink=%d\n", failedLink);
//	G->AddIntAttrDatE(failedLink, INF,WEIGHTATTR);
//	getSP(G, src, dest, sp);
//	printsp(sp);

	return 0;
}
