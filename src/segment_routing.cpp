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

#define DEBUG 1
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
//// what type of graph do you want to use?
typedef PNEANet MyGraph; // undirected graph
//typedef PNGraph PGraph;  //   directed graph
//typedef PNEGraph PGraph;  //   directed multigraph
//typedef TPt<TNodeNet<TInt> > PGraph;
//typedef TPt<TNodeEdgeNet<TInt, TInt> > PGraph;

struct Segment {
	bool isNode; // otherwise link
	int id;
};
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

void printsp(Path& sp) {
	for (int i = 0; i < sp.size() - 1; ++i) {
		printf("%d,", sp[i]);
	}
	if (sp.size()) {
		printf("%d\n", sp[sp.size() - 1]);
	}
}

bool inPSpace(MyGraph& G, MyNode node, MyNode S, MyNode F) {
	Path sp_;
	getSP(G, S, node, sp_);
	for (int i = 1; i < sp_.size(); ++i) {
		if (sp_[i - 1] == S && sp_[i] == F) {
			printf("%d not in PSpace; ", node);
			printf("SP(S, %d)=", node);
			printsp(sp_);
			return false;
		}
	}
	return true;
}

bool inQSpace(MyGraph& G, MyNode node, MyNode S, MyNode F, MyNode D) {
	Path sp_;
	getSP(G, node, D, sp_);
	for (int i = sp_.size() - 2; i > -1; --i) {
		if (sp_[i] == S && sp_[i + 1] == F) {
			printf("%d not in QSpace; ", node);
			printf("SP(%d,D)=", node);
			printsp(sp_);
			return false;
		}
	}
	return true;
}

int main() {
	/* initialize random seed: */
	srand(time(NULL));
//	ifstream topgen;
//	topgen.open("/Users/mahmoud/Downloads/weights-dist/1221/latencies.intra");

//	if (!topgen.is_open()) {
//		mylog("file not opened.");
//		return nodes;
//	}

//	G= TSnap::GenRndGnm<MyGraph>(10,20);

	const TStr inFile =
			"/Users/mahmoud/Downloads/weights-dist/3967/latencies.intra",
			inFile1 = "in.txt";
// load graph from file
	TStrHash<TInt> InStrToNIdH;
	MyGraph G = LoadWeightedEdgeListStr<MyGraph>(inFile, 0, 1, 2, WEIGHTATTR,
			InStrToNIdH);
	int nofLinks = G->GetEdges();
	int nofNodes = G->GetNodes();
	printf("Loaded graph: %d nodes, %d links\n", nofNodes, nofLinks);

	TIntStrH NodeLabelH;
	int dest;
	for (auto NI = G->BegNI(); NI < G->EndNI(); NI++) {
		NodeLabelH.AddDat(NI.GetId(), InStrToNIdH.GetKey(NI.GetId()));
		dest = NI.GetId();
	}

	DrawGViz(G, TGVizLayout::gvlDot, "mygraph.png", "Loaded Graph", true);
	DrawGViz(G, TGVizLayout::gvlDot, "mygraph_labeled.png", "Loaded Graph",
			NodeLabelH);

	Dest_Link_Table table;
	// TODO iterate over all destinations
	BackupPaths bp; // for a fixed destination
	table[dest] = bp;

	for (TNEANet::TEdgeI EI = G->BegEI(); EI < G->EndEI(); EI++) {
		int S = EI.GetSrcNId();
		int F = EI.GetDstNId();
		int D = dest; // iterate on all possible destinations
		printf("link (%d, %d), D=%d\n", S, F, D);

		Path sp;
		int len = getSP(G, S, D, sp);
		if (len >= INF) {
			mylog << "D is not reachable from S\n";
			continue;
		}
		printf("first shortest path weight=%d: ", len);
		printsp(sp);
		if (sp[1] != F) {
			mylog << "SP not using the link\n";
			continue;
		}

		const int w = G->GetIntAttrDatE(EI.GetId(), WEIGHTATTR);
		G->AddIntAttrDatE(EI.GetId(), INF, WEIGHTATTR); // cripple the link (S,F)
		len = getSP(G, S, D, sp); // shortest path without (S,F) => post convergence path
		G->AddIntAttrDatE(EI.GetId(), w, WEIGHTATTR);
		printf("second shortest path weight=%d: ", len);
		printsp(sp);
		if (len >= INF) {
			mylog << "no backup path\n";
			continue;
		}

		int Pspace_sup = S, Qspace_inf = D;
		for (int i = 0; i < sp.size() && inPSpace(G, sp[i], S, F); ++i) {
			printf("%d in PSpace\n", sp[i]);
			Pspace_sup = i;
		}
		for (int i = sp.size() - 1;
				i >= Pspace_sup && inQSpace(G, sp[i], S, F, D); --i) { // search backward for the infimum of QSpace, stop when it touches the PSpace
			printf("%d in QSpace\n", sp[i]);
			Qspace_inf = i;
		}

		Assert(Qspace_inf == Pspace_sup + 1 || Pspace_sup == Qspace_inf);
		MyNode p = sp[Pspace_sup], q = sp[Qspace_inf];
		printf("p=%d q=%d\n", p, q);
		// extract the segments from p and q

		Segment seg = { true, p };
		SegmentStack stack = { seg };
		table[dest][EI.GetId()] = stack;
		if (Qspace_inf > Pspace_sup) { // also push a "tunnel link"
			table[dest][EI.GetId()].push_back( { false, G->GetEId(p, q) });
		}
	}

	printf("\nmap size for dest node=%u\n", table[dest].size());
	printf("table size=%u\n", table.size());

//	int randNode = sp[rand() % (sp.size() - 1)];
//	printf("randNode=%d", randNode);
//	int failedLink = G->GetEI(sp[randNode], sp[randNode + 1]).GetId();
//	printf("failedLink=%d\n", failedLink);
//	G->AddIntAttrDatE(failedLink, INF,WEIGHTATTR);
//	getSP(G, src, dest, sp);
//	printsp(sp);

	return 0;
}
