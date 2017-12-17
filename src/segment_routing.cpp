#include <iostream>
#include "Snap.h"
#include <fstream>
#include <unordered_map>

using namespace std;
using namespace TSnap;

#define DEBUG 1

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

int main() {

//	ifstream topgen;
//	topgen.open("/Users/mahmoud/Downloads/weights-dist/1221/latencies.intra");

//	if (!topgen.is_open()) {
//		mylog("file not opened.");
//		return nodes;
//	}

//// what type of graph do you want to use?
	typedef PNEANet MyGraph; // undirected graph
	//typedef PNGraph PGraph;  //   directed graph
	//typedef PNEGraph PGraph;  //   directed multigraph
	//typedef TPt<TNodeNet<TInt> > PGraph;
	//typedef TPt<TNodeEdgeNet<TInt, TInt> > PGraph;

//	G= TSnap::GenRndGnm<MyGraph>(10,20);

	const TStr weightAttr = "weight";
	TStrHash<TInt> InStrToNIdH;
	MyGraph G = LoadWeightedEdgeListStr<MyGraph>(
			"/Users/mahmoud/Downloads/weights-dist/1221/latencies.intra", 0, 1,
			2, weightAttr, InStrToNIdH);

	TIntStrH NodeLabelH;
	for (auto NI = G->BegNI(); NI < G->EndNI(); NI++) {
		NodeLabelH.AddDat(NI.GetId(), InStrToNIdH.GetKey(NI.GetId()));
	}

	DrawGViz(G, TGVizLayout::gvlDot, "mygraph.png", "Loaded Graph", true);
	DrawGViz(G, TGVizLayout::gvlDot, "mygraph_labeled.png", "Loaded Graph",
			NodeLabelH);

	TIntH SP_tree;
	//GetWeightedShortestPath(G, 0, SP_tree, weightAttr);
	return 0;
}
