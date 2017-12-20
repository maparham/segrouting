#include <iostream>
#include "Snap.h"
#include <fstream>
#include <unordered_map>
#include <vector>


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

int GetWeightedShortestPath11(const PNEANet Graph, const int& SrcNId,
		TIntFltH& NIdDistH, const TFltV& Attr) {
	TIntV frontier;

	NIdDistH.Clr(false);
	NIdDistH.AddDat(SrcNId, 0);
	frontier.Add(SrcNId);
	while (!frontier.Empty()) {
		const int NId = findMinimum(frontier, NIdDistH);
		printf("\nNId=%d", NId);
		const PNEANet::TObj::TNodeI NodeI = Graph->GetNI(NId);
		for (int v = 0; v < NodeI.GetOutDeg(); v++) {
			int DstNId = NodeI.GetOutNId(v);
			int EId = NodeI.GetOutEId(v);

			float elen = Attr[EId];
			if (!NIdDistH.IsKey(DstNId)) {
				NIdDistH.AddDat(DstNId, NIdDistH.GetDat(NId) + Attr[EId]);
				frontier.Add(DstNId);
				printf("\nNIdDistH.AddDat(%d,  Attr[%d]=%f);\n", DstNId, EId,
						NIdDistH.GetDat(NId) + Attr[EId]);

			} else {
				printf("\nNIdDistH.GetDat(%d)=%f", DstNId,
						NIdDistH.GetDat(DstNId));
				printf("\nNIdDistH[%d]=%f \n", DstNId,
						NIdDistH[DstNId]);
				printf(
						"for node %d, testing the parent %d, new distance=%f, the current dist=%f\n",
						DstNId, NId, NIdDistH.GetDat(NId) + Attr[EId],
						NIdDistH.GetDat(DstNId));

				if (NIdDistH[DstNId] > NIdDistH.GetDat(NId) + Attr[EId]) {
					NIdDistH[DstNId] = NIdDistH.GetDat(NId) + Attr[EId];
					printf("DstNId=%d NIdDistH[DstNId]=%f", DstNId,
							NIdDistH[DstNId]);
				}
			}
		}
	}
	return 0;
}

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

	const TStr weightAttr = "weight", inFile =
			"/Users/mahmoud/Downloads/weights-dist/1221/latencies.intra",
			inFile1 = "in.txt";
	// load graph from file
	TStrHash<TInt> InStrToNIdH;
	MyGraph G = LoadWeightedEdgeListStr<MyGraph>(inFile, 0, 1, 2, weightAttr,
			InStrToNIdH);

	TIntStrH NodeLabelH;
	int dest, src = 0;
	for (auto NI = G->BegNI(); NI < G->EndNI(); NI++) {
		NodeLabelH.AddDat(NI.GetId(), InStrToNIdH.GetKey(NI.GetId()));
		dest = NI.GetId();
	}

	DrawGViz(G, TGVizLayout::gvlDot, "mygraph.png", "Loaded Graph", true);
	DrawGViz(G, TGVizLayout::gvlDot, "mygraph_labeled.png", "Loaded Graph",
			NodeLabelH);

//	int EId = G->GetEI(src, 1).GetId();
//	 G->AddIntAttrDatE(EId, 9999999, weightAttr);

	TIntH SP_tree;
	GetWeightedShortestPathTree(G, src, SP_tree, weightAttr);

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
	vector<int> sp;
	while (SP_tree.IsKey(N) && N != src) {
		sp.push_back(N); // save the shortest path
		mylog << N << ",";
		int parent = SP_tree.GetDat(N);

		int EID = G->GetEI(parent, N).GetId();
		spLen += G->GetIntAttrDatE(EID, weightAttr);
		N = parent;
	}
	mylog << "\nspLen=" << spLen;
	return 0;
}
