#include <TILFA.hpp>
#include <time.h>

/*struct RandomSelector {
 int n;
 const int choose;
 const int nofLinks;

 RandomSelector(const int nofLinks, const int choose = 2) :
 nofLinks(nofLinks), choose(choose) {
 srand(time(NULL));
 }
 bool next(vector<int>& id) {
 id.resize(choose);
 for (int i = 0; i < choose; ++i) {
 id[i] = rand() % nofLinks;
 for (int j = 0; j < i; ++j) {
 if (id[j] == id[i]) {
 --i; // repeat the same iteration for i until distinct
 break;
 }
 }
 }
 return true;
 }
 };*/

class Klaus: public TILFA {

	Dest_Link_Table_3D table;

	/*	bool inTable(int k1, int k2, int k3) {
	 return table.find(k1) != table.end()
	 && table[k1].find(k2) != table[k1].end()
	 && table[k1][k2].find(k3) != table[k1][k2].end();
	 }
	 SegmentStack& newTableEntry(const int k0, const int k1, const int k2) {
	 table.insert(make_pair(k0, SSMap_2D()));
	 table[k0].insert(make_pair(k1, SSMap()));
	 table[k0][k1].insert(make_pair(k2, SegmentStack()));
	 }*/

	pair<int, int> eval_double_failure() {
		int fail = 0, success = 0;
		for (TNEANet::TNodeI NI = G->BegNI(); NI < G->EndNI(); NI++) {
			int D = NI.GetId();

			vector<int> EIds;
			int N = 100;
			while (N-- > 0) { // TODO
				do {
					EIds[0] = G->GetRndEId();
					EIds[1] = G->GetRndEId();
				} while (table[D][EIds[0]][EIds[1]].size() > 0);
			}
			SegmentStack stack = table[D][EIds[0]][EIds[1]];
			for (TNEANet::TEdgeI EI1 = G->BegEI(); EI1 < G->EndEI(); EI1++) {
				int S = EI1.GetSrcNId();
				int F = EI1.GetDstNId();
				Link L1 = Link(S, F, G);
				Path bp1;
				vector<int> destinationMap1;
				PRINTF("L1=(%d,%d), D=%d: ", S, F, D);
				double len = backup_path(L1, stack, bp1, destinationMap1);
				print_path(bp1);
				Assert(L1.inPath(bp1) == false);
				// all second failures on the first backup path
				for (int i = 1; i < bp1.size(); ++i) {
					Link L2 = Link(bp1[i - 1], bp1[i], G);
					Path bp2;
					vector<int> destinationMap2;
					PRINTF("L2=(%d,%d), D=%d: ", L2.tail, L2.head,
							destinationMap1[S]);
					//					len += backup_path(L2, destinationMap1[S], table, bp2,
					//							destinationMap2); // case1: keep the current stack
					len += backup_path(L2, stack, bp2, destinationMap2);// case2: flush the stack down to D
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
		return pair<int, int>(fail, success);
	}
}
;

