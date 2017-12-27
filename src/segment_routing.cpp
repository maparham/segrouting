#include <execinfo.h>
#include <signal.h>
#include <iostream>
#include <fstream>
#include <algorithm>    // std::reverse
#include <vector>
#include <map>
#include <time.h>
#include <TILFA.hpp>

using namespace std;
using namespace TSnap;

#define INF numeric_limits<int>::max()/3
#define WEIGHTATTR "weight"

#define __DEBUG__ 1
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
	TILFA tilfa;
	while (filelist >> path) {
		PRINTF(" Loading %s\n", path.c_str());
		tilfa.load_graph(path);

		Dest_Link_Table table;
		tilfa.compute_SR_table(table);
		PRINTF("table size=%u\n", table.size());

		tilfa.eval_double_failure(table);
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
