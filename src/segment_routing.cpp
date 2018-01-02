#include <execinfo.h>
#include <signal.h>
#include <iostream>
#include <fstream>
#include <algorithm>    // std::reverse
#include <vector>
#include <map>
#include <time.h>
#include <TILFA.hpp>
#include <Klaus.hpp>

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
		if (path[0] == '#') {
			continue;
		}
		PRINTF(" Loading %s\n", path.c_str());
		tilfa.load_graph(path);
		//tilfa.drawGraph(path);
		tilfa.printInfo();

		tilfa.compute_SR_table();

		pair<int, int> res = tilfa.eval_double_failure();
		PRINTF("fail=%d, success=%d, ratio=%f\n", res.first, res.second,
				(double) res.first / (res.first + res.second));

	}

//	int randNode = sp[rand() % (sp.size() - 1)];
//	PRINTF("randNode=%d", randNode);
//	int failedLink = G->GetEI(sp[randNode], sp[randNode + 1]).GetId();
//	PRINTF("failedLink=%d\n", failedLink);
//	G->AddIntAttrDatE(failedLink, INF,WEIGHTATTR);
//	getSP(G, src, dest, sp);
//	printsp(sp);

	return 0;
}
