#include <execinfo.h>
#include <signal.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <time.h>

//#define doubleTILFA
//#define FLUSH_STACK
#include <TILFA.hpp>

#define __DEBUG__ 1
#include <utils.h>

using namespace std;
using namespace TSnap;

#define INF numeric_limits<int>::max()/3
#define WEIGHTATTR "weight"

int main() {
	init();
	clock_t begin = clock();
#ifdef doubleTILFA
	printf("double-fail TILFA\n");
#else
	printf("single-fail TILFA\n");
#endif
#ifdef FLUSH_STACK
	printf("with flush\n");
#else
	printf("no flush\n");
#endif
//	G= TSnap::GenRndGnm<MyGraph>(10,20);

	std::ifstream filelist("files.txt");
	string path;
	while (filelist >> path) {
		if (path[0] == '#') {
			continue;
		}
		printf("\nLoading %s\n", path.c_str());
		TILFA tilfa;
		tilfa.load_graph(path);
		tilfa.printInfo();

		Reporter reporter(path);

		//tilfa.drawGraph(path);

		Result res = tilfa.eval_double_failure(reporter);

		PRINTF("fail=%d, success=%d, ratio=%f; maxSS=%d\n", res.fail, res.success,
				(double) res.fail / (res.fail + res.success), res.maxStackSize);

		clock_t now = clock();
		clock_t elapsed_secs = double(now - begin) / CLOCKS_PER_SEC;
		if (elapsed_secs < 60) {
			PRINTF("elapsed time=%d sec\n", elapsed_secs);
		} else {
			PRINTF("elapsed time=%d min\n", elapsed_secs / 60);
		}
	}

	clock_t end = clock();
	clock_t elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	PRINTF("total time=%d min\n", elapsed_secs / 60);
	return 0;
}
