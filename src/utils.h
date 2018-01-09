/*
 * utils.h
 *
 *  Created on: 3 Jan 2018
 *      Author: mahmoud
 */

//#define LOG2FILE 1
#ifdef __DEBUG__
#ifdef Log2File
#define PRINTF(...) fprintf(f,__VA_ARGS__)
#else
#define PRINTF printf
#endif
#else
#define PRINTF(format, args...) ((void)0)
#endif

#ifndef UTILS_H_
#define UTILS_H_

#define RESULT_DIR string("results")

#include <string>
#include <sstream>
#include <vector>
#include <iterator>
#include <utils.h>

using namespace std;

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

FILE *f;
void init() {
	signal(SIGSEGV, handler); // install our handler
	/* initialize random seed: */
	srand(time(NULL));
	f = fopen("log.txt", "w");
	if (f == NULL) {
		printf("Error opening file!\n");
		exit(1);
	}
}

template<typename Out>
void split(const std::string &s, char delim, Out result) {
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim)) {
		*(result++) = item;
	}
}

std::vector<std::string> split(const std::string &s, char delim) {
	std::vector<std::string> elems;
	split(s, delim, std::back_inserter(elems));
	return elems;
}

struct Reporter {
	string name;
	ofstream of;
	ofstream all;
	size_t escapes = 0;

	static string resultDir() {
		string dir = RESULT_DIR + '/';

#ifdef FLUSH_STACK
		dir+="withflush_";
#else
		dir += "noflush_";
#endif

#ifdef doubleTILFA
		dir+="doubleTILFA";
#else
		dir += "TILFA";
#endif
		dir += '/';
		return dir;
	}

	// constructor
	Reporter(const string& path) {
		// prepare result file
		vector<string> tokens = split(path, '/');
		reverse(tokens.begin(), tokens.end());
		name = tokens[1] + '_' + tokens[0];
		string dir = resultDir();
		string resultPath = dir + name + ".txt";

		printf("Result Path: %s\n", resultPath.c_str());

		of.open(resultPath);
		all.open(dir + "all.txt", std::ios_base::app);
		if (!of.is_open() || !all.is_open()) {
			printf("file not opened: %s", resultPath.c_str());
			exit(1);
		}
	}
	template<typename T>
	Reporter& operator<<(T input) {
		of << input;
		all << input;
		return *this;
	}
	~Reporter() {
		all.close();
		of.close();
	}
};

#endif /* UTILS_H_ */
