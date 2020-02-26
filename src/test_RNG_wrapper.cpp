#include <iostream>
#include <list>
#include <fstream>
#include <random>
#include <vector>
#include <string>
#include <limits.h>
#include <float.h>
#include <assert.h>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <map>
#include <queue>
#include <ctime>
#include <boost/algorithm/string.hpp>
#include <boost/tokenizer.hpp>
#include <signal.h>


// forward class declaration and global variables
#include "ForwardClassDeclarations.h"
#include "GlobalVariables.h"

// RNG
#include "RNG_wrapper.h"
#include "RNG_wrapper.cpp"



int main()
{
	RNG_wrapper rngw(15);
	size_t nbZeros = 0;
	size_t nbOnes = 0;
	for (int i = 0 ; i<500000; ++i)
		rngw.get_1b() ? nbOnes++: nbZeros++;
	std::cout << nbZeros << ":" << nbOnes << "\n";

	{
		size_t N = 8;
		std::vector<size_t> counts(N);
		for (int i = 0 ; i<500000000; ++i)
			counts[rngw.uniform_int_distribution(N-1)]++;

		for (int i = 0 ; i < counts.size(); ++i)
			std::cout << counts[i] << ":";
		std::cout << "\n";
	}
	
	{
		size_t N = 8;
		std::uniform_int_distribution<size_t> d(0,N-1);
		std::vector<size_t> counts(N);
		for (int i = 0 ; i<500000000; ++i)
			counts[d(rngw.getRNG())]++;

		for (int i = 0 ; i < counts.size(); ++i)
			std::cout << counts[i] << ":";
		std::cout << "\n";
	}	
}

