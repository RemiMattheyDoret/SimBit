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

#include "RNG_wrapper.h"
#include "RNG_wrapper.cpp"
#include "Walker.h"
#include "Walker.cpp"


int main ()
{
	RNG_wrapper rngw(12);

	std::vector<double> cumSumFits = {0.8, 1.6, 1.7, 1.8, 5};
	Walker walker(cumSumFits);

	std::vector<size_t> counts(5,0);
	for (int i = 0 ; i < 1e5; ++i)
	{
		counts[walker(rngw)]++;
	}
	
	for (auto& elem : counts)
		std::cout << elem << ":";
	std::cout << "\n";

}
