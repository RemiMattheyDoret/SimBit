#include <iostream>
#include <fstream>
#include <vector>
#include <limits.h>
#include <assert.h>
#include <sstream>
#include <algorithm>

#include "CompressedSortedDeque.h"
#include "CompressedSortedDeque.cpp"


int main()
{
	std::cout << "Start testing...\n";

	std::vector<uint32_t> stdv = {65, 90, 9380, 9384, 27685, 35000, 35040, 85040, 125040, 185040, 355360, 1553000, 1853000, 1853002};

	std::cout << "print stdv\n";
	for (auto& elem : stdv)
	{
		std::cout << elem << " "  << std::flush;	
	}
	std::cout << "\n"  << std::flush;

	std::cout << "Create CompressedSortedDeque\n"  << std::flush;
	CompressedSortedDeque v(6553000);

	std::cout << "push_backs\n";
	for (auto& elem : stdv)
	{
		v.push_back(elem);	
	}


	std::cout << "INSERTIONS\n";
	auto it0 = v.insert(303);
	std::cout << "I just inserted the value " << *it0 << "\n";
	std::cout << "The value that follows it is " << *(it0+1) << "\n";
	std::cout << "The value that precedes it is " << *(it0-1) << "\n";
	v.insert(it0,302);
	v.insert(301);
	v.insert(300);
	
	std::cout << "1 Let's iterate through this vector\n";
	for (auto it = v.begin() ; it < v.end() ; ++it)
	{
		std::cout << *it << " ";
	}
	std::cout << "\n";
	std::cout << "iteration done\n";
	
		

	std::cout << "2 Let's iterate through this vector again\n";
	for (auto it = v.begin() ; it < v.end() ; ++it)
	{
		std::cout << *it << " ";
	}
	std::cout << "\n";
	std::cout << "iteration done\n";

	
	std::cout << "search for 27685\n";
	auto it = v.lower_bound(27685);
	std::cout << "Value found is "<< *it <<". Erase it now\n";
	v.erase(it);
	std::cout << "I just erased this element\n";

	std::cout << "Search element 9380 (who exists). I get " << *v.lower_bound(9380) << "\n";
	std::cout << "Search element 9381 (who does not exist). I get " << *v.lower_bound(9381) << "\n";
	std::cout << "Search element 27685 (who does not exist anymore). I get " << *v.lower_bound(27685) << "\n";

	std::cout << "4 Let's iterate through this vector\n";
	for (auto it = v.begin() ; it < v.end() ; ++it)
	{
		std::cout << *it << " ";
	}
	std::cout << "\n";

	std::cout << "Let's create v2\n";
	CompressedSortedDeque v2(6553000);
	std::cout << "v2 created. Now let's copy over from 302 to 126000\n";
	v2.extend(302, 126000, v);

	std::cout << "Let's iterate through v2\n";
	for (auto it = v2.begin() ; it < v2.end() ; ++it)
	{
		std::cout << *it << " ";
	}
	std::cout << "\n";

	return 0;
}
