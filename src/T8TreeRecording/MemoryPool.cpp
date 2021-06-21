template<typename T>
MemoryPool<T>::MemoryPool()
{
	assert(multiplier > 1.0);
}

template<typename T>
MemoryPool<T>::~MemoryPool()
{
	//std::cout << "When deleting the memory pool, there is " << available.size() << " / " << objects.size() << " elements available. Throughout the simulation, " << nbNEWs << " elements have been required and " << nbDELETESs << " elements have been returned\n";
}

template<typename T>
void MemoryPool<T>::setMultiplier(float m)
{
	this->multiplier = m;
	assert(this->multiplier > 1.0);
}

template<typename T>
void MemoryPool<T>::increaseSize(size_t newSize)
{
	std::cout << "Increasing size from " << objects.size() << " to " << newSize << " (multiplier = "<<multiplier <<")\n";
	assert(objects.size() < newSize);

	auto from = objects.size();
	for (size_t i = from ; i < newSize ; ++i)
	{
		objects.push_back(T()); // no need to reserve, they are deques
		small_available.push_front(&objects.back()); // no need to reserve, they are deques
	}
}

template<typename T>
void MemoryPool<T>::reallocate()
{
	assert(objects.size()); // The user should do the first resize
	size_t oldSize = objects.size();
	size_t addedSize = (size_t) ceil(oldSize * (multiplier - 1.0));
	size_t newSize = oldSize + addedSize;
	assert(newSize > oldSize);

	increaseSize(newSize);
}


template<typename T>
void MemoryPool<T>::resize(size_t newSize)
{
	int addedSize = newSize - objects.size();

	if (addedSize == 0) return;
	if (addedSize > 0)
	{
		increaseSize(newSize);
	} else
	{
		std::cout << "In function 'template<typename T> void MemoryPool<T>::resize(size_t newSize)', the newsize ("<<newSize<<") is smaller than the old size ("<<objects.size()<<"). Sorry the memoryPool does not allow to reduce size (because it does not know what objects to remove).\n";
	}
}

template<typename T>
size_t MemoryPool<T>::size() const
{
	return objects.size();
}


template<typename T>
T* MemoryPool<T>::NEW(bool isSmall)
{
	++nbNEWs;

	/*
	if (GP->CurrentGeneration > 9990)
	{
		if (isSmall) std::cout << "asking small\n"; else std::cout << "asking big\n";
		std::cout << "big_available.size() = " << big_available.size() << "\n";
		std::cout << "small_available.size() = " << small_available.size() << "\n--------";
	}
	*/

	if (!isSmall)
	{
		if (!big_available.empty())
		{
			auto r = big_available.back();
			big_available.pop_back();
			return r;
		}
	}

	
	// If needs small or if there is no big available
	
	if (small_available.empty())
	{
		if (big_available.empty())
		{
			reallocate();
			assert(!small_available.empty());
			auto r = small_available.back();
			small_available.pop_back();
			return r;
		} else
		{
			auto r = big_available.back();
			big_available.pop_back();
			return r;
		}		
	} else
	{
		auto r = small_available.back();
		small_available.pop_back();
		return r;
	}
}


template<typename T>
void MemoryPool<T>::DELETE(T* ptr, bool isSmall)
{
	//std::cout << "deleting " << ptr << "\n";

	/*
	if (GP->CurrentGeneration > 9990)
	{
		if (isSmall) std::cout << "deleting small\n"; else std::cout << "deleting big\n";
		std::cout << "big_available.size() = " << big_available.size() << "\n";
		std::cout << "small_available.size() = " << small_available.size() << "\n--------";
	}
	*/

	ptr->clear();

	if (isSmall)
	{
		small_available.push_back(ptr);
	} else
	{
		if (ptr->geneticData.capacity())
			big_available.push_back(ptr);
			//std::cout << "adding to big\n";
		else
			small_available.push_back(ptr);
	}
		
	++nbDELETESs;
}


template<typename T>
void MemoryPool<T>::makeBig(T* o)
{
	if (!big_available.empty())
	{
		auto big = big_available.back();
		if (big->geneticData.capacity() > o->geneticData.capacity())
		{
			big_available.pop_back();
			assert(big->geneticData.size() == 0);

			big->geneticData.insert(big->geneticData.begin(), o->geneticData.begin(), o->geneticData.end());
			big->geneticData.swap(o->geneticData);
			big->geneticData.clear();

			small_available.push_back(big);
		}
	}
}



template<typename T>
void MemoryPool<T>::shrink_all_to_fit()
{
	for (T& obj : objects)
		if (obj.geneticData.size() != obj.geneticData.capacity()) obj.shrink_to_fit();
}

template<typename T>
void MemoryPool<T>::computeSomeStats()
{
	std::cout << "objects.size() = "<<objects.size()<<"\n";
	std::cout << "small_available.size() = "<<small_available.size()<<"\n";
	std::cout << "big_available.size() = "<<big_available.size()<<"\n";

	size_t totalCapacity = 0;
	size_t totalSize = 0;
	size_t capacityUnusedInSmall = 0;
	size_t capacityUnusedInBig = 0;

	for (T* ptr : small_available)
	{
		capacityUnusedInSmall += ptr->geneticData.capacity();
	}

	for (T* ptr : big_available)
	{
		capacityUnusedInBig += ptr->geneticData.capacity();
	}

	for (auto& obj : objects)
	{
		totalCapacity += obj.geneticData.capacity();
		totalSize += obj.geneticData.size();
	}

	//std::cout << "capacityUnusedInSmall = "<<capacityUnusedInSmall<<"\n";
	//std::cout << "capacityUnusedInBig = "<<capacityUnusedInBig<<"\n";
	//std::cout << "sizeUsed = "<<totalCapacity - capacityUnusedInSmall - capacityUnusedInBig<<"\n";
	
	std::cout << "-> fraction of total capacity that is used = "<< (double)totalSize / (double)totalCapacity * (double)100.0<<"%\n";

	std::cout << "-> fraction of total capacity that is available in small = "<< (double)capacityUnusedInSmall / (double)totalCapacity * (double)100.0<<"%\n";

	std::cout << "-> fraction of total capacity that is available in big = "<< (double)capacityUnusedInBig / (double)totalCapacity * (double)100.0<<"%\n";

	std::cout << "-> fraction of total capacity that is distributed and used = "<< (totalSize) / (double)(totalCapacity - capacityUnusedInBig - capacityUnusedInSmall) * (double)100.0<<"%\n";
	
}
	

template<typename T>
void MemoryPool<T>::moveBigToSmall()
{

	while (!big_available.empty())
	{
		small_available.push_back(big_available.back());
		big_available.pop_back();
	}
}

template<typename T>
void MemoryPool<T>::emptySmallMemory()
{
	for (T* ptr : small_available)
		std::vector<uint32_t>().swap(ptr->geneticData);
}


