

/*

	  []       []          [][][][][][][][][][][]       []
	biggest  biggest-1                                smallest

	biggest - size - 1 = smallest

*/



template<typename node_t>
template<typename nodeID_type>
nodeID_type T4TreeRec::HaplotypesContainer<node_t>::getIndex(nodeID_type& i) const
{
	return biggestID - i;
}

template<typename node_t>
T4TreeRec::HaplotypesContainer<node_t>::HaplotypesContainer()
:biggestID(0)
{}

template<typename node_t>
const std::deque<node_t*>& T4TreeRec::HaplotypesContainer<node_t>::getData() const
{
	return data;
}

template<typename node_t>
template<typename nodeID_type>
node_t* T4TreeRec::HaplotypesContainer<node_t>::getHaploP(nodeID_type& nodeID) const
{
	//std::cout << "getting oldNodeID = " << nodeID << "\n";
	auto index = getIndex(nodeID);
	assert(index >= 0 && index < data.size());
	return data[index];
}

template<typename node_t>
template<typename nodeID_type>
bool T4TreeRec::HaplotypesContainer<node_t>::doesAlreadyExist(nodeID_type& nodeID) const
{
	auto index = getIndex(nodeID);
	if (index < data.size() && index >= 0)
	{
		return data[index] != nullptr;
	}
	return false;
}

template<typename node_t>
template<typename nodeID_type>
node_t& T4TreeRec::HaplotypesContainer<node_t>::getHaplo(nodeID_type& nodeID) const
{
	return *getHaploP(nodeID);
}

template<typename node_t>
template<typename nodeID_type>
void T4TreeRec::HaplotypesContainer<node_t>::deleteHaplo(nodeID_type& nodeID)
{
	auto index = getIndex(nodeID);
	assert(data[index] != nullptr);
	delete data[index];
	data[index] = nullptr;
	if (biggestID == nodeID)
	{
		while (data.front() == nullptr)
		{
			data.pop_front();
			--biggestID;
		}
	} else
	{
		assert(data.front() != nullptr);
	}
}

template<typename node_t>
void T4TreeRec::HaplotypesContainer<node_t>::deleteAllHaplos()
{
	for (uint32_t index = 0 ; index < data.size() ; ++index)
	{
		if (data[index] != nullptr)
		{
			delete data[index];
		}
	}
	std::deque<node_t*>().swap(data);
}


template<typename node_t>
template<typename nodeID_type>
void T4TreeRec::HaplotypesContainer<node_t>::insertHaploP(nodeID_type nodeID, node_t* H)
{
	/*
	std::cout << "biggestID = " << biggestID << "\n";
	std::cout << "nodeID = " << nodeID << "\n";
	std::cout << "data.size() = " << data.size() << "\n";
	std::cout << "getIndex(nodeID) = " << getIndex(nodeID) << "\n";
	std::cout << "---------------------\n";
	*/

	//std::cout << "inserting oldNodeID = " << nodeID << "\n";

	bool hasBeebInsertedYet = false;

	if (nodeID > biggestID)
	{
		for (; biggestID != nodeID ;++biggestID)
		{
			data.push_front(nullptr);
		}
		assert(getIndex(nodeID) == 0);
		data.front() = H;
		hasBeebInsertedYet = true;
	}

	auto index = getIndex(nodeID);
	assert(index >= 0);
	if (index >= data.size())
	{
		while ( !(index < data.size()) )
		{
			data.push_back(nullptr);
		}
		data.back() = H;
		hasBeebInsertedYet = true;
	}

	assert(index < data.size());
	if ( !hasBeebInsertedYet )
	{
		data[index] = H;
	}
}

template<typename node_t>
void T4TreeRec::HaplotypesContainer<node_t>::shrink_to_fit()
{
	data.shrink_to_fit();
}

template<typename node_t>
void T4TreeRec::HaplotypesContainer<node_t>::iterator_restart()
{
	iterator_index = 0;
	while (data[iterator_index] == nullptr)
	{
		++iterator_index;
	}	
}

template<typename node_t>
typename T4TreeRec::HaplotypesContainer<node_t>::ObjForIteration T4TreeRec::HaplotypesContainer<node_t>::iterator_next()
{
	auto oldID = biggestID - iterator_index;
	assert(oldID >= 0);
	T4TreeRec::HaplotypesContainer<node_t>::ObjForIteration r(data[iterator_index], oldID, data[iterator_index]->getNewNodeID());

	++iterator_index;
	while (iterator_index < data.size() && data[iterator_index] == nullptr)
	{
		++iterator_index;
	}

	return r;
}


template<typename node_t>
bool T4TreeRec::HaplotypesContainer<node_t>::iterator_isMore()
{
	return iterator_index != data.size();
}

template<typename node_t>
uint32_t T4TreeRec::HaplotypesContainer<node_t>::size1() const
{
	uint32_t size = 0;
	for (uint32_t i = 0 ; i < data.size() ; ++i)
	{
		if (data[i] != nullptr) ++size;
	}
	return size;
}

template<typename node_t>
uint32_t T4TreeRec::HaplotypesContainer<node_t>::size2() const
{
	return data.size();
}



