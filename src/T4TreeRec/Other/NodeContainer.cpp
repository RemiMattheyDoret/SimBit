typedef NodeContainerSearchHelper int

NodeContainer
{
private:
	int firstGeneration; // int in case of negative generations
	int currentGeneration;
	size_t nbNodes;
	std::vector<size_t> data;
	// data[generation - firstGeneration] contains the first node ID that starts at generation 'generation'

public:
	template<typename generation_type>
	
	NodeContainer(const generation_type& gen, const size_t nbNodesGenerationsToReserveFor )
	:firstGeneration(gen), currentGeneration(gen), nbNodes(0)
	{
		data.reserve(nbNodesGenerationsToReserveFor);
	}

	template<typename generation_type>
	void addNode(const generation_type& generation)
	{
		if (currentGeneration != generation)
		{
			data.push_back(nbNodes);
			currentGeneration = generation;
		}
		++nbNodes;
	}

	template<typename generation_type>
	size_t getFirstIDOfGeneration(const generation_type& generation) const
	{
		return data[generation - firstGeneration];
	}

	template<typename INT>
	void reverseIDs(INT& totalNbNodes)
	{
		for (auto& elem : data)
			elem = totalNbNodes - elem;
	}

	template<typename nodeID_type>
	size_t findGenerationOfNodeID(nodeID_type& nodeID)
	{
		return firstGeneration + (std::lower_bound(data.begin(), data.end(), nodeID) - data.begin());
	}

	template<typename nodeID_type>
	size_t findGenerationOfNodeID(nodeID_type& nodeID, NodeContainerSearchHelper& searchHelper, bool goingForwardInTime)
	{
		// Binary search is probably not ideal here. Better use a searchHelper to look through table generation per generation
		// return firstGeneration + (std::lower_bound(data.begin() + knownAtLeastGeneration, data.end(), nodeID) - data.begin());
		if (goingForwardInTime)
		{
			for (int generation_index = searchHelper ; generation_index < data.size() ; ++generation_index)
			{
				if (data[generation_index] < nodeID)
				{
					searchHelper = generation_index;
					return firstGeneration + generation_index;
				}
			}
			assert(currentGeneration == firstGeneration + data.size() );
			searchHelper = currentGeneration;
			return currentGeneration;
		} else
		{
			for (int generation_index = searchHelper ; generation_index >= 0 ; --generation_index)
			{
				if (data[generation_index] < nodeID)
				{
					searchHelper = generation_index;
					return firstGeneration + generation_index;
				}
			}
			assert(currentGeneration == firstGeneration + data.size() );
			searchHelper = currentGeneration;
			return currentGeneration;
		}
	}
};