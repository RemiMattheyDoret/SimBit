
class NodeGeneticsContainer
{
private:
	size_t smallestID;
	std::deque<NodeGenetics*> data; // vectors of pointers so to not allocate memory about size of container and how much memory is reserved
	// first element of data is the last ID and it goes decreasing in IDs

	template<typename nodeID_type>
	nodeID_type getIndex(nodeID_type& i)
	{
		return i - smallestID;
	}

public:
	NodeGeneticsContainer():smallestID(std::numeric_limits<size_t>::max()){}

	const std::deque<NodeGenetics*>& getData() const
	{
		return data;
	}

	template<typename nodeID_type>
	NodeGenetics* getHaploP(nodeID_type& nodeID) const
	{
		return data[getIndex(nodeID)];
	}

	template<typename nodeID_type>
	bool doesAlreadyExist(nodeID_type& nodeID) const
	{
		auto index = getIndex(nodeID)
		if (index < data.size())
		{
			return data[getIndex(nodeID)] != nullptr;
		}
		return false;
	}

	template<typename nodeID_type>
	NodeGenetics& getHaplo(nodeID_type& nodeID) const
	{
		return *data[getIndex(nodeID)];
	}

	template<typename nodeID_type>
	void deleteHaplo(nodeID_type& nodeID)
	{
		auto index = getIndex(nodeID);
		assert(data[index] != nullptr);
		delete data[index];
		data[index] = nullptr;
		if (smallestID == nodeID)
		{
			while (data.front() == nullptr)
			{
				data.pop_front();
				++smallestID;
			}
		} else
		{
			assert(data.front() != nullptr);
		}
	}

	void deleteAllHaplo()
	{
		for (size_t index = 0 ; index < data.size() ; ++index)
		{
			if (data[index] != nullptr)
			{
				delete data[index];
			}
		}
		std::deque<NodeGenetics*>().swap(data);
	}


	template<typename nodeID_type>
	void insertHaploP(nodeID_type& nodeID, NodeGenetics* H)
	{
		if (nodeID < smallestID) 
		{
			// It needs more memory at front?
			// Ideally that should almost never happen

			auto toAddAtFront = smallestID - nodeID - 1 ;
			for ( size_t i = 0 ; i < toAddAtFront  ; ++i )
			{
				data.push_front(nullptr);
			}
			data.push_front(H);
			smallestID = nodeID;
		} else
		{
			auto index = getIndex(nodeID);
			if (index > data.size())
			{
				// It needs more memory at back?
				auto toAddAtFront = index - data.size() - 1;
				for ( size_t i = 0 ; i < toAddAtFront  ; ++i )
				{
					data.push_back(nullptr);
				}
				assert(index == data.size());
				data.push_back(H);
			}

			// insert
			data[index] = H;
		}
	}

	void deleteAllHaplos()
	{
		for (auto p : data) delete p;
		data.clear();
	}

	void shrink_to_fit()
	{
		data.shrink_to_fit();
	}
};
