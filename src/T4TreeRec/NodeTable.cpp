
void T4TreeRec::NodeTable::swap(NodeTable& other)
{
	this->data.swap(other.data);
}

void T4TreeRec::NodeTable::clear()
{
	data.clear();
}

void T4TreeRec::NodeTable::reserve(uint32_t n)
{
	data.reserve(n);
}

void T4TreeRec::NodeTable::reverse()
{
	std::reverse(data.begin(), data.end());
}

uint32_t T4TreeRec::NodeTable::addNode(Node& node)
{
	//if (data.size() == data.capacity())
	//	std::cout << "Increase size of NodeTable!\n";
	data.push_back(node);
	return data.size()-1;
}

uint32_t T4TreeRec::NodeTable::addNode(Node&& node)
{
	//if (data.size() == data.capacity())
	//	std::cout << "Increase size of NodeTable!\n";
	data.push_back(node);	
	return data.size()-1;
}


uint32_t T4TreeRec::NodeTable::size() const
{
	return data.size();
}

void T4TreeRec::NodeTable::print() const
{
	std::cout << "NodeTable:\n";
	for (uint32_t nodeID = 0 ; nodeID < data.size(); ++nodeID)
	{
		std::cout << "\t{" << data[nodeID].generation << ", " << data[nodeID].patch << "}\n";
	}
	std::cout << "\n";
}

template<typename INT>
const T4TreeRec::Node& T4TreeRec::NodeTable::operator[](const INT i) const
{
	return data[i];
}

template<typename INT>
T4TreeRec::Node& T4TreeRec::NodeTable::operator[](const INT i)
{
	return data[i];
}

T4TreeRec::Node& T4TreeRec::NodeTable::back()
{
	return data.back();
}
