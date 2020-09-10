
void T4TreeRec::EdgeTable::swap(T4TreeRec::EdgeTable& other)
{
	this->data.swap(other.data);
}

void T4TreeRec::EdgeTable::clear()
{
	data.clear();
}

/*void T4TreeRec::EdgeTable::addEdge(const T4TreeRec::Edge& input)
{
	data.push_back(input);
}*/

void T4TreeRec::EdgeTable::addEdge(const T4TreeRec::Edge input)
{
	//if (data.size() == data.capacity()) std::cout << "must reallocate edges!\n";
	data.push_back(input);
}

void T4TreeRec::EdgeTable::addEdge(uint32_t l, uint32_t r, int p, int c)
{
	//if (data.size() == data.capacity()) std::cout << "must reallocate edges!\n";
	data.push_back({l, r, p, c});
}

uint32_t T4TreeRec::EdgeTable::size() const
{
	return data.size();
}

void T4TreeRec::EdgeTable::reverseAndMerge(const uint32_t nbNodes)
{
	// Must reverse both the ordering and the IDs from newID to newnewID
	auto first = data.begin();
	auto last = data.end() - 1;

	// First change the IDs and reverse the vector
	while (first < last)
	{
		// Change IDs
		first->parent = nbNodes - first->parent - 1;
		first->child  = nbNodes - first->child - 1;
		last->parent  = nbNodes - last->parent - 1;
		last->child   = nbNodes - last->child - 1;

		// Reverse ordering
		std::swap(first->parent, last->parent);
		std::swap(first->child, last->child);
		std::swap(first->left, last->left);
		std::swap(first->right, last->right);
		++first;
		--last;
	}

	if (first == last)
	{
		first->parent = nbNodes - first->parent - 1;
		first->child  = nbNodes - first->child - 1;
	}

	// Second merge adjacent edges that need to be merged
	#ifdef DEBUG
	std::cout << "before mergeing\n";
	print();
	#endif

	/*
	std::sort(data.begin(), data.end(), [](const Edge& lhs, const Edge& rhs) -> bool
	{
		if (lhs.parent == rhs.parent)
			return lhs.child < rhs.child;
		else 
			return lhs.parent < rhs.parent;
	}
	);
	*/
	
	

	for (uint32_t i = 0 ; i < data.size()-1 ; ++i)
	{
		auto l = i;
		while (data.size() < l+1 && data[l].parent == data[l+1].parent && data[l].child == data[l+1].child && data[l].left == data[l+1].right)
		{
			++l;
		}
		if (l != i)
		{
			auto newLeft = data[l].left;
			assert(newLeft < data[i].left);
			data[i].left = newLeft;
			#ifdef DEBUG
			std::cout << "i = " << i << "\n";
			std::cout << "l = " << l << "\n";
			#endif
			data.erase(data.begin() + i + 1, data.begin() + l + 1);
		}
	}
	#ifdef DEBUG
	std::cout << "after mergeing\n";
	print();
	#endif
	/*
	data[0].print();
	data[1].print();
	data[2].print();
	data[3].print();
	*/
	
	/*
	// change IDs from newID to newnewID
	for (uint32_t i = 0 ; i < data.size() ; ++i)
	{
		data[i].parent = nbNodes - data[i].parent - 1;
		data[i].child  = nbNodes - data[i].child - 1;
	}

	std::sort(
		data.begin(),
		data.end(),
		[](const Edge& lhs, const Edge& rhs)
		{
			lhs.parent < rhs.parent
		}
	);
	*/
}

void T4TreeRec::EdgeTable::print() const
{
	std::cout << "EdgeTable:\n\tleft\tright\tparent\tchild\n";
	for (auto& edge : data)
	{
		std::cout << "\t" << edge.left << "\t" << edge.right << "\t" << edge.parent << "\t" << edge.child << "\n";
	}
	std::cout << "\n";
}

template<typename INT>
const T4TreeRec::Edge& T4TreeRec::EdgeTable::operator[](INT i) const
{
	return data[i];
}

template<typename INT>
T4TreeRec::Edge& T4TreeRec::EdgeTable::operator[](INT i)
{
	return data[i];
}

T4TreeRec::Edge& T4TreeRec::EdgeTable::back()
{
	return data.back();
}

T4TreeRec::EdgeTable::EdgeTable(std::vector<T4TreeRec::Edge> d):data(d){}
T4TreeRec::EdgeTable::EdgeTable(){}

void T4TreeRec::EdgeTable::reserve(uint32_t size)
{
	data.reserve(size);
}
