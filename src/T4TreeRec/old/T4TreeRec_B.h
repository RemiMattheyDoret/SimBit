

/*class T4TreeRecEntry
{
	friend class T4TreeRec;

	size_t generation;
	std::vector<size_t> RecPoints;
	std::pair<size_t, size_t> parents;
	std::vector<size_t> children;

	T4TreeRecEntry(size_t& g, std::vector<size_t>& RP, std::pair<size_t, size_t>& p);
	T4TreeRecEntry(size_t& g, std::vector<size_t>& RP, std::pair<size_t, size_t>& p, std::vector<size_t>& c);
	void addChild(size_t child);
};*/




class LastGenerationIDmap
{
private:
	size_t min;
	std::vector<size_t> newIDsForOldIDs;

public:
	template<typename INT1, typename INT2>
	LastGenerationIDmap(INT1 nbNodes, INT2 totalPatchSize)
	{
		assert(nbNodes >= totalPatchSize*2);
		min = nbNodes - totalPatchSize*2;
		newIDsForOldIDs.resize(totalPatchSize*2);
	}

	void add(const size_t oldID, const size_t newID)
	{
		auto i = oldID - min;
		assert(i < newIDsForOldIDs.size());
		newIDsForOldIDs[i] = newID;
	}

	size_t getNewID(const size_t oldID) const
	{
		auto i = oldID - min;
		assert(i >= 0 && i < newIDsForOldIDs.size());
		return newIDsForOldIDs[i];
	}
};


class Edge
{
public:
	size_t left;
	size_t right;
	int parent;
	int child;

	Edge(){}
	Edge(size_t l, size_t r, int p, int c):left(l), right(r), parent(p), child(c){}
	Edge(const Edge& e):left(e.left), right(e.right), parent(e.parent), child(e.child){}
	Edge(const Edge&& e):left(e.left), right(e.right), parent(e.parent), child(e.child){}
	Edge operator=(Edge& other)
	{
		left = other.left;
		right = other.right;
		parent = other.parent;
		child = other.child;
		return *this;
	}

	Edge operator=(const Edge&& other)
	{
		left = std::move(other.left);
		right = std::move(other.right);
		parent = std::move(other.parent);
		child = std::move(other.child);
		return *this;
	}

	void swap(Edge& other)
	{
		auto l = left;
		auto r = right;
		auto p = parent;
		auto c = child;

		left = other.left;
		right = other.right;
		parent = other.parent;
		child = other.child;

		other.left = l;
		other.right = r;
		other.parent = p;
		other.child = c;
	}


	bool operator>(const Edge& other) const
	{
		if (this->parent == other.parent)
		{
			if (this->child == other.child)
			{
				return this->left < other.left;
			} else
			{
				return this->child > other.child;
			}
		} else
		{
			return this->parent > other.parent;
		}
	}

	bool operator<(const Edge& other) const
	{
		if (this->parent == other.parent)
		{
			if (this->child == other.child)
			{
				return this->left > other.left;
			} else
			{
				return this->child < other.child;
			}
		} else
		{
			return this->parent < other.parent;
		}
	}
};


class Node
{
public:
	int birth;
	Node(int b):birth(b){}
	Node(const Node& n):birth(n.birth){}
	Node(const Node&& n):birth(n.birth){}
	Node():birth(-1){}
	Node operator=(const Node& other)
	{
		birth = other.birth;
		return *this;
	}
};


/////////////
// Segment //
/////////////


class Segment
{
public:
	size_t left;
	size_t right;
	int newID;
	Segment(size_t l, size_t r, int n):left(l), right(r), newID(n){assert(left<=right);};
	Segment(){};
	Segment operator=(const Segment& other)
	{
		std::cout << ".h line 167\n";
		left = other.left;
		std::cout << ".h line 169\n";
		right = other.right;
		newID = other.newID;
		return *this;
	}

	bool operator>(const Segment& other)const {return left > other.left;}
	bool operator<(const Segment& other)const {return left < other.left;}
};


/*
class SegmentContainer
{
	std::priority_queue<Segment, std::vector<Segment>, std::greater<Segment>> data;

	const Segment& top()
	{
		return data.top();
	}

	void pop()
	{
		data.pop();
	}

	void push(Segment& seg)
	{
		data.push(seg);
	}

	void push(Segment& seg)
	{
		data.push(seg);
	}
};
*/	

class NodeGenetics
{

public:
	std::vector<size_t> mutations;  // T5 style representation
	int birth;
	//std::priority_queue<Segment, std::vector<Segment>, std::less<Segment>> whatHasBeenSetYet;

	NodeGenetics(){}
	NodeGenetics(int g):birth(g){}
	NodeGenetics(const NodeGenetics& n):mutations(n.mutations), birth(n.birth){}
	NodeGenetics(const NodeGenetics&& n):mutations(n.mutations), birth(n.birth){}
	NodeGenetics(std::vector<size_t> m, int b):mutations(m), birth(b){}
	//void assertWasFullySet() const;
	NodeGenetics operator=(const NodeGenetics& other)
	{
		mutations = other.mutations;
		birth = other.birth;
		return *this;
	}

	void clear()
	{
		mutations.clear();
		//whatHasBeenSetYet = std::priority_queue<Segment, std::vector<Segment>, std::less<Segment>>();
	}
	void shrink_to_fit()
	{
		mutations.shrink_to_fit();
	}

	void propagateMutationsForSegment(const NodeGenetics& parent, const size_t left, const size_t right);

	template<typename INT>
	void mutateLocus(INT MutPosition);
};




template<typename T>
class ReversibleDeque
{
private:
	std::deque<T> data;
	bool isEndAtEnd;

public:
	ReversibleDeque():isEndAtEnd(true){}
	ReversibleDeque(const ReversibleDeque& other)
	:data(other.data), isEndAtEnd(other.isEndAtEnd)
	{}

	ReversibleDeque(const ReversibleDeque&& other)
	:data(std::move(other.data)), isEndAtEnd(std::move(other.isEndAtEnd))
	{}

	ReversibleDeque& operator=(const ReversibleDeque& other)
	{
		data = other.data;
		isEndAtEnd = other.isEndAtEnd;
		return *this;
	}

	ReversibleDeque& operator=(const ReversibleDeque&& other)
	{
		data = std::move(other.data);
		isEndAtEnd = std::move(other.isEndAtEnd);
		return *this;
	}

	
	void push_back(const T& elem)
	{
		if (isEndAtEnd)
		{
			data.push_back(elem);
		} else
		{
			data.push_front(elem);
		}
	}

	void push_front(const T& elem)
	{
		if (isEndAtEnd)
		{
			data.push_front(elem);
		} else
		{
			data.push_back(elem);
		}
	}

	template<typename INT>
	T& operator[](const INT i)
	{
		if (isEndAtEnd)
		{
			return data[i];
		} else
		{
			return data[data.size()-i-1];
		}	
	}

	void reverse()
	{
		isEndAtEnd = !isEndAtEnd;
	}

	size_t size(){return data.size();} const

	void swap(ReversibleDeque& other)
	{
		this->data.swap(other.data);
	}

	void clear()
	{
		data.clear();
	}

	T& back() 
	{
		if (isEndAtEnd)
		{
			return data.back();
		} else
		{
			return data.front();
		}	
	}

	T& front() 
	{
		if (isEndAtEnd)
		{
			return data.front();
		} else
		{
			return data.back();
		}	
	}
};


template<typename T>
class ReversibleVector
{
private:
	std::vector<T> data_back;
	std::vector<T> data_front;
	bool isEndAtEnd;

public:
	ReversibleVector():isEndAtEnd(true){}
	ReversibleVector(const ReversibleVector& other)
	:data_back(other.data_back), data_front(other.data_front), isEndAtEnd(other.isEndAtEnd)
	{}

	ReversibleVector(const ReversibleVector&& other)
	:data_back(std::move(other.data_back)), data_front(std::move(other.data_front)), isEndAtEnd(std::move(other.isEndAtEnd))
	{}

	ReversibleVector& operator=(const ReversibleVector& other)
	{
		data_back = other.data_back;
		data_front = other.data_front;
		isEndAtEnd = other.isEndAtEnd;
		return *this;
	}

	ReversibleVector& operator=(const ReversibleVector&& other)
	{
		data_back = std::move(other.data_back);
		data_front = std::move(other.data_front);
		isEndAtEnd = std::move(other.isEndAtEnd);
		return *this;
	}

	
	void push_back(const T& elem)
	{
		if (isEndAtEnd)
		{
			data_back.push_back(elem);
		} else
		{
			data_front.push_back(elem);
		}
	}

	void push_front(const T& elem)
	{
		if (isEndAtEnd)
		{
			data_front.push_back(elem);
		} else
		{
			data_back.push_back(elem);
		}
	}

	template<typename INT>
	T& operator[](const INT i)
	{
		/*
			[.....|..]
			[..|.....]
		*/
		if (isEndAtEnd)
		{
			if (i >= data_front.size())
			{
				return data_back[i - data_front.size()];
			} else
			{
				return data_front[data_front.size() - i - 1];
			}
		} else
		{
			if (i >= data_back.size())
			{
				return data_front[i - data_back.size()];
			} else
			{
				return data_back[data_back.size() - i - 1];
			}
		}	
	}

	void reverse()
	{
		isEndAtEnd = !isEndAtEnd;
	}

	size_t size(){return data.size();} const

	void swap(ReversibleVector& other)
	{
		this->data_back.swap(other.data_back);
		this->data_front.swap(other.data_front);
		this->isEndAtEnd.swap(other.isEndAtEnd);
	}

	void clear()
	{
		data_back.clear();
		data_front.clear();
	}

	T& back() 
	{
		/*
			[.....|..]
			[..|.....]
		*/
		if (isEndAtEnd)
		{
			return data_back.back();
		} else
		{
			return data_front.back();
		}	
	}

	T& front() 
	{
		/*
			[.....|..]
			[..|.....]
		*/
		if (isEndAtEnd)
		{
			return data_front.back();
		} else
		{
			return data_back.back();
		}	
	}
};


class NodeTable
{
public:
	ReversibleVector<Node> data;

	void swap(NodeTable& other)
	{
		this->data.swap(other.data);
	}

	void clear()
	{
		data.clear();
	}

	void reverseData()
	{
		data.reverse();
	}

	int addNode(Node input)
	{
		return addNode(input.birth);
	}

	int addNode(int gen)
	{
		data.push_back({gen});
		//std::cout << "Added node birth = " << r.birth << " ID = " << r.ID << "\n";
		return data.size()-1;
	}

	size_t size()
	{
		return data.size();
	}

	void print()
	{
		std::cout << "NodeTable:\n";
		for (size_t nodeID = 0 ; nodeID < data.size(); ++nodeID)
		{
			std::cout << "\t" << data[nodeID].birth << "\n";
		}
		std::cout << "\n";
	}

	template<typename INT>
	Node& operator[](INT i)
	{
		return data[i];
	}

	Node& back()
	{
		return data.back();
	}

	//NodeTable(std::vector<Node> d):data(d){}
	NodeTable(){}

};

class EdgeTable
{
public:
	std::vector<Edge> data;

	void swap(EdgeTable& other)
	{
		this->data.swap(other.data);
	}

	void clear()
	{
		data.clear();
	}

	void addEdge(const Edge& input)
	{
		data.push_back(input);
	}

	void addEdge(const Edge input)
	{
		data.push_back(input);
	}

	void addEdge(size_t l, size_t r, int p, int c)
	{
		data.push_back({l, r, p, c});
	}
	
	size_t size()
	{
		return data.size();
	}

	void print()
	{
		std::cout << "EdgeTable:\n\tleft\tright\tparent\tchild\n";
		for (auto& edge : data)
		{
			std::cout << "\t" << edge.left << "\t" << edge.right << "\t" << edge.parent << "\t" << edge.child << "\n";
		}
		std::cout << "\n";
	}

	template<typename INT>
	Edge& operator[](INT i)
	{
		return data[i];
	}

	Edge& back()
	{
		return data.back();
	}

	EdgeTable(std::vector<Edge> d):data(d){}
	EdgeTable(){}
};

class T4TreeRec
{
private:
	EdgeTable edges;
	NodeTable nodes;

	std::vector<NodeGenetics> ancestorsGenetics;
	std::vector<size_t> ancestorsID;
	std::vector<bool> isAncestor;
	int lastGenerationSimplified = -1;
	//std::vector<size_t> ancestorGeneticsMapToCurrentID; // Watch out if I put this back on because it si not initliaze anymore. ancestorGeneticsMapToCurrentID[currentID] = ancestralID;
	
	void setPopToUniqueID(Pop& pop, size_t T4ID);
	int searchForClonesAmongAncestors(const NodeGenetics& input) const ;
	void simplify_mergeSegmentsGoingThroughNode(std::vector<Segment>& segments);
	void simplify_transferChildSegments(std::vector<Segment>& parentSegments, std::vector<Segment>& childSegments, size_t left, size_t right);
	void simplify(Pop& pop);
	void printInfoForDebug(Pop& pop)  ;

public:
	T4TreeRec(){}
	T4TreeRec(EdgeTable ET, NodeTable NT):edges(ET), nodes(NT){}
	void initialize(Pop& pop);
	size_t addHaplotype(std::vector<int>& RP, std::pair<size_t, size_t> p);
	size_t addAncestor(NodeGenetics& newAncestor);
	void simplify_ifNeeded(Pop& pop);
	void print();
	std::vector<std::vector<std::vector<size_t>>> placeMutations(Pop& pop, bool isNeedSimplify = true);
	
};






