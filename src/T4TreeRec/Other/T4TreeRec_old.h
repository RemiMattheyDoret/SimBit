

/*class T4TreeRecEntry
{
	friend class T4TreeRec;

	uint32_t generation;
	std::vector<uint32_t> RecPoints;
	std::pair<uint32_t, uint32_t> parents;
	std::vector<uint32_t> children;

	T4TreeRecEntry(uint32_t& g, std::vector<uint32_t>& RP, std::pair<uint32_t, uint32_t>& p);
	T4TreeRecEntry(uint32_t& g, std::vector<uint32_t>& RP, std::pair<uint32_t, uint32_t>& p, std::vector<uint32_t>& c);
	void addChild(uint32_t child);
};*/

class Edge
{
public:
	uint32_t left;
	uint32_t right;
	int parent;
	int child;

	Edge(){}
	Edge(uint32_t l, uint32_t r, int p, int c):left(l), right(r), parent(p), child(c){}
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
	uint32_t left;
	uint32_t right;
	int node;
	Segment(uint32_t l, uint32_t r, int n):left(l), right(r), node(n){assert(left<=right);};
	Segment(){};
	Segment operator=(const Segment& other)
	{
		left = other.left;
		right = other.right;
		node = other.node;
		return *this;
	}

	bool operator>(const Segment& other)const {return left > other.left;}
	bool operator<(const Segment& other)const {return left < other.left;}
};

class NodeGenetics
{

public:
	std::vector<uint32_t> mutations;  // T5 style representation
	int birth;
	//std::priority_queue<Segment, std::vector<Segment>, std::less<Segment>> whatHasBeenSetYet;

	NodeGenetics(){}
	NodeGenetics(int g):birth(g){}
	NodeGenetics(const NodeGenetics& n):mutations(n.mutations), birth(n.birth){}
	NodeGenetics(const NodeGenetics&& n):mutations(n.mutations), birth(n.birth){}
	NodeGenetics(std::vector<uint32_t> m, int b):mutations(m), birth(b){}
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

	void propagateMutationsForSegment(const NodeGenetics& parent, const uint32_t left, const uint32_t right);

	template<typename INT>
	void mutateLocus(INT MutPosition);
};


class NodeTable
{
public:
	std::vector<Node> data;

	void swap(NodeTable& other)
	{
		this->data.swap(other.data);
	}

	void clear()
	{
		data.clear();
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

	uint32_t size()
	{
		return data.size();
	}

	void print()
	{
		std::cout << "NodeTable:\n";
		for (auto& node : data)
		{
			std::cout << "\t" << node.birth << "\n";
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

	NodeTable(std::vector<Node> d):data(d){}
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

	void addEdge(uint32_t l, uint32_t r, int p, int c)
	{
		data.push_back({l, r, p, c});
	}
	
	uint32_t size()
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
	std::vector<uint32_t> ancestorsID;
	std::vector<bool> isAncestor;
	int lastGenerationSimplified = -1;
	//std::vector<uint32_t> ancestorGeneticsMapToCurrentID; // Watch out if I put this back on because it si not initliaze anymore. ancestorGeneticsMapToCurrentID[currentID] = ancestralID;
	
	void setPopToUniqueID(Pop& pop, uint32_t T4ID);
	int searchForClonesAmongAncestors(const NodeGenetics& input) const ;
	void simplify(Pop& pop);
	void printInfoForDebug(Pop& pop)  ;

public:
	T4TreeRec(){}
	T4TreeRec(EdgeTable ET, NodeTable NT):edges(ET), nodes(NT){}
	void initialize(Pop& pop);
	uint32_t addHaplotype(std::vector<uint32_t>& RP, std::pair<uint32_t, uint32_t> p);
	uint32_t addAncestor(NodeGenetics& newAncestor);
	void simplify_ifNeeded(Pop& pop);
	void print();
	std::vector<std::vector<std::vector<uint32_t>>> placeMutations(Pop& pop, bool isNeedSimplify = true, bool gatherAndOutputData = true);

	std::vector<uint32_t> getMutationsOfID(const uint32_t id) const;
	//template<typename INT>
	//bool getAlleleOfID(const uint32_t id, const INT locus) const; Thatt would be too slow

	std::string getPaintedHaplotype(){return "Painting chromosome is not available in current version";}
	
};

