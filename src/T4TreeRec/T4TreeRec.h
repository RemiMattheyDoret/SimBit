
class T4TreeRec
{
public:
	/*
		#####################
		### Inner classes ###
		#####################
	*/

	class Segment
	{
		public:
			uint32_t left;
			uint32_t right;
			int child;
			Segment(uint32_t l, uint32_t r, int c);
			friend bool operator<(const Segment& lhs, const Segment& rhs)
			{
				return lhs.left < rhs.left;
			}
			void print() const;
	};

	struct PointForFindingOverlap
	{
	    uint32_t location; 
	    bool overlap;    // does this point start/close a new segment
	    int ID;
	};

	struct PaintedSegment
	{
		Segment segment; // segment.child is T4ID of focal child in current generation
		int patch_index; // patch of ancestor
		int ind_ID;      // T4ID of ancestor

		friend bool operator<(const PaintedSegment& lhs, const PaintedSegment& rhs)
		{
			return lhs.segment < rhs.segment;
		}
	};



	/*
	struct PaintedSegmentFrequency
	{
		uint32_t left;
		uint32_t right;
		double freq;
	};
	*/

	struct PaintedSegmentDiversity
	{
		uint32_t left;
		uint32_t right;
		size_t nbColors;
		int colorHighestFrequency;
		size_t freqOfColorHighestFrequency;
		double heterozygosity;

		PaintedSegmentDiversity(){};
		PaintedSegmentDiversity(uint32_t l, uint32_t r, size_t c, int mainColor, size_t freqMainColor, double h);

		std::string toString() const;
	};

	struct PaintedSegmentOrigin
	{
		size_t color;
		double contribution; 
		int originPatch;

		PaintedSegmentOrigin(){};
		PaintedSegmentOrigin(size_t a, double b, int c)
		:color(a), contribution(b), originPatch(c)
		{
			assert(contribution > 0.0);
		};

		std::string toString() const;
	};

	struct PaintedSegmentSummary
	{
		PaintedSegmentDiversity paintedSegmentDiversity;
		PaintedSegmentOrigin paintedSegmentOrigin;

		PaintedSegmentSummary(PaintedSegmentDiversity a, PaintedSegmentOrigin b)
		:paintedSegmentDiversity(a), paintedSegmentOrigin(b)
		{}

		std::string toString(unsigned char outputType) const;
	};


	class LastGenerationIDmap
	{
		private:
			uint32_t min;
			std::vector<uint32_t> newIDsForOldIDs;

		public:
			template<typename INT1, typename INT2> LastGenerationIDmap(INT1 nbNodes, INT2 totalPatchSize);
			void add(const uint32_t oldID, const uint32_t newID);
			uint32_t getNewID(const uint32_t oldID) const;
	};

	template<typename T>
	class ReversibleVector
	{
		private:
			std::vector<T> data_back;
			std::vector<T> data_front;
			bool isEndAtEnd;

		public:
			ReversibleVector<T>();
			ReversibleVector<T>(const ReversibleVector<T>& other);
			ReversibleVector<T>(const ReversibleVector<T>&& other);
			ReversibleVector<T>& operator=(const ReversibleVector<T>& other);
			ReversibleVector<T>& operator=(const ReversibleVector<T>&& other);
			void reserve(uint32_t n);
			void push_back(const T& elem);
			void push_front(const T& elem);
			template<typename int_type> const T& operator[](const int_type i) const;
			template<typename int_type> T& operator[](const int_type i);
			void reverse();
			void shrink_to_fit();
			uint32_t size() const;
			void swap(ReversibleVector<T>& other);
			void clear();
			T& back();
			T& front();
	};

	class NodeGenetics
	{
		private:
			std::vector<uint32_t> mutations;
			int birth;
			int ID;

		public:
			friend bool operator<(NodeGenetics lhs, NodeGenetics rhs)
			{
				return lhs.ID < rhs.ID;
			}
			NodeGenetics();
			NodeGenetics(int id);
			NodeGenetics(std::vector<uint32_t> m, int b, int id);
			NodeGenetics(std::vector<uint32_t>& m, int b, int id);
			void resetID(int newID);
			void clear();
			void shrink_to_fit();
			int getID() const;
			int getGeneration() const;
			void setGeneration(const int g);
			void print();

			void push_backMutations(std::vector<uint32_t>& m);
			void push_backMutations(std::vector<uint32_t>&& m);
			std::vector<uint32_t>& getMutations();
			const std::vector<uint32_t>& getMutations() const;
			std::vector<uint32_t>&& moveMutations();
			template<typename INT> void mutateLocus(INT MutPosition);
			template<typename INT> NodeGenetics getSubsetMutations(INT from, INT to) const;

			void PrintBinaryFile(OutputFile& file) const;
			void readFromBinaryFile(BinaryFileToRead& binfile);
	};

	class Edge
	{
		public:
			uint32_t left;
			uint32_t right;
			int parent;
			int child;

			void swap(Edge& other);
			bool operator>(const Edge& other) const;
			bool operator<(const Edge& other) const;
			void print() const;
			void PrintBinaryFile(OutputFile& file) const;
			void readFromBinaryFile(BinaryFileToRead& binfile);
	};

	class EdgeTable
	{
		private:
			std::vector<Edge> data;
		public:

			void swap(EdgeTable& other);
			void clear();
			//void addEdge(const Edge& input);
			void addEdge(const Edge input);
			void addEdge(uint32_t l, uint32_t r, int p, int c);
			uint32_t size() const;
			void reverseAndMerge(const uint32_t nbNodes);
			void print() const;
			template<typename INT> const Edge& operator[](INT i) const;
			template<typename INT> Edge& operator[](INT i);
			Edge& back();
			EdgeTable(std::vector<Edge> d);
			EdgeTable();
			void reserve(uint32_t size);
			void deleteAll();
			void PrintBinaryFile(OutputFile& file) const;
			void readFromBinaryFile(BinaryFileToRead& binfile);
	};

	struct Node
	{
		int generation;
		int patch;

		void PrintBinaryFile(OutputFile& file) const
		{
			file.writeBinary(generation);
			file.writeBinary(patch);
		}

		void readFromBinaryFile(BinaryFileToRead& binfile)
		{
			binfile.read(generation);
			binfile.read(patch);
		}
	};

	class NodeTable
	{
		private:
			std::vector<Node> data;
		public:

			NodeTable();
			void swap(NodeTable& other);
			void clear();
			void reserve(uint32_t n);
			void reverse();
			uint32_t addNode(Node& node);
			uint32_t addNode(Node&& node);
			uint32_t size() const;
			void print() const;
			template<typename INT> const Node& operator[](const INT i) const;
			template<typename INT>  Node& operator[](const INT i);
			Node& back();
			void deleteAll();
			void PrintBinaryFile(OutputFile& file) const;
			void readFromBinaryFile(BinaryFileToRead& binfile);
	};

	class HaplotypeOfSegments
	{
		private:
			std::vector<Segment> segments;
			Node myNode;
			int newNodeID = -1;
			bool alreadySegmentized = false;

			void assertOrderingAndNoOverlap();
			void removeMovingSequence(Segment& movingSegment, uint32_t& childSegmentIndex);
			Segment getMovingSegment(const Edge& edge, Segment& childSegment);
			
		public:

			void segmentize(NodeTable& nodeTable, EdgeTable& edgeTable, bool shouldNodeBeKeptIfItHasASegment);
			void print() const;
			std::vector<Segment>::iterator begin();
			std::vector<Segment>::iterator end();
			uint32_t size() const;
			template<typename int_type> Segment& operator[](int_type& i);
			int getNewNodeID() const;
			int getGeneration() const;
			int getPatchIndex() const;
			bool getAlreadySegmentized() const;
			void reset_alreadySegmentized(bool b);
			void set_newNodeID(int id);
			void sortAndMerge();

			HaplotypeOfSegments(std::vector<Segment>& segs, Node&& n);
			HaplotypeOfSegments(std::vector<Segment>& segs, Node& n);
			HaplotypeOfSegments(std::vector<Segment> segs, Node n);
			HaplotypeOfSegments(Node n);

			std::vector<Segment>& getSegmentsRef();
			
			void transferSegmentsFromChild(HaplotypeOfSegments& childSegments, const Edge& edge, EdgeTable& newEdgeTable, NodeTable& newNodeTable, bool shouldNodeBeKeptIfItHasASegment);
			void simpleTransferSegmentsFromChild(HaplotypeOfSegments& childSegments, const Edge& edge, bool isOrdered = true);
	};


	template<typename node_t>
	class HaplotypesContainer
	{
		struct ObjForIteration
		{
			node_t* haploP;
			int oldID;
			int newID;
			ObjForIteration(node_t* a, int b, int c):haploP(a), oldID(b), newID(c){}
		};

		private:
			uint32_t iterator_index;
			std::deque<node_t*> data; // vectors of pointers so to not allocate memory about size of container and how much memory is reserved
			int biggestID;
			// first element of data is the last ID and it goes decreasing in IDs

			template<typename nodeID_type>
			nodeID_type getIndex(nodeID_type& i) const;

		public:
			HaplotypesContainer();
			const std::deque<node_t*>& getData() const;
			template<typename nodeID_type> node_t* getHaploP(nodeID_type& nodeID) const;
			template<typename nodeID_type> bool doesAlreadyExist(nodeID_type& nodeID) const;
			template<typename nodeID_type> node_t& getHaplo(nodeID_type& nodeID) const;
			template<typename nodeID_type> void deleteHaplo(nodeID_type& nodeID);
			template<typename nodeID_type> void insertHaploP(nodeID_type nodeID, node_t* H);
			void deleteAllHaplos();
			void shrink_to_fit();
			void iterator_restart();
			ObjForIteration iterator_next();
			bool iterator_isMore();
			uint32_t size1() const;
			uint32_t size2() const;
	};


	class PatchCountData
	{
		public:
			std::vector<uint32_t> vecFreqs;
			std::map<uint32_t, uint32_t> mapFreqs;
			bool isVec;

			template<typename INT> uint32_t& operator[](INT i)
			{
				if (isVec) return vecFreqs[i]; else return mapFreqs[i];
			}

			template<typename INT> void resize(INT n)
			{
				if (isVec) vecFreqs.resize(n);
			}

			template<typename INT> void resize(INT n, double value)
			{
				if (isVec) vecFreqs.resize(n, value);
			}

			template<typename INT> PatchCountData(bool isMuchDiversity, INT nbLoci)
			{
				isVec = isMuchDiversity;
				this->resize(nbLoci, 0.0);
			}

			void clear()
			{
				if (isVec)
				{
					std::vector<uint32_t>().swap(vecFreqs);
				} else
				{
					std::map<uint32_t, uint32_t>().swap(mapFreqs);
				}
			}


	};

	/*
		End of inner classes
	*/


private:
	EdgeTable edges;
	NodeTable nodes;
	int ancestralGeneration;

	std::vector<NodeGenetics> ancestorsGenetics;
	int lastGenerationSimplified = -1;
	unsigned char haveAllLociCoallesced_info = 'u';
	//std::vector<uint32_t> ancestorGeneticsMapToCurrentID; // Watch out if I put this back on because it si not initliaze anymore. ancestorGeneticsMapToCurrentID[currentID] = ancestralID;
	
	void setPopToUniqueID(Pop& pop, uint32_t T4ID);
	int searchForClonesAmongAncestors(const NodeGenetics& input) const ;
	void simplify_mergeSegmentsGoingThroughNode(std::vector<Segment>& segments);
	void simplify_transferChildSegments(std::vector<Segment>& parentSegments, std::vector<Segment>& childSegments, uint32_t left, uint32_t right);
	


	void doStuffWithAncestors(HaplotypesContainer<HaplotypeOfSegments>& A, NodeTable& No, EdgeTable& Eo);
	//void resetAncestorsEdgesAndNodesAfterPlacingMutations(std::vector<NodeGenetics>& allNodes, bool shouldCompleteDeleteAllTree);

	static void propagateMutationsForSegment(NodeGenetics& child, const NodeGenetics& parent, const uint32_t left, const uint32_t right);
	//void printInfoForDebug(Pop& pop)  ;

	std::vector<std::vector<std::vector<T4TreeRec::PaintedSegment>>> computePaintedHaplotypes(const int paintedGeneration, const std::vector<uint32_t>& focalT4IDs, const std::vector<uint32_t>& focalT4IDs_patches) const; // computePaintedHaplotypes()[patch_index][T4ID_index][segment_index]
	std::vector<T4TreeRec::PaintedSegmentSummary> computeSmallSegmentsSummaryFromPaintedHaplotypes(const std::vector<std::vector<T4TreeRec::PaintedSegment>>& allGatheredSegments, const unsigned char outputType) const;
	/*std::vector<T4TreeRec::PaintedSegmentFrequency> computeLargeSegmentsFrequenciesFromPaintedHaplotypes(const std::vector<std::vector<T4TreeRec::PaintedSegment>>& allGatheredSegments) const;*/

	void computePaintedHaplotypes_exploreTree(HaplotypesContainer<HaplotypeOfSegments>& allSegments, const int paintedGeneration) const;

	std::vector<std::vector<std::vector<PaintedSegment>>> computePaintedHaplotypes_gatherSegments(HaplotypesContainer<HaplotypeOfSegments>& allSegments, const std::vector<uint32_t>& focalT4IDs, const std::vector<uint32_t>& focalT4IDs_patches) const;

	int killOnDemand_T4Locus = -1;
	bool killOnDemand_isT4LocusFixed = false;
	bool isAlreadyAsSimplifiedAsPossible = true;

	std::vector<std::vector<std::vector<uint32_t>>> getMutations(Pop& pop);
	template<typename INT> const std::vector<uint32_t>& getMutationsOfID(INT ID);

	//void setPartOfNewAncestor(std::vector<NodeGenetics>& allNodes, std::vector<NodeGenetics>& newAncestorsGenetics);

public:
	static std::vector<int> generationsToKeepInTheTree;

	T4TreeRec(){edges.reserve(100); nodes.reserve(100);}
	T4TreeRec(EdgeTable ET, NodeTable NT):edges(ET), nodes(NT){}
	void initialize(Pop& pop);
	uint32_t addHaplotype(std::vector<uint32_t>& RP, std::pair<uint32_t, uint32_t> p, int patch_index);
	uint32_t addAncestor(NodeGenetics& newAncestor, int patch_index, int generation);
	void simplify_ifNeeded(Pop& pop);
	void simplify(Pop& pop);
	void print() const;
	bool haveAllLociCoallesced_ifKnown() const;
	
	//template<typename INT> const std::vector<uint32_t>& getMutationsOfID(std::vector<std::vector<std::vector<uint32_t>>>& mutations, Pop& pop, INT ID);
	std::pair< std::vector<std::vector<std::vector<uint32_t>>>, std::vector<PatchCountData> > placeMutations(Pop& pop, bool shouldDeleteTree, bool shouldCompleteDeleteAllTree, unsigned char outputType);
	void writePaintedHaplotypes(const std::vector<uint32_t>& focalT4IDs, const std::vector<uint32_t>& focalT4IDs_patches, const int paintedGeneration, const int observedGeneration, OutputFile& file) const;
	void writePaintedHaplotypesSummary(const std::vector<uint32_t>& focalT4IDs, const std::vector<uint32_t>& focalT4IDs_patches, const int paintedGeneration, const int observedGeneration, OutputFile& file, const unsigned char outputType ) const;
	void shift_generations_after_burn_in();

	void setLocusForWhichFixationMustBeComputedAtTheNextSimplify(int locus);
	bool isLocusForWhichFixationHadToBeComputedFixed();
	void assertIsFullySimplified() const;
	void deleteAll();
	void PrintBinaryFile(OutputFile& file) const;
	void readFromBinaryFile(BinaryFileToRead& binfile);
};

