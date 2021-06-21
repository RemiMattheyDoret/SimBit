typedef uint32_t T8ID_type;

class T8TreeRecording
{
public:
	static void printT8IDs(Pop& pop);
	void printT8IDs();
private:

	friend class T8Segment;

	std::vector<std::vector<T8Segment*>> lastGenerationSegments;
	std::vector<std::vector<T8Segment*>> currentGenerationSegments;
	size_t currentGenerationIDcounter = 0;
	//size_t nbGenerationsSinceLastSimplify = 0;
	size_t lastGenerationMemPoolsWereShrunkToFit;
	MemoryPool<T8Segment> memPool_segments;
	int lastTimeMutationsWerePlaced;
	static std::vector<uint32_t> staticGeneticData;


	void prune(std::vector<std::vector<T8Segment*>>& childrenSegments);
	void propagateMutations(T8Segment* focal);
	//std::vector<uint32_t> computeGeneticDataForRecombiningSegment(T8Segment* mother, T8Segment* father, bool& isMother, size_t& beginningOfSegment, std::vector<int>& recPositions, size_t& recIndex);
	void placeMutations();
	void moveGeneticDataFromTo(T8Segment* giver, T8Segment* receiver);
	void copyGeneticDataFromTo(T8Segment* giver, T8Segment* receiver);
	void moveSortedGeneticDataFromTo(T8Segment* giver, T8Segment* receiver);
	void copySortedGeneticDataFromTo(T8Segment* giver, T8Segment* receiver);
	void deleteNode(T8Segment* node);

	//void deleteParentInProfitOfChild(T8Segment*& parent, T8Segment*& child);

	// Statics-like - just to avoid reallocating
	std::stack<T8Segment*> pointersToChildren; 
	std::vector<uint32_t> geneticDataAtRecPlaces;
	//std::vector<size_t> hash;

	/*
	//size_t SL_childID;
	//double SL_haplotypeFitness;
	bool SL_isMother;
	size_t SL_recIndex;
	uint32_t SL_beginningOfSegment;
	uint32_t SL_endOfSegment;
	uint32_t SL_from;
	uint32_t SL_to;
	size_t SL_segmentIndex;
	*/


public:
	T8TreeRecording();
	~T8TreeRecording();
	void initialize(Pop& pop);
	void announceStartOfGeneration();
	void announceEndOfGeneration();
	void deleteTree();
	//void simplify(); // made public just in case but should be called from 'announceEndOfGeneration'
	std::pair<T8ID_type, fitnesstype> addOffspring(T8ID_type motherID, T8ID_type fatherID, std::vector<uint32_t>& recPositions);
	std::vector<std::vector<std::vector<uint32_t>>> getMutationsData(Pop& pop, bool shouldDeleteTree);
	void computeFrequenciesFromRawData(std::vector<std::vector<std::vector<uint32_t>>>& data, std::vector<std::map<uint32_t, double>>& freqs, bool shouldClearData);
};
