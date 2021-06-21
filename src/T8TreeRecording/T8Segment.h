class T8Segment
{
	friend class T8TreeRecording;
	template<class T> friend class MemoryPool;

	T8Segment* parent = nullptr;
	std::vector<uint32_t> geneticData;
	fitnesstype w = 1.0;
	uint16_t nbChildren = 0;

	// Statics for performance
	//static size_t static_mutate_position;
	//static size_t static_mutate_newsize;
	//static size_t static_mutate_mut;

	static T8Segment* buildChild(T8Segment* parent, T8Segment* child);
	static T8Segment* buildChild(std::vector<uint32_t>& genData, T8Segment* child);
	template<typename INT> void mutate(const INT segmentIndex, bool mustBeSorted);
	void forceComputeFitness();
	void clear_shrink();
	void shrink_to_fit();
	void clear();
	void printGeneticData();
};