class SampleSequenceDataContainer
{
	size_t SSD_index;
	std::vector<SampleSequenceData> SSDs;

public:
	void start();
	SampleSequenceData& getNextSequence();
	bool areThereMoreSequences() const;
	void readInput(InputReader& input);
};