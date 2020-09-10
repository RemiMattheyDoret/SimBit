class SampleSequenceData
{
	uint32_t from;
	uint32_t to;
	uint32_t nextL;
	std::vector<uint32_t> lociList; 

public:
	size_t patch;
	size_t ind;
	unsigned char haplo;

	size_t nbLoci() const;
	bool isMoreLoci();
	uint32_t nextLocus();
	void start();

	SampleSequenceData(InputReader& input);
};