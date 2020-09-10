
class GeneticMap_ChangeDescription
{
public:
	bool willIChange;
	uint32_t TxLocus;
};

class T56LocusGender
{
public:
	bool isNtrl;
	uint32_t locusInGender;

	T56LocusGender(bool b, uint32_t i):isNtrl(b), locusInGender(i){}
};

class LocusDescription
{
public:
	unsigned char type;
	uint32_t locusInType;

	LocusDescription(unsigned char t, uint32_t l):type(t), locusInType(l){}
};


class GeneticMap
{
private:

	std::vector<uint32_t> _FromT1LocusToLocus;
	std::vector<uint32_t> _FromT2LocusToLocus;
	std::vector<uint32_t> _FromT3LocusToLocus;
	std::vector<uint32_t> _FromT4LocusToLocus;
	std::vector<uint32_t> _FromT56LocusToLocus;
	std::vector<uint32_t> _FromT56ntrlLocusToLocus;
	std::vector<uint32_t> _FromT56selLocusToLocus;
	std::vector<uint32_t> _FromLocusToNextT1Locus;
	std::vector<uint32_t> _FromLocusToNextT2Locus;
	std::vector<uint32_t> _FromLocusToNextT3Locus;
	std::vector<uint32_t> _FromLocusToNextT4Locus;
	std::vector<uint32_t> _FromLocusToNextT56Locus;
	std::vector<uint32_t> _FromLocusToNextT56ntrlLocus;
	std::vector<uint32_t> _FromLocusToNextT56selLocus;
	std::vector<T56LocusGender>  _FromT56LocusToT56genderLocus;

	unsigned char onlyTypeUsed;

	bool isT1Used = false;
	bool isT2Used = false;
	bool isT3Used = false;
	bool isT4Used = false;
	bool isT56Used = false;
	bool isT56ntrlUsed = false;
	bool isT56selUsed = false;

	uint32_t nbT56LociGotFromFitnessValues = 0; // Only used to check inputs


	void shrink_to_fit();

	void print();
	template<typename T>
	void print(std::string name, std::vector<T>& v);
	void print(std::string name, std::vector<T56LocusGender>& v);
	template<typename T>
	void print(std::string name, T& v);

	int readType(const std::string& s_type) const;

	template<typename INT>
	uint32_t FromLocusToNextTxLocus(const INT i, const bool isTxUsed, const std::vector<uint32_t>& map, const unsigned char typeValue) const;
	



public:

	uint32_t T1_nbLoci;
	uint32_t T1_nbChars;
	uint32_t T1_nbLociLastByte;
	uint32_t T2_nbLoci;
	uint32_t T3_nbLoci;
	uint32_t T4_nbLoci;
	uint32_t T5_nbLoci;
	uint32_t T5sel_nbLoci;
	uint32_t T5ntrl_nbLoci;
	uint32_t T6_nbLoci;
	uint32_t T6sel_nbLoci;
	uint32_t T6ntrl_nbLoci;
	uint32_t T56_nbLoci;
	uint32_t T56sel_nbLoci;
	uint32_t T56ntrl_nbLoci;
	uint32_t TotalNbLoci;

	bool isT56selCompress;
	bool isT56ntrlCompress;


	template<typename INT>
	uint32_t FromT1LocusToLocus(const INT i) const;
	template<typename INT>
	uint32_t FromT2LocusToLocus(const INT i) const;
	template<typename INT>
	uint32_t FromT3LocusToLocus(const INT i) const;
	template<typename INT>
	uint32_t FromT4LocusToLocus(const INT i) const;
	template<typename INT>
	uint32_t FromT56LocusToLocus(const INT i) const;
	template<typename INT>
	uint32_t FromT56ntrlLocusToLocus(const INT i) const;
	template<typename INT>
	uint32_t FromT56selLocusToLocus(const INT i) const;
	template<typename INT>
	uint32_t FromLocusToNextT1Locus(const INT i) const;
	template<typename INT>
	uint32_t FromLocusToNextT2Locus(const INT i) const;
	template<typename INT>
	uint32_t FromLocusToNextT3Locus(const INT i) const;
	template<typename INT>
	uint32_t FromLocusToNextT4Locus(const INT i) const;
	template<typename INT>
	uint32_t FromLocusToNextT56Locus(const INT i) const;
	template<typename INT>
	uint32_t FromLocusToNextT56ntrlLocus(const INT i) const;
	template<typename INT>
	uint32_t FromLocusToNextT56selLocus(const INT i) const;
	template<typename INT>
	bool isT56neutral(const INT i) const;
	template<typename INT>
	T56LocusGender FromT56LocusToT56genderLocus(const INT i) const;
	template<typename INT>
	unsigned char getLocusType(const INT i) const;
	template<typename INT>
	LocusDescription getLocusTypeAndItsIndex(const INT i) const;



	void setT56GenderLoci(std::vector<std::vector<double>>& T56_fit, bool alreadyKnowItIsAllNtrl = false, bool alreadyKnowItIsAllSel = false);
	void readLoci(InputReader& input);
	void readT56Compression(InputReader& input);


};



