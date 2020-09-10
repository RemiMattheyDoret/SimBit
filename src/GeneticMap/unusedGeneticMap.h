
class GeneticMap_ChangeDescription
{
public:
	bool willIChange;
	size_t TxLocus;
};

class T56LocusGender
{
public:
	bool isNtrl;
	size_t locusInGender;

	T56LocusGender(bool b, size_t i):isNtrl(b), locusInGender(i){}

};


class GeneticMap
{
private:

	//void expandMapToChangeIndexWithRepeatedElement(std::vector<size_t>& v, size_t e, size_t n);
    //void expandMapToChangeIndexWithIncreasingElement(std::vector<size_t>& v, size_t from, size_t to);
    void expandMapToChangeIndex(std::vector<size_t>& v, std::vector<size_t>& Tx_atChange, unsigned char type, size_t change_index);
    void push_back_repeated_value(std::vector<size_t>& v, size_t value, size_t nbRepeats);

	class PreviousLocusType // This is a little bit of a useless class but it helped me write this code
    {
        private:
        unsigned char type; // some non-sense number when undefined
        
        public:

        PreviousLocusType():type(250){}

        bool isSameType(unsigned char& otherType)
        {
            return type == otherType;
        }

        void setType(unsigned char& newType){assert(type != 5);type = newType;}
    };


	bool isOnlyOneType = false; // Note btw that T56ntrl and T56sel are different
	unsigned char onlyType = 250;

    unsigned char FT1Type = 1;
    unsigned char FT2Type = 2;
    unsigned char FT3Type = 3;
    unsigned char FT4Type = 4;
    unsigned char FT56Type = 5;
    unsigned char FT56ntrlType = 50;
    unsigned char FT56selType = 51;

    unsigned char typeFormat(std::string& s, bool allowSubTypes) const;
    std::string typeFormat(unsigned char& s) const;



	// Maps to last change index
	std::vector<size_t> fromLocusToLastChangeIndex;
	std::vector<size_t> fromT1LocusToLastChangeIndex;
	std::vector<size_t> fromT2LocusToLastChangeIndex;
	std::vector<size_t> fromT3LocusToLastChangeIndex;
	std::vector<size_t> fromT4LocusToLastChangeIndex;
	std::vector<size_t> fromT56ntrlLocusToLastChangeIndex;
	std::vector<size_t> fromT56selLocusToLastChangeIndex;
	std::vector<size_t> fromT56LocusToLastChangeIndex;

	// Size of the number of changes (a change is a switch from using one typpe to using another)
	std::vector<unsigned char> whoWillIncreaseFromChange;

	// Size zero or of size of the number of changes (a change is a switch from using one typpe to using another)
	std::vector<size_t> TallLocus_atChange;
	std::vector<size_t> T1Locus_atChange;
	std::vector<size_t> T2Locus_atChange;
	std::vector<size_t> T3Locus_atChange;
	std::vector<size_t> T4Locus_atChange;
	std::vector<size_t> T56ntrlLocus_atChange;
	std::vector<size_t> T56selLocus_atChange;
	std::vector<size_t> T56Locus_atChange;



	

	template<typename INT>
	size_t getLocusFromTxLocus(INT TxLocus, const std::vector<size_t>& fromTxLocusToLastChangeIndex, const std::vector<size_t>& TxLocus_atChange) const;

	template<typename INT>
	size_t getTxLocusFromLocus(INT Locus, const std::vector<size_t>& TxLocus_atChange, const unsigned char type) const;

	void readLoci_SetNewChangeOfType(unsigned char& type, PreviousLocusType& previousLocus);

	void shrink_to_fit();

	void print();
	template<typename T>
	void print(std::string name, std::vector<T>& v);
	template<typename T>
	void print(std::string name, T& v);
	



public:

	size_t T1_nbLoci;
	size_t T1_nbChars;
	size_t T1_nbLociLastByte;
	size_t T2_nbLoci;
	size_t T3_nbLoci;
	size_t T4_nbLoci;
	size_t T5_nbLoci;
	size_t T5sel_nbLoci;
	size_t T5ntrl_nbLoci;
	size_t T6_nbLoci;
	size_t T6sel_nbLoci;
	size_t T6ntrl_nbLoci;
	size_t T56_nbLoci;
	size_t T56sel_nbLoci;
	size_t T56ntrl_nbLoci;
	size_t TotalNbLoci;

	bool isT5selCompress;
	bool isT5ntrlCompress;



	GeneticMap():
    isOnlyOneType(false),
    onlyType(250),
    FT1Type (1),
	FT2Type (2),
	FT3Type (3),
	FT4Type (4),
	FT56Type (5),
	FT56ntrlType (50),
	FT56selType (51)
    {}



	void readLoci(InputReader& input, std::vector<bool>& isT56LocusUnderSelection, bool isT5ntrlCompress, bool isT5selCompress);
	
	template<typename INT>
	T56LocusGender getT56GLocusFromT56Locus(INT T56Locus) const;

	template<typename INT>
	bool isT56neutral(INT T56Locus) const;
	

	template<typename INT>
	size_t getLocusFromT1Locus(INT T1Locus) const;

	template<typename INT>
	size_t getLocusFromT2Locus(INT T2Locus) const;

	template<typename INT>
	size_t getLocusFromT3Locus(INT T3Locus) const;

	template<typename INT>
	size_t getLocusFromT4Locus(INT T4Locus) const;

	template<typename INT>
	size_t getLocusFromT56ntrlLocus(INT T5ntrlLocus) const;

	template<typename INT>
	size_t getLocusFromT56selLocus(INT T5selLocus) const;

	template<typename INT>
	size_t getLocusFromT56Locus(INT T5Locus) const;

	template<typename INT>
	size_t getT1LocusFromLocus(INT locus) const;

	template<typename INT>
	size_t getT2LocusFromLocus(INT locus) const;	

	template<typename INT>
	size_t getT3LocusFromLocus(INT locus) const;
	
	template<typename INT>
	size_t getT4LocusFromLocus(INT locus) const;

	template<typename INT>
	size_t getT56ntrlLocusFromLocus(INT locus) const;
	
	template<typename INT>
	size_t getT56selLocusFromLocus(INT locus) const;

	template<typename INT>
	unsigned char getLocusSubtypeFromLocus(INT locus) const;

	//template<typename INT> Not implemented!
	//size_t getT56LocusFromLocus(INT locus) const;
};



