/*

 Author: Remi Matthey-Doret

    MIT License

    Copyright (c) 2017 Remi Matthey-Doret

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

 */


class Haplotype // a Haplotype is a whole halotype but can correspond to several independent chromosomes.
{
friend class T56_memoryManager;
private:
    std::vector<unsigned char> T1_Alleles;   // Type 1 Each bit is a binary site
    std::vector<unsigned char> T2_Alleles;   // Type 2 Each byte is a locus for which number of mutations are counted.
    std::vector<int16_t> T3_Alleles;          // Type 3 QTLs

    // T4 types are tracked with T4Tree

    std::vector<uint32_t> T5sel_Alleles;          // Type 5 SLimM style. Just like T1 except that only mutations are keeping tracked of
    std::vector<uint32_t> T5ntrl_Alleles;

    CompressedSortedDeque T6sel_Alleles;       // Type 6 is like Type 5 except that it tries to reduce RAM by separating values into prefix and suffix
    CompressedSortedDeque T6ntrl_Alleles;

    std::vector<T7Gene> T7_Alleles; // Type 7 use for ENTWINE type of development and produce a phenotype on the same space than T3

    std::vector<double> W_T1;
    std::vector<double> W_T2;
    std::vector<double> W_T56;
    // No W_T3 as the fitness makes sense for the individual only

    
    std::vector<uint32_t> getT56ntrlTrueHaplotype();
    std::vector<uint32_t> getT56selTrueHaplotype();

public:
    // Used for T4Tree
    // ID is a little bit of a weird attribute because Haplotype does not manage its value but T4Tree and LifeCycle do. It will remain uninitialized unless T4Tree wants to do something with it
    uint32_t T4ID;


    ZipIterator<std::vector<uint32_t>, std::vector<uint32_t>::iterator> T5ntrl_AllelesBegin();
    ZipIterator<std::vector<uint32_t>, std::vector<uint32_t>::iterator> T5ntrl_AllelesEnd();
    std::vector<uint32_t>::iterator T5sel_AllelesBegin();
    std::vector<uint32_t>::iterator T5sel_AllelesEnd();


    ZipIterator<CompressedSortedDeque, CompressedSortedDeque::iterator> T6ntrl_AllelesBegin();
    ZipIterator<CompressedSortedDeque, CompressedSortedDeque::iterator> T6ntrl_AllelesEnd();
    CompressedSortedDeque::iterator T6sel_AllelesBegin();
    CompressedSortedDeque::iterator T6sel_AllelesEnd();


    ZipIterator<std::vector<uint32_t>, std::vector<uint32_t>::iterator> T56ntrl_AllelesBegin(ZipIterator<std::vector<uint32_t>, std::vector<uint32_t>::iterator> nothing);
    ZipIterator<std::vector<uint32_t>, std::vector<uint32_t>::iterator> T56ntrl_AllelesEnd(ZipIterator<std::vector<uint32_t>, std::vector<uint32_t>::iterator> nothing);
    std::vector<uint32_t>::iterator T56sel_AllelesBegin(std::vector<uint32_t>::iterator& nothing);
    std::vector<uint32_t>::iterator T56sel_AllelesEnd(std::vector<uint32_t>::iterator& nothing);

    ZipIterator<CompressedSortedDeque, CompressedSortedDeque::iterator> T56ntrl_AllelesBegin(ZipIterator<CompressedSortedDeque, CompressedSortedDeque::iterator> nothing);
    ZipIterator<CompressedSortedDeque, CompressedSortedDeque::iterator> T56ntrl_AllelesEnd(ZipIterator<CompressedSortedDeque, CompressedSortedDeque::iterator> nothing);
    CompressedSortedDeque::iterator T56sel_AllelesBegin(CompressedSortedDeque::iterator& nothing);
    CompressedSortedDeque::iterator T56sel_AllelesEnd(CompressedSortedDeque::iterator& nothing);

    ZipIterator<std::vector<uint32_t>, std::vector<uint32_t>::iterator> T56ntrl_AllelesIterator(ZipIterator<std::vector<uint32_t>, std::vector<uint32_t>::iterator> nothing, unsigned value);
    std::vector<uint32_t>::iterator T56sel_AllelesIterator(std::vector<uint32_t>::iterator& nothing, unsigned value);

    ZipIterator<CompressedSortedDeque, CompressedSortedDeque::iterator> T56ntrl_AllelesIterator(ZipIterator<CompressedSortedDeque, CompressedSortedDeque::iterator> nothing, unsigned value);
    CompressedSortedDeque::iterator T56sel_AllelesIterator(CompressedSortedDeque::iterator& nothing, unsigned value);


    double getW_T1(int fitnessMapIndex);
    double getW_T2(int fitnessMapIndex);
    double getW_T56(int fitnessMapIndex);
    void setAllW_T1(double w);
    void setAllW_T2(double w);
    void setAllW_T56(double w);
    void setW_T1(double w, int fitnessMapIndex);
    void setW_T2(double w, int fitnessMapIndex);
    void setW_T56(double w, int fitnessMapIndex);
    size_t getW_T56_size();

    template<typename INT = uint32_t> void removeT7Gene(INT index);
    template<typename INT = uint32_t> void duplicateT7Gene(INT index);
    void clearT7Genes();

    template<typename INT = uint32_t>
    unsigned char getT1_char(INT T1_char_index);
    template<typename INT = uint32_t>
    bool getT1_Allele(const INT T1Locus);
    template<typename INT1 = uint32_t, typename INT2 = uint32_t>
    bool getT1_Allele(const INT1 char_index, const INT2 bit_index);
    template<typename INT = uint32_t>
    unsigned char getT2_Allele(const INT char_index);
    template<typename INT = uint32_t>
    int16_t getT3_Allele(const INT index);
    template<typename INT = uint32_t>
    T7Gene& getT7_Allele(const INT index);
    size_t nbT7Genes();
    //bool getT56_Allele(const int Locus);
    //bool getT5ntrl_Allele(const int Locus);
    //bool getT5sel_Allele(const int Locus);
    //int getT5ntrl_nthMutation(const int n);
    //int getT5sel_nthMutation(const int n);

    template<typename valueType, typename INT = uint32_t>
    void setT1_Allele(const INT char_index, const INT bit_index, const valueType& value);
    template<typename INT = uint32_t>
    void setT1_AlleleToOne(INT char_index, INT bit_index);
    template<typename INT = uint32_t>
    void setT1_AlleleToZero(INT char_index, INT bit_index);
    template<typename INT = uint32_t>
    void setT1_char(INT T1_char_index, unsigned char& c);
    template<typename INT = uint32_t>
    void setT1_char(INT T1_char_index, unsigned char&& c);
    template<typename INT = uint32_t>
    void setT2_Allele(const INT char_index, const unsigned char value);
    template<typename INT = uint32_t>
    void setT3_Allele(const INT index, const double value);
    
    template<typename INT = uint32_t>
    void setEntireT5_Allele(std::vector<INT>& t5a);
    template<typename INT = uint32_t>
    void setT5_Allele(const INT locus, const bool& value);
    template<typename INT = uint32_t>
    void setT5ntrl_Allele(const INT locus, const bool& value);
    template<typename INT = uint32_t>
    void setT5sel_Allele(const INT locus, const bool& value);

    template<typename INT = uint32_t>
    void setT5sel_AlleleToOne(INT locus);
    template<typename INT = uint32_t>
    void setT5sel_AlleleToZero(INT locus);
    template<typename INT = uint32_t>
    void setT5sel_AlleleToOne_JustPushBack(INT locus);
    template<typename INT = uint32_t>
    void setT5ntrl_AlleleToOne(INT locus);
    template<typename INT = uint32_t>
    void setT5ntrl_AlleleToZero(INT locus);
    template<typename INT = uint32_t>
    void setT5ntrl_AlleleToOne_JustPushBack(INT locus);

    template<typename INT = uint32_t>
    void setEntireT6_Allele(std::vector<INT>& t6a);
    template<typename INT = uint32_t>
    void setT6_Allele(const INT locus, const bool& value);
    template<typename INT = uint32_t>
    void setT6ntrl_Allele(const INT locus, const bool& value);
    template<typename INT = uint32_t>
    void setT6sel_Allele(const INT locus, const bool& value);

    template<typename INT = uint32_t>
    void setT6sel_AlleleToOne(INT locus);
    template<typename INT = uint32_t>
    void setT6sel_AlleleToZero(INT locus);
    template<typename INT = uint32_t>
    void setT6sel_AlleleToOne_JustPushBack(INT locus);
    template<typename INT = uint32_t>
    void setT6ntrl_AlleleToOne(INT locus);
    template<typename INT = uint32_t>
    void setT6ntrl_AlleleToZero(INT locus);
    template<typename INT = uint32_t>
    void setT6ntrl_AlleleToOne_JustPushBack(INT locus);

    template<typename INT = uint32_t>
    void mutateT1_Allele(INT MutPosition, int& Habitat);
    template<typename INT = uint32_t>
    void toggleT1_Allele(INT byte_index, INT bit_index);
    template<typename INT = uint32_t>
    void AddMutT2_Allele(INT char_index);
    template<typename INT = uint32_t>
    void AddMutT2_Allele(INT char_index, int Habitat);
    template<typename INT = uint32_t>
    void mutateT3_Allele(INT index);
    
    template<typename INT = uint32_t>
    void mutateT56ntrl_Allele(INT MutPosition);
    template<typename INT = uint32_t>
    void mutateT56sel_Allele(INT MutPosition, int Habitat);


    template<typename INT = uint32_t>
    void mutateT5ntrl_Allele(std::vector<INT>& MutPositions);
    template<typename INT = uint32_t>
    void mutateT5ntrl_Allele(INT MutPosition);
    template<typename INT = uint32_t>
    void mutateT5sel_Allele(INT MutPosition, int Habitat);
    template<typename INT = uint32_t>
    void toggleT5ntrl_Allele(INT MutPosition);
    template<typename INT = uint32_t>
    void toggleT5sel_Allele(INT MutPosition);


    template<typename INT = uint32_t>
    void mutateT6ntrl_Allele(INT MutPosition);
    template<typename INT = uint32_t>
    void mutateT6sel_Allele(INT MutPosition, int Habitat);
    template<typename INT = uint32_t>
    void toggleT6ntrl_Allele(INT MutPosition);
    template<typename INT = uint32_t>
    void toggleT6sel_Allele(INT MutPosition);
    

    template<typename INT = uint32_t>
    void copyIntoT1(INT from, INT to, Haplotype& SourceChromo);
    template<typename INT = uint32_t>
    void copyIntoT2(INT from, INT to, Haplotype& SourceChromo);
    template<typename INT = uint32_t>
    void copyIntoT3(INT from, INT to, Haplotype& SourceChromo);
    void clearT56Alleles();
    template<typename INT = uint32_t>
    void copyIntoT56ntrl(INT from, INT to, Haplotype& SourceChromo);
    template<typename INT = uint32_t>
    void copyIntoT56sel(INT from, INT to, Haplotype& SourceChromo);
    
    void print(bool WithRecDist, std::string& prefix);
    void AssertBitSetSize(int T1_nbChars);
    Haplotype(std::vector<unsigned char> T1_info, std::vector<unsigned char> T2_info, std::vector<int16_t> T3_info, uint32_t t4id, std::vector<uint32_t> T56_info);
    Haplotype(const std::vector<unsigned char>& T1_Allel);
    Haplotype(const int patch_index,char Abiogenesis, int indHaplo_index);
    Haplotype(bool ShouldReadPopFromBinary);
    Haplotype();                              // default constructor
    //Haplotype( Haplotype&& other);             // MOVE constructor
    Haplotype(const Haplotype& other);             // copy constructor
    Haplotype& operator=(const Haplotype& other);  // copy assignment operator
    //Haplotype& operator=(Haplotype&& other); // move assignment operator
    void swap(Haplotype& other);
    void PrintBinaryFile(OutputFile& file);

    int getT3_AllelesSize() //For debug purposes
    {
        return T3_Alleles.size();
    }; 

    bool isFreeFromMutations();
    bool isFreeFromMutations(int T1_locusFrom, int T1_locusTo);
/*
    std::vector<uint32_t>::const_iterator T5sel_AllelesCBegin();
    std::vector<uint32_t>::const_iterator T5sel_AllelesCEnd();
    std::vector<uint32_t>::const_iterator T5sel_AllelesCiterator(int locus, std::vector<uint32_t>::const_iterator from);
    std::vector<uint32_t>::const_iterator T5sel_AllelesCiterator(int locus);


    std::vector<uint32_t>::const_iterator T5ntrl_AllelesCBegin();
    std::vector<uint32_t>::const_iterator T5ntrl_AllelesCEnd();
    std::vector<uint32_t>::const_iterator T5ntrl_AllelesCiterator(int locus, std::vector<uint32_t>::const_iterator from);
    std::vector<uint32_t>::const_iterator T5ntrl_AllelesCiterator(int locus);*/

    //int T5_AllelesPosition(int locus, int from);
    //int T5_AllelesPosition(int locus);
    //int T5_howManyMutations();
    void toggleT56ntrlLoci(std::vector<int>& lociToToggle);
    //void toggleT56selLoci(std::vector<int>& lociToToggle, int Habitat);
    void toggleFromT5ntrl_Allele(int& MutPosition, uint32_t& from);
    void toggleFromT6ntrl_Allele(int& MutPosition, uint32_t& blockIndexFrom, uint32_t& fromInBlock);
    void toggleFromT5sel_Allele(int& MutPosition, uint32_t& from, int Habitat);
    void toggleFromT6sel_Allele(int& MutPosition, uint32_t& blockIndexFrom, uint32_t& fromInBlock, int Habitat);

    int getNbT5ntrl();
    int getNbT6ntrl();
    int getNbT5sel();
    int getNbT6sel();

    void printT5sel_Alleles();
    void printT5ntrl_Alleles();
    bool isT5ntrlMutation(int locus);
    bool isT5selMutation(int locus);

    void printT6sel_Alleles();
    void printT6ntrl_Alleles();
    bool isT6ntrlMutation(int locus);
    bool isT6selMutation(int locus);

    void assertT5orderAndUniqueness();
    template<typename INT = uint32_t>
    void updateFitnessAfterT56Mutation(INT MutPosition, bool isNowFoundInAlleles, int Habitat);
    template<typename ITERATOR, typename CONTAINER, typename INT = uint32_t>
    char T56finishMutation(ITERATOR& haploP, CONTAINER& container, INT MutPosition);
    uint32_t nbT56muts();
    uint32_t nbT56muts(int fitnessMapIndex);


    // Multiplicity fitness calculator
    
    double CalculateT1FitnessMultiplicity(const int& Habitat);

    double CalculateT2Fitness(const int& Habitat);
    
    double CalculateT56FitnessMultiplicity(const int& Habitat);
    double CalculateT5FitnessMultiplicity(const int& Habitat);
    double CalculateT6FitnessMultiplicity(const int& Habitat);

    double CalculateT1FitnessMultiplicityOnSubsetOfLoci(const int& Habitat, const std::vector<int>& LociSet);
    double CalculateT2FitnessOnSubsetOfLoci(const int& Habitat, const std::vector<int>& LociSet);
    template<typename Iterator>
    double CalculateT56FitnessMultiplicityOnSubsetOfLoci(const int& Habitat, const std::vector<int>& LociSet, Iterator it, Iterator itEnd);
    double CalculateT56FitnessMultiplicityOnSubsetOfLoci(const int& Habitat, const std::vector<int>& LociSet);
    void freeT56Memory();
    void shrink_to_fitT56();
};

