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
private:
    std::vector<unsigned char> T1_Alleles;   // Type 1 Each bit is a binary site
    std::vector<unsigned char> T2_Alleles;   // Type 2 Each byte is a locus for which number of mutations are counted.
    std::vector<char> T3_Alleles;            // Type 3 Each byte is the allelic effect along a single phenotypic dimension

    // T4 types are tracked with T4Tree

    std::vector<unsigned int> T5sel_Alleles;          // Type 5 SLimM style. Just like T1 except that only mutations are keeping tracked of
    std::vector<unsigned int> T5ntrl_Alleles;

    CompressedSortedDeque T6sel_Alleles;       // Tyoe 6 is like Type 5 except that it tries to reduce RAM by separating values into prefix and suffix
    CompressedSortedDeque T6ntrl_Alleles;

    std::vector<double> W_T1;
    std::vector<double> W_T2;
    std::vector<double> W_T56;
    // No W_T3 as the fitness makes sense for the individual only

    
    std::vector<unsigned int> getT56ntrlTrueHaplotype();
    std::vector<unsigned int> getT56selTrueHaplotype();
public:


    ZipIterator<std::vector<unsigned int>, std::vector<unsigned int>::iterator> T5ntrl_AllelesBegin();
    ZipIterator<std::vector<unsigned int>, std::vector<unsigned int>::iterator> T5ntrl_AllelesEnd();
    std::vector<unsigned int>::iterator T5sel_AllelesBegin();
    std::vector<unsigned int>::iterator T5sel_AllelesEnd();


    ZipIterator<CompressedSortedDeque, CompressedSortedDeque::iterator> T6ntrl_AllelesBegin();
    ZipIterator<CompressedSortedDeque, CompressedSortedDeque::iterator> T6ntrl_AllelesEnd();
    CompressedSortedDeque::iterator T6sel_AllelesBegin();
    CompressedSortedDeque::iterator T6sel_AllelesEnd();


    ZipIterator<std::vector<unsigned int>, std::vector<unsigned int>::iterator> T56ntrl_AllelesBegin(ZipIterator<std::vector<unsigned int>, std::vector<unsigned int>::iterator> nothing);
    ZipIterator<std::vector<unsigned int>, std::vector<unsigned int>::iterator> T56ntrl_AllelesEnd(ZipIterator<std::vector<unsigned int>, std::vector<unsigned int>::iterator> nothing);
    std::vector<unsigned int>::iterator T56sel_AllelesBegin(std::vector<unsigned int>::iterator& nothing);
    std::vector<unsigned int>::iterator T56sel_AllelesEnd(std::vector<unsigned int>::iterator& nothing);

    ZipIterator<CompressedSortedDeque, CompressedSortedDeque::iterator> T56ntrl_AllelesBegin(ZipIterator<CompressedSortedDeque, CompressedSortedDeque::iterator> nothing);
    ZipIterator<CompressedSortedDeque, CompressedSortedDeque::iterator> T56ntrl_AllelesEnd(ZipIterator<CompressedSortedDeque, CompressedSortedDeque::iterator> nothing);
    CompressedSortedDeque::iterator T56sel_AllelesBegin(CompressedSortedDeque::iterator& nothing);
    CompressedSortedDeque::iterator T56sel_AllelesEnd(CompressedSortedDeque::iterator& nothing);

    ZipIterator<std::vector<unsigned int>, std::vector<unsigned int>::iterator> T56ntrl_AllelesIterator(ZipIterator<std::vector<unsigned int>, std::vector<unsigned int>::iterator> nothing, unsigned value);
    std::vector<unsigned int>::iterator T56sel_AllelesIterator(std::vector<unsigned int>::iterator& nothing, unsigned value);

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


    unsigned char getT1_char(int& T1_char_index);
    bool getT1_Allele(const int T1Locus);
    bool getT1_Allele(const int char_index, const int bit_index);
    unsigned char getT2_Allele(const int char_index);
    char getT3_Allele(const int char_index);
    //bool getT5ntrl_Allele(const int Locus);
    //bool getT5sel_Allele(const int Locus);
    //int getT5ntrl_nthMutation(const int n);
    //int getT5sel_nthMutation(const int n);

    void setT1_Allele(const int& char_index, const int& bit_index, const int& value);
    void setT1_AlleleToOne(int& char_index, int& bit_index);
    void setT1_AlleleToZero(int& char_index, int& bit_index);
    void setT1_char(int& T1_char_index, unsigned char& c);
    void setT1_char(int& T1_char_index, unsigned char&& c);
    void setT2_Allele(const int char_index, const unsigned char value);
    void setT3_Allele(const int char_index, const char value);
    
    void setEntireT5_Allele(std::vector<unsigned int>& t5a);
    void setT5_Allele(const int& locus, const bool& value);
    void setT5ntrl_Allele(const int& locus, const bool& value);
    void setT5sel_Allele(const int& locus, const bool& value);

    void setT5sel_AlleleToOne(int& locus);
    void setT5sel_AlleleToZero(int& locus);
    void setT5sel_AlleleToOne_JustPushBack(unsigned& locus);
    void setT5ntrl_AlleleToOne(int& locus);
    void setT5ntrl_AlleleToZero(int& locus);
    void setT5ntrl_AlleleToOne_JustPushBack(unsigned& locus);

    void setEntireT6_Allele(std::vector<unsigned int>& t6a);
    void setT6_Allele(const int& locus, const bool& value);
    void setT6ntrl_Allele(const int& locus, const bool& value);
    void setT6sel_Allele(const int& locus, const bool& value);

    void setT6sel_AlleleToOne(int& locus);
    void setT6sel_AlleleToZero(int& locus);
    void setT6sel_AlleleToOne_JustPushBack(unsigned& locus);
    void setT6ntrl_AlleleToOne(int& locus);
    void setT6ntrl_AlleleToZero(int& locus);
    void setT6ntrl_AlleleToOne_JustPushBack(unsigned& locus);

    void mutateT1_Allele(int& MutPosition, int& Habitat);
    void toggleT1_Allele(int& byte_index, int& bit_index);
    void AddMutT2_Allele(int& char_index);
    void AddMutT2_Allele(int& char_index, int Habitat);
    void AddMutT3_Allele(int& char_index);
    
    template<typename INT>
    void mutateT56ntrl_Allele(INT MutPosition);
    template<typename INT>
    void mutateT56sel_Allele(INT MutPosition, int Habitat);


    template<typename INT>
    void mutateT5ntrl_Allele(std::vector<INT>& MutPositions);
    template<typename INT>
    void mutateT5ntrl_Allele(INT MutPosition);
    template<typename INT>
    void mutateT5sel_Allele(INT MutPosition, int Habitat);
    template<typename INT>
    void toggleT5ntrl_Allele(INT MutPosition);
    template<typename INT>
    void toggleT5sel_Allele(INT MutPosition);


    template<typename INT>
    void mutateT6ntrl_Allele(INT MutPosition);
    template<typename INT>
    void mutateT6sel_Allele(INT MutPosition, int Habitat);
    template<typename INT>
    void toggleT6ntrl_Allele(INT MutPosition);
    template<typename INT>
    void toggleT6sel_Allele(INT MutPosition);
    

    void copyIntoT1(int from, int to, Haplotype& SourceChromo);
    void copyIntoT2(int from, int to, Haplotype& SourceChromo);
    void copyIntoT3(int from, int to, Haplotype& SourceChromo);
    void clearT56Alleles();
    void copyIntoT56ntrl(int from, int to, Haplotype& SourceChromo);
    void copyIntoT56sel(int from, int to, Haplotype& SourceChromo);
    
    void print(bool WithRecDist, std::string& prefix);
    void AssertBitSetSize(int T1_nbChars);
    Haplotype(const std::vector<unsigned char>& T1_Allel);
    Haplotype(const int patch_index,char Abiogenesis, int indHaplo_index);
    Haplotype(bool ShouldReadPopFromBinary);
    Haplotype();                              // default constructor
    Haplotype(const Haplotype& other);             // copy constructor
    Haplotype& operator=(const Haplotype& other);  // copy assignment operator
    Haplotype& operator=(Haplotype&& other); // move assignment operator
    void swap(Haplotype& other);
    void PrintBinaryFile(OutputFile& file);

    int getT3_AllelesSize() //For debug purposes
    {
        return T3_Alleles.size();
    }; 

    bool isFreeFromMutations();
    bool isFreeFromMutations(int T1_locusFrom, int T1_locusTo);
/*
    std::vector<unsigned int>::const_iterator T5sel_AllelesCBegin();
    std::vector<unsigned int>::const_iterator T5sel_AllelesCEnd();
    std::vector<unsigned int>::const_iterator T5sel_AllelesCiterator(int locus, std::vector<unsigned int>::const_iterator from);
    std::vector<unsigned int>::const_iterator T5sel_AllelesCiterator(int locus);


    std::vector<unsigned int>::const_iterator T5ntrl_AllelesCBegin();
    std::vector<unsigned int>::const_iterator T5ntrl_AllelesCEnd();
    std::vector<unsigned int>::const_iterator T5ntrl_AllelesCiterator(int locus, std::vector<unsigned int>::const_iterator from);
    std::vector<unsigned int>::const_iterator T5ntrl_AllelesCiterator(int locus);*/

    //int T5_AllelesPosition(int locus, int from);
    //int T5_AllelesPosition(int locus);
    //int T5_howManyMutations();
    void toggleT56ntrlLoci(std::vector<int>& lociToToggle);
    //void toggleT56selLoci(std::vector<int>& lociToToggle, int Habitat);
    void toggleFromT5ntrl_Allele(int& MutPosition, unsigned int& from);
    void toggleFromT6ntrl_Allele(int& MutPosition, unsigned int& blockIndexFrom, unsigned int& fromInBlock);
    void toggleFromT5sel_Allele(int& MutPosition, unsigned int& from, int Habitat);
    void toggleFromT6sel_Allele(int& MutPosition, unsigned int& blockIndexFrom, unsigned int& fromInBlock, int Habitat);

    int getNbT5ntrl();
    int getNbT6ntrl();

    void printT5sel_Alleles();
    void printT5ntrl_Alleles();
    bool isT5ntrlMutation(int locus);
    bool isT5selMutation(int locus);

    void printT6sel_Alleles();
    void printT6ntrl_Alleles();
    bool isT6ntrlMutation(int locus);
    bool isT6selMutation(int locus);

    void assertT5orderAndUniqueness();
    template<typename INT>
    void updateFitnessAfterT56Mutation(INT MutPosition, bool isNowFoundInAlleles, int Habitat);
    template<typename ITERATOR, typename CONTAINER, typename INT>
    bool T56finishMutation(ITERATOR& haploP, CONTAINER& container, INT MutPosition);
    size_t nbT56muts();
    size_t nbT56muts(int fitnessMapIndex);


    // Multiplicity fitness calculator
    
    double CalculateT1FitnessMultiplicity(const int& Habitat);

    double CalculateT2Fitness(const int& Habitat);
    
    template<typename Iterator>
    double CalculateT56FitnessMultiplicity(const int& Habitat, Iterator it, Iterator itEnd);
    double CalculateT56FitnessMultiplicity(const int& Habitat);

    double CalculateT1FitnessMultiplicityOnSubsetOfLoci(const int& Habitat, const std::vector<int>& LociSet);
    double CalculateT2FitnessOnSubsetOfLoci(const int& Habitat, const std::vector<int>& LociSet);
    template<typename Iterator>
    double CalculateT56FitnessMultiplicityOnSubsetOfLoci(const int& Habitat, const std::vector<int>& LociSet, Iterator it, Iterator itEnd);
    double CalculateT56FitnessMultiplicityOnSubsetOfLoci(const int& Habitat, const std::vector<int>& LociSet);

};

