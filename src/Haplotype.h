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

    std::vector<unsigned int> T5_Alleles;          // Type 5 Just like T1 except that only 

    std::vector<double> W_T1;
    std::vector<double> W_T2;
    std::vector<double> W_T5;
    // No W_T3 as the fitness makes sense for the individual only
public:


    double getW_T1(int fitnessMapIndex);
    double getW_T2(int fitnessMapIndex);
    double getW_T5(int fitnessMapIndex);
    void setAllW_T1(double w);
    void setAllW_T2(double w);
    void setAllW_T5(double w);
    void setW_T1(double w, int fitnessMapIndex);
    void setW_T2(double w, int fitnessMapIndex);
    void setW_T5(double w, int fitnessMapIndex);


    unsigned char getT1_char(int& T1_char_index);
    bool getT1_Allele(const int T1Locus);
    bool getT1_Allele(const int char_index, const int bit_index);
    unsigned char getT2_Allele(const int char_index);
    char getT3_Allele(const int char_index);
    bool getT5_Allele(const int Locus);
    size_t getT5_nthMutation(const int n);

    void setT1_Allele(const int& char_index, const int& bit_index, const int& value);
    void setT1_AlleleToOne(int& char_index, int& bit_index);
    void setT1_AlleleToZero(int& char_index, int& bit_index);
    void setT1_char(int& T1_char_index, unsigned char& c);
    void setT1_char(int& T1_char_index, unsigned char&& c);
    void setT2_Allele(const int char_index, const unsigned char value);
    void setT3_Allele(const int char_index, const char value);
    void setT5_Allele(const int& locus, const bool& value);
    void setEntireT5_Allele(std::vector<unsigned int>& t5a);
    void setT5_AlleleToOne(int& locus);
    void setT5_AlleleToZero(int& locus);
    void setT5_AlleleToOne_JustPushBack(int& locus);

    void toggleT1_Allele(int& MutPosition);
    void toggleT1_Allele(int& MutPosition, int Habitat);
    void AddMutT2_Allele(int& char_index);
    void AddMutT2_Allele(int& char_index, int Habitat);
    void AddMutT3_Allele(int& char_index);
    void toggleT5_Allele(int& MutPosition);
    void toggleT5_Allele(int& MutPosition, int Habitat);
    

    void copyIntoT1(int from, int to, Haplotype& SourceChromo);
    void copyIntoT2(int from, int to, Haplotype& SourceChromo);
    void copyIntoT3(int from, int to, Haplotype& SourceChromo);
    void clearT5Alleles();
    void copyIntoT5(int from, int to, Haplotype& SourceChromo);
    
    void print(bool WithRecDist, std::string& prefix);
    void AssertBitSetSize(int T1_nbChars);
    Haplotype(const std::vector<unsigned char>& T1_Allel);
    Haplotype(const int patch_index,char Abiogenesis, double indHaplo_index);
    Haplotype(bool ShouldReadPopFromBinary);
    Haplotype();                             // default constructor
    Haplotype(const Haplotype& other);             // copy constructor
    Haplotype& operator=(const Haplotype& other);  // copy assignment operator
    Haplotype& operator=(Haplotype&& other); // move assignment operator
    void PrintBinaryFile(OutputFile& file);

    int getT3_AllelesSize() //For debug purposes
    {
        return T3_Alleles.size();
    }; 

    bool isFreeFromMutations();
    bool isFreeFromMutations(int T1_locusFrom, int T1_locusTo);

    std::vector<unsigned int>::const_iterator T5_AllelesCBegin();
    std::vector<unsigned int>::const_iterator T5_AllelesCEnd();
    bool isT5Mutation(unsigned int locus);
    std::vector<unsigned int>::const_iterator T5_AllelesCiterator(int locus, std::vector<unsigned int>::const_iterator from);
    std::vector<unsigned int>::const_iterator T5_AllelesCiterator(int locus);
    int T5_AllelesPosition(int locus, int from);
    int T5_AllelesPosition(int locus);
    int T5_howManyMutations();

};

