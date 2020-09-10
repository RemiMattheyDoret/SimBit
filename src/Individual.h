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


class Individual
{
private:
    // Individuals are necessarily diploid for the moment
    Haplotype haplo0;
    Haplotype haplo1;
    static std::vector<double> fitnessComponents;
    
public:
    static std::vector<double> T3_IndPhenotype; // This variable beeing a static would cause trouble for multithreaded version of SimBit

    Haplotype& getHaplo(const int& haplo_index);

    const std::vector<double>& CalculateFitnessComponents(const int& Habitat);
    
    double CalculateFitness(const int& patch_index);
    double CalculateT1FitnessNoMultiplicity(const int& Habitat);
    double CalculateT1EpistaticFitness(const int& Habitat);
    double CalculateT2Fitness(const int& Habitat, int fitnessMapIndex, int T2_locusFrom, int T2_locusTo);
    void CalculateT3Phenotype(const int& Habitat);
    static double CalculateT3Fitness(const int& Habitat);
    

    template<typename Iterator>
    double CalculateT56FitnessNoMultiplicity(const int& Habitat, Iterator itHaplo0, Iterator itHaplo1, Iterator itHaplo0End, Iterator itHaplo1End);

    template<typename Iterator>
    double CalculateT56FitnessNoMultiplicityOnSubsetOfLoci(const int& Habitat, const std::vector<int>& LociSet, Iterator itHaplo0, Iterator itHaplo1, Iterator itHaplo0End, Iterator itHaplo1End);


    static bool isLocusIsInSet(const int locus, const std::vector<int>& LociSet);
    std::vector<double> CalculateFitnessComponentsOnSubsetOfLoci(const int& Habitat, const int lociSetIndex);
    double CalculateT1FitnessNoMultiplicityOnSubsetOfLoci(const int& Habitat, const std::vector<int>& LociSet);
    double CalculateT1EpistaticFitnessOnSubsetOfLoci(const int& Habitat, const std::vector<int>& LociSet);
    double CalculateT2FitnessOnSubsetOfLoci(const int& Habitat, const std::vector<int>& LociSet);
    void CalculateT3PhenotypeOnSubsetOfLoci(const int& Habitat, const std::vector<int>& LociSet);

    void SetHaplo(int haplo_index, Haplotype& chrom);
    Individual(Haplotype& h0, Haplotype& h1);
    Individual(const int patch_index,char Abiogenesis, int ind_index);
    Individual(bool ShouldReadPopFromBinary);
    Individual(Haplotype& knownHaplotype, char Abiogenesis);
    Individual(const Haplotype& knownHaplotype);   // copy constructor
    Individual();
    /*Individual(const Individual& I);
    Individual(const Individual&& I);
    Individual operator=(const Individual& I);
    Individual operator=(const Individual&& I);*/

    void PrintBinaryFile(OutputFile& file);
    bool isFreeFromMutations();
    bool isFreeFromMutations(int T1_locusFrom, int T1_locusTo);
    //Individual(Individual&&) = default;
    //Individual& operator=(Individual&& ) = default;

    void toggleT56LociFromHaplotypes(std::vector<int>& T5ntrlLociToToggle, std::vector<int>& T5selLociToToggle, int Habitat);

    void freeT56Memory();
};
