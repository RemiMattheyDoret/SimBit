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


Note for Remi of things to do:
    Test how much slower it is to have N=1e5, L=10 vs N=10, L=1e5 to estimate the cost of having Individuals not contiguous in memory

    When using several environment, fitnessMap should take migration rate into account. This migration rate can vary through time and therefore fitnessMap should be redefined for faster simulations

 */


class Individual
{
private:
    // Individuals are necessarily diploid for the moment
    Haplotype haplo0;
    Haplotype haplo1;
    
public:
    static std::vector<double> T3_IndPhenotype; // This variable beeing a static would cause trouble for multithreaded version of SimBit

    Haplotype& getHaplo(const int& haplo_index);

    std::vector<double> CalculateFitnessComponents(const int& Habitat);
    double CalculateFitness(const int& patch_index);
    double CalculateT1FitnessMultiplicity(const int& Habitat, int fitnessMapIndex, int T1_locusFrom, int T1_locusTo);
    double CalculateT1FitnessNoMultiplicity(const int& Habitat);
    double CalculateT1EpistaticFitness(const int& Habitat);
    double CalculateT2Fitness(const int& Habitat, int fitnessMapIndex, int T2_locusFrom, int T2_locusTo);
    void CalculateT3Phenotype(const int& Habitat);
    static double CalculateT3Fitness(const int& Habitat);

    void SetHaplo(int haplo_index, Haplotype& chrom);
    Individual(Haplotype& matChrom, Haplotype& patChrom);
    Individual(const int patch_index,char Abiogenesis);
    Individual(bool ShouldReadPopFromBinary);
    Individual(Haplotype& knownHaplotype, char Abiogenesis);
    Individual(const Haplotype& knownHaplotype);   // copy constructor
    void PrintBinaryFile(OutputFile& file);
    bool isFreeFromMutations();
    bool isFreeFromMutations(int T1_locusFrom, int T1_locusTo);
    //Individual(Individual&&) = default;
    //Individual& operator=(Individual&& ) = default;
};

