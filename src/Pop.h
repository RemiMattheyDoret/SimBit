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



class Pop
{
private:
    std::vector<Patch> patches;
    static std::vector<int>  T2LociToCorrect;
    std::vector<int>         indexFirstMale;
    std::vector<std::vector<Walker>>      walkers; // Walkers[patch_index][sex]


    std::vector<int> findWhatMustBeToggledAndUpdateFlipped(const std::vector<unsigned>& popFreqs, std::vector<uint32_t>& flipped, const int& nbLoci,  const double freqThreshold);
    void toggleT56FixedMutations();
    void checkIfCumSumFitsIsNotTooSmall(int patch_index);

    
public:
    std::vector<std::vector<std::vector<double>>> CumSumFits; // each patch, each gender, each individual fitness (cumulated)

    Patch& getPatch(const int& patch_index);
    int getNbPatches();
    void AddPatch(Patch& newPatch);
    void AddPatch();
    void RemoveLastPatch();
    int SelectionParent(int patch_from, int sex);
    void CalculateFitnesses();
    //double CalculateFitnessForNextGeneration(Individual& Offspring, int patch_index, int ind_index);
    //void prepareNextGenerationAndIndexFirstMale(const int patch_index, const std::vector<int>& patchSizeNextGeneration);
    
    int SelectionOriginPatch(uint32_t patch_to, double rndIfNotStochastic = -1.0);
    int patchSizeNextGeneration(int patch_index);
    void toggleT56MutationsIfNeeded();
    void toggleT56LociFromEveryone(std::vector<int>& T5ntrlLociToToggle, std::vector<int>& T5selLociToToggle);
    
    Pop(bool ShouldReadPopFromBinary);
    /*Pop();
    Pop(const Pop&& p);
    Pop(const Pop& p);
    Pop operator=(const Pop&& p);
    Pop operator=(const Pop& p);*/

    void PrintBinaryFile();

    static void addT2LocusToCorrect(int T2Locus);
    int correctT2Loci();

    static void updatePops(Pop& pop1, Pop& pop2, int speciesIndex, int oldNbPatches, std::vector<int>& previousPatchSizes);

    std::vector<unsigned> computeT56ntrlFrequencies();
    std::vector<unsigned> computeT56selFrequencies();
    std::vector<unsigned> computeT56Frequencies();
    std::vector<double> computeT56RelativeFrequencies();
    std::vector<std::vector<double>> computeT1RelativeFrequenciesPerPatch();
    std::vector<std::vector<unsigned>> computePatchSpecificT56ntrlFrequencies();
    std::vector<std::vector<unsigned>> computePatchSpecificT56selFrequencies();
    std::vector<std::vector<unsigned>> computePatchSpecificT56Frequencies();
    void freeMemory();
    void freeT56Memory(); // reset all T56 to zero
    void shrink_to_fitT56();

    std::vector<T1_locusDescription> listT1PolymorphicLoci();

};

