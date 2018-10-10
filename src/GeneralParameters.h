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


class GeneralParameters // is a singleton
{
    // attribute starting with '__' means that they store info about sequential parameters. Attributes being used have the same name without the '__' and the are reference of a single "line" of the element starting with '__'. '__GenerationChange' indicate when the parameters should be updated. All demographic parameters must have an entry for the same generations.
public:
    // patchSize
    std::vector<std::vector<int>>                   allSpeciesPatchSizes;  // allSpeciesPatchSizes[patch_index][speciesIndex];
    std::vector<std::string>                        speciesNames;

    std::vector<int>                                __GenerationChange;
    
    // Technic
    std::mt19937                                    mt;
    int                                             random_seed;
    int                                             nbThreads;
    int                                             OverwriteMode;
    bool                                            DryRun;


    // BinaryFileToRead The path will be species specific though and are stored in SSP
    std::ifstream                                   BinaryFileToRead;

    // sequencing error rate
    double                                          sequencingErrorRate;
    
    // nbGenerations
    int                                             nbGenerations;
    int                                             startAtGeneration;
    int                                             CurrentGeneration;

    // SpeciesStuff
    int                                             nbSpecies;
    std::vector<std::vector<double>>                speciesEcoRel_effect;
    std::vector<std::vector<char>>                  speciesEcoRel_type;

    // nbPatches
    int                                             PatchNumber;
    int                                             maxEverPatchNumber;
    std::vector<int>                                __PatchNumber;    
    
    // Random Distributions. If change number of distributions, don't forget to change the incremenet of 'nbParamSet' in 'setRandomDistributions'
    std::uniform_real_distribution<double>           random_0and1;
    std::uniform_int_distribution<int>               random_0or1;

    // Manage outputs
    int                                             nbWarningsSentFrom_Haplotype_AddMutT2_Allele = 0;
    std::vector<int>                                output_FST_nbPatchesToConsider;

    // methods
    void readTemporalChanges(InputReader input);
    void readPatchNumber(InputReader input);
    void readT1_FST_info(InputReader input);
    void readSeed(InputReader input);
    void saveSSPPatchSize_toGP();
    //void saveSSPPatchSize_toGP_lowerSecurity();
    //void UpdateParametersallSpeciesPatchSizes();
    void UpdateParametersPatchNumber(int generation_index);
    
    void initializeAllSpeciesPatchSizes();
};
