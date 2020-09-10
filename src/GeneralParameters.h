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

class oneSpInteraction
{
public:
    char type; 
    // type == 'A': just magnitude
    // type == 'B': magnitude multiplied by causal species patch size
    // type == 'C': magnitude multiplied by recipient species patch size
    // type == 'D': magnitude multiplied by both species patch size
    double magnitude;

    oneSpInteraction():type('A'), magnitude(0.0){}
    oneSpInteraction(char t, double m):type(t), magnitude(m)
    {
        assert(type == '0' || type == 'A' || type == 'B' || type == 'C' || type == 'D');
    }
};


class GeneralParameters // is a singleton
{
    // attribute starting with '__' means that they store info about sequential parameters. Attributes being used have the same name without the '__' and the are reference of a single "line" of the element starting with '__'. '__GenerationChange' indicate when the parameters should be updated. All demographic parameters must have an entry for the same generations.
public:
    // patchSize
    //std::vector<std::vector<int>>                   allSpeciesPatchSizes;  // allSpeciesPatchSizes[patch_index][speciesIndex];
    std::vector<std::vector<int>>                   allSpeciesPatchSizePreviousGeneration; // allSpeciesPatchSizePreviousGeneration[patch_index][speciesIndex];
    std::vector<std::string>                        speciesNames;

    std::vector<int>                                __GenerationChange;
    
    // Technic
    RNG_wrapper                                     rngw;
    int                                             nbThreads;
    int                                             OverwriteMode;
    bool                                            DryRun;
    bool                                            printProgress;


    // BinaryFileToRead The path will be species specific though and are stored in SSP
    std::ifstream                                   BinaryFileToRead;

    // sequencing error rate
    double                                          sequencingErrorRate;
    
    // nbGenerations
    int                                             nbGenerations;
    int                                             startAtGeneration;
    int                                             CurrentGeneration;
    int                                             burnInUntilT4Coal_check_every_N_generations;

    // SpeciesStuff
    int                                             nbSpecies;
    //std::vector<char>                               __typeOfSpeciesInteraction;
    std::vector<std::vector<std::vector<oneSpInteraction>>>   __speciesInteraction;
    std::vector<std::vector<std::vector<double>>>   __speciesCompetition;
    //char                                            typeOfSpeciesInteraction;
    std::vector<std::vector<oneSpInteraction>>                speciesInteraction; // speciesInteraction[toSpIndex][fromSpIndex]
    std::vector<std::vector<double>>                speciesCompetition; // speciesCompetition[toSpIndex][fromSpIndex]

    // nbPatches
    int                                             PatchNumber;
    int                                             maxEverPatchNumber;
    std::vector<int>                                __PatchNumber;    

    // Manage outputs
    int                                             nbWarningsSentFrom_Haplotype_AddMutT2_Allele = 0;
    std::vector<int>                                output_FST_nbPatchesToConsider;

    // methods
    //void readTemporalChanges(InputReader& input); // not in use anymore
    void readTemporalChanges(std::vector<int>& T);
    void readPatchNumber(InputReader& input);
    void readSpeciesEcologicalRelationships(InputReader& input);
    void readT1_FST_info(InputReader& input);
    void readSeed(InputReader& input);
    void readBurnInUntilT4Coal(InputReader& input);
    void testIfEndOfT4BurnIn(Pop& pop_Offspring, std::vector<bool>& neutralBurnIn_hasSpeciesCoalesced);
    //void saveSSPPatchSize_toGP();
    //void saveSSPPatchSize_toGP_lowerSecurity();
    //void UpdateParametersallSpeciesPatchSizes();
    void update(int generation_index);
    void setAllPatchSizePreviousGenerationIfNeeded();
    
    //void initializeAllSpeciesPatchSizes();
    GeneralParameters();
};
