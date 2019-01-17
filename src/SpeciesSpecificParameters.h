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

class SpeciesSpecificParameters
{
       // attribute starting with '__' means that they store info about sequential parameters. Attributes being used have the same name without the '__' and the are reference of a single "line" of the element starting with '__'. '__GenerationChange' indicate when the parameters should be updated. All demographic parameters must have an entry for the same generations.
public:

    SimulationTracker                               simTracker;
    ResetGenetics                                   resetGenetics;
    
    // species ID
    std::string                                     speciesName;
    int                                             speciesIndex;
    bool                                            isSpeciesExtinct = false;

    //Initialization
    std::string                                     readPopFromBinaryPath;
    bool                                            readPopFromBinary;

    // Basic Demography
    std::vector<int>                                patchCapacity;
    std::vector<int>                                maxEverpatchCapacity;
    std::vector<std::vector<int>>                   __patchCapacity;
    std::vector<int>                                patchSize;
    int                                             TotalpatchCapacity;
    int                                             nbSubGenerationsPerGeneration;

    double                                          cloningRate;
    double                                          selfingRate;
    bool                                            malesAndFemales;
    double                                          sexRatio;


    // Dispersal
    DispersalData                                   dispersalData;
    bool                                            DispWeightByFitness;
    bool                                            gameteDispersal;
    
    // Genetic Map
    std::vector<int>                                ChromosomeBoundaries;
    std::vector<FromLocusToTXLocusElement>          FromLocusToTXLocus; // points to the current or previous TXLocus from Locus
    std::vector<int>                                FromT1LocusToLocus;
    std::vector<int>                                FromT2LocusToLocus;
    std::vector<int>                                FromT3LocusToLocus;
    std::vector<int>                                FromT4LocusToLocus;
    std::vector<int>                                FromT5LocusToLocus;

    // Fitness Map
    std::vector<int>                                FromLocusToFitnessMapIndex;
    std::vector<FromLocusToTXLocusElement>          FromLocusToFitnessMapBoundaries;
    int                                             NbElementsInFitnessMap;
    double                                          FitnessMapProbOfEvent;
    double                                          FitnessMapCoefficient;
    int                                             FitnessMapMinimNbLoci;

    // Other performance related things
    int                recomputeLociOverWhichFitnessMustBeComputedEveryHowManyGenerations;
    //bool allowToCorrectRecomputeLociOverWhichFitnessMustBeComputedEveryHowManyGenerations;

    // Genetics and selection Both T1, T2, T3, T4 and T5
    int                                             ploidy;
    std::vector<double>                             RecombinationRate;
    double                                          TotalRecombinationRate;

    int                                             TotalNbLoci;
    bool                                            T1_isMultiplicitySelection;
    bool                                            T5_isMultiplicitySelection;
    bool                                            recRateOnMismatch_bool;
    int                                             recRateOnMismatch_halfWindow;
    double                                          recRateOnMismatch_factor;
    double                                          fecundityForFitnessOfOne;

    std::vector<std::vector<double>>                __growthK;
    std::vector<double>                             growthK; // It is a double to avoid recasting but values must be integers.
    
    // Genetics and selection T1
    std::vector<std::vector<double>>                T1_Initial_AlleleFreqs;
    bool                                            T1_Initial_AlleleFreqs_AllZeros;
    bool                                            T1_Initial_AlleleFreqs_AllOnes;
    int                                             T1_nbChars;
    int                                             T1_nbBits;
    std::vector<std::vector<double>>                T1_FitnessEffects;
    bool                                            T1_isSelection;
    bool                                            T1_isEpistasis;
    std::vector<double>                             T1_MutationRate;  // cumulative
    double                                          T1_Total_Mutation_rate;
    
    std::vector<std::vector<std::vector<T1_locusDescription>>>   T1_Epistasis_LociIndices;
    std::vector<std::vector<std::vector<double>>>                T1_Epistasis_FitnessEffects; //get fitness value with this->FitnessEffects_ForASingleHabitat[habitat][groupOfLoci][fitnessValueIndex], where the fitnessValueIndex is function of the genotype
    int                                             T1_nbBitsLastByte;

    
    // Genetics and selection T2
    int                                             T2_nbChars;
    std::vector<std::vector<double>>                T2_FitnessEffects;
    bool                                            T2_isSelection;
    std::vector<double>                             T2_MutationRate;  // cumulative
    double                                          T2_Total_Mutation_rate;    

    // Genetics and selection T3
    int                                             T3_nbChars;
    std::vector<std::vector<double>>                T3_PhenotypicEffects; // even if it is not selection directly, it is habitat-specific to allow plasticity
    char                                            T3_fitnessLandscapeType;
    //std::vector<double>                             T3_fitnessLandscapeGaussMean;
    //std::vector<double>                             T3_fitnessLandscapeGaussVar;
    std::vector<std::vector<double>>                T3_fitnessLandscapeOptimum;
    std::vector<std::vector<double>>                T3_fitnessLandscapeLinearGradient;
    std::vector<std::vector<double>>                T3_fitnessLandscapeGaussStrength;
    bool                                            T3_isSelection;
    std::vector<double>                             T3_MutationRate;  // cumulative
    double                                          T3_Total_Mutation_rate;
    int                                             T3_PhenoNbDimensions;
    std::vector<std::vector<double>>                T3_DevelopmentalNoiseStandardDeviation;

    // Genetics T4
    int                                             T4_nbBits;
    std::vector<double>                             T4_MutationRate; // cumulative
    Tree                                            T4Tree;
    double                                          T4_maxAverageNbNodesPerHaplotypeBeforeRecalculation;

    // Genetics T5
    int                                             T5_nbBits;
    std::vector<double>                             T5_MutationRate;  // cumulative
    double                                          T5_Total_Mutation_rate;
    std::vector<std::vector<double>>                T5_FitnessEffects;
    bool                                            T5_isSelection;
    bool                                            T5_isMuliplicitySelection;

    // Ecology
    std::vector<int>                                Habitats;
    std::vector<std::vector<int>>                   __Habitats;
    int                                             MaxHabitat;
    std::vector<int>                                __MaxHabitat;
    
    // Random Distributions. If change number of distributions, don't forget to change the incremenet of 'nbParamSet' in 'setRandomDistributions'
    std::poisson_distribution<int>                   rpois_nbRecombination;
    std::poisson_distribution<int>                   T1_rpois_nbMut;
    std::poisson_distribution<int>                   T2_rpois_nbMut;
    std::poisson_distribution<int>                   T3_rpois_nbMut;
    std::poisson_distribution<int>                   T5_rpois_nbMut;
    
    std::uniform_int_distribution<int>               runiform_int_ForRecPos;    // Used when const recombination rate
    std::uniform_real_distribution<double>           runiform_double_ForRecPos;  // Used when variation in recombination rate
    
    std::uniform_int_distribution<int>               T1_runiform_int_ForMutPos;    // Used when const mutation rate
    std::uniform_int_distribution<int>               T2_runiform_int_ForMutPos;    // Used when const mutation rate
    std::uniform_int_distribution<int>               T3_runiform_int_ForMutPos;    // Used when const mutation rate
    std::uniform_int_distribution<int>               T5_runiform_int_ForMutPos;    // Used when const mutation rate
    std::uniform_real_distribution<double>           T1_runiform_double_ForMutPos;  // Used when variation in mutation rate
    std::uniform_real_distribution<double>           T2_runiform_double_ForMutPos;  // Used when variation in mutation rate
    std::uniform_real_distribution<double>           T3_runiform_double_ForMutPos;  // Used when variation in mutation rate
    std::uniform_real_distribution<double>           T5_runiform_double_ForMutPos;  // Used when variation in mutation rate
    
    
    // Other
    int                                              centralT1LocusForExtraGeneticInfo;
    std::vector<double>                              outputSFSbinSizes;
    
    std::vector<std::vector<int>>                    subsetT1LociForfitnessSubsetLoci_file;
    std::vector<std::vector<int>>                    subsetT2LociForfitnessSubsetLoci_file;
    std::vector<std::vector<int>>                    subsetT3LociForfitnessSubsetLoci_file;
    std::vector<std::vector<int>>                    subsetT5LociForfitnessSubsetLoci_file;
    std::vector<std::vector<int>>                    subsetT1epistasisLociForfitnessSubsetLoci_file;
    
    // methods
    void readLoci(InputReader& input);
    void readT1_Initial_AlleleFreqs(InputReader& input);
    void readnbSubGenerations(InputReader& input);
    void readRecombinationRate(InputReader& input);
    void readRecRateOnMismatch(InputReader& input);
    void readT1_MutationRate(InputReader& input);
    void readT2_MutationRate(InputReader& input);
    void readT3_MutationRate(InputReader& input);
    void readT4_MutationRate(InputReader& input);
    void readT5_MutationRate(InputReader& input);
    void readadditiveEffectAmongLoci(InputReader& input);
    void readT3_PhenotypicEffects(InputReader& input);
    void readT1_EpistaticFitnessEffects(InputReader& input);
    void readT1_FitnessEffects(InputReader& input);
    void readT2_FitnessEffects(InputReader& input);
    void readT3_FitnessLandscape(InputReader& input);
    void readT3_DevelopmentalNoise(InputReader& input);
    void readT4_maxAverageNbNodesPerHaplotype(InputReader& input);
    void readT5_FitnessEffects(InputReader& input);
    void readResetGenetics(InputReader& input);
    void readHabitats(InputReader& input);
    void readGrowthK(InputReader& input);
    void readResetTrackedT1Muts(InputReader& input);
    void readGameteDispersal(InputReader& input);
    void readpatchCapacity(InputReader& input);
    //void readT1_vcfOutput_sequence(InputReader& input);
    void readDispMat(InputReader& input);
    void readCentralT1LocusForExtraGeneticInfo(InputReader& input);
    void readInitialpatchSize(InputReader& input);
    void readCloningRate(InputReader& input);
    void readSelfingRate(InputReader& input);
    void readDispWeightByFitness(InputReader& input);
    void readPloidy(InputReader& input);
    void readFitnessMapInfo(InputReader& input);
    void readfecundityForFitnessOfOne(InputReader& input);
    void readMatingSystem(InputReader& input);
    void readReadPopFromBinary(InputReader& input);
    void readSubsetLociForfitnessSubsetLoci_file(InputReader& input);
    void readOutputSFSbinSize(InputReader& input);
    void readT4_printTree(InputReader& input);
    int selectNonEmptyPatch(int firstPatchToLookAt, std::vector<int>& PSs, bool increasingOrder);

    void IsThereSelection();
    //void setInitialT1_AlleleFreqTo(const int uniqueFreq);
    //void ClearT1_Initial_AlleleFreqs();
    bool setFromLocusToFitnessMapIndex_DecidingFunction(double sumOfProb, int nbLoci);
    void setFromLocusToFitnessMapIndex();

    void setRandomDistributions();
    std::vector<int> UpdateParameters(int generation_index);

    SpeciesSpecificParameters(std::string sN, int sI);
    ~SpeciesSpecificParameters();
};
