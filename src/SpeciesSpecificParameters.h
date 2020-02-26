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
private:
    void resetT56_freqThresholToDefault();
    bool T5_freqThresholdWasSetToDefault;

    bool isT56LocusUnderSelection(size_t locus);


public:

    Genealogy                                       genealogy;
    ResetGenetics                                   resetGenetics;
    
    // species ID
    std::string                                     speciesName;
    int                                             speciesIndex;
    int                                             whenDidExtinctionOccur = -1; // -1 means not extinct (yet)

    //Initialization
    std::string                                     readPopFromBinaryPath;
    bool                                            readPopFromBinary;
    std::vector<std::vector<int>>                   funkyMathForQuasiRandomT1AllFreqInitialization;
    std::vector<std::vector<int>>                   funkyMathForQuasiRandomT56AllFreqInitialization;

    // Basic Demography
    std::vector<int>                                patchCapacity;
    std::vector<int>                                maxEverpatchCapacity;
    std::vector<std::vector<int>>                   __patchCapacity;
    std::vector<int>                                patchSize;
    int                                             TotalpatchCapacity;
    int                                             TotalpatchSize;
    int                                             nbSubGenerationsPerGeneration;

    double                                          cloningRate;
    double                                          selfingRate;
    bool                                            malesAndFemales;
    double                                          sexRatio;


    // basic selection
    char                                            selectionOn;
    bool                                            T1_isLocalSelection;
    bool                                            T2_isLocalSelection;
    bool                                            T3_isLocalSelection;
    bool                                            T56_isLocalSelection;
    
    // selectionOn = 0; "fertility"
    // selectionOn = 1; "viability"
    // selectionOn = 2; "both"


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
    std::vector<int>                                FromT56selLocusToLocus;
    std::vector<int>                                FromT56ntrlLocusToLocus;
    std::vector<std::pair<bool, unsigned>>          FromT56LocusToT56genderLocus; // true means ntrl.
    std::vector<bool>                               isT56neutral;

    // Fitness Map
    std::vector<int>                                FromLocusToFitnessMapIndex;
    std::vector<int>                                FitnessMapInfo_wholeDescription;
    std::vector<FromLocusToTXLocusElement>          FromFitnessMapIndexToTXLocus;
    int                                             NbElementsInFitnessMap;
    double                                          FitnessMapProbOfEvent;
    double                                          FitnessMapT5WeightProbOfEvent;
    double                                          FitnessMapCoefficient;
    int                                             FitnessMapMinimNbLoci;

    // Other performance related things
    //int                 recomputeLociOverWhichFitnessMustBeComputedEveryHowManyGenerations;
    bool                                            SwapInLifeCycle;
    //bool allowToCorrectRecomputeLociOverWhichFitnessMustBeComputedEveryHowManyGenerations;

    // Genetics and selection Both T1, T2, T3, T4 and T5
    int                                             ploidy;
    std::vector<double>                             RecombinationRate;
    double                                          TotalRecombinationRate;

    int                                             TotalNbLoci;
    bool                                            T1_isMultiplicitySelection;
    bool                                            T56_isMultiplicitySelection;
    bool                                            recRateOnMismatch_bool;
    int                                             recRateOnMismatch_halfWindow;
    double                                          recRateOnMismatch_factor;
    double                                          fecundityForFitnessOfOne;

    std::vector<std::vector<double>>                __growthK;
    std::vector<double>                             growthK; // It is a double to avoid recasting but values must be integers.
    bool                                            stochasticGrowth;
    
    // Genetics and selection T1
    std::vector<std::vector<double>>                T1_Initial_AlleleFreqs;
    bool                                            T1_Initial_AlleleFreqs_AllZeros;
    bool                                            T1_Initial_AlleleFreqs_AllOnes;
    
    std::map<std::string, Individual>               individualTypes;
    std::vector<std::vector<std::string>>           individualTypesForInitialization;

    bool                                            isIndividualInitialization;
    int                                             T1_nbChars;
    int                                             T1_nbLoci;
    std::vector<std::vector<double>>                T1_FitnessEffects;
    bool                                            T1_isSelection;
    bool                                            T1_isEpistasis;
    std::vector<double>                             T1_MutationRate;  // cumulative
    double                                          T1_Total_Mutation_rate;
    
    std::vector<std::vector<std::vector<T1_locusDescription>>>   T1_Epistasis_LociIndices;
    std::vector<std::vector<std::vector<double>>>                T1_Epistasis_FitnessEffects; //get fitness value with this->FitnessEffects_ForASingleHabitat[habitat][groupOfLoci][fitnessValueIndex], where the fitnessValueIndex is function of the genotype
    int                                             T1_nbLociLastByte;
    

    
    // Genetics and selection T2
    int                                             T2_nbLoci;
    std::vector<std::vector<double>>                T2_FitnessEffects;
    bool                                            T2_isSelection;
    std::vector<double>                             T2_MutationRate;  // cumulative
    double                                          T2_Total_Mutation_rate;    

    // Genetics and selection T3
    int                                             T3_nbLoci;
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
    int                                             T5_nbLoci;
    int                                             T5ntrl_nbLoci;
    int                                             T5sel_nbLoci;
    std::vector<unsigned int>                       T5ntrl_flipped;
    //std::vector<unsigned int>                       T5sel_flipped;


    // T6 (compressed T5)
    int                                             T6_nbLoci;
    int                                             T6ntrl_nbLoci;
    int                                             T6sel_nbLoci;
    CompressedSortedDeque                           T6ntrl_flipped;
    //CompressedSortedDeque                           T6sel_flipped;

    // T56 (both T5 and T6)
    int                                             T56_nbLoci;
    int                                             T56ntrl_nbLoci;
    int                                             T56sel_nbLoci;
    bool                                            T56_isSelection;
    bool                                            T56_isMuliplicitySelection;
    std::vector<std::vector<double>>                T56_Initial_AlleleFreqs;
    bool                                            T56_Initial_AlleleFreqs_AllZeros;
    bool                                            T56_Initial_AlleleFreqs_AllOnes;
    int                                             T56_toggleMutsEveryNGeneration;
    int                                             T56_toggleMutsEveryNGeneration_nextGeneration;
    double                                          T56ntrl_frequencyThresholdForFlippingMeaning;
    double                                          T56sel_frequencyThresholdForFlippingMeaning;
    double                                          T56_approximationForNtrl;
    bool                                            T56ntrl_compress;
    bool                                            T56sel_compress;
    int                                             quickScreenAtL_T56_nbLoci;
    std::vector<double>                             T56_MutationRate;  // cumulative
    double                                          T56_Total_Mutation_rate;
    std::vector<std::vector<double>>                T56_FitnessEffects;
    
    
    // general fitness
    bool                                            isAnySelection;


    // Ecology
    std::vector<int>                                Habitats;
    std::vector<std::vector<int>>                   __Habitats;
    int                                             MaxHabitat;
    std::vector<int>                                __MaxHabitat;
    int                                             MaxEverHabitat;

    // Random shit
    GeneticSampler                                  geneticSampler;
    bool                                            geneticSampling_withWalker;
    bool                                            individualSampling_withWalker;
    
    
    
    // Other
    int                                              centralT1LocusForExtraGeneticInfo;
    std::vector<double>                              outputSFSbinSizes;
    KillOnDemand                                     killOnDemand;
    std::string                                      killOnDemand_FunctionName;
    std::string                                      KillOnDemand_homemade_parameters;
    std::string                                      ifFixedKillOnDemand_loci;
    
    std::vector<std::vector<int>>                    subsetT1LociForfitnessSubsetLoci_file;
    std::vector<std::vector<int>>                    subsetT2LociForfitnessSubsetLoci_file;
    std::vector<std::vector<int>>                    subsetT3LociForfitnessSubsetLoci_file;
    std::vector<std::vector<int>>                    subsetT5LociForfitnessSubsetLoci_file;
    std::vector<std::vector<int>>                    subsetT1epistasisLociForfitnessSubsetLoci_file;
    
    // methods
    void readLoci(InputReader& input);
    void readKillOnDemand(InputReader& input);
    void readT1_Initial_AlleleFreqs(InputReader& input);
    void readT56_Initial_AlleleFreqs(InputReader& input);
    void readnbSubGenerations(InputReader& input);
    void readRecombinationRate(InputReader& input);
    void readRecRateOnMismatch(InputReader& input);
    void readT1_MutationRate(InputReader& input);
    void readT2_MutationRate(InputReader& input);
    void readT3_MutationRate(InputReader& input);
    void readT4_MutationRate(InputReader& input);
    void readT56_MutationRate(InputReader& input);
    void readadditiveEffectAmongLoci(InputReader& input);
    void readT3_PhenotypicEffects(InputReader& input);
    void readT1_EpistaticFitnessEffects(InputReader& input);
    void readT1_FitnessEffects(InputReader& input);
    void readT2_FitnessEffects(InputReader& input);
    void readT3_FitnessLandscape(InputReader& input);
    void readT3_DevelopmentalNoise(InputReader& input);
    void readT4_maxAverageNbNodesPerHaplotype(InputReader& input);
    void readT56_FitnessEffects(InputReader& input);
    void readT56_compress(InputReader& input);
    void readT56_approximationForNtrl(InputReader& input);
    void readT56_freqThreshold(InputReader& input);
    //void readreverseFixedT5selMuts(InputReader& input);
    void readT1mutsDirectional(InputReader& input);
    //void readT5mutsDirectional(InputReader& input);
    void readResetGenetics(InputReader& input);
    ResetGeneticsEvent_A readResetGenetics_readEventA(InputReader& input, int& eventIndex);
    ResetGeneticsEvent_B readResetGenetics_readEventB(InputReader& input, int& eventIndex);
    void readHabitats(InputReader& input);
    void readSelectionOn(InputReader& input);
    void readPopGrowthModel(InputReader& input);
    //void readResetTrackedT1Muts(InputReader& input);
    void readGameteDispersal(InputReader& input);
    void readpatchCapacity(InputReader& input);
    //void readT1_vcfOutput_sequence(InputReader& input);
    void readSwapInLifeCycle(InputReader& input);
    void readDispMat(InputReader& input);
    void readCentralT1LocusForExtraGeneticInfo(InputReader& input);
    Haplotype getHaplotypeForIndividualType(InputReader& input, bool haploIndex, std::string& IndividualTypeName);
    void readIndividualTypes(InputReader& input);
    void readIndividualInitialization(InputReader& input);
    void readInitialpatchSize(InputReader& input);
    void readCloningRate(InputReader& input);
    void readSelfingRate(InputReader& input);
    void readStochasticGrowth(InputReader& input);
    void readDispWeightByFitness(InputReader& input);
    void readPloidy(InputReader& input);
    void readFitnessMapInfo(InputReader& input);
    void readT56_toggleMutsEveryNGeneration(InputReader& input);
    void readfecundityForFitnessOfOne(InputReader& input);
    void readMatingSystem(InputReader& input);
    void readReadPopFromBinary(InputReader& input);
    void readSubsetLociForfitnessSubsetLoci_file(InputReader& input);
    void readOutputSFSbinSize(InputReader& input);
    void readT4_printTree(InputReader& input);
    void readGeneticSampling_withWalker(InputReader& input);
    void readIndividualSampling_withWalker(InputReader& input);
    int selectNonEmptyPatch(int firstPatchToLookAt, std::vector<int>& PSs, bool increasingOrder);

    void IsThereSelection();
    //void setInitialT1_AlleleFreqTo(const int uniqueFreq);
    //void ClearT1_Initial_AlleleFreqs();
    bool setFromLocusToFitnessMapIndex_DecidingFunction(double sumOfProb, int nbLoci);
    void setFromLocusToFitnessMapIndex();

    void setRandomDistributions();
    std::vector<int> UpdateParameters(int generation_index);

    SpeciesSpecificParameters(std::string sN, int sI);
    SpeciesSpecificParameters();
    ~SpeciesSpecificParameters();


    void readQuickScreenOfOptionL(InputReader& input);

};
