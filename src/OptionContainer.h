
class OptionContainer
{
public:
	std::vector<Option> options =
    {
        {{"seed", "random_seed"},                       {}, false},
        {{"printProgress"},                             {}, false},
        {{"removeInputFileAfterReading"},               {}, false},
        
    	// Basics simulations
        {{"nbGens","nbGenerations"},                    {}, false},
        {{"startAtGeneration"},                         {"nbGens"}, false},
        {{"nbThreads"},                                 {}, false},
        {{"S", "species"},                              {}, false},
        {{"nbSubGens","nbSubGenerations"},              {"nbGens", "S"}, false},
        //{{"T", "TemporalChanges"}, 		  				{"nbGens"}, false},   // all species must have temporal change at the same time THIS IS NOT AN OPTION ANYMORE
        {{"PN","PatchNumber"},                          {}, false},
        {{"N", "patchCapacity"},                        {"PN","S",}, false},
        {{"H", "Habitats"},                           	{"S","PN",}, false}, // not all species need to have same habitat distinction
    
        // T5_fit
        {{"T5_approximationForNtrl","T56_approximationForNtrl"},                   {"H","S", "seed"}, false},
        {{"T5_fit", "T5_FitnessEffects"},               {"H","S", "seed", "T5_approximationForNtrl"}, false}, // First quickScreenOfL, then read T5_fit and finally read --L. In --L the T56_FitnessEffects objects is reduced to remove what is not under selection.
        {{"T5_compressData", "T5_compress"},            {"H","S", "seed", "T5_fit"}, false},
        
        // Genetics and Selection Both T1, T2 and T3 (first part)
        {{"L","Loci"},                           		{"S", "nbSubGens", "T5_fit", "T5_approximationForNtrl", "T5_compress", "PN", "N"}, false},
        {{"ploidy"},                      				{"S"}, false},
        {{"fec","fecundityForFitnessOfOne"},            {"S"}, false},
        {{"fecundityDependentOfFitness"},               {"S"}, false},

        // Dispersal
        {{"m", "DispMat"},                              {"PN","S"}, false},
        {{"DispWeightByFitness"},                       {"S","fec"}, false},
        {{"gameteDispersal"},                           {"S"}, false},
        {{"forcedMigration"},                           {"PN", "N", "S"}, false},
        {{"stochasticMigration"},                       {"S"}, false},        


        // Basics Demography
        {{"InitialpatchSize"},                          {"PN","S","N","fec"}, false},
        {{"cloningRate"},                               {"S"}, false},
        {{"selfingRate"},                               {"S"}, false},
        {{"matingSystem"},                              {"S","fec","N","PN"}, false}, // if fec != -1.0, then I must make sure that there will be at least one male and one female in every patch at all times. This assumption might be relaxed later
        {{"killIndividuals"},                           {"S","fec","N","PN"}, false},

        // Other stuff about outputs
        {{"LogfileType"},                               {}, false},
        {{"sequencingErrorRate"},                       {}, false},
        
        // Selection
        {{"additiveEffectAmongLoci"},                   {"S"}, false},
        {{"selectionOn"},                               {"L","S"}, false},


        // Genetics and Selection T1
        {{"T1_mu", "T1_MutationRate"},                  {"S", "L"}, false},
        {{"T1_epistasis","T1_EpistaticFitnessEffects"}, {"H","S", "L"}, false},
        {{"T1_mutDirection"},                           {"L","S"}, false},
        //{{"T1mutsDirectional"},                         {"S", "L"}, false},

        // Genetics and Selection T2
        {{"T2_mu","T2_MutationRate"},                   {"S", "L"}, false},

        // Genetics and Selection T4
        {{"T4_mu","T4_MutationRate"},                   {"S", "L"}, false},
        {{"T4_mutDirection"},                           {"L","S"}, false},
        {{"T4_nbRunsToPlaceMutations"},                 {"L","S"}, false},

        // Genetics and Selection T5
        {{"T5_mu","T5_MutationRate"},                   {"S", "L"}, false},
        {{"T56_mutDirection"},                           {"L","S"}, false},
        // T5_fit appears before --L
        {{"T5_toggleMutsEveryNGeneration"},                   {"S", "L"}, false},
        {{"T5_freqThreshold", "T5_frequencyThresholdForFlippingMeaning"},              {"S", "L","PN","N"}, false},
        // {{"reverseFixedT5selMuts"},              {"S", "L"}, false},
        //{{"T5mutsDirectional"},                         {"S", "L", "H","T5_fit"}, false},

        // Genetics and Selection T3
        {{"T3_mu","T3_MutationRate"},                   {"S", "L"}, false},
        {{"T3_mutationalEffect"},                       {"S", "L"}, false},
        {{"T3_pheno", "T3_PhenotypicEffects"},          {"S", "H", "L"}, false},
        {{"T3_fit","T3_FitnessLandscape"},              {"S", "H", "L", "T3_pheno"}, false},
        {{"T3_DN","T3_DevelopmentalNoise"},             {"S", "H", "L", "T3_pheno"}, false},

        {{"T1_fit", "T1_FitnessEffects"},               {"H","S", "seed", "L"}, false},
        {{"T2_fit", "T2_FitnessEffects"},               {"S","H", "L"}, false},

        // Genetics and Selection for T1, T2, T3 and T4 (second part)
        {{"r","RecombinationRate"},                     {"S", "L"}, false},
        {{"recRateOnMismatch"},                         {"S", "L"}, false},
        {{"FitnessMapInfo"},           {"S", "H", "L","T1_mu","T2_mu","T3_mu","r","m","T1_fit","T2_fit"}, false},
        //{{"resetTrackedT1Muts"},               {"S","L","T1_mu","T1_fit","N","PN"}, false},


        {{"T4_simplifyEveryNGenerations"},               {"S", "L", "N", "r"}, false},

        // Initializer
        {{"indTypes","individualTypes"},                {"L", "T5_fit", "T1_fit", "T2_fit", "T3_fit", "FitnessMapInfo", "InitialpatchSize", "nbGenerations"}, false},
        {{"redefIndTypes"},                             {"indTypes"}, false},
        {{"resetGenetics"},                             {"indTypes","L","nbGenerations","S","PN","N", "T5_MutationRate"}, false},
        {{"indIni","individualInitialization"},         {"indTypes", "L", "T5_fit", "T1_fit", "T2_fit", "T3_fit", "FitnessMapInfo", "InitialpatchSize"}, false},
        {{"T1_ini", "T1_Initial_AlleleFreqs"},          {"S", "L","PN","indIni"}, false},
        {{"T5_ini", "T5_Initial_AlleleFreqs"},          {"S", "L","PN", "indIni"}, false},
        {{"burnInUntilT4Coal", "burnIn"},                         {"S","L"}, false},


        // Outputs
        {{"SegDiversityFile_includeMainColor"},     {"L"}, true},
        {{"SNPfreqCalculationAssumption"},          {"L", "T4_mu", "N", "T5_mu", "T1_mu"}, false},
        {{"GP", "GeneralPath"},                     {}, false},
        {{"T1_vcf_file","T1_VCF_file"},             {"GP", "S", "startAtGeneration", "L"}, true},
        {{"T1_LargeOutput_file"},                   {"GP", "S", "startAtGeneration", "L"}, true},
        {{"T1_AlleleFreq_file"},                    {"GP", "S", "startAtGeneration", "L"}, true},
        {{"Log","Logfile", "Logfile_file"},         {"GP", "S", "startAtGeneration", "L","LogfileType"}, false},
        {{"T1_MeanLD_file"},                        {"GP", "S", "startAtGeneration", "L"}, true},
        {{"T1_LongestRun_file"},                    {"GP", "S", "startAtGeneration", "L"}, true},
        {{"T1_HybridIndex_file"},                   {"GP", "S", "startAtGeneration", "L"}, true},
        {{"T1_AverageHybridIndex_file"},            {"GP", "S", "startAtGeneration", "L"}, true},
        {{"T1_ExpectiMinRec_file"},                 {"GP", "S", "startAtGeneration", "L"}, true},
        {{"T2_LargeOutput_file"},                   {"GP", "S", "startAtGeneration", "L"}, false},
        {{"SaveBinary_file"},                       {"GP", "S", "startAtGeneration", "L"}, false},
        {{"T3_LargeOutput_file"},                   {"GP", "S", "startAtGeneration", "L"}, false},
        {{"T3_MeanVar_file"},                       {"GP", "S", "startAtGeneration", "L"}, false},
        {{"fitness_file"},                          {"GP", "S", "startAtGeneration", "L"}, false},
        {{"fitnessSubsetLoci_file"},                {"GP", "S", "startAtGeneration", "L"}, false},
        {{"fitnessStats_file"},                     {"GP", "S", "startAtGeneration", "L"}, false},
        {{"T1_FST_file"},                           {"GP", "S", "startAtGeneration", "L"}, true},
        {{"T1_FST_info"},                           {"GP", "S", "startAtGeneration", "T1_FST_file","PN", "L"}, false},
        {{"extraGeneticInfo_file"},                 {"GP", "S", "startAtGeneration", "L"}, true},
        {{"patchSize_file"},                        {"GP", "S", "startAtGeneration", "L"}, false},
        {{"extinction_file"},                       {"GP", "S", "startAtGeneration", "L"}, false},
        {{"genealogy_file"},                        {"GP", "S", "startAtGeneration", "L"}, false},
        {{"coalesce","shouldGenealogyBeCoalesced"}, {"GP", "S", "startAtGeneration", "L"}, false},
        {{"T4_LargeOutput_file"},                   {"GP", "S", "startAtGeneration", "L"}, true},
        {{"T4_vcf_file","T4_VCF_file"},             {"GP", "S", "startAtGeneration", "L"}, true},
        {{"T4_SFS_file"},                           {"GP", "S", "startAtGeneration", "L"}, true},
        {{"T4_SNPfreq_file"},                       {"GP", "S", "startAtGeneration", "L"}, true},
        {{"Tx_SNPfreq_file"},                       {"GP", "S", "startAtGeneration", "L"}, true},
        {{"T1_SFS_file"},                           {"GP", "S", "startAtGeneration", "L"}, true},
        {{"T1_haplotypeFreqs_file"},                {"GP", "S", "startAtGeneration", "L"}, true},
        {{"burnInLength_file"},                     {"GP", "S", "startAtGeneration", "L"}, true},
        //{{"T4_printTree"},                          {"GP", "S", "startAtGeneration", "L"}, true},
        //{{"T4_coalescenceFst_file"},                {"GP", "S", "startAtGeneration", "L"}, true},

        {{"T5_vcf_file","T5_VCF_file"},             {"GP", "S", "startAtGeneration", "L"}, true},
        {{"T5_SFS_file"},                                {"GP", "S", "startAtGeneration", "L"}, true},
        {{"T5_AlleleFreq_file"},                    {"GP", "S", "startAtGeneration", "L"}, true},
        {{"T5_LargeOutput_file"},                   {"GP", "S", "startAtGeneration", "L"}, true},

        {{"outputSFSbinSize"},                      {"T4_SFS_file", "T1_SFS_file", "GP", "S", "startAtGeneration", "L"}, true},
        {{"sampleSeq_file"},                        {"GP", "S", "PN", "L"}, true},

        {{"T4_paintedHaplo_ignorePatchSizeSecurityChecks"},                  {}, true},
        {{"T4_paintedHaplo_file"},                  {"GP", "S", "startAtGeneration", "L", "T4_vcf_file", "T4_SFS_file", "sampleSeq_file", "T4_paintedHaplo_ignorePatchSizeSecurityChecks"}, true},
        {{"T4_paintedHaploSegmentsDiversity_file"},{"GP", "S", "startAtGeneration", "L", "T4_vcf_file", "T4_SFS_file", "sampleSeq_file", "T4_paintedHaplo_ignorePatchSizeSecurityChecks"}, true},

        
        {{"T4_nbMutationPlacingsPerOutput"},                        {"S", "L"}, true},
        {{"T4_respectPreviousMutationsOutputs"},                    {"S", "L", "T4_nbMutationPlacingsPerOutput", "T4_paintedHaplo_file","T4_paintedHaploSegmentsDiversity_file"}, true},
        

        // Species interaction
        {{"eco","speciesEcologicalRelationships"},	{"S","seed"}, false},
        {{"popGrowthModel"},                        {"eco", "S", "PN", "N", "fec"}, false},
        {{"stochasticGrowth"},                        {"popGrowthModel"}, false},
        {{"swapInLifeCycle"},                       {"S","r","L"}, false},
        

        // other technical parameters
        {{"Overwrite"},                   		            {}, false},
        {{"readPopFromBinary"},           		            {"S"}, false},
        {{"DryRun"},                         		        {}, false},
        {{"centralT1LocusForExtraGeneticInfo"},             {}, false},
        {{"killOnDemand"},                                  {}, false},
        {{"geneticSampling_withWalker"},                    {}, false},
        {{"individualSampling_withWalker"},                 {"PN", "N"}, false},
        {{"shrinkT56EveryNGeneration"},                      {}, false}
    };

    //void received(std::string& optionNameReceived);
    //bool hasReceived(std::string name);
    void listOptions();
    int getNbOptions();
    Option& getOption(std::string name);
    int HowManyOptionsWereInitiated();
    bool wasInitiated(std::string name);
    std::string getOptionFirstName(int optionIndex);
    std::string renameFlag(std::string flag);


    uint32_t minVector(const std::vector<uint32_t> v);
    double mean_levenshtein_distance(const std::string& s1, const std::string& s2);
};


