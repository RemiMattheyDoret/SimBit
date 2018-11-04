
class OptionContainer
{
public:
	std::vector<Option> options =
    {
        {{"seed", "random_seed"},                       {}, false},

    	// Basics simulations
        {{"nbGens","nbGenerations"},                    {}, false},
        {{"startAtGeneration"},                         {"nbGens"}, false},
        {{"nbThreads"},                                 {}, false},
        {{"S", "species"},                              {}, false},
        {{"nbSubGens","nbSubGenerations"},              {"nbGens", "S"}, false},
        {{"T", "TemporalChanges"}, 		  				{"nbGens"}, false},   // all species must have temporal change at the same time
        {{"PN","PatchNumber"},                          {"T"}, false},
        {{"H", "Habitats"},                           	{"S","PN","T"}, false}, // not all species need to have same habitat distinction
    
        
        // Genetics and Selection Both T1, T2 and T3 (first part)
        {{"L","Loci"},                           		{"S"}, false},
        {{"ploidy"},                      				{"S"}, false},
        {{"fec","fecundityForFitnessOfOne"},            {"S"}, false},

        // Dispersal
        {{"m", "DispMat"},                              {"PN","S"}, false},
        {{"DispWeightByFitness"},                       {"S","fec"}, false},
        {{"gameteDispersal"},                           {"S"}, false},

        // Basics Demography
        {{"N", "patchCapacity"},                        {"PN","S","T"}, false},
        {{"InitialpatchSize"},                          {"PN","S","N","fec"}, false},
        {{"cloningRate"},                               {"S"}, false},
        {{"selfingRate"},                               {"S"}, false},
        {{"matingSystem"},                              {"S","T","fec","N","PN"}, false}, // if fec != -1.0, then I must make sure that there will be at least one male and one female in every patch at all times. This assumption might be relaxed later

        // Other stuff about outputs
        {{"LogfileType"},                               {}, false},
        {{"sequencingErrorRate"},                       {}, false},
        
        // Selection
        {{"additiveEffectAmongLoci"},                   {"S"}, false},

        // Genetics and Selection T1
        {{"T1_mu", "T1_MutationRate"},                  {"S", "L"}, false},
        {{"T1_fit", "T1_FitnessEffects"},               {"H","S", "seed", "L"}, false},
        {{"T1_ini", "T1_Initial_AlleleFreqs"},          {"S", "L"}, false},
        {{"T1_epistasis","T1_EpistaticFitnessEffects"}, {"H","S", "L"}, false},
        
        // Genetics and Selection T2
        {{"T2_mu","T2_MutationRate"},                   {"S", "L"}, false},
        {{"T2_fit", "T2_FitnessEffects"},               {"S","H", "L"}, false},

        // Genetics and Selection T3
        {{"T3_mu","T3_MutationRate"},                   {"S", "L"}, false},
        {{"T3_pheno", "T3_PhenotypicEffects"},          {"S", "H", "L"}, false},
        {{"T3_fit","T3_FitnessLandscape"},              {"S", "H", "L", "T3_pheno"}, false},
        {{"T3_DN","T3_DevelopmentalNoise"},             {"S", "H", "L", "T3_pheno"}, false},

        // Genetics and Selection Both T1, T2 and T3 (second part)
        {{"r","RecombinationRate"},                     {"S", "L"}, false},
        {{"recRateOnMismatch"},                         {"S", "L"}, false},
        {{"FitnessMapInfo"},           {"S", "T", "H", "L","T1_mu","T2_mu","T3_mu","r","m","T1_fit","T2_fit"}, false},
        {{"resetTrackedT1Muts"},               {"S","L","T1_mu","T1_fit","N","PN"}, false},

        // Outputs
        {{"GP", "GeneralPath"},                     {}, false},
        {{"T1_vcf_file","T1_VCF_file"},             {"GP", "S", "startAtGeneration", "L"}, true},
        {{"T1_LargeOutput_file"},                   {"GP", "S", "startAtGeneration", "L"}, true},
        {{"T1_AlleleFreq_file"},                    {"GP", "S", "startAtGeneration", "L"}, true},
        {{"Log","Logfile", "Logfile_file"},         {"GP", "S", "startAtGeneration", "L"}, false},
        {{"T1_MeanLD_file"},                        {"GP", "S", "startAtGeneration", "L"}, true},
        {{"T1_LongestRun_file"},                    {"GP", "S", "startAtGeneration", "L"}, true},
        {{"T1_HybridIndex_file"},                   {"GP", "S", "startAtGeneration", "L"}, true},
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

        // Species interaction
        {{"eco", "ecoRelation","SpeciesEcologicalRelationships"},	{"S","seed"}, false},
        {{"growthK", "growthCarryingCapacity"},         {"eco", "S", "PN", "N"}, false},

        // other technical parameters
        {{"Overwrite"},                   		            {}, false},
        {{"readPopFromBinary"},           		            {"S"}, false},
        {{"DryRun"},                         		        {}, false},
        {{"centralT1LocusForExtraGeneticInfo"},             {}, false},
        {{"resetGenetics"},                                 {"L","nbGenerations","T","S","PN","N"}, false}
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


    unsigned int minVector(const std::vector<unsigned int> v);
    double mean_levenshtein_distance(const std::string& s1, const std::string& s2);
};


