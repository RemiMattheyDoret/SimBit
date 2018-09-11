
class OptionContainer
{
public:
	std::vector<Option> options =
    {
        {{"seed", "random_seed"},                       {}},

    	// Basics simulations
        {{"nbGens","nbGenerations"},                    {}},
        {{"startAtGeneration"},                         {"nbGens"}},
        {{"nbThreads"},                                 {}},
        {{"S", "species"},                              {}},
        {{"nbSubGens","nbSubGenerations"},              {"nbGens", "S"}},
        {{"T", "TemporalChanges"}, 		  				{"nbGens"}},   // all species must have temporal change at the same time
        {{"PN","PatchNumber"},                          {"T"}},
        {{"H", "Habitats"},                           	{"S","PN","T"}}, // not all species need to have same habitat distinction

        // Paths and times
        {{"GP", "GeneralPath"},         		    {}},
        {{"T1_vcf_file","T1_VCF_file"},             {"GP", "S", "startAtGeneration"}},
        {{"T1_LargeOutput_file"},           		{"GP", "S", "startAtGeneration"}},
        {{"T1_AlleleFreq_file"},            		{"GP", "S", "startAtGeneration"}},
        {{"Log","Logfile", "Logfile_file"},         {"GP", "S", "startAtGeneration"}},
        {{"T1_MeanLD_file"},             			{"GP", "S", "startAtGeneration"}},
        {{"T1_LongestRun_file"},           			{"GP", "S", "startAtGeneration"}},
        {{"T1_HybridIndex_file"},          			{"GP", "S", "startAtGeneration"}},
        {{"T1_ExpectiMinRec_file"},         		{"GP", "S", "startAtGeneration"}},
        {{"T2_LargeOutput_file"},          			{"GP", "S", "startAtGeneration"}},
        {{"SaveBinary_file"},           			{"GP", "S", "startAtGeneration"}},
        {{"T3_LargeOutput_file"},          			{"GP", "S", "startAtGeneration"}},
        {{"T3_MeanVar_file"},        	  			{"GP", "S", "startAtGeneration"}},
        {{"fitness_file"},         		 			{"GP", "S", "startAtGeneration"}},
        {{"fitnessStats_file"},          			{"GP", "S", "startAtGeneration"}},
        {{"T1_allTypesFST_file"},          			{"GP", "S", "startAtGeneration"}},
        {{"extraGeneticInfo_file"},         		{"GP", "S", "startAtGeneration"}},
        {{"patchSize_file"},          				{"GP", "S", "startAtGeneration"}},
        {{"extinction_file"},                       {"GP", "S", "startAtGeneration"}},
        {{"genealogy_file"},                        {"GP", "S", "startAtGeneration"}},
        {{"coalesce","shouldGenealogyBeCoalesced"}, {"GP", "S", "startAtGeneration"}},
    
        
        // Dispersal
        {{"m", "DispMat"},                          	{"PN","S"}},
        {{"DispWeightByFitness"},         				{"S"}},
        {{"gameteDispersal"},                           {"S"}},
        
        // Genetics and Selection Both T1, T2 and T3 (first part)
        {{"L","Loci"},                           		{"S"}},
        {{"ploidy"},                      				{"S"}},
        {{"fec","fecundityForFitnessOfOne"},            {"S"}},

        // Basics Demography
        {{"N", "patchCapacity"},                        {"PN","S","T"}},
        {{"InitialpatchSize"},                          {"PN","S","N","fec"}},
        {{"cloningRate"},                               {"S"}},
        {{"selfingRate"},                               {"S"}},
        {{"matingSystem"},                              {"S","T","fec","N","PN"}}, // if fec != -1.0, then I must make sure that there will be at least one male and one female in every patch at all times. This assumption might be relaxed later

        // Other stuff about outputs
        {{"LogfileType"},                               {}},
        {{"T1_SubsetOut","T1_output_sequence"},         {"L"}}, // specify [from, to]
        {{"sequencingErrorRate"},                       {}},
        
        // Selection
        {{"additiveEffectAmongLoci"},                   {"S"}},

        // Genetics and Selection T1
        {{"T1_mu", "T1_MutationRate"},                  {"S", "L"}},
        {{"T1_fit", "T1_FitnessEffects"},               {"H","S", "seed", "L"}},
        {{"T1_ini", "T1_Initial_AlleleFreqs"},          {"S", "L"}},
        {{"T1_epistasis","T1_EpistaticFitnessEffects"}, {"H","S", "L"}},
        
        // Genetics and Selection T2
        {{"T2_mu","T2_MutationRate"},                   {"S", "L"}},
        {{"T2_fit", "T2_FitnessEffects"},               {"S","H", "L"}},

        // Genetics and Selection T3
        {{"T3_mu","T3_MutationRate"},                   {"S", "L"}},
        {{"T3_pheno", "T3_PhenotypicEffects"},          {"S", "H", "L"}},
        {{"T3_fit","T3_FitnessLandscape"},              {"S", "H", "L", "T3_pheno"}},
        {{"T3_DN","T3_DevelopmentalNoise"},             {"S", "H", "L", "T3_pheno"}},

        // Genetics and Selection Both T1, T2 and T3 (second part)
        {{"r","RecombinationRate"},                     {"S", "L"}},
        {{"recRateOnMismatch"},                         {"S", "L"}},
        {{"FitnessMapInfo"},           {"S", "T", "H", "L","T1_mu","T2_mu","T3_mu","r","m","T1_fit","T2_fit"}},
        {{"resetTrackedT1Muts"},               {"S","L","T1_mu","T1_fit","N","PN"}},

        // Species interaction
        {{"eco", "ecoRelation","SpeciesEcologicalRelationships"},	{"S","seed"}},
        {{"growthK", "growthCarryingCapacity"},         {"eco", "S", "PN", "N"}},

        // other technical parameters
        {{"Overwrite"},                   		            {}},
        {{"readPopFromBinary"},           		            {"S"}},
        {{"DryRun"},                         		        {}},
        {{"centralT1LocusForExtraGeneticInfo"},             {}},
        {{"resetGenetics"},                                 {"L","nbGenerations","T","S","PN","N"}}
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


