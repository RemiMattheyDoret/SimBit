struct T7DevParameters
{
    size_t maxAge;
    size_t maxDeltaT;
    double basal_transcription_rate;
    double EPSILON;
    bool stochasticDevelopment;
    bool fitnessOverTime;
    size_t nbCisSites;
    double basic_signal_trans_effect;
    size_t basic_signal_conc;

    double protein_decayRate;
    double mRNA_decayRate;

    double translationRate;
};


struct T7MutParameters
{
	double duplication=0.0;
	double deletion=0.0;
	double cisEffectMu=0.0;
	double cisEffectSD=0.0;
	double transEffectMu=0.0;
	double transEffectSD=0.0;
	double addingInteraction=0.0;
	double changeTargetMu=0.0;
	double kmu=0.0;
    
    double totalMutationRatePerGene=0.0;
};

struct T7PhenParameters
{
    std::vector<size_t> agesAtwhichPhenotypeIsSampled;
    size_t nbDimensions;

    char fitnessLandscapeType;
    double initialPhenotype;
    std::vector<std::vector<std::vector<double>>> fitnessLandscapeOptimum;
    std::vector<std::vector<double>> fitnessLandscapeLinearGradient;
    std::vector<std::vector<double>> fitnessLandscapeGaussStrength;
};
