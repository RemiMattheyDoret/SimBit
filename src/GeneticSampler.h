class GeneticSampler
{
public:
	GeneticSampler();

	size_t 			get_nbRecombinations();
    size_t 			get_T1_nbMuts();
    size_t 			get_T2_nbMuts();
    size_t 			get_T3_nbMuts();
    size_t 			get_T4_nbMuts(const size_t nbGenerationsInBetween, const size_t left, const size_t right);
    size_t 			get_T56_nbMuts();
    size_t 			get_recombinationPosition();
    size_t 			get_T1_mutationPosition();
    size_t 			get_T2_mutationPosition();
    size_t 			get_T3_mutationPosition();
    size_t 			get_T4_mutationPosition(const size_t left, const size_t right);
    size_t 			get_T56_mutationPosition();


    void 			set_nbRecombinations(const double& rate);
    void 			set_T1_nbMuts(const double& rate);
    void 			set_T2_nbMuts(const double& rate);
    void 			set_T3_nbMuts(const double& rate);
    void 			set_T4_nbMuts(const double& rate);
    void 			set_T56_nbMuts(const double& rate);
    void 			set_recombinationPosition(const std::vector<double>& cumSumRates);
    void 			set_T1_mutationPosition(const std::vector<double>& cumSumRates);
    void 			set_T2_mutationPosition(const std::vector<double>& cumSumRates);
    void 			set_T3_mutationPosition(const std::vector<double>& cumSumRates);
    void 			set_T4_mutationPosition(const std::vector<double>& cumSumRates);
    void 			set_T56_mutationPosition(const std::vector<double>& cumSumRates);


private:

	/* Set the poisson distributions first as I use the rates to assert that the cumulative sum of rate seem correct*/

	std::poisson_distribution<size_t> poisson_nbRecombinations;
	std::poisson_distribution<size_t> poisson_T1_nbMuts;
	std::poisson_distribution<size_t> poisson_T2_nbMuts;
	std::poisson_distribution<size_t> poisson_T3_nbMuts;
	std::poisson_distribution<size_t> poisson_T56_nbMuts;
	double recombinationtotalRate;
	double T1_mut_totalRate;
	double T2_mut_totalRate;
	double T3_mut_totalRate;
	double T4_mut_totalRate;
	double T56_mut_totalRate;

	
	// Non_cst rate that needs walker (no need for T4 walker)
	Walker walker_recombinationPosition;
	Walker walker_T1_mutationPosition;
	Walker walker_T2_mutationPosition;
	Walker walker_T3_mutationPosition; 
	Walker walker_T56_mutationPosition;

	// Non_cst rate that does not need walker
	std::vector<double> cumSumProbs_recombinationPosition;
	std::vector<double> cumSumProbs_T1_mutationPosition;
	std::vector<double> cumSumProbs_T2_mutationPosition;
	std::vector<double> cumSumProbs_T3_mutationPosition;
	std::vector<double> cumSumProbs_T4_mutationPosition;
	std::vector<double> cumSumProbs_T56_mutationPosition;

	// Cst rate that does not need walker
	bool isCstRate_recombinationPosition = false;
	bool isCstRate_T1_mutationPosition = false;
	bool isCstRate_T2_mutationPosition = false;
	bool isCstRate_T3_mutationPosition = false;
	bool isCstRate_T4_mutationPosition = false;
	bool isCstRate_T56_mutationPosition = false;

	// T4 special
	double T4_lastComputedTotalRate = -1.0; // set to -1 to create an error if I mess things up

};
