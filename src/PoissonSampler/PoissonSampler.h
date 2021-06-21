
typedef __uint128_t PSuint_t;

class PoissonSampler
{
	std::vector<std::pair<size_t,double>> probs;
	bool usingSTL = false;
	.. STLProbDistribution;

	PoissonSampler(double lambda);
	size_t sample();
	bool isMultiplicationOverflow(PSuint_t x, PSuint_t a, PSuint_t b);
};