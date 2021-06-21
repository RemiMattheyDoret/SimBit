

bool PoissonSampler::isMultiplicationOverflow(PSuint_t x, PSuint_t a, PSuint_t b)
{
	if (a != 0 && x / a != b)
	    return true;
	return false;
}



PoissonSampler::PoissonSampler(double lambda)
{
	size_t upTo = 0;
	double sumOfProb = 0.0;
	double desiredSumOfProb = 1.0 - 1e-10;

	PSuint_t k = 0;
	PSuint_t factorialK = 1;
	while (sumOfProb < desiredSumOfProb)
	{
		double x = ((double)pow(lambda, k) * (double)exp(-lambda)) / (double)factorialK;
		probs.push_back({k, x});

		++k;
		auto newFactorialK = factorialK * k
		if (isMultiplicationOverflow(newFactorialK, factorialK, k))
		{
			usingSTL = true;
			STLProbDistribution = ..;

		}
		factorialK = newFactorialK;
	}

	sort and make it cumulative.
}


size_t PoissonSampler::sample()
{
	if (usingSTL)
	{

	} else
	{

	}
}
