
// I might want to implement a Poisson distribution with this formula (https://stats.stackexchange.com/questions/35658/simple-approximation-of-poisson-cumulative-distribution-in-long-tail/376393) and read bits by bits. It'll be fast for low rates.


RNG_wrapper::RNG_wrapper(int seed)
:
initial_seed(seed),
rng(seed),
cache_index(NB_BITS_IN_CACHE)
//nbDrawsLeft8b(0),
//nbDrawsLeft16b(0),
//nbDrawsLeft32b(0)
{
	assert(RNG_type::max() == maxValue_64b);
	assert(RNG_type::min() == 0);
}


RNG_wrapper::RNG_wrapper()
{
	assert(RNG_type::max() == maxValue_64b);
	assert(RNG_type::min() == 0);
};

RNG_type& RNG_wrapper::getRNG()
{
	return rng;
}

void RNG_wrapper::setRNG(RNG_type& in)
{
	rng = in;
}

std::uint32_t RNG_wrapper::operator()()
{
	return rng();
};


RNG_wrapper RNG_wrapper::operator=(const RNG_wrapper other)
{
	this->initial_seed = other.initial_seed;
	this->rng = other.rng;
	this->cache = other.cache;
	this->cache_index = other.cache_index;
/*	this->randomValue8b = other.randomValue8b;
	this->randomValue16b = other.randomValue16b;
	this->randomValue32b = other.randomValue32b;
	this->nbDrawsLeft1b = other.nbDrawsLeft1b;
	this->nbDrawsLeft8b = other.nbDrawsLeft8b;
	this->nbDrawsLeft16b = other.nbDrawsLeft16b;
	this->nbDrawsLeft32b = other.nbDrawsLeft32b;*/
	return *this;
}



int RNG_wrapper::getInitialSeed() const
{
	return initial_seed;
}


RNG_type RNG_wrapper::getRNGState() const
{
	return rng;
}


void RNG_wrapper::set_seed(int seed)
{
	initial_seed = seed;
	rng = RNG_type(seed);
}


void RNG_wrapper::set_seed(std::ifstream& file)
{
	file >> rng;
    initial_seed = std::numeric_limits<int>::quiet_NaN();
}

	


bool RNG_wrapper::get_1b()
{
	if (cache_index == 0)
	{
		cache = rng();
		cache_index = NB_BITS_IN_CACHE;
	}
	bool r = cache & 1;
	cache >>= 1;
	cache_index--;
	//std::cout << r << " ";
	return r;
}


/*
std::uint8_t RNG_wrapper::get_8b()
{
	if (nbDrawsLeft8b == 0)
	{
		randomValue8b = rng();
		nbDrawsLeft8b = 8;
	}
	std::uint8_t r = randomValue8b & 255;
	randomValue8b >>= 8;
	nbDrawsLeft8b--;
	return r;
}


std::uint16_t RNG_wrapper::get_16b()
{
	if (nbDrawsLeft16b == 0)
	{
		randomValue16b = rng();
		nbDrawsLeft16b = 4;
	}
	std::uint16_t r = randomValue16b & 65535;
	randomValue16b >>= 16;
	nbDrawsLeft16b--;
	return r;
}

std::uint32_t RNG_wrapper::get_32b()
{
	if (nbDrawsLeft32b == 0)
	{
		randomValue32b = rng();
		nbDrawsLeft32b = 2;
	}
	std::uint16_t r = randomValue32b & 4294967295;
	randomValue32b >>= 32;
	nbDrawsLeft32b--;
	return r;
}*/


template<class int_type>
int_type RNG_wrapper::uniform_int_distribution(int_type diff) // diff is not included unlike in the random header logic
{
	assert(diff);
	if (diff == 1)
	{
		return 0;
	} else if (diff == 2)
	{
		return get_1b();
	} else
	{
		//std::cout << "diff = " << diff << "\n";
		std::uniform_int_distribution<int_type> d(0,diff-1);
		//int_type x = d(rng);
		//std::cout << "x = " << x << "\n";
		return d(rng);
		//return rng() % (diff);
	}
	//return (int_type) (uniform_real_distribution((double) to) + 0.5);
}

double RNG_wrapper::uniform_real_distribution(double diff)
{
	/*double rr = rng();
	double r = (double) rr / (maxValue_64b+0.1) * diff;
	std::cout << "-----\n";
	std::cout << "rr = " << rr << "\n";
	std::cout << "diff = " << diff << "\n";
	std::cout << "maxValue_64b = " << maxValue_64b << "\n";
	std::cout << "r = " << r << "\n";
	std::cout << "-----\n";*/

	std::uniform_real_distribution<double> d(0,diff);
	return d(rng);

	//return (double) rng() / (maxValue_64b+0.1) * diff;
}

template<class int_type>
int_type RNG_wrapper::uniform_int_distribution(int_type from, int_type to) // to is not included unlike in the random header filethis logic
{
	std::uniform_int_distribution<int_type> d(from,to-1);
	return d(rng);
	//return from + uniform_int_distribution(to - from);
}

double RNG_wrapper::uniform_real_distribution(double from, double to)
{
	std::uniform_real_distribution<double> d(from,to);
	return d(rng);
}


template<class T>
double RNG_wrapper::normal(T mean, T sd)
{
	std::normal_distribution<double> d(mean,sd);
	return d(rng);
}

template<class T>
size_t RNG_wrapper::poisson(T mean)
{
	std::poisson_distribution<size_t> d(mean);
	return d(rng);
}


