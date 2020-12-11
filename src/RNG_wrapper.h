
using RNG_type = std::mt19937_64; // The code below assumes that the min and max value are 0 and 2^64-1 (which is the case for mersenne twister). So be careful if I use another RNG (there is an assertion in the constructors for security)
constexpr unsigned char NB_BITS_IN_CACHE = sizeof(RNG_type::result_type) * 8;

class RNG_wrapper
{
private:
	int initial_seed;
	RNG_type rng;

	RNG_type::result_type cache;
	unsigned char cache_index = 0;
	
	/*std::uint64_t randomValue8b;
	std::uint64_t randomValue16b;
	std::uint64_t randomValue32b;
	unsigned char nbDrawsLeft8b;
	unsigned char nbDrawsLeft16b;
	unsigned char nbDrawsLeft32b;*/

	//static constexpr double maxValue_8b = 255.0;
	//static constexpr double maxValue_16b = 65535.0;
	//static constexpr double maxValue_32b = 4294967295.0;
	static constexpr double maxValue_64b = 18446744073709551616.0;


public:
	RNG_wrapper(int seed);
	RNG_wrapper();

	std::uint32_t operator()();
	RNG_wrapper operator=(const RNG_wrapper other);

	void set_seed(int seed);
	void set_seed(std::ifstream& file);
	void setRNG(RNG_type& in);
	int getInitialSeed() const;
	RNG_type getRNGState() const;

	RNG_type& getRNG();

	bool get_1b();
	/*std::uint8_t get_8b();
	std::uint16_t get_16b();
	std::uint32_t get_32b();
	std::uint64_t get_64b();*/


	template<class int_type>
	int_type uniform_int_distribution(int_type diff);
	double uniform_real_distribution(double diff);

	template<class int_type>
	int_type uniform_int_distribution(int_type from, int_type to);
	double uniform_real_distribution(double from, double to);

	template<class T>
	double normal(T mean, T sd);

	template<class T>
	size_t poisson(T mean);
};

