
struct ContainerReturnedFromT7Develop
{
    std::vector<std::vector<double>> phenotypeOverTime; // a vector for each dimension
    std::vector<double> timesAtwhichPhenotypeIsSampled;
    ContainerReturnedFromT7Develop(std::vector<std::vector<double>> a, std::vector<double> b):phenotypeOverTime(a), timesAtwhichPhenotypeIsSampled(b) {}
};

struct OneProtEffect
{
	unsigned indexCausal;
	bool isEnhancer;
	double magnitude;
	double K_val;
    std::vector<PhenotypicEffect> phenotypicEffects;

    /*    
	OneProtEffect(unsigned a, bool b, double c, double d, std::vector<PhenotypicEffect>& e):indexCausal(a), isEnhancer(b), magnitude(c), K_val(d), phenotypicEffects(e){}

    OneProtEffect(){}

    OneProtEffect(OneProtEffect& other)
    :
        indexCausal(other.indexCausal),
        isEnhancer(other.isEnhancer),
        magnitude(other.magnitude),
        K_val(other.K_val),
        phenotypicEffects(other.phenotypicEffects)
    {}

    OneProtEffect(OneProtEffect&& other)
    :
        indexCausal(std::move(other.indexCausal)),
        isEnhancer(std::move(other.isEnhancer)),
        magnitude(std::move(other.magnitude)),
        K_val(std::move(other.K_val)),
        phenotypicEffects(std::move(other.phenotypicEffects))
    {}

    OneProtEffect operator=(OneProtEffect other)
    {
        OneProtEffect me(other);
        return me;
    }

    OneProtEffect operator=(OneProtEffect& other)
    {
        OneProtEffect me(other);
        return me;
    }
    */
};


