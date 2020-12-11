struct CIS
{
	double cisEffect;
	unsigned char k;
};


struct PhenotypicEffect
{
    unsigned char phenotypeDimension;
    double magnitude;

    PhenotypicEffect(unsigned a, double b):phenotypeDimension(a), magnitude(b){}
};

struct T7Gene
{
	static std::vector<double> K_values_map;
	
	uint64_t ID;
	std::vector<uint64_t> cisSites;
	std::vector<CIS> cis; // howOthersAffectMe[site].first is cisEffect and .second is k_val
	unsigned char targetCisSite;
	double trans;
	std::vector<PhenotypicEffect> phenotypicEffects;

	uint16_t T7Locus;

	void attemptMutation();

	void mutatecisEffectSD(size_t index);
	void mutateCisTarget(size_t index);
	void mutatetransEffectSD();
	void mutateaddingInteraction();
	void mutateremovingInteraction();
	void mutatek(size_t index);
};