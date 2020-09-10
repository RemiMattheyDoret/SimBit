
class TraitSpecificFitnessMap
{
private:
	std::vector<std::vector<double>> fitnessEffects;
	std::vector<double>* fitnessEffectsInHabitat;

public:
	bool isSelection;
	bool isMuliplicitySelection;
	bool isAllSameEffect;
};

class FitnessMap
{
private:

	std::vector<std::vector<double>>   	T1_FitnessEffects;
	std::vector<std::vector<double>>   	T2_FitnessEffects;
	std::vector<std::vector<double>>   	T56_FitnessEffects;

	std::vector<double>*                T1_FitnessEffects_inHabitat;
	std::vector<double>*                T2_FitnessEffects_inHabitat;
	std::vector<double>*                T56_FitnessEffects_inHabitat;

	
    bool                                            T1_isSelection;
    bool                                            T2_isSelection;
    bool                                            T56_isSelection;

    bool                                            T1_isAllSameFitEffect;
    bool                                            T2_isAllSameFitEffect;
    bool                                            T56_isAllSameFitEffect;

    bool                                            T1_isMuliplicitySelection;
    bool                                            T56_isMuliplicitySelection;



    bool                                            T1_isEpistasis;
    std::vector<std::vector<std::vector<T1_locusDescription>>>   T1_Epistasis_LociIndices;
    std::vector<std::vector<std::vector<double>>>                T1_Epistasis_FitnessEffects; //get fitness value with this->FitnessEffects_ForASingleHabitat[habitat][groupOfLoci][fitnessValueIndex], where the fitnessValueIndex is function of the genotype


public:

	template<typename INT>
	void setHabitat(const INT habitat) const;

	template<typename INT>
	uint32_t get(const INT i) const;
};



