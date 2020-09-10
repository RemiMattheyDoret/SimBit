template<typename INT>
void FitnessMap::setHabitat(const INT habitat) const
{
    if (T1_isSelection)
    {
        #ifdef DEBUG
        assert(T1_FitnessEffects.size() > habitat);
        #endif

        T1_FitnessEffects_inHabitat = &T1_FitnessEffects[habitat];

        #ifdef DEBUG
        if (T1_isAllSameFitEffect)
        {
            if (T1_isMuliplicitySelection)
                assert(T1_FitnessEffects_inHabitat.size() == 1);
            else
                assert(T1_FitnessEffects_inHabitat.size() == 3);
        } else
        {
            if (T1_isMuliplicitySelection)
                assert(T1_FitnessEffects_inHabitat.size() == SSP->Gmap.T1_nbLoci);
            else
                assert(T1_FitnessEffects_inHabitat.size() == 3 * SSP->Gmap.T1_nbLoci);
        }
        #endif
    }

    if (T2_isSelection)
    {
        #ifdef DEBUG
        assert(T2_FitnessEffects.size() > habitat);
        #endif

        T2_FitnessEffects_inHabitat = &T2_FitnessEffects[habitat];

        #ifdef DEBUG
        if (T2_isAllSameFitEffect)
        {
            assert(T2_FitnessEffects_inHabitat.size() == 1);
        } else
        {
            assert(T2_FitnessEffects_inHabitat.size() == SSP->Gmap.T2_nbLoci);            
        }
        #endif
    }

    if (T56_isSelection)
    {
        #ifdef DEBUG
        assert(T56_FitnessEffects.size() > habitat);
        #endif

        T56_FitnessEffects_inHabitat = &T56_FitnessEffects[habitat];

        #ifdef DEBUG
        if (T56_isAllSameFitEffect)
        {
            if (T56_isMuliplicitySelection)
                assert(T56_FitnessEffects_inHabitat.size() == 1);
            else
                assert(T56_FitnessEffects_inHabitat.size() == 2);
        } else
        {
            if (T56_isMuliplicitySelection)
                assert(T56_FitnessEffects_inHabitat.size() == SSP->Gmap.T1_nbLoci);
            else
                assert(T56_FitnessEffects_inHabitat.size() == 2 * SSP->Gmap.T1_nbLoci);
        }
        #endif
    }
}

template<typename locus_t, typename genotype_t>
uint32_t FitnessMap::getT1Fit(const locus_t locus, const genotype_t genotype) const
{
#ifdef DEBUG
    assert(T1_isSelection);
#endif
    if (T1_isAllSameFitEffect)
    {
        if (T1_isMuliplicitySelection)
        {
            return (*T1_FitnessEffects_inHabitat)[0];
        } else
        {
            return (*T1_FitnessEffects_inHabitat)[genotype];
        }
    } else
    {
        if (T1_isMuliplicitySelection)
        {
            return (*T1_FitnessEffects_inHabitat)[locus];
        } else
        {
            return (*T1_FitnessEffects_inHabitat)[locus * 3 + genotype];
        }
    }
}
