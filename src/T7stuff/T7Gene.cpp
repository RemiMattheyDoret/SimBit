

void T7Gene::mutatecisEffectSD(size_t index)
{
    assert(index < cis.size());
    auto x = GP->rngw.normal(0.0, SSP->T7mutpars.cisEffectSD);
    cis[index].cisEffect += x;
}

void T7Gene::mutatetransEffectSD()
{
    auto x = GP->rngw.normal(0.0, SSP->T7mutpars.transEffectSD);
    trans += x;
}

void T7Gene::mutateCisTarget(size_t index)
{
    if (GP->rngw.get_1b())
    {
        ++targetCisSite;
    } else
    {
        --targetCisSite;
    }
}

void T7Gene::mutatek(size_t index)
{
    assert(K_values_map.size() == 5);
    if (GP->rngw.get_1b())
    {
        if (cis[index].k < 4) ++cis[index].k;
    } else
    {
        if (cis[index].k > 0) --cis[index].k;
    }
}



void T7Gene::attemptMutation()
{
    if (SSP->T7mutpars.cisEffectMu > 0)
    {
        auto nbMuts = GP->rngw.poisson(SSP->T7mutpars.transEffectMu);
        for (size_t muti = 0 ; muti < nbMuts ; ++muti)
            mutatecisEffectSD(GP->rngw.uniform_int_distribution(SSP->T7devpars.nbCisSites));
    }

    if (SSP->T7mutpars.transEffectMu > 0 && GP->rngw.uniform_real_distribution(1.0) < SSP->T7mutpars.transEffectMu)
    {
        mutatetransEffectSD();
    }


    if (SSP->T7mutpars.changeTargetMu > 0)
    {
        auto nbMuts = GP->rngw.poisson(SSP->T7mutpars.changeTargetMu);
        for (size_t muti = 0 ; muti < nbMuts ; ++muti)
            mutateCisTarget(GP->rngw.uniform_int_distribution(cis.size()));
    }

    if (SSP->T7mutpars.kmu > 0)
    {
        auto nbMuts = GP->rngw.poisson(SSP->T7mutpars.kmu);
        for (size_t muti = 0 ; muti < nbMuts ; ++muti)
            mutatek(GP->rngw.uniform_int_distribution(cis.size()));
    }
}

/*
void T7Gene::removeInexistentT7GeneIDIfNeeded(Pop& pop)
{
    if (SSP->Gmap.T7_nbLoci == 0) return;
    if (!(GP->CurrentGeneration %% SSP->lookForInexistentT7GeneIDEveryNGenerations == 0)) return;


    // Clear whoAffectMe from inexistent gene ID
    for (size_t patch_index = 0 ; patch_index < GP->patchNumber ; ++patch_index)
    {
        for (size_t ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ++ind_index)
        {
            for (size_t haplo_index = 0 ; haplo_index < 2 ; ++haplo_index)
            {
                auto& genes = pop.getPatch(patch_index).getInd(ind_index).getHaplo(haplo_index).getT7Genes();
                for (auto it = gene.whoAffectsMe.end() ; it != gene.whoAffectsMe.begin() ; --it)
                {
                    auto itf = std::lower_bund(listGeneIDs.begin(), listGeneIDs.end(), *it);
                    if (itf == geneIDs.end() || *itf != *it)
                    {
                        it = gene.whoAffectsMe.erase(it);
                    }
                }
            }
        }   
    }
}
*/
