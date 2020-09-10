void T56_memoryManager::doStuff(Haplotype& haplo)
{
    auto& ntrl = haplo.T5ntrl_Alleles;
    auto& sel = haplo.T5sel_Alleles;

    /////////////////
    // Record info //
    /////////////////
    if (shouldRecordInfo)
    {
        if (SSP->Gmap.T5ntrl_nbLoci && ntrl.size() > maxNtrlSizeLastGeneration)
            maxNtrlSizeLastGeneration = ntrl.size();

        if (SSP->Gmap.T5sel_nbLoci && sel.size() > maxSelSizeLastGeneration)
            maxSelSizeLastGeneration = sel.size();
    }


    ////////////
    // Shrink //
    ////////////
    if (shouldShrink)
    {
        if (SSP->Gmap.T5ntrl_nbLoci)
        {
            if (ntrl.capacity() > 1.1 * (double)ntrl.size()) // if there is some non-negligible loss
            {
                auto oldCapacity = ntrl.capacity();
                ntrl.shrink_to_fit();
                if (oldCapacity > maxNtrlSizeLastGeneration && ntrl.size() < maxNtrlSizeLastGeneration && ntrl.size() > (double)maxNtrlSizeLastGeneration / 1.5)
                {
                    ntrl.reserve(maxNtrlSizeLastGeneration);
                } else
                {
                    ntrl.reserve(1.1 * (double)ntrl.size());
                }
            }
        }

        if (SSP->Gmap.T5sel_nbLoci)
        {
            if (sel.capacity() > 1.1 * (double)sel.size()) // if there is some non-negligible loss
            {
                auto oldCapacity = sel.capacity();
                sel.shrink_to_fit();
                if (oldCapacity > maxNtrlSizeLastGeneration && sel.size() < maxSelSizeLastGeneration && sel.size() > (double)maxSelSizeLastGeneration / 1.5)
                {
                    sel.reserve(maxSelSizeLastGeneration);
                } else
                {
                    sel.reserve(1.1 * (double)sel.size());
                }
            }
        }
    }
}


void T56_memoryManager::setGenerationInfo()
{
    if (attempt_shrinking_every_N_generation != -1)
    {
        // If info was recorded last generation
        if (shouldRecordInfo)
            shouldShrink = true;
        else
            shouldShrink = false;
        
        // Check if it is time to record information
        if (GP->CurrentGeneration % attempt_shrinking_every_N_generation == 0)
            shouldRecordInfo = true;
        else
            shouldRecordInfo = false;
    }
}

template<typename INT>
void T56_memoryManager::set_attempt_shrinking_every_N_generation(INT i)
{
    if (i != -1 && i < 2)
    {
        std::cout << "For option --shrinkT56EveryNGeneration, received a value that is lower than two and differentt from -1 (-1 means no shrinking)\n";
        abort();
    }
    attempt_shrinking_every_N_generation = i;
}
