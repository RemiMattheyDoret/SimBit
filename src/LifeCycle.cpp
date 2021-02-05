/*

 Author: Remi Matthey-Doret

    MIT License

    Copyright (c) 2017 Remi Matthey-Doret

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

 */




void LifeCycle::BREEDING_SELECTION_DISPERSAL(Pop& Offspring_pop, Pop& Parent_pop)
{   
    //if (GP->CurrentGeneration >= 0) std::cout << "Enters in 'BREEDING_SELECTION_DISPERSAL'\n";
    #ifdef DEBUG
    std::cout << "Enters in 'BREEDING_SELECTION_DISPERSAL'\n";
    #endif
    /*
    std::cout << "patchSize when entering LifeCycle::BREEDING_SELECTION_DISPERSAL: ";
    for (auto& e : SSP->patchSize) std::cout << e << " ";
    std::cout << "\n";   
    */
    SSP->T56_memManager.setGenerationInfo();

    //uint32_t nbSwaps = 0;

    //std::cout << "\nBEGIN: Offspring_pop.getPatch(0).getInd(2).getHaplo(0).T4ID = " << Offspring_pop.getPatch(0).getInd(2).getHaplo(0).T4ID << "\n";
    //std::cout << "BEGIN: Parent_pop.getPatch(0).getInd(2).getHaplo(0).T4ID = " << Parent_pop.getPatch(0).getInd(2).getHaplo(0).T4ID << "\n";

    // Calculate fitness
    //SSP->simTracker.prepareT1SitesForFitness(Parent_pop);
    if (SSP->isAnySelection && SSP->selectionOn != 1) Parent_pop.CalculateFitnesses(); // also set Index First Male; Will compute fitness only if it is needed

    // 0. ReSet Dispersal info if patch size may vary and // 2. Compute patch size for next generation
    const std::vector<int>& patchSizeNextGeneration = SSP->dispersalData.setBackwardMigrationIfNeededAndGetNextGenerationPatchSizes(Parent_pop.CumSumFits);
    assert(patchSizeNextGeneration.size() == GP->PatchNumber);
    //std::cout << "patchSizeNextGeneration[0] = " << patchSizeNextGeneration[0] << "\n";
    //std::cout << "SSP->patchSize[0] = " << SSP->patchSize[0] << "\n";
    //SSP->dispersalData.print();


    // 4.5 Genealogy
    if (SSP->genealogy.isTime())
        SSP->genealogy.startNewGeneration();

    // Loop through each patch and through each individual offspring to be created
    //#pragma omp parallel for
    
    if (SSP->SwapInLifeCycle && SSP->selectionOn == 0) // no selection on viability
    {
        findAllParents(Parent_pop, patchSizeNextGeneration); // set PD (of class ParentsData)
    } else
    {
        PD.resizeCloneInfo(patchSizeNextGeneration);
    }


    // Create offsprings
    SSP->TotalpatchSize = 0;
    for (int patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
    {
        SSP->TotalpatchSize += patchSizeNextGeneration[patch_index];
        // Prepare Next Generation CumSumFits and Index First Male.
        //Offspring_pop.prepareNextGenerationAndIndexFirstMale(patch_index, patchSizeNextGeneration);


        for (int offspring_index = 0 ; offspring_index < patchSizeNextGeneration[patch_index] ; ++offspring_index)
        {
            //std::cout << "Creating offspring " << offspring_index << " of patch " << patch_index << "...\n";
            //std::cout << "patchSizeNextGeneration["<<patch_index<<"] = " << patchSizeNextGeneration[patch_index] << "\n";
            
            //std::cout << "LifeCycle line 74 offspring_index = "<<offspring_index<<"\n";
            //std::cout << " patchSizeNextGeneration["<<patch_index<<"] = " <<  patchSizeNextGeneration[patch_index] << "   GP->PatchNumber = "<<GP->PatchNumber<<"\n";

            // Get offspring haplotypes
            Individual& Offspring = Offspring_pop.getPatch(patch_index).getInd(offspring_index);
            Haplotype& Offspring_mHaplo = Offspring.getHaplo(0); // maternally inherited haplotype
            Haplotype& Offspring_fHaplo = Offspring.getHaplo(1); // paternally inherited haplotype

            //std::cout << "Before: Offspring_mHaplo.T4ID = " << Offspring_mHaplo.T4ID << "\n";
            //std::cout << "Offspring_fHaplo.T4ID = " << Offspring_fHaplo.T4ID << "\n";

            //std::cout << "\n\n\n\n\n\nBegin Offspring_mHaplo.nbT56muts() = "<< Offspring_mHaplo.nbT56muts() << "\n";
            //std::cout << "Begin Offspring_fHaplo.nbT56muts() = "<< Offspring_fHaplo.nbT56muts() << "\n";


            redoOffspring:

            // Get parents
            auto&& CD = SSP->SwapInLifeCycle ? PD.couples[patch_index][offspring_index] : findCouple(Parent_pop,patch_index, offspring_index, PD, patchSizeNextGeneration[patch_index]);
            
            auto& mother = Parent_pop.getPatch(CD.mother.patch).getInd(CD.mother.ind);
            auto& father = Parent_pop.getPatch(CD.father.patch).getInd(CD.father.ind);

            //std::cout << "mother.getHaplo0(0).T4ID = " << mother.getHaplo(0).T4ID << "\n";
            //std::cout << "mother.getHaplo0(1).T4ID = " << mother.getHaplo(1).T4ID << "\n";

            //std::cout << "Begin mother.getHaplo(CD.mother.segregationIndex).nbT56muts() = "<< mother.getHaplo(CD.mother.segregationIndex).nbT56muts() << "\n";
            //std::cout << "Begin father.getHaplo(CD.father.segregationIndex).nbT56muts() = "<< father.getHaplo(CD.father.segregationIndex).nbT56muts() << "\n";

            //std::cout << "M: " << CD.mother.patch << " " << CD.mother.ind  << " " << CD.mother.segregationIndex << " | " << CD.father.patch << " " << CD.father.ind << " " << CD.father.segregationIndex << "\n";
            

            // Clear T56 vectors
            if (SSP->Gmap.T56sel_nbLoci || SSP->Gmap.T56ntrl_nbLoci)
            {
                Offspring_mHaplo.clearT56Alleles();
                Offspring_fHaplo.clearT56Alleles();
            }
            

            if (PD.shouldIClone(patch_index, offspring_index))
            {
                assert(CD.mother.patch == CD.father.patch);
                assert(CD.mother.ind == CD.father.ind);
                assert(CD.mother.segregationIndex != CD.father.segregationIndex);
                assert(CD.mother.nbRecs == 0);
                assert(CD.father.nbRecs == 0);

                std::vector<uint32_t> breakpoints = {std::numeric_limits<uint32_t>::max()};
                
                if (SSP->SwapInLifeCycle && CD.mother.nbRecs == 0 && PD.isLastOffspring(CD.mother, patch_index, offspring_index,0))
                {
                    //nbSwaps+=2;
                    //std::cout << "swaps while cloning\n";
                    assert(PD.isLastOffspring(CD.father, patch_index, offspring_index,1));
                    reproduceThroughSwap(mother, Offspring_mHaplo, CD.mother, patch_index);
                    reproduceThroughSwap(father, Offspring_fHaplo, CD.father, patch_index);
                } else
                {
                    reproduceThroughCopy(
                       mother,
                       Offspring_mHaplo,
                       CD.mother,
                       patch_index
                    );

                    reproduceThroughCopy(
                       father,
                       Offspring_fHaplo,
                       CD.father,
                       patch_index
                    );
                }
            } else
            {
                // Mother
                if (SSP->SwapInLifeCycle && CD.mother.nbRecs == 0 && PD.isLastOffspring(CD.mother, patch_index, offspring_index, 0))
                {
                    //nbSwaps++;
                    //std::cout << "swaps mother\n";
                    reproduceThroughSwap(mother, Offspring_mHaplo, CD.mother, patch_index);

                } else
                {
                    //std::cout << "copies mother\n";
                    reproduceThroughCopy(
                       mother,
                       Offspring_mHaplo,
                       CD.mother,
                       patch_index
                    );
                }

                
                // Father
                if (SSP->SwapInLifeCycle && CD.father.nbRecs == 0 && PD.isLastOffspring(CD.father, patch_index, offspring_index,1))
                {
                    //nbSwaps++;
                    //std::cout << "swaps father\n";
                    reproduceThroughSwap(father, Offspring_fHaplo, CD.father, patch_index);

                } else
                {
                    //std::cout << "copies father\n";
                    reproduceThroughCopy(
                       father,
                       Offspring_fHaplo,
                       CD.father,
                       patch_index
                    );
                }
            } // end of if else should I clone

            // Fitness for next generation. Set fitness to zero if the habitat has changed. Could do better here as not all of them will necessarily have habitat specific selection
            if (SSP->Habitats[patch_index] != SSP->Habitats[CD.mother.patch])
            {
                if (SSP->T1_isLocalSelection && SSP->T1_isMultiplicitySelection)
                {
                    Offspring_mHaplo.setAllW_T1(-1.0);
                }
                if (SSP->T2_isLocalSelection && SSP->T2_isSelection)
                {
                    Offspring_mHaplo.setAllW_T2(-1.0);
                }
                if (SSP->T56_isLocalSelection && SSP->T56_isMultiplicitySelection)
                {
                    Offspring_mHaplo.setAllW_T56(-1.0);
                }
            }

            if (SSP->Habitats[patch_index] != SSP->Habitats[CD.father.patch])
            {
                if (SSP->T1_isLocalSelection && SSP->T1_isMultiplicitySelection)
                {
                    Offspring_fHaplo.setAllW_T1(-1.0);
                }
                if (SSP->T2_isLocalSelection && SSP->T2_isSelection)
                {
                    Offspring_fHaplo.setAllW_T2(-1.0);
                }
                if (SSP->T56_isLocalSelection && SSP->T56_isMultiplicitySelection)
                {
                    Offspring_fHaplo.setAllW_T56(-1.0);
                }
            }

            //std::cout << "just before mutation Offspring_mHaplo.nbT56muts() = "<< Offspring_mHaplo.nbT56muts() << "\n";
            //std::cout << "just before mutation Offspring_fHaplo.nbT56muts() = "<< Offspring_fHaplo.nbT56muts() << "\n";
            Mutate(
                Offspring_mHaplo,
                SSP->Habitats[patch_index] // needs the habitat to adjust the fitness
            );

            Mutate(
                Offspring_fHaplo,
                SSP->Habitats[patch_index] // needs the habitat to adjust the fitness
            );
            //std::cout << "line 235: Offspring_mHaplo.T4ID = " << Offspring_mHaplo.T4ID << "\n";
            //std::cout << "just after mutation Offspring_mHaplo.nbT56muts() = "<< Offspring_mHaplo.nbT56muts() << "\n";
            //std::cout << "just after mutation Offspring_fHaplo.nbT56muts() = "<< Offspring_fHaplo.nbT56muts() << "\n";
            

            // Selection on viability
            if (SSP->selectionOn != 0)
            {
                double fitness = Offspring_pop.getPatch(patch_index).getInd(offspring_index).CalculateFitness(patch_index);
                if (GP->rngw.uniform_real_distribution(1.0) > fitness)
                {
                    CD = findCouple(Parent_pop, patch_index, offspring_index, PD, patchSizeNextGeneration[patch_index]); // Will modify cloneInfo too. 
                    if (SSP->SwapInLifeCycle)
                    {
                        PD.couples[patch_index][offspring_index] = CD;
                    } 
                    
                    goto redoOffspring;
                }
            }

            SSP->T56_memManager.doStuff(Offspring_mHaplo);
            SSP->T56_memManager.doStuff(Offspring_fHaplo);


            // 4.5 Genealogy. Will do something only if the last isTime was true
            SSP->genealogy.addOffspringIfIsTime(patch_index, offspring_index, CD.mother.patch, CD.mother.ind, CD.father.patch, CD.father.ind);

            //std::cout << "End Offspring_mHaplo.nbT56muts() = "<< Offspring_mHaplo.nbT56muts() << "\n";
            //std::cout << "End Offspring_fHaplo.nbT56muts() = "<< Offspring_fHaplo.nbT56muts() << "\n";

            //std::cout << "Begin mother.getHaplo(CD.mother.segregationIndex).nbT56muts() = "<< mother.getHaplo(CD.mother.segregationIndex).nbT56muts() << "\n";
            //std::cout << "Begin father.getHaplo(CD.father.segregationIndex).nbT56muts() = "<< father.getHaplo(CD.father.segregationIndex).nbT56muts() << "\n\n\n\n\n";
        }
    }
    
    
    

    // 4.5 Genealogy. Will do something only if the last isTime was true
    //std::cout << "line 184\n";
    SSP->genealogy.printBufferIfIsTime();
    //std::cout << "line 186\n";

    // Check if species is extinct and set patch size
    if (SSP->fecundityForFitnessOfOne != -1.0)
    {
        assert(SSP->whenDidExtinctionOccur == -1); // assert it has nt gone extinct yet
        bool didGetExtinct = true;
        for (int patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
        {
            SSP->patchSize[patch_index] = patchSizeNextGeneration[patch_index];
            if (patchSizeNextGeneration[patch_index] > 0)
            {
                didGetExtinct = false;
            }
        }
        if (didGetExtinct)
        {
            SSP->whenDidExtinctionOccur = GP->CurrentGeneration;
        }
    }

    //std::cout << "END0: Offspring_pop.getPatch(0).getInd(2).getHaplo(0).T4ID = " << Offspring_pop.getPatch(0).getInd(2).getHaplo(0).T4ID << "\n";
    //std::cout << "END0: Parent_pop.getPatch(0).getInd(2).getHaplo(0).T4ID = " << Parent_pop.getPatch(0).getInd(2).getHaplo(0).T4ID << "\n";
    

    if (SSP->whenDidExtinctionOccur == -1 && SSP->Gmap.T4_nbLoci > 0)
    {
        SSP->T4Tree.simplify_ifNeeded(Offspring_pop);
    }

    if (SSP->whenDidExtinctionOccur == -1 && SSP->Gmap.T56_nbLoci)
    {
        Offspring_pop.toggleT56MutationsIfNeeded();
    }

    /*
    std::cout << "patchSize when Exiting LifeCycle::BREEDING_SELECTION_DISPERSAL: ";
    for (auto& e : SSP->patchSize) std::cout << e << " ";
    std::cout << "\n";  
    */

    //std::cout << "nbSwaps = " << nbSwaps << "\n";
}


void LifeCycle::reproduceThroughSwap(Individual& parent, Haplotype& offspringHaplotype, HaplotypeData& parentHaploData, int& patch_index)
{
    auto& parentalHaplo = parent.getHaplo(parentHaploData.segregationIndex);

    auto parentalT4ID = parentalHaplo.T4ID; // only used for T4
    parentalHaplo.swap(offspringHaplotype);

    if (SSP->Gmap.T4_nbLoci > 0)
    {
        //std::cout <<"\tLIFECYCLE - CLONING!!! b = NA: SSP->Gmap.T4_nbLoci = "<<SSP->Gmap.T4_nbLoci<<"\n";
        std::vector<uint32_t> RecPos = {std::numeric_limits<uint32_t>::max()};
        offspringHaplotype.T4ID = SSP->T4Tree.addHaplotype(RecPos, {parentalT4ID, std::numeric_limits<std::uint32_t>::max()}, patch_index);
    }
}


/*int LifeCycle::recombination_nbRecs()
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'recombination_nbRecs'\n";
#endif       
    int nbRecs;
    
    if (SSP->TotalRecombinationRate != 0)
    {
        nbRecs = SSP->rpois_nbRecombination(GP->rngw.getRNG()); // For the moment does not allow variation in recombination rate.
    } else
    {
        nbRecs = 0;
    }
    return nbRecs;
}*/

void LifeCycle::recombination_RecPositions(int& nbRecs, Individual& parent)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'recombination_RecPositions'\n";
#endif  
    // Find positions (here a position is a separation between bit or between a bit and a byte or anything of interest) where recombination occurs

    if (nbRecs == 0)
    {
        recombination_breakpoints.resize(1);
        recombination_breakpoints[0] = std::numeric_limits<uint32_t>::max();
    } else
    {
        recombination_breakpoints.resize(0);
        //recombination_breakpoints.reserve(nbRecs+1);

        //std::cout << "nbRecs = "<< nbRecs << std::endl;
        std::vector<int> RecPositions;
        RecPositions.reserve(nbRecs);
        
        for (int i = 0 ; i < nbRecs ; i++)
        {
            RecPositions.push_back(SSP->geneticSampler.get_recombinationPosition());
        }

            /*
            while(true)
            {
                if (SSP->RecombinationRate.size() == 1) // Constant value
                {
                    RecPosition = SSP->runiform_int_ForRecPos(GP->rngw.getRNG()); // is bounded between 0 and size - 1.
                } else // Not constant
                {
                    assert(SSP->RecombinationRate.size() == SSP->Gmap.TotalNbLoci - 1);

                    rnd = SSP->runiform_double_ForRecPos(GP->rngw.getRNG()); // is bounded between 0 and total_recombinationRate.
                    
                    // binary search SSP->RecombinationRate must be cumulative
                   RecPosition = distance(
                                    SSP->RecombinationRate.begin(),
                                    std::upper_bound(
                                        SSP->RecombinationRate.begin(),
                                        SSP->RecombinationRate.end(),
                                        rnd
                                    )
                                );

                }
                assert(RecPosition >= 0 && RecPosition < SSP->Gmap.TotalNbLoci - 1);

                
                if (SSP->recRateOnMismatch_bool)
                {
                    int fromT1Locus = std::max(0,RecPosition - SSP->recRateOnMismatch_halfWindow);
                    int toT1Locus   = std::min(SSP->Gmap.T1_nbLoci,RecPosition + SSP->recRateOnMismatch_halfWindow);
                    //std::cout << fromT1Locus << " " << toT1Locus << " | ";
                    double sumAtExponentTerm = 0.0;
                    for (int T1Locus = fromT1Locus ; T1Locus < toT1Locus ; T1Locus++)
                    {
                        if (parent.getHaplo(0).getT1_Allele(T1Locus) != parent.getHaplo(1).getT1_Allele(T1Locus))
                        {
                            sumAtExponentTerm += 1 / abs(T1Locus - RecPosition);
                        }
                    }
                    double probabilityToKeepRecEvent = exp(SSP->recRateOnMismatch_factor * sumAtExponentTerm);
                    assert(probabilityToKeepRecEvent >= 0.0 & probabilityToKeepRecEvent <= 1.0);
                    //std::cout << "probabilityToKeepRecEvent = " << probabilityToKeepRecEvent << std::endl;

                    if (probabilityToKeepRecEvent > GP->rngw.uniform_real_distribution(1.0));
                    {
                        RecPositions.push_back(RecPosition);
                        break;
                    } // else keep looping over the while(true) loop
                } else
                {
                    RecPositions.push_back(RecPosition);
                    break;
                }
            }*/
        

        // sort positions
        sort( RecPositions.begin(), RecPositions.end() );
        
        // Keep only odd number of recombinations
        assert(RecPositions.size() == nbRecs);
        int nbRepeats = 1;
        
        for (int i = 1 ; i < RecPositions.size() ; i++)
        {
            if (RecPositions[i-1] == RecPositions[i])
            {
                nbRepeats++;
            } else
            {
                if (nbRepeats % 2)
                {
                    recombination_breakpoints.push_back(RecPositions[i-1]);
                }
                nbRepeats = 1;
            }
        }
        
        if (nbRepeats % 2)
        {
            //std::cout << RecPositions.back() << "\n";
            assert(RecPositions.back() >= 0 && RecPositions.back() < SSP->Gmap.TotalNbLoci);
            recombination_breakpoints.push_back(RecPositions.back());
        }
        
        // Add segregation   
        //if (SSP->ChromosomeBoundaries.size() != 0) assert(SSP->ChromosomeBoundaries.back() != std::numeric_limits<uint32_t>::max());
        //if (breakpoints.size()!=0) assert(breakpoints.back() != std::numeric_limits<uint32_t>::max());
        /*
        for (int b : SSP->ChromosomeBoundaries)
        {
            if (GP->rngw.get_1b())
            {
                recombination_breakpoints.insert(
                    std::upper_bound(recombination_breakpoints.begin(), recombination_breakpoints.end(), b), 
                    b
                );
            }
        }
        */

        // add upper bound
        recombination_breakpoints.push_back(std::numeric_limits<uint32_t>::max());
    }
}

void LifeCycle::copyOver(Individual& parent, Haplotype& TransmittedChrom, int segregationIndex)
{
    assert(segregationIndex == 0 || segregationIndex == 1 );
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'copyOver'\n";
#endif
    if (recombination_breakpoints.size()==1)
    {
        
        //if (GP->CurrentGeneration >= 50) std::cout << "parent.getHaplo(segregationIndex).getW_T1(0) = " << parent.getHaplo(segregationIndex).getW_T1(0) << "\n";
        TransmittedChrom = parent.getHaplo(segregationIndex); // ID should be correctly incremented here
        //if (GP->CurrentGeneration >= 50) std::cout << "TransmittedChrom.getW_T1(0) = " << TransmittedChrom.getW_T1(0) << "\n";
        
    } else if (recombination_breakpoints.size()>1)
    {
        
        //std::cout << "breakpoints.size() " << breakpoints.size() << "\n";
        #ifdef DEBUG
        long long Safety_T1_Absolutefrom = std::numeric_limits<uint32_t>::max();
        long long Safety_T1_Absoluteto   = -1;
        
        long long Safety_T2_Absolutefrom = std::numeric_limits<uint32_t>::max();
        long long Safety_T2_Absoluteto   = -1;
        
        long long Safety_T3_Absolutefrom = std::numeric_limits<uint32_t>::max();
        long long Safety_T3_Absoluteto   = -1;

        long long Safety_T4_Absolutefrom = std::numeric_limits<uint32_t>::max();
        long long Safety_T4_Absoluteto   = -1;
        
        long long Safety_T56ntrl_Absolutefrom = std::numeric_limits<uint32_t>::max();
        long long Safety_T56ntrl_Absoluteto   = -1;

        long long Safety_T56sel_Absolutefrom = std::numeric_limits<uint32_t>::max();
        long long Safety_T56sel_Absoluteto   = -1;
        #endif
        
        uint32_t T1_from = 0;
        uint32_t T2_from = 0;
        uint32_t T3_from = 0;
        uint32_t T4_from = 0;
        uint32_t T56ntrl_from = 0;
        uint32_t T56sel_from = 0;
        uint32_t T1_to = 0;
        uint32_t T2_to = 0;
        uint32_t T3_to = 0;
        uint32_t T4_to = 0;
        uint32_t T56ntrl_to = 0;
        uint32_t T56sel_to = 0;
        uint32_t fitnessMapIndexFrom = 0;

        int haplo_index = segregationIndex; // There is really no reason for copying segregationIndex. I could just use segregationIndex all the way through...but he.... whatever!

        
        /*std::cout << "\t";
        for (auto& b : breakpoints)
            std::cout << b << " ";
        std::cout << std::endl;*/
        
        

        for (auto b : recombination_breakpoints)
        {
            
            // Get positions for the two traits
            if (b == std::numeric_limits<uint32_t>::max())
            {
                
                // Copy is over the range [from, to)
                T1_to = SSP->Gmap.T1_nbLoci;
                T2_to = SSP->Gmap.T2_nbLoci;
                T3_to = SSP->Gmap.T3_nbLoci;
                T4_to = SSP->Gmap.T4_nbLoci;
                T56ntrl_to = SSP->Gmap.T56ntrl_nbLoci;
                T56sel_to = SSP->Gmap.T56sel_nbLoci;
                b = SSP->Gmap.TotalNbLoci - 1;
                
            } else
            {
                assert(b >= 0);
                // because copy is over the range [from, to) and if recombination happens at position 4, then the 4th element belongs to the previous stretch of DNA while the 5th element belongs to the next stretch of DNA.
                
                T1_to = SSP->Gmap.FromLocusToNextT1Locus(b);
                T2_to = SSP->Gmap.FromLocusToNextT2Locus(b);
                T3_to = SSP->Gmap.FromLocusToNextT3Locus(b);
                T4_to = SSP->Gmap.FromLocusToNextT4Locus(b);
                T56ntrl_to = SSP->Gmap.FromLocusToNextT56ntrlLocus(b);
                T56sel_to = SSP->Gmap.FromLocusToNextT56selLocus(b);
                
            }

            #ifdef DEBUG
            
            assert(T1_to <= SSP->Gmap.T1_nbLoci);
            assert(T2_to <= SSP->Gmap.T2_nbLoci);
            assert(T3_to <= SSP->Gmap.T3_nbLoci);
            assert(T4_to <= SSP->Gmap.T4_nbLoci);
            assert(T56ntrl_to <= SSP->Gmap.T56ntrl_nbLoci);
            assert(T56sel_to <= SSP->Gmap.T56sel_nbLoci);
            
            #endif

            //////////////////
            // Set Fitnesses//
            //////////////////


            if (SSP->T1_isMultiplicitySelection || SSP->T2_isSelection || SSP->T56_isMultiplicitySelection)
            {  

                //// Find from to and what to recompute if anything

                int fitnessMapIndexToRecompute;
                int fitnessMapIndexTo;
            
                // Before
                int fitnessMapIndexBefore = SSP->FromLocusToFitnessMapIndex[b];

                // After
                int fitnessMapIndexAfter;
                if (b == SSP->Gmap.TotalNbLoci-1)
                {
                    fitnessMapIndexAfter = SSP->NbElementsInFitnessMap;
                    assert(SSP->FromLocusToFitnessMapIndex[b]+1 == fitnessMapIndexAfter);
                } else
                {
                    fitnessMapIndexAfter = SSP->FromLocusToFitnessMapIndex[b+1];
                }

                // Figure out if there is a need to recompute a fitnessMap value (if the crossover happened in between two fitnessMapIndex, then there is no need to recompute anything)
                if (fitnessMapIndexBefore == fitnessMapIndexAfter && fitnessMapIndexBefore >= fitnessMapIndexFrom) // ' && fitnessMapIndexBefore >= fitnessMapIndexFrom' is actually not needed normally
                {
                    fitnessMapIndexToRecompute = fitnessMapIndexBefore; // Must force recalculation of fitness
                    assert(b < SSP->Gmap.TotalNbLoci-1);
                } else
                {
                    fitnessMapIndexToRecompute = -1; // no need to recompute fitness
                }

                // Figure out end of things to copy over (fitnessMapIndexTo, excluded)
                fitnessMapIndexTo = fitnessMapIndexBefore + 1; //  excluded. (fitnessMapIndexBefore + 1 will rarely equal fitnessMapIndexAfter but it can.) 


                // Figure out beginning of things to copy over (fitnessMapIndexFrom, included)
                assert(fitnessMapIndexFrom <= fitnessMapIndexTo);
                assert(fitnessMapIndexToRecompute == -1 || (fitnessMapIndexToRecompute >= fitnessMapIndexFrom && fitnessMapIndexToRecompute < fitnessMapIndexTo));
                /*if (fitnessMapIndexFrom > fitnessMapIndexTo)
                {
                    assert(fitnessMapIndexFrom == fitnessMapIndexTo + 1);
                    fitnessMapIndexFrom = fitnessMapIndexAfter;
                }*/


                            
                
                //// Set fitnesses for each type of trait
                if (fitnessMapIndexFrom < fitnessMapIndexTo)
                {
                    if (SSP->T1_isMultiplicitySelection)
                    {
                        
                        /*
                        std::cout << "fitnessMapIndexToRecompute = " << fitnessMapIndexToRecompute << "\n";
                        std::cout << "fitnessMapIndexFrom = " << fitnessMapIndexFrom << "\n";
                        std::cout << "fitnessMapIndexTo = " << fitnessMapIndexTo << "\n";
                        */

                        // Copy fitnesses that can be copied
                        for (int fitnessMapIndex = fitnessMapIndexFrom; fitnessMapIndex < fitnessMapIndexTo ; fitnessMapIndex++)
                        {
                            //std::cout << "Assgning FMI "<<fitnessMapIndex<< " to ";
                            if (fitnessMapIndex != fitnessMapIndexToRecompute)
                            {
                                //std::cout << parent.getHaplo(haplo_index).getW_T1(fitnessMapIndex) << "\n";
                                TransmittedChrom.setW_T1(
                                    parent.getHaplo(haplo_index).getW_T1(fitnessMapIndex),
                                    fitnessMapIndex
                                );
                            } else
                            {
                                //std::cout << "-1.0" << "\n";
                                TransmittedChrom.setW_T1(-1.0,fitnessMapIndexToRecompute);
                            }
                        }
                    }

                    if (SSP->T2_isSelection)
                    {
                        // Copy fitnesses that can be copied
                        for (int fitnessMapIndex = fitnessMapIndexFrom; fitnessMapIndex < fitnessMapIndexTo ; fitnessMapIndex++)
                        {
                            if (fitnessMapIndex != fitnessMapIndexToRecompute)
                            {
                                //std::cout << "Assgning FMI = " << fitnessMapIndex << "\n";
                                TransmittedChrom.setW_T2(
                                    parent.getHaplo(haplo_index).getW_T2(fitnessMapIndex),
                                    fitnessMapIndex
                                );
                            } else
                            {
                                TransmittedChrom.setW_T2(-1.0, fitnessMapIndexToRecompute);
                            }
                        }
                    }

                    if (SSP->T56_isMultiplicitySelection)
                    {   

                        // Copy fitnesses that can be copied
                        for (int fitnessMapIndex = fitnessMapIndexFrom; fitnessMapIndex < fitnessMapIndexTo ; fitnessMapIndex++)
                        {
                            if (fitnessMapIndex != fitnessMapIndexToRecompute)
                            {
                                //std::cout << "Assigning FMI = " << fitnessMapIndex << "\n";
                                TransmittedChrom.setW_T56(
                                    parent.getHaplo(haplo_index).getW_T56(fitnessMapIndex),
                                    fitnessMapIndex
                                );
                            } else
                            {
                                TransmittedChrom.setW_T56(-1.0, fitnessMapIndexToRecompute);
                            }
                        }
                    }
                }
                    

                //// Set fitnessMapIndexFrom for next breakpoint
                fitnessMapIndexFrom = fitnessMapIndexTo;
                
            }


            //////////
            // Copy //
            //////////

            //std::cout << "In copyOver: T1_from = "<< T1_from <<" T1_to = "<< T1_to << std::endl;
            //std::cout << "In copyOver: T2_from = "<< T2_from <<" T2_to = "<< T2_to << std::endl;

            // Copy for T1. from included, to excluded
            if (T1_from < T1_to)
            {
                
                assert(T1_to > 0 && T1_to <= SSP->Gmap.T1_nbLoci + 1);
                TransmittedChrom.copyIntoT1(T1_from, T1_to, parent.getHaplo(haplo_index));
                #ifdef DEBUG
                Safety_T1_Absolutefrom = myMin(Safety_T1_Absolutefrom, T1_from);
                Safety_T1_Absoluteto = myMax(Safety_T1_Absoluteto, T1_to);
                #endif
            }       
            
            // Copy for T2. from included, to excluded
            if (T2_from < T2_to)
            {
                
                assert(T2_to > 0 && T2_to <= SSP->Gmap.T2_nbLoci + 1);
                TransmittedChrom.copyIntoT2(T2_from, T2_to, parent.getHaplo(haplo_index));
                #ifdef DEBUG
                Safety_T2_Absolutefrom = myMin(Safety_T2_Absolutefrom, T2_from);
                Safety_T2_Absoluteto = myMax(Safety_T2_Absoluteto, T2_to);
                #endif
            }

            // Copy for T3. from included, to excluded
            if (T3_from < T3_to)
            {
                assert(T3_to > 0 && T3_to <= SSP->Gmap.T3_nbLoci + 1);
                TransmittedChrom.copyIntoT3(T3_from, T3_to, parent.getHaplo(haplo_index));
                #ifdef DEBUG
                Safety_T3_Absolutefrom = myMin(Safety_T3_Absolutefrom, T3_from);
                Safety_T3_Absoluteto = myMax(Safety_T3_Absoluteto, T3_to);
                #endif
            }

            // Copy for T4. from included, to excluded
            if (T4_from < T4_to)
            {
                // nothing to do
                
                assert(T4_to > 0 && T4_to <= SSP->Gmap.T4_nbLoci + 1);
                //std::cout <<"\tLIFECYCLE!! b = "<<b<<" T4_from = "<<T4_from<<" T4_to = "<<T4_to<<"\n";
                #ifdef DEBUG
                Safety_T4_Absolutefrom = myMin(Safety_T4_Absolutefrom, T4_from);
                Safety_T4_Absoluteto = myMax(Safety_T4_Absoluteto, T4_to);
                #endif
            }

            // Copy for T56ntrl. from included, to excluded
            if (T56ntrl_from < T56ntrl_to)
            {
                assert(T56ntrl_to > 0 && T56ntrl_to <= SSP->Gmap.T56ntrl_nbLoci + 1);
                
                
                TransmittedChrom.copyIntoT56ntrl(T56ntrl_from, T56ntrl_to, parent.getHaplo(haplo_index));
                #ifdef DEBUG
                Safety_T56ntrl_Absolutefrom = myMin(Safety_T56ntrl_Absolutefrom, T56ntrl_from);
                Safety_T56ntrl_Absoluteto = myMax(Safety_T56ntrl_Absoluteto, T56ntrl_to);
                #endif
            }

            // Copy for T56sel. from included, to excluded
            if (T56sel_from < T56sel_to)
            {
                
                assert(T56sel_to > 0 && T56sel_to <= SSP->Gmap.T56sel_nbLoci + 1);
                
                TransmittedChrom.copyIntoT56sel(T56sel_from, T56sel_to, parent.getHaplo(haplo_index));
                #ifdef DEBUG
                Safety_T56sel_Absolutefrom = myMin(Safety_T56sel_Absolutefrom, T56sel_from);
                Safety_T56sel_Absoluteto = myMax(Safety_T56sel_Absoluteto, T56sel_to);
                #endif
            }
            
            // Set new 'from'
            T1_from = T1_to;
            T2_from = T2_to;
            T3_from = T3_to;
            T4_from = T4_to;
            T56ntrl_from = T56ntrl_to;
            T56sel_from = T56sel_to;
            
            // Switch parent chromosome
            haplo_index = !haplo_index;

            
        }

#ifdef DEBUG
        
        // Security Checks
        if (SSP->Gmap.T1_nbChars)
        {
            assert(T1_to == SSP->Gmap.T1_nbLoci);
            assert(Safety_T1_Absoluteto == SSP->Gmap.T1_nbLoci);
            assert(Safety_T1_Absolutefrom == 0);
        }

        if (SSP->Gmap.T2_nbLoci)
        {
            assert(T2_to == SSP->Gmap.T2_nbLoci);
            assert(Safety_T2_Absoluteto == SSP->Gmap.T2_nbLoci);
            assert(Safety_T2_Absolutefrom == 0);
        }

        if (SSP->Gmap.T3_nbLoci)
        {
            assert(T3_to == SSP->Gmap.T3_nbLoci);
            assert(Safety_T3_Absoluteto == SSP->Gmap.T3_nbLoci);
            assert(Safety_T3_Absolutefrom == 0);
        }

        if (SSP->Gmap.T4_nbLoci)
        {
            assert(T4_to == SSP->Gmap.T4_nbLoci);
            assert(Safety_T4_Absoluteto == SSP->Gmap.T4_nbLoci);
            assert(Safety_T4_Absolutefrom == 0);
        }

        if (SSP->Gmap.T56ntrl_nbLoci)
        {
            assert(T56ntrl_to == SSP->Gmap.T56ntrl_nbLoci);
            assert(Safety_T56ntrl_Absoluteto == SSP->Gmap.T56ntrl_nbLoci);
            assert(Safety_T56ntrl_Absolutefrom == 0);
        }

        if (SSP->Gmap.T56sel_nbLoci)
        {
            assert(T56sel_to == SSP->Gmap.T56sel_nbLoci);
            assert(Safety_T56sel_Absoluteto == SSP->Gmap.T56sel_nbLoci);
            assert(Safety_T56sel_Absolutefrom == 0);
        }
        
#endif

        
    } else
    {
        std::cout << "Internal error. 'breakpoints' has size of zero. It should always contain at least the element std::numeric_limits<uint32_t>::max()\n";
        abort();
    }
    
#ifdef DEBUG
    TransmittedChrom.assertT5orderAndUniqueness();
#endif
    
}

void LifeCycle::reproduceThroughCopy(Individual& parent, Haplotype& TransmittedChrom, HaplotypeData& parentData, int& patch_index)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'reproduceThroughCopy'\n";
#endif     

    //std::cout << "TransmittedChrom.T4ID = " << TransmittedChrom.T4ID << "\n";
    //std::cout << "parent.getHaplo(parentData.segregationIndex).T4ID = " << parent.getHaplo(parentData.segregationIndex).T4ID<< "\n";
    //std::cout << "parent.getHaplo(!parentData.segregationIndex).T4ID = " << parent.getHaplo(!parentData.segregationIndex).T4ID<< "\n";
    
    recombination_RecPositions(parentData.nbRecs, parent); // parent is used if recRateOnMismatch_bool

    // copy from parents chromosomes to TransmittedChrom. The function also set the new values for W_T1 and W_T2
    copyOver(parent, TransmittedChrom, parentData.segregationIndex); // This will copy data from parents to offspring so it must always run

    if (SSP->Gmap.T4_nbLoci > 0)
    {
        TransmittedChrom.T4ID = SSP->T4Tree.addHaplotype(recombination_breakpoints, {parent.getHaplo(parentData.segregationIndex).T4ID, parent.getHaplo(!parentData.segregationIndex).T4ID}, patch_index);
    }
}



void LifeCycle::Mutate(Haplotype& TransmittedChrom, int Habitat)
{
    if (SSP->T1_Total_Mutation_rate>0) Mutate_T1(TransmittedChrom, Habitat);
    if (SSP->T2_Total_Mutation_rate>0) Mutate_T2(TransmittedChrom, Habitat);
    if (SSP->T3_Total_Mutation_rate>0) Mutate_T3(TransmittedChrom);
    if (SSP->T56_Total_Mutation_rate>0) Mutate_T56(TransmittedChrom, Habitat);
    if (SSP->T7mutpars.totalMutationRatePerGene>0) Mutate_T7(TransmittedChrom);
}

void LifeCycle::Mutate_T7(Haplotype& TransmittedChrom)
{
    // deletion
    if (SSP->T7mutpars.deletion > 0 && TransmittedChrom.nbT7Genes())
    {
        auto nbMuts = GP->rngw.poisson(SSP->T7mutpars.deletion * TransmittedChrom.nbT7Genes());
        if (nbMuts >= TransmittedChrom.nbT7Genes())
        {
            TransmittedChrom.clearT7Genes();
        } else
        {
            for (size_t muti = 0 ; muti < nbMuts ; ++muti)
            {
                auto pos = GP->rngw.uniform_int_distribution(TransmittedChrom.nbT7Genes());
                TransmittedChrom.removeT7Gene(pos);
            }
        }
    }

    // Duplications
    if (SSP->T7mutpars.duplication > 0 && TransmittedChrom.nbT7Genes() > 0 && TransmittedChrom.nbT7Genes() < SSP->Gmap.T7_nbLoci)
    {
        auto nbMuts = GP->rngw.poisson(SSP->T7mutpars.duplication * TransmittedChrom.nbT7Genes());
        for (size_t muti = 0 ; muti < nbMuts ; ++muti)
        {
            auto pos = GP->rngw.uniform_int_distribution(TransmittedChrom.nbT7Genes());
            TransmittedChrom.duplicateT7Gene(pos);
        }
    }

    // Normal mutations
    for (size_t geneIndex = 0 ; geneIndex < TransmittedChrom.nbT7Genes() ; ++geneIndex)
    {
        TransmittedChrom.getT7_Allele(geneIndex).attemptMutation();
    }
}


void LifeCycle::Mutate_T1(Haplotype& TransmittedChrom, int Habitat)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'Mutate_T1'\n";
#endif     
    int nbMuts = SSP->geneticSampler.get_T1_nbMuts();
    //std::cout << "nbMuts = " << nbMuts << "\n";
    
    for (int i = 0 ; i < nbMuts ; i++)
    {
        auto MutPosition = SSP->geneticSampler.get_T1_mutationPosition();
        // Make the mutation
        TransmittedChrom.mutateT1_Allele(MutPosition, Habitat);         // toggle bit
        // TransmittedChrom.setT1_AlleleToOne(byte_index,bit_index,);   // set bit to one
    }
}


void LifeCycle::Mutate_T2(Haplotype& TransmittedChrom, int Habitat)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'Mutate_T2'\n";
#endif       
    int nbMuts = SSP->geneticSampler.get_T2_nbMuts();
    
    for (int i = 0 ; i < nbMuts ; i++)
    {
        auto MutPosition = SSP->geneticSampler.get_T2_mutationPosition();
        
        // Make the mutation
        //std::cout << "MutPosition = " <<  MutPosition << "    SSP->Gmap.T2_nbLoci = " << SSP->Gmap.T2_nbLoci << "\n";
        assert(MutPosition < SSP->Gmap.T2_nbLoci);
        TransmittedChrom.AddMutT2_Allele(MutPosition, Habitat);  // add 1 and changes T2_W
    }
}

void LifeCycle::Mutate_T3(Haplotype& TransmittedChrom)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'Mutate_T3'\n";
#endif       
    int nbMuts = SSP->geneticSampler.get_T3_nbMuts();
    //std::cout << "nbMuts = " << nbMuts << "\n";
    for (int i = 0 ; i < nbMuts ; i++)
    {
        auto MutPosition = SSP->geneticSampler.get_T3_mutationPosition();

        // Make the mutation
        //std::cout << "MutPosition = " <<  MutPosition << "    SSP->Gmap.T3_nbLoci = " << SSP->Gmap.T3_nbLoci << "\n";
        assert(MutPosition < SSP->Gmap.T3_nbLoci);   
        TransmittedChrom.mutateT3_Allele(MutPosition);  // add 1 and changes T2_W   
    }  
}

/*void LifeCycle::Mutate_T56(Haplotype& TransmittedChrom, int Habitat)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'Mutate_T56'\n";
#endif     
    int nbMuts = SSP->T56_rpois_nbMut(GP->rngw.getRNG());
    //std::cout << "nbMuts = " << nbMuts << "\n";
    std::vector<int> T5ntrlMutations;
    if (SSP->Gmap.T5ntrl_nbLoci)
        T5ntrlMutations.reserve(nbMuts);
    
    for (int i = 0 ; i < nbMuts ; i++)
    {
        //std::cout << "nbMuts = " << nbMuts << "\n";
        int MutPosition;
        double rnd;
        // Find Position
        //std::cout << "SSP->T56_MutationRate.size() = " << SSP->T56_MutationRate.size() << "\n";
        if (SSP->T56_MutationRate.size() == 1)
        {
            MutPosition = SSP->T56_runiform_int_ForMutPos(GP->rngw.getRNG());
        } else
        {
            rnd = SSP->T56_runiform_double_ForMutPos(GP->rngw.getRNG());
            
            // binary search
            MutPosition = distance(SSP->T56_MutationRate.begin(),
                                   std::upper_bound(SSP->T56_MutationRate.begin(), SSP->T56_MutationRate.end(), rnd)
                                   );
        }
        
        // Make the mutation - toggle bit
        //std::cout << "MutPosition = " << MutPosition << "\n";

        auto& locusGender = SSP->FromT56LocusToT56genderLocus[MutPosition];
        //std::cout << "locusGender.second = " << locusGender.second << "\n";
        if (locusGender.first)
        {
            if (SSP->Gmap.isT56ntrlCompress)
            {
                TransmittedChrom.mutateT56ntrl_Allele(locusGender.second);
            } else
            {
                auto it = std::lower_bound(T5ntrlMutations.begin(), T5ntrlMutations.end(), locusGender.second);
                if (it == T5ntrlMutations.end())
                {
                    T5ntrlMutations.push_back(locusGender.second);
                } else
                {
                    if (*it == locusGender.second)
                    {
                        T5ntrlMutations.erase(it);
                    } else
                    {
                        T5ntrlMutations.insert(it, locusGender.second);
                    }
                }
            }
        } else
        {
            TransmittedChrom.mutateT56sel_Allele(locusGender.second, Habitat);
        }
    }

    // T5ntrlMutations and selMutations are already sorted an without duplicates

    if (T5ntrlMutations.size())
        TransmittedChrom.mutateT5ntrl_Allele(T5ntrlMutations);
}*/

void LifeCycle::Mutate_T56(Haplotype& TransmittedChrom, int Habitat)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'Mutate_T56'\n";
#endif     
    int nbMuts = SSP->geneticSampler.get_T56_nbMuts();;
    //std::cout << "nbMuts = " << nbMuts << "\n";
    
    for (int i = 0 ; i < nbMuts ; i++)
    {
        auto MutPosition = SSP->geneticSampler.get_T56_mutationPosition();
        
        // Make the mutation - toggle bit
        //std::cout << "MutPosition = " << MutPosition << "\n";

        auto locusGender = SSP->Gmap.FromT56LocusToT56genderLocus(MutPosition);
        //std::cout << "locusGender.second = " << locusGender.second << "\n";
        if (locusGender.isNtrl)
        {
            TransmittedChrom.mutateT56ntrl_Allele(locusGender.locusInGender);
        } else
        {
            TransmittedChrom.mutateT56sel_Allele(locusGender.locusInGender, Habitat);
        }
    }
}


void LifeCycle::findAllParents(Pop& pop, const std::vector<int>& patchSizeNextGeneration)
{
    PD.resizeForNewGeneration(patchSizeNextGeneration);

    assert(patchSizeNextGeneration.size() == GP->PatchNumber);
    for (int patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
    {
        for (int offspring_index = 0 ; offspring_index < patchSizeNextGeneration[patch_index] ; ++offspring_index)
        {
            CoupleData CD = findCouple(pop, patch_index, offspring_index, PD, patchSizeNextGeneration[patch_index]); // Give PD for cloneInfo

            
            PD.lastOffspring[CD.mother.segregationIndex][CD.mother.patch][CD.mother.ind] = HaplotypeData(patch_index, offspring_index, 0, CD.mother.nbRecs);
            if (CD.mother.nbRecs > 0)
            {
                PD.lastOffspring[!CD.mother.segregationIndex][CD.mother.patch][CD.mother.ind] = HaplotypeData(patch_index, offspring_index, 0, CD.mother.nbRecs);
            }
            
            PD.lastOffspring[CD.father.segregationIndex][CD.father.patch][CD.father.ind] = HaplotypeData(patch_index, offspring_index, 1, CD.father.nbRecs);
            if (CD.father.nbRecs > 0)
            {
                PD.lastOffspring[!CD.father.segregationIndex][CD.father.patch][CD.father.ind] = HaplotypeData(patch_index, offspring_index, 1, CD.father.nbRecs);
            }
            
                
            PD.couples[patch_index][offspring_index] = CD;
        }
    }

    
    /*{
        assert(GP->PatchNumber == 1);

        std::vector<uint32_t> count(2*SSP->TotalpatchSize,0);
        for (uint32_t ind_index = 0 ; ind_index < patchSizeNextGeneration[0] ; ++ind_index)
        {
            auto& CD = PD.couples[0][ind_index];
            ++(count[2 * CD.mother.ind + CD.mother.segregationIndex]);
            ++(count[2 * CD.father.ind + CD.father.segregationIndex]);
        }

        std::vector<uint32_t> dist(2*SSP->TotalpatchSize,0);
        for (uint32_t i = 0 ; i < count.size(); ++i)
        {
            ++(dist[count[i]]);
        }

        std::cout << "\n";
        uint32_t nbChildren = 0;
        for (uint32_t i = 0 ; i < dist.size(); ++i)
        {
            nbChildren += i * dist[i];
            if (dist[i]) std::cout << i << ": " << dist[i] << "\n";
            //for (uint32_t j = 0 ; j < dist[i] ; ++j)
            //{
             //   std::cout << "*";
            //}
            //std::cout << std::endl;
        }
        assert(nbChildren == 2*patchSizeNextGeneration[0]);
    }*/
    
}



LifeCycle::CoupleData LifeCycle::findCouple(Pop& pop, int& patch_index, int& offspring_index, ParentsData& PD, const int& nextGenPatchSize)
{
    int mother_patch_from;
    if (SSP->isStochasticMigration)
        mother_patch_from = pop.SelectionOriginPatch(patch_index);
    else
        mother_patch_from = pop.SelectionOriginPatch(patch_index, nextGenPatchSize==1 ? 0 : (double) offspring_index / (double) (nextGenPatchSize-1) );

    #ifdef DEBUG
    if (mother_patch_from != patch_index) NBMIGRANTS++;
    #endif
    // Select the mother
    int mother_index = pop.SelectionParent(mother_patch_from, 0); // Selection on fertility (if applicable) in there

    if (SSP->cloningRate != 0.0 && (SSP->cloningRate == 1.0 || (GP->rngw.uniform_real_distribution(1.0) < SSP->cloningRate)))
    {
        // Cloning
        PD.cloneInfo[patch_index][offspring_index] = true;
        int segregationIndex = GP->rngw.get_1b();
        return CoupleData(
            HaplotypeData(mother_patch_from, mother_index, segregationIndex, 0), // mother
            HaplotypeData(mother_patch_from, mother_index, !segregationIndex, 0)  // father
        );
        
    } else
    {
        // Not cloning
        if (SSP->cloningRate != 0.0)
        {  
            PD.cloneInfo[patch_index][offspring_index] = false;
        }
        

        int father_patch_from;
        int father_index;
        if (SSP->selfingRate > 0.0 && GP->rngw.uniform_real_distribution(1.0) < SSP->selfingRate)
        {
            // selfing
            father_patch_from = mother_patch_from;
            father_index = mother_index;
        } else
        {
            // not selfing
            if (SSP->selfingRate == 0.0)
            {
                // No selfing at all -> Not even following Wright-Fisher model
                father_index = -1;
                father_patch_from = -1;
                while (father_index == mother_index && father_patch_from == mother_patch_from)
                {
                    if (SSP->gameteDispersal)
                    {
                        if (SSP->isStochasticMigration)
                            father_patch_from = pop.SelectionOriginPatch(patch_index);
                        else
                            father_patch_from = pop.SelectionOriginPatch(patch_index, nextGenPatchSize==1 ? 0 : (double) offspring_index / (double) (nextGenPatchSize-1) );
                    } else
                    {
                        father_patch_from = mother_patch_from;
                    }
                    father_index = pop.SelectionParent(father_patch_from, SSP->malesAndFemales); // selection on fertility in there
                }
                assert(father_index >= 0);
                assert(father_patch_from >= 0);
            } else
            {
                // Wright-Fisher model. slefing with prob 1/2N (assuming neutrality) -> Default
                if (SSP->gameteDispersal)
                {
                    if (SSP->isStochasticMigration)
                        father_patch_from = pop.SelectionOriginPatch(patch_index);
                    else
                        father_patch_from = pop.SelectionOriginPatch(patch_index, nextGenPatchSize==1 ? 0 : (double) offspring_index / (double) (nextGenPatchSize-1) );
                } else
                {
                    father_patch_from = mother_patch_from;
                }
                father_index = pop.SelectionParent(father_patch_from, SSP->malesAndFemales); // selection on fertility in there
            }
        }

        //std::cout << "FC: " << mother_patch_from << " " << mother_index << " | " << father_patch_from << " " << father_index << "\n";

        auto segregationIndexA = GP->rngw.get_1b();
        auto segregationIndexB = GP->rngw.get_1b();
        auto nbRecsA = SSP->geneticSampler.get_nbRecombinations();
        auto nbRecsB = SSP->geneticSampler.get_nbRecombinations();
        // I set these values before because the order of evaluation is compiler dependent otherwise which would make simulations not reproducible
        return CoupleData(
            HaplotypeData(mother_patch_from, mother_index, segregationIndexA, nbRecsA), // mother
            HaplotypeData(father_patch_from, father_index, segregationIndexB, nbRecsB)  // father
        );
    }
}
