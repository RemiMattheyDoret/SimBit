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
    /* 
    #ifdef DEBUG
    std::cout << "Enters in 'BREEDING_SELECTION_DISPERSAL'\n";
    #endif
    */   
    //size_t nbSwaps = 0;

    if (SSP->T4_nbBits > 0)
    {
        SSP->T4Tree.newGeneration();
    }


    // Calculate fitness
    //SSP->simTracker.prepareT1SitesForFitness(Parent_pop);
    if (SSP->isAnySelection && SSP->selectionOn != 1) Parent_pop.CalculateFitnesses(); // also set Index First Male; Will compute fitness only if it is needed

    // 0. ReSet Dispersal info if patch size may vary and // 2. Compute patch size for next generation
    const std::vector<int>& patchSizeNextGeneration = SSP->dispersalData.setBackwardMigrationIfNeededAndGetNextGenerationPatchSizes(Parent_pop.CumSumFits);
    assert(patchSizeNextGeneration.size() == GP->PatchNumber);


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
            //std::cout << "LifeCycle line 74 offspring_index = "<<offspring_index<<"\n";
            //std::cout << " patchSizeNextGeneration["<<patch_index<<"] = " <<  patchSizeNextGeneration[patch_index] << "   GP->PatchNumber = "<<GP->PatchNumber<<"\n";

            // Get offspring haplotypes
            Individual& Offspring = Offspring_pop.getPatch(patch_index).getInd(offspring_index);
            Haplotype& Offspring_mHaplo = Offspring.getHaplo(0); // maternally inherited haplotype
            Haplotype& Offspring_fHaplo = Offspring.getHaplo(1); // paternally inherited haplotype

            //std::cout << "\n\n\n\n\n\nBegin Offspring_mHaplo.nbT56muts() = "<< Offspring_mHaplo.nbT56muts() << "\n";
            //std::cout << "Begin Offspring_fHaplo.nbT56muts() = "<< Offspring_fHaplo.nbT56muts() << "\n";


            redoOffspring:

            // Get parents
            auto&& CD = SSP->SwapInLifeCycle ? PD.couples[patch_index][offspring_index] : findCouple(Parent_pop,patch_index, offspring_index, PD);
            
            auto& mother = Parent_pop.getPatch(CD.mother.patch).getInd(CD.mother.ind);
            auto& father = Parent_pop.getPatch(CD.father.patch).getInd(CD.father.ind);

            

            //std::cout << "Begin mother.getHaplo(CD.mother.segregationIndex).nbT56muts() = "<< mother.getHaplo(CD.mother.segregationIndex).nbT56muts() << "\n";
            //std::cout << "Begin father.getHaplo(CD.father.segregationIndex).nbT56muts() = "<< father.getHaplo(CD.father.segregationIndex).nbT56muts() << "\n";

            //std::cout << "M: " << CD.mother.patch << " " << CD.mother.ind  << " " << CD.mother.segregationIndex << " | " << CD.father.patch << " " << CD.father.ind << " " << CD.father.segregationIndex << "\n";


            // Clear T56 vectors
            if (SSP->T56sel_nbLoci || SSP->T56ntrl_nbLoci)
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

                std::vector<int> breakpoints = {INT_MAX};
                
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

            //std::cout << "just after mutation Offspring_mHaplo.nbT56muts() = "<< Offspring_mHaplo.nbT56muts() << "\n";
            //std::cout << "just after mutation Offspring_fHaplo.nbT56muts() = "<< Offspring_fHaplo.nbT56muts() << "\n";
            

            // Selection on viability
            if (SSP->selectionOn != 0)
            {
                double fitness = Offspring_pop.getPatch(patch_index).getInd(offspring_index).CalculateFitness(patch_index);
                if (GP->rngw.uniform_real_distribution(1.0) > fitness)
                {
                    CD = findCouple(Parent_pop, patch_index, offspring_index, PD); // Will modify cloneInfo too. 
                    if (SSP->SwapInLifeCycle)
                    {
                        PD.couples[patch_index][offspring_index] = CD;
                    } 
                    
                    goto redoOffspring;
                }
            }


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
    

    if (SSP->T4_nbBits > 0)
    {
        SSP->T4Tree.pruneDeadLineages();
    }

    if (SSP->T56_nbLoci)
    {
        Offspring_pop.toggleT56MutationsIfNeeded();
    }

    //std::cout << "nbSwaps = " << nbSwaps << "\n";
}


void LifeCycle::reproduceThroughSwap(Individual& parent, Haplotype& offspringHaplotype, HaplotypeData& parentHaploData, int& patch_index)
{
    
    parent.getHaplo(parentHaploData.segregationIndex).swap(offspringHaplotype);

    if (SSP->T4_nbBits > 0)
    {
        //std::cout <<"\tLIFECYCLE - CLONING!!! b = NA: SSP->T4_nbBits = "<<SSP->T4_nbBits<<"\n";
        SSP->T4Tree.addChildHaplotype_setParentInfo(parentHaploData.patch, parentHaploData.ind);
        SSP->T4Tree.addChildHaplotype_addNode(parentHaploData.segregationIndex,0,SSP->T4_nbBits);
        SSP->T4Tree.addChildHaplotype_finished(patch_index);
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
        recombination_breakpoints[0] = INT_MAX;
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
                    assert(SSP->RecombinationRate.size() == SSP->TotalNbLoci - 1);

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
                assert(RecPosition >= 0 && RecPosition < SSP->TotalNbLoci - 1);

                
                if (SSP->recRateOnMismatch_bool)
                {
                    int fromT1Locus = std::max(0,RecPosition - SSP->recRateOnMismatch_halfWindow);
                    int toT1Locus   = std::min(SSP->T1_nbLoci,RecPosition + SSP->recRateOnMismatch_halfWindow);
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
            assert(RecPositions.back() >= 0 && RecPositions.back() < SSP->TotalNbLoci);
            recombination_breakpoints.push_back(RecPositions.back());
        }
        
        // Add segregation   
        //if (SSP->ChromosomeBoundaries.size() != 0) assert(SSP->ChromosomeBoundaries.back() != INT_MAX);
        //if (breakpoints.size()!=0) assert(breakpoints.back() != INT_MAX);
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

        // add upper bound
        recombination_breakpoints.push_back(INT_MAX);
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
        TransmittedChrom = parent.getHaplo(segregationIndex);
        //if (GP->CurrentGeneration >= 50) std::cout << "TransmittedChrom.getW_T1(0) = " << TransmittedChrom.getW_T1(0) << "\n";
        if (SSP->T4_nbBits > 0)
        {
            //std::cout <<"\tLIFECYCLE!! b = NA: SSP->T4_nbBits = "<<SSP->T4_nbBits<<"\n";
            SSP->T4Tree.addChildHaplotype_addNode(segregationIndex, 0, SSP->T4_nbBits);
        }
        
    } else if (recombination_breakpoints.size()>1)
    {
        //std::cout << "breakpoints.size() " << breakpoints.size() << "\n";
        #ifdef DEBUG
        int Safety_T1_Absolutefrom = INT_MAX;
        int Safety_T1_Absoluteto   = -1;
        
        int Safety_T2_Absolutefrom = INT_MAX;
        int Safety_T2_Absoluteto   = -1;
        
        int Safety_T3_Absolutefrom = INT_MAX;
        int Safety_T3_Absoluteto   = -1;

        int Safety_T4_Absolutefrom = INT_MAX;
        int Safety_T4_Absoluteto   = -1;
        
        int Safety_T56ntrl_Absolutefrom = INT_MAX;
        int Safety_T56ntrl_Absoluteto   = -1;

        int Safety_T56sel_Absolutefrom = INT_MAX;
        int Safety_T56sel_Absoluteto   = -1;
        #endif
        
        int T1_from = 0;
        int T2_from = 0;
        int T3_from = 0;
        int T4_from = 0;
        int T56ntrl_from = 0;
        int T56sel_from = 0;
        int T1_to = 0;
        int T2_to = 0;
        int T3_to = 0;
        int T4_to = 0;
        int T56ntrl_to = 0;
        int T56sel_to = 0;
        int fitnessMapIndexFrom = 0;

        int haplo_index = segregationIndex; // There is really no reason for copying segregationIndex. I could just use segregationIndex all the way through...but he.... whatever!

        
        /*std::cout << "\t";
        for (auto& b : breakpoints)
            std::cout << b << " ";
        std::cout << std::endl;*/
        
        

        for (auto& b : recombination_breakpoints)
        {
            // Get positions for the two traits
            if (b == INT_MAX)
            {
                // Copy is over the range [from, to)
                T1_to = SSP->T1_nbLoci;
                T2_to = SSP->T2_nbLoci;
                T3_to = SSP->T3_nbLoci;
                T4_to = SSP->T4_nbBits;
                T56ntrl_to = SSP->T56ntrl_nbLoci;
                T56sel_to = SSP->T56sel_nbLoci;
                b = SSP->TotalNbLoci - 1;
            } else
            {
                assert(b >= 0);
                // because copy is over the range [from, to) and if recombination happens at position 4, then the 4th element belongs to the previous stretch of DNA while the 5th element belongs to the next stretch of DNA.
                T1_to = SSP->FromLocusToTXLocus[b].T1;
                T2_to = SSP->FromLocusToTXLocus[b].T2;
                T3_to = SSP->FromLocusToTXLocus[b].T3;
                T4_to = SSP->FromLocusToTXLocus[b].T4;
                T56ntrl_to = SSP->FromLocusToTXLocus[b].T56ntrl;
                T56sel_to = SSP->FromLocusToTXLocus[b].T56sel;
            }

            #ifdef DEBUG
            
            assert(T1_to <= SSP->T1_nbLoci);
            assert(T2_to <= SSP->T2_nbLoci);
            assert(T3_to <= SSP->T3_nbLoci);
            assert(T4_to <= SSP->T4_nbBits);
            assert(T56ntrl_to <= SSP->T56ntrl_nbLoci);
            assert(T56sel_to <= SSP->T56sel_nbLoci);
            
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
                if (b == SSP->TotalNbLoci-1)
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
                    assert(b < SSP->TotalNbLoci-1);
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
                assert(T1_to > 0 && T1_to <= SSP->T1_nbLoci + 1);
                TransmittedChrom.copyIntoT1(T1_from, T1_to, parent.getHaplo(haplo_index));
                #ifdef DEBUG
                Safety_T1_Absolutefrom = std::min(Safety_T1_Absolutefrom, T1_from);
                Safety_T1_Absoluteto = std::max(Safety_T1_Absoluteto, T1_to);
                #endif
            }       
            
            // Copy for T2. from included, to excluded
            if (T2_from < T2_to)
            {
                assert(T2_to > 0 && T2_to <= SSP->T2_nbLoci + 1);
                TransmittedChrom.copyIntoT2(T2_from, T2_to, parent.getHaplo(haplo_index));
                #ifdef DEBUG
                Safety_T2_Absolutefrom = std::min(Safety_T2_Absolutefrom, T2_from);
                Safety_T2_Absoluteto = std::max(Safety_T2_Absoluteto, T2_to);
                #endif
            }

            // Copy for T3. from included, to excluded
            if (T3_from < T3_to)
            {
                assert(T3_to > 0 && T3_to <= SSP->T3_nbLoci + 1);
                TransmittedChrom.copyIntoT3(T3_from, T3_to, parent.getHaplo(haplo_index));
                #ifdef DEBUG
                Safety_T3_Absolutefrom = std::min(Safety_T3_Absolutefrom, T3_from);
                Safety_T3_Absoluteto = std::max(Safety_T3_Absoluteto, T3_to);
                #endif
            }

            // Copy for T4. from included, to excluded
            if (T4_from < T4_to)
            {
                assert(T4_to > 0 && T4_to <= SSP->T4_nbBits + 1);
                //std::cout <<"\tLIFECYCLE!! b = "<<b<<" T4_from = "<<T4_from<<" T4_to = "<<T4_to<<"\n";
                SSP->T4Tree.addChildHaplotype_addNode(haplo_index, T4_from, T4_to);
                #ifdef DEBUG
                Safety_T4_Absolutefrom = std::min(Safety_T4_Absolutefrom, T4_from);
                Safety_T4_Absoluteto = std::max(Safety_T4_Absoluteto, T4_to);
                #endif
            }

            // Copy for T56ntrl. from included, to excluded
            if (T56ntrl_from < T56ntrl_to)
            {
                assert(T56ntrl_to > 0 && T56ntrl_to <= SSP->T56ntrl_nbLoci + 1);

                
                TransmittedChrom.copyIntoT56ntrl(T56ntrl_from, T56ntrl_to, parent.getHaplo(haplo_index));
                #ifdef DEBUG
                Safety_T56ntrl_Absolutefrom = std::min(Safety_T56ntrl_Absolutefrom, T56ntrl_from);
                Safety_T56ntrl_Absoluteto = std::max(Safety_T56ntrl_Absoluteto, T56ntrl_to);
                #endif
            }

            // Copy for T56sel. from included, to excluded
            if (T56sel_from < T56sel_to)
            {
                assert(T56sel_to > 0 && T56sel_to <= SSP->T56sel_nbLoci + 1);
                
                TransmittedChrom.copyIntoT56sel(T56sel_from, T56sel_to, parent.getHaplo(haplo_index));
                #ifdef DEBUG
                Safety_T56sel_Absolutefrom = std::min(Safety_T56sel_Absolutefrom, T56sel_from);
                Safety_T56sel_Absoluteto = std::max(Safety_T56sel_Absoluteto, T56sel_to);
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
        if (SSP->T1_nbChars)
        {
            assert(T1_to == SSP->T1_nbLoci);
            assert(Safety_T1_Absoluteto == SSP->T1_nbLoci);
            assert(Safety_T1_Absolutefrom == 0);
        }

        if (SSP->T2_nbLoci)
        {
            assert(T2_to == SSP->T2_nbLoci);
            assert(Safety_T2_Absoluteto == SSP->T2_nbLoci);
            assert(Safety_T2_Absolutefrom == 0);
        }

        if (SSP->T3_nbLoci)
        {
            assert(T3_to == SSP->T3_nbLoci);
            assert(Safety_T3_Absoluteto == SSP->T3_nbLoci);
            assert(Safety_T3_Absolutefrom == 0);
        }

        if (SSP->T4_nbBits)
        {
            assert(T4_to == SSP->T4_nbBits);
            assert(Safety_T4_Absoluteto == SSP->T4_nbBits);
            assert(Safety_T4_Absolutefrom == 0);
        }

        if (SSP->T56ntrl_nbLoci)
        {
            assert(T56ntrl_to == SSP->T56ntrl_nbLoci);
            assert(Safety_T56ntrl_Absoluteto == SSP->T56ntrl_nbLoci);
            assert(Safety_T56ntrl_Absolutefrom == 0);
        }

        if (SSP->T56sel_nbLoci)
        {
            assert(T56sel_to == SSP->T56sel_nbLoci);
            assert(Safety_T56sel_Absoluteto == SSP->T56sel_nbLoci);
            assert(Safety_T56sel_Absolutefrom == 0);
        }
#endif

        
    } else
    {
        std::cout << "Internal error. 'breakpoints' has size of zero. It should always contain at least the element INT_MAX\n";
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

    if (SSP->T4_nbBits > 0)
    {
        SSP->T4Tree.addChildHaplotype_setParentInfo(parentData.patch, parentData.ind);
    }

    recombination_RecPositions(parentData.nbRecs, parent); // parent is used if recRateOnMismatch_bool

    // copy from parents chromosomes to TransmittedChrom. The function also set the new values for W_T1 and W_T2
    copyOver(parent, TransmittedChrom, parentData.segregationIndex); // This will copy data from parents to offspring so it must always run

     if (SSP->T4_nbBits > 0)
    {
        SSP->T4Tree.addChildHaplotype_finished(patch_index);
    }
}



void LifeCycle::Mutate(Haplotype& TransmittedChrom, int Habitat)
{
    if (SSP->T1_Total_Mutation_rate>0) Mutate_T1(TransmittedChrom, Habitat);
    if (SSP->T2_Total_Mutation_rate>0) Mutate_T2(TransmittedChrom, Habitat);
    if (SSP->T3_Total_Mutation_rate>0) Mutate_T3(TransmittedChrom);
    if (SSP->T56_Total_Mutation_rate>0) Mutate_T56(TransmittedChrom, Habitat);
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
        //std::cout << "MutPosition = " <<  MutPosition << "    SSP->T2_nbLoci = " << SSP->T2_nbLoci << "\n";
        assert(MutPosition < SSP->T2_nbLoci);
        TransmittedChrom.AddMutT2_Allele(MutPosition, Habitat);  // add 1 and changes T2_W
    }
}

void LifeCycle::Mutate_T3(Haplotype& TransmittedChrom)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'Mutate_T3'\n";
#endif       
    int nbMuts = SSP->geneticSampler.get_T3_nbMuts();
    for (int i = 0 ; i < nbMuts ; i++)
    {
        auto MutPosition = SSP->geneticSampler.get_T3_mutationPosition();

        // Make the mutation
        //std::cout << "MutPosition = " <<  MutPosition << "    SSP->T3_nbLoci = " << SSP->T3_nbLoci << "\n";
        assert(MutPosition < SSP->T3_nbLoci);   
        TransmittedChrom.AddMutT3_Allele(MutPosition);  // add 1 and changes T2_W   
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
    if (SSP->T5ntrl_nbLoci)
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
            if (SSP->T56ntrl_compress)
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

        auto& locusGender = SSP->FromT56LocusToT56genderLocus[MutPosition];
        //std::cout << "locusGender.second = " << locusGender.second << "\n";
        if (locusGender.first)
        {
            TransmittedChrom.mutateT56ntrl_Allele(locusGender.second);
        } else
        {
            TransmittedChrom.mutateT56sel_Allele(locusGender.second, Habitat);
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
            CoupleData CD = findCouple(pop, patch_index, offspring_index, PD); // Give PD for cloneInfo

            
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

        std::vector<size_t> count(2*SSP->TotalpatchSize,0);
        for (size_t ind_index = 0 ; ind_index < patchSizeNextGeneration[0] ; ++ind_index)
        {
            auto& CD = PD.couples[0][ind_index];
            ++(count[2 * CD.mother.ind + CD.mother.segregationIndex]);
            ++(count[2 * CD.father.ind + CD.father.segregationIndex]);
        }

        std::vector<size_t> dist(2*SSP->TotalpatchSize,0);
        for (size_t i = 0 ; i < count.size(); ++i)
        {
            ++(dist[count[i]]);
        }

        std::cout << "\n";
        size_t nbChildren = 0;
        for (size_t i = 0 ; i < dist.size(); ++i)
        {
            nbChildren += i * dist[i];
            if (dist[i]) std::cout << i << ": " << dist[i] << "\n";
            //for (size_t j = 0 ; j < dist[i] ; ++j)
            //{
             //   std::cout << "*";
            //}
            //std::cout << std::endl;
        }
        assert(nbChildren == 2*patchSizeNextGeneration[0]);
    }*/
    
}



LifeCycle::CoupleData LifeCycle::findCouple(Pop& pop, int& patch_index, int& offspring_index, ParentsData& PD)
{
    int mother_patch_from = pop.SelectionOriginPatch(patch_index);
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
                        father_patch_from = pop.SelectionOriginPatch(patch_index);
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
                    father_patch_from = pop.SelectionOriginPatch(patch_index);
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
