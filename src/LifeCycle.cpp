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


Note for Remi of things to do:
    Test how much slower it is to have N=1e5, L=10 vs N=10, L=1e5 to estimate the cost of having Individuals not contiguous in memory

    When using several environment, fitnessMap should take migration rate into account. This migration rate can vary through time and therefore fitnessMap should be redefined for faster simulations
 */




void LifeCycle::BREEDING_SELECTION_DISPERSAL(Pop& Offspring_pop, Pop& Parent_pop)
{   
    /* 
    #ifdef DEBUG
    std::cout << "Enters in 'BREEDING_SELECTION_DISPERSAL'\n";
    #endif
    */



    // Calculate fitness
    SSP->simTracker.prepareT1SitesForFitness(Parent_pop);
    Parent_pop.CalculateFitnesses(); // also set Index First Male;

    // 0. ReSet Dispersal info if patch size may vary and // 2. Compute patch size for next generation
    std::vector<int> patchSizeNextGeneration;
    if (SSP->fecundityForFitnessOfOne != -1.0 || SSP->DispWeightByFitness)
    {
        patchSizeNextGeneration = SSP->dispersalData.SetBackwardMigrationAndGetNextGenerationPatchSizes(
                Parent_pop.CumSumFits, 
                false
            );
        assert(patchSizeNextGeneration.size() == GP->PatchNumber);
    } else
    {
        patchSizeNextGeneration = SSP->patchCapacity;
    }

    // 4.5 Genealogy
    if (SSP->simTracker.genealogy.isTime())
        SSP->simTracker.genealogy.startNewGeneration();

    // Loop through each patch and through each individual offspring to be created
    //#pragma omp parallel for
    for (int patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
    {
        //std::cout << " patchSizeNextGeneration["<<patch_index<<"] = " <<  patchSizeNextGeneration[patch_index] << "\n";
         // Loop over all offspring
        for (int offspring_index = 0 ; offspring_index < patchSizeNextGeneration[patch_index] ; ++offspring_index)
        {
            //std::cout << "LifeCycle line 74 offspring_index = "<<offspring_index<<"\n";
            //std::cout << " patchSizeNextGeneration["<<patch_index<<"] = " <<  patchSizeNextGeneration[patch_index] << "   GP->PatchNumber = "<<GP->PatchNumber<<"\n";
            // Get offspring haplotypes
            Haplotype& Offspring_mHaplo = Offspring_pop.getPatch(patch_index).getInd(offspring_index).getHaplo(0); // maternally inherited haplotype
            Haplotype& Offspring_fHaplo = Offspring_pop.getPatch(patch_index).getInd(offspring_index).getHaplo(1); // paternally inherited haplotype

            // 3: select where the mother came from
            int mother_patch_from = Parent_pop.SelectionOriginPatch(patch_index);
            #ifdef DEBUG
            if (mother_patch_from != patch_index) NBMIGRANTS++;
            #endif          
            // 4: Select parents.
            // Select the mother
            int mother_index = Parent_pop.SelectionParent(mother_patch_from,0);

            bool shouldIClone = false;
            if (SSP->cloningRate != 0.0)
            {
                if (GP->random_0and1(GP->mt) < SSP->cloningRate)
                {
                    shouldIClone = true;
                }
            }
            Individual& mother = Parent_pop.getPatch(mother_patch_from).getInd(mother_index);
            
            // 4: Select parents.
            // Select the father
            int father_patch_from;
            int father_index;

            if (shouldIClone)
            {
                Offspring_mHaplo = mother.getHaplo(0);
                Offspring_fHaplo = mother.getHaplo(1);
                father_index = mother_index;
                father_patch_from = mother_patch_from;
            } else
            {
                if (SSP->selfingRate > 0.0 && GP->random_0and1(GP->mt) < SSP->selfingRate)
                {
                    father_index = mother_index;
                    father_patch_from = mother_patch_from;
                } else
                {
                    if (SSP->gameteDispersal)
                    {
                        father_patch_from = Parent_pop.SelectionOriginPatch(patch_index);
                    } else
                    {
                        father_patch_from = mother_patch_from;
                    }
                        
                    father_index = Parent_pop.SelectionParent(father_patch_from, SSP->malesAndFemales);
                }
                Individual& father = Parent_pop.getPatch(father_patch_from).getInd(father_index);

                // 6. The function 'recombination' direclty write the offsprings (one chromosome per call of 'recombination' in 'Offsrping_pop'   
                recombination(
                   mother,
                   Offspring_mHaplo
                );

                recombination(
                   father,
                   Offspring_fHaplo
                );

            } // end of if ShouldIClone         
            // Mutate
            Mutate(
                Offspring_mHaplo,
                SSP->Habitats[patch_index] // needs the habitat to adjust the fitness
            );

            Mutate(
                Offspring_fHaplo,
                SSP->Habitats[patch_index] // needs the habitat to adjust the fitness
            );
            // 1. Fitness for next generation. Set fitness to zero if the habitat has changed
            if (SSP->Habitats[patch_index] != SSP->Habitats[mother_patch_from])
            {
                if (SSP->T1_isSelection)
                {
                    Offspring_mHaplo.setAllW_T1(-1.0);
                }
                if (SSP->T2_isSelection)
                {
                    Offspring_mHaplo.setAllW_T2(-1.0);
                }
            }
            if (SSP->Habitats[patch_index] != SSP->Habitats[father_patch_from])
            {
                if (SSP->T1_isSelection)
                {
                    Offspring_fHaplo.setAllW_T1(-1.0);
                }
                if (SSP->T2_isSelection)
                {
                    Offspring_fHaplo.setAllW_T2(-1.0);
                }
            }
            // 4.5 Genealogy. Will do something only if the last isTime was true
            SSP->simTracker.genealogy.addOffspringIfIsTime(patch_index, offspring_index, mother_patch_from, mother_index, father_patch_from, father_index);         
        }
    }
    // set patch size
    SSP->patchSize = patchSizeNextGeneration;
    GP->saveSSPPatchSize_toGP();

    // 4.5 Genealogy. Will do something only if the last isTime was true
    //std::cout << "line 184\n";
    SSP->simTracker.genealogy.printBufferIfIsTime();
    //std::cout << "line 186\n";

/*    if (SSP->fecundityForFitnessOfOne != -1.0)
    {
        assert(!SSP->isSpeciesExtinct);
        SSP->isSpeciesExtinct = true;
        for (int patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
        {
            SSP->patchSize[patch_index] = patchSizeNextGeneration[patch_index];
            assert(Offspring_pop.CumSumFits[patch_index].size() == patchSizeNextGeneration[patch_index]);
            if (patchSizeNextGeneration[patch_index] > 0)
            {
                SSP->isSpeciesExtinct = false;
            }
        }
        if (SSP->isSpeciesExtinct)
        {
            SSP->simTracker.gotExtinct();
        }
    }
    */

}


int LifeCycle::recombination_nbRecs()
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'recombination_nbRecs'\n";
#endif       
    int nbRecs;
    
    if (SSP->TotalRecombinationRate != 0)
    {
        nbRecs = SSP->rpois_nbRecombination(GP->mt); // For the moment does not allow variation in recombination rate.
    } else
    {
        nbRecs = 0;
    }
    return nbRecs;
}

void LifeCycle::recombination_RecPositions(int& nbRecs, std::vector<int>& breakpoints, Individual& parent)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'recombination_RecPositions'\n";
#endif  
    // Find positions (here a position is a separation between bit or between a bit and a byte or anything of interest) where recombination occurs
    if (nbRecs)
    {
        //std::cout << "nbRecs = "<< nbRecs << std::endl;
        std::vector<int> RecPositions;
        RecPositions.reserve(nbRecs);
        int RecPosition;
        double rnd;
        
        for (int i = 0 ; i < nbRecs ; i++)
        {
            while(true)
            {
                if (SSP->RecombinationRate.size() == 1) // Constant value
                {
                    RecPosition = SSP->runiform_int_ForRecPos(GP->mt); // is bounded between 0 and size - 1.
                } else // Not constant
                {
                    assert(SSP->RecombinationRate.size() == SSP->TotalNbLoci - 1);

                    rnd = SSP->runiform_double_ForRecPos(GP->mt); // is bounded between 0 and total_recombinationRate.
                    
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
                    int toT1Locus   = std::min(SSP->T1_nbBits,RecPosition + SSP->recRateOnMismatch_halfWindow);
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

                    if (probabilityToKeepRecEvent > GP->random_0and1(GP->mt))
                    {
                        RecPositions.push_back(RecPosition);
                        break;
                    } // else keep looping over the while(true) loop
                } else
                {
                    RecPositions.push_back(RecPosition);
                    break;
                }
            }
        }
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
                    breakpoints.push_back(RecPositions[i-1]);
                }
                nbRepeats = 1;
            }
        }
        
        if (nbRepeats % 2)
        {
            //std::cout << RecPositions.back() << "\n";
            assert(RecPositions.back() >= 0 && RecPositions.back() < SSP->TotalNbLoci);
            breakpoints.push_back(RecPositions.back());
        }
    }
    
    // Add segregation   
    if (SSP->ChromosomeBoundaries.size() != 0) assert(SSP->ChromosomeBoundaries.back() != INT_MAX);
    if (breakpoints.size()!=0) assert(breakpoints.back() != INT_MAX);
    for (int b : SSP->ChromosomeBoundaries)
    {
        if (GP->random_0or1(GP->mt)) breakpoints.push_back(b);
    }

    // Sort and add upper bound
    sort(breakpoints.begin(), breakpoints.end());
    breakpoints.push_back(INT_MAX);
}

void LifeCycle::recombination_copyOver(Individual& parent,Haplotype& TransmittedChrom, std::vector<int>& breakpoints, int segregationInfo)
{    
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'recombination_copyOver'\n";
#endif  
    if (breakpoints.size()==1)
    {
        int SegregationIndex = GP->random_0or1(GP->mt);
        TransmittedChrom = parent.getHaplo(SegregationIndex);
        
    } else if (breakpoints.size()>1)
    {
        int Safety_T1_Absolutefrom = INT_MAX;
        int Safety_T1_Absoluteto   = -1;
        
        int Safety_T2_Absolutefrom = INT_MAX;
        int Safety_T2_Absoluteto   = -1;
        
        int Safety_T3_Absolutefrom = INT_MAX;
        int Safety_T3_Absoluteto   = -1;
        
        
        int T1_from = 0;
        int T2_from = 0;
        int T3_from = 0;
        int T1_to = 0;
        int T2_to = 0;
        int T3_to = 0;
        int fitnessMapIndexFrom = 0;

        int haplo_index;
        if (segregationInfo == 0 || segregationInfo == 1)
        {
            haplo_index = segregationInfo;
        } else
        {
            haplo_index = GP->random_0or1(GP->mt); // Here is first chromosome segregation happening.
        }
        
        
        /*std::cout << "\t";
        for (auto& b : breakpoints)
            std::cout << b << " ";
        std::cout << std::endl;*/
        
        

        for (auto& b : breakpoints)
        {
            // Get positions for the two traits
            if (b == INT_MAX)
            {
                // Copy is over the range [from, to)
                T1_to = SSP->T1_nbBits;
                T2_to = SSP->T2_nbChars;
                T3_to = SSP->T3_nbChars;
                b = SSP->TotalNbLoci - 1;
            } else
            {
                assert(b >= 0);
                // because copy is over the range [from, to) and if recombination happens at position 4, then the 4th element belongs to the previous stretch of DNA while the 5th element belongs to the next stretch of DNA.
                T1_to = SSP->FromLocusToTXLocus[b].T1;
                T2_to = SSP->FromLocusToTXLocus[b].T2;
                T3_to = SSP->FromLocusToTXLocus[b].T3;
            }

            #ifdef DEBUG
            
            assert(T1_to <= SSP->T1_nbBits);
            assert(T2_to <= SSP->T2_nbChars);
            assert(T3_to <= SSP->T3_nbChars);
            
            #endif

            // Set Fitnesses
            if (SSP->T1_isSelection || SSP->T2_isSelection)
            {  
                //std::cout << "b = " << b << "  SSP->FromLocusToFitnessMapIndex.size() = " << SSP->FromLocusToFitnessMapIndex.size() << " INT_MAX = " << INT_MAX << "\n";

                //std::cout << "b = "<< b << "\n";
                //assert(SSP->TotalNbLoci == SSP->FromLocusToFitnessMapIndex.size());

                // There is a fitnessMapIndex juste before and just after the point where the recombination occurred. Check the fitnessMap just before and just after the interlocus where the crossOver happened

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
                int fitnessMapIndexToRecompute;
                if (fitnessMapIndexBefore == fitnessMapIndexAfter) 
                {
                    fitnessMapIndexToRecompute = fitnessMapIndexBefore; // Must force recalculation of fitness
                    assert(b < SSP->TotalNbLoci-1);
                } else
                {
                    fitnessMapIndexToRecompute = -1; // no need to recompute fitness
                }

                // Figure out end of things to copy over (fitnessMapIndexTo, excluded)
                int fitnessMapIndexTo = fitnessMapIndexAfter; // excluded


                // Figure out beginning of things to copy over (fitnessMapIndexFrom, included)
                if (fitnessMapIndexFrom > fitnessMapIndexTo)
                {
                    assert(fitnessMapIndexFrom == fitnessMapIndexTo + 1);
                    fitnessMapIndexFrom = fitnessMapIndexAfter;
                }

                /*std::cout << "\n";
                std::cout << "T1_from = " << T1_from << " SSP->FromT1LocusToLocus[T1_from] = " << SSP->FromT1LocusToLocus[T1_from] << " SSP->FromLocusToFitnessMapIndex[SSP->FromT1LocusToLocus[T1_from]] = " << SSP->FromLocusToFitnessMapIndex[SSP->FromT1LocusToLocus[T1_from]] << "\n";
                if (T1_to == SSP->TotalNbLoci)
                {
                    std::cout << "T1_to == SSP->TotalNbLoci\n";
                } else
                {
                    std::cout << "T1_to = " << T1_to << " SSP->FromT1LocusToLocus[T1_to] = " << SSP->FromT1LocusToLocus[T1_to] << " SSP->FromLocusToFitnessMapIndex[SSP->FromT1LocusToLocus[T1_to]] = " << SSP->FromLocusToFitnessMapIndex[SSP->FromT1LocusToLocus[T1_to]] << "\n";
                }
                
                std::cout << "fitnessMapIndexToRecompute = " << fitnessMapIndexToRecompute << "\n"; 
                std::cout << "fitnessMapIndexFrom = " << fitnessMapIndexFrom << "\n";
                std::cout << "fitnessMapIndexTo = " << fitnessMapIndexTo << "\n";*/

                /*assert(SSP->FromLocusToFitnessMapIndex[SSP->FromT1LocusToLocus[T1_from]] <= fitnessMapIndexFrom);
                if (T1_to != SSP->TotalNbLoci)
                {
                    assert(SSP->FromLocusToFitnessMapIndex[SSP->FromT1LocusToLocus[T1_to]] >= fitnessMapIndexTo);
                }*/


                /*if (!(fitnessMapIndexFrom <= fitnessMapIndexTo))
                {
                    std::cout << "SSP->FromLocusToFitnessMapIndex:\n";
                    for (auto& e : SSP->FromLocusToFitnessMapIndex)
                    {
                        std::cout << e << " ";
                    }
                    std::cout << "\n";
                }*/
                
                

                if (SSP->T1_isSelection && SSP->FitModel_T1_isMultiplicity && SSP->T2_isSelection)
                {
                    // Copy fitnesses that can be copied
                    for (int fitnessMapIndex = fitnessMapIndexFrom; fitnessMapIndex < fitnessMapIndexTo ; fitnessMapIndex++)
                    {
                        if (fitnessMapIndex != fitnessMapIndexToRecompute)
                        {
                            //std::cout << "Assgning FMI = " << fitnessMapIndex << "\n";
                            TransmittedChrom.setW_T1(
                                parent.getHaplo(haplo_index).getW_T1(fitnessMapIndex),
                                fitnessMapIndex
                            );

                            TransmittedChrom.setW_T2(
                                parent.getHaplo(haplo_index).getW_T2(fitnessMapIndex),
                                fitnessMapIndex
                            );
                        }
                    }

                    if (fitnessMapIndexToRecompute != -1)
                    {
                        // Set fitnesses of the affected fitnessMap part to -1.0. Will need to be recalculated
                        TransmittedChrom.setW_T1(-1.0, fitnessMapIndexToRecompute);
                        TransmittedChrom.setW_T2(-1.0, fitnessMapIndexToRecompute);
                    }
                } else
                {
                    if (SSP->T1_isSelection && SSP->FitModel_T1_isMultiplicity)
                    {
                        
                        /*
                        std::cout << "fitnessMapIndexToRecompute = " << fitnessMapIndexToRecompute << "\n";
                        std::cout << "fitnessMapIndexFrom = " << fitnessMapIndexFrom << "\n";
                        std::cout << "fitnessMapIndexTo = " << fitnessMapIndexTo << "\n";
                        */

                        // Copy fitnesses that can be copied
                        for (int fitnessMapIndex = fitnessMapIndexFrom; fitnessMapIndex < fitnessMapIndexTo ; fitnessMapIndex++)
                        {
                            if (fitnessMapIndex != fitnessMapIndexToRecompute)
                            {
                                //std::cout << "Assgning FMI = " << fitnessMapIndex << "\n";
                                TransmittedChrom.setW_T1(
                                    parent.getHaplo(haplo_index).getW_T1(fitnessMapIndex),
                                    fitnessMapIndex
                                );
                            }
                        }

                        if (fitnessMapIndexToRecompute != -1)
                        {
                            // Set fitnesses of the affected fitnessMap part to -1.0. Will need to be recalculated
                            TransmittedChrom.setW_T1(-1.0, fitnessMapIndexToRecompute);
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
                            }
                        }

                        // Set fitness of the affected fitnessMap part to -1.0. Will need to be recalculated
                        if (fitnessMapIndexToRecompute != -1)
                        {
                            // Set fitnesses of the affected fitnessMap part to -1.0. Will need to be recalculated
                            TransmittedChrom.setW_T2(-1.0, fitnessMapIndexToRecompute);
                        }
                    }
                }
                if (fitnessMapIndexToRecompute == -1)
                {
                    fitnessMapIndexFrom = fitnessMapIndexTo;
                } else
                {
                    fitnessMapIndexFrom = fitnessMapIndexTo + 1;
                }
            }

            //std::cout << "In recombination_CopyOver: T1_from = "<< T1_from <<" T1_to = "<< T1_to << std::endl;
            //std::cout << "In recombination_CopyOver: T2_from = "<< T2_from <<" T2_to = "<< T2_to << std::endl;

            // Copy for T1. from included, to excluded
            if (T1_from < T1_to)
            {
                assert(T1_to > 0 && T1_to <= SSP->T1_nbBits + 1);
                TransmittedChrom.copyIntoT1(T1_from, T1_to, parent.getHaplo(haplo_index));
                Safety_T1_Absolutefrom = std::min(Safety_T1_Absolutefrom, T1_from);
                Safety_T1_Absoluteto = std::max(Safety_T1_Absoluteto, T1_to);
            }       
            
            // Copy for T2. from included, to excluded
            if (T2_from < T2_to)
            {
                assert(T2_to > 0 && T2_to <= SSP->T2_nbChars + 1);
                TransmittedChrom.copyIntoT2(T2_from, T2_to, parent.getHaplo(haplo_index));
                Safety_T2_Absolutefrom = std::min(Safety_T2_Absolutefrom, T2_from);
                Safety_T2_Absoluteto = std::max(Safety_T2_Absoluteto, T2_to);
            }

            // Copy for T3. from included, to excluded
            if (T3_from < T3_to)
            {
                assert(T3_to > 0 && T3_to <= SSP->T3_nbChars + 1);
                TransmittedChrom.copyIntoT3(T3_from, T3_to, parent.getHaplo(haplo_index));
                Safety_T3_Absolutefrom = std::min(Safety_T3_Absolutefrom, T3_from);
                Safety_T3_Absoluteto = std::max(Safety_T3_Absoluteto, T3_to);
            }
            
            // Set new 'from'
            T1_from = T1_to;
            T2_from = T2_to;
            T3_from = T3_to;
            
            // Switch parent chromosome
            haplo_index = !haplo_index;
        }

        // Security Checks
        if (SSP->T1_nbChars)
        {
            assert(T1_to == SSP->T1_nbBits);
            assert(Safety_T1_Absoluteto == SSP->T1_nbBits);
            assert(Safety_T1_Absolutefrom == 0);
        }

        if (SSP->T2_nbChars)
        {
            assert(T2_to == SSP->T2_nbChars);
            assert(Safety_T2_Absoluteto == SSP->T2_nbChars);
            assert(Safety_T2_Absolutefrom == 0);
        }

        if (SSP->T3_nbChars)
        {
            assert(T3_to == SSP->T3_nbChars);
            assert(Safety_T3_Absoluteto == SSP->T3_nbChars);
            assert(Safety_T3_Absolutefrom == 0);
        }

        
    } else
    {
        std::cout << "Internal error. 'breakpoints' has size of zero. It should always contain at least the element INT_MAX\n";
        abort();
    }
}

void LifeCycle::recombination(Individual& parent,Haplotype& TransmittedChrom)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'recombination'\n";
#endif     
    int nbRecs = recombination_nbRecs();

    std::vector<int> breakpoints;
    recombination_RecPositions(nbRecs, breakpoints, parent); // parent is used if recRateOnMismatch_bool

    // copy from parents chromosomes to TransmittedChrom. The function also set the new values for W_T1 and W_T2
    recombination_copyOver(parent, TransmittedChrom, breakpoints, -1); // This will copy data from parents to offspring so it must always run
}



void LifeCycle::Mutate(Haplotype& TransmittedChrom, int Habitat)
{
    if (SSP->T1_Total_Mutation_rate>0) Mutate_T1(TransmittedChrom, Habitat);
    if (SSP->T2_Total_Mutation_rate>0) Mutate_T2(TransmittedChrom, Habitat);
    if (SSP->T3_Total_Mutation_rate>0) Mutate_T3(TransmittedChrom);
}


void LifeCycle::Mutate_T1(Haplotype& TransmittedChrom, int Habitat)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'Mutate_T1'\n";
#endif     
    int nbMuts = SSP->T1_rpois_nbMut(GP->mt);
    
    for (int i = 0 ; i < nbMuts ; i++)
    {
        //std::cout << "nbMuts = " << nbMuts << "\n";
        int MutPosition;
        double rnd;
        // Find Position
        if (SSP->T1_MutationRate.size() == 1)
        {
            MutPosition = SSP->T1_runiform_int_ForMutPos(GP->mt);
        } else
        {
            rnd = SSP->T1_runiform_double_ForMutPos(GP->mt);
            
            // binary search
            MutPosition = distance(SSP->T1_MutationRate.begin(),
                                   std::upper_bound(SSP->T1_MutationRate.begin(), SSP->T1_MutationRate.end(), rnd)
                                   );
        }
        
        // Make the mutation
        TransmittedChrom.toggleT1_Allele(MutPosition, Habitat);         // toggle bit
        // TransmittedChrom.setT1_AlleleToOne(byte_index,bit_index,);   // set bit to one
    }
}


void LifeCycle::Mutate_T2(Haplotype& TransmittedChrom, int Habitat)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'Mutate_T2'\n";
#endif       
    int nbMuts = SSP->T2_rpois_nbMut(GP->mt);
    
    for (int i = 0 ; i < nbMuts ; i++)
    {
        double rnd;
        int MutPosition;

        // Find Position
        if (SSP->T2_MutationRate.size() == 1)
        {
            MutPosition = SSP->T2_runiform_int_ForMutPos(GP->mt);
        } else
        {
            rnd = SSP->T2_runiform_double_ForMutPos(GP->mt);
            
            // binary search
            MutPosition = distance(SSP->T2_MutationRate.begin(),
                                   std::upper_bound(SSP->T2_MutationRate.begin(), SSP->T2_MutationRate.end(), rnd)
                                   );
        }
        
        // Make the mutation
        //std::cout << "MutPosition = " <<  MutPosition << "    SSP->T2_nbChars = " << SSP->T2_nbChars << "\n";
        assert(MutPosition < SSP->T2_nbChars);
        TransmittedChrom.AddMutT2_Allele(MutPosition, Habitat);  // add 1 and changes T2_W
    }
}

void LifeCycle::Mutate_T3(Haplotype& TransmittedChrom)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'Mutate_T3'\n";
#endif       
    int nbMuts = SSP->T3_rpois_nbMut(GP->mt);  
    for (int i = 0 ; i < nbMuts ; i++)
    {
        int MutPosition;
        double rnd; 

        // Find Position
        if (SSP->T3_MutationRate.size() == 1)
        {
            MutPosition = SSP->T3_runiform_int_ForMutPos(GP->mt);
        } else
        {
            rnd = SSP->T3_runiform_double_ForMutPos(GP->mt);
            
            // binary search
            MutPosition = distance(SSP->T3_MutationRate.begin(),
                   std::upper_bound(SSP->T3_MutationRate.begin(), SSP->T3_MutationRate.end(), rnd)
                   );
        }   
        // Make the mutation
        //std::cout << "MutPosition = " <<  MutPosition << "    SSP->T3_nbChars = " << SSP->T3_nbChars << "\n";
        assert(MutPosition < SSP->T3_nbChars);   
        TransmittedChrom.AddMutT3_Allele(MutPosition);  // add 1 and changes T2_W   
    }  
}

