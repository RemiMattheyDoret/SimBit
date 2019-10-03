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

void SimulationTracker::gotExtinct()
{
    assert(this->whenDidExtinctionOccur == -1);
    this->whenDidExtinctionOccur = GP->CurrentGeneration;
}    

void SimulationTracker::addMutation(int byte_index, int bit_index, int MutPosition)
{
    //std::cout << "mut at " << MutPosition << "\n";
    if (SSP->recomputeLociOverWhichFitnessMustBeComputedEveryHowManyGenerations != -1 && SSP->T1_isSelection)
    {
        T1_locusDescription TrackedMutation(byte_index, bit_index, MutPosition);
        if (SSP->T1_isMultiplicitySelection)
        {
            //std::cout << MutPosition << " is in fitness block " << SSP->FromLocusToFitnessMapIndex[MutPosition] << "\n";
            T1SitesForFitnessMultiplicityCalculationMutationsForNextGeneration[ SSP->FromLocusToFitnessMapIndex[SSP->FromT1LocusToLocus[MutPosition]]].push_back(TrackedMutation);
        } else
        {
            T1SitesForFitnessNoMultiplicityCalculationMutationsForNextGeneration.push_back(TrackedMutation);
        }
    }
}


void SimulationTracker::Initialization(Pop& pop)
{
    #ifdef CALLENTRANCEFUNCTIONS
        std::cout << "Enters in 'SimulationTracker::Initialization(Pop pop)'" << std::endl;
    #endif
    // Just an assertion
    int x = 0;
    assert(SSP->patchCapacity.size() == GP->PatchNumber);
    for (auto& OnepatchCapacity : SSP->patchCapacity)
    {
        x += OnepatchCapacity;
    }
    assert(x == SSP->TotalpatchCapacity);

    //this->MinimalSizeForReducingPolymorphicT1Sites = std::min(SSP->T1_nbBits / 500, SSP->TotalpatchCapacity * (SSP->T1_nbChars + SSP->T2_nbChars)); // 'this->MinimalSizeForReducingPolymorphicT1Sites' is initiated so that 'PolymorphicT1Sites' should not take more than a twelfth of the RAM taken by all genetic data.
    
    if (SSP->T1_isSelection)
    {
        if (SSP->T1_isMultiplicitySelection)
        {
            prepareFitnessMapIndexForT1SitesForFitnessMultiplicity();
            if (!SSP->T1_Initial_AlleleFreqs_AllZeros || SSP->recomputeLociOverWhichFitnessMustBeComputedEveryHowManyGenerations == -1)
            {
                recomputeT1SitesForFitnessMultiplicityCalculation(pop);
            }
        } else // if there is possibly any genetic variance (which include the case where data is read form a binary file)
        {
            this->recomputeT1SitesForFitnessNoMultiplicityCalculation(pop);

            if (SSP->recomputeLociOverWhichFitnessMustBeComputedEveryHowManyGenerations == -1)
                assert(T1SitesForFitnessNoMultiplicityCalculation.size() == SSP->T1_nbBits);
        }
        assert(!(
            (T1SitesForFitnessNoMultiplicityCalculation.size() != 0)
            &&
            (T1SitesForFitnessMultiplicityCalculation.size() != 0)
        ));
        //std::cout << "In initialization: T1SitesForFitnessNoMultiplicityCalculation.size() = " << T1SitesForFitnessNoMultiplicityCalculation.size() << "\n";
    }
    forceToRecomputeNextTime = false;
}

void SimulationTracker::prepareT1SitesForFitness(Pop& pop)
{
    
    /*std::cout << "\nT1SitesForFitnessNoMultiplicityCalculationMutationsForNextGeneration.size() = " << T1SitesForFitnessNoMultiplicityCalculationMutationsForNextGeneration.size() << "\n";
    std::cout << "T1SitesForFitnessNoMultiplicityCalculation.size() = " << T1SitesForFitnessNoMultiplicityCalculation.size() << "\n";
    std::cout << "\n\nforceToRecomputeNextTime = " << forceToRecomputeNextTime << "\n\n\n";*/
    
    if (forceToRecomputeNextTime)
    {
        //std::cout << "being forced to recompute\n";
        if (SSP->T1_isMultiplicitySelection)
        {
            recomputeT1SitesForFitnessMultiplicityCalculation(pop);
        } else
        {
            recomputeT1SitesForFitnessNoMultiplicityCalculation(pop);
        }

    } else if (SSP->recomputeLociOverWhichFitnessMustBeComputedEveryHowManyGenerations != -1 && SSP->T1_isSelection)
    {
        if (SSP->T1_isMultiplicitySelection)
        {
            // add the new mutations
            addMutsSortAndRemoveDuplicateOfT1SitesForFitnessMultiplicity();

            if ( GP->CurrentGeneration!= 0 && GP->CurrentGeneration % SSP->recomputeLociOverWhichFitnessMustBeComputedEveryHowManyGenerations == 0 )
            {
                RemoveNoFitnessEffectT1SitesForFitnessMultiplicity(pop);
            }
         
        } else
        {
            addMutsSortAndRemoveDuplicateOfT1SitesForFitnessNoMultiplicity();

            if ( GP->CurrentGeneration!= 0 && GP->CurrentGeneration % SSP->recomputeLociOverWhichFitnessMustBeComputedEveryHowManyGenerations == 0 )
            {
                RemoveNoFitnessEffectT1SitesForFitnessNoMultiplicity(pop);
            }

            
            //std::cout << "T1SitesForFitnessNoMultiplicityCalculation.size() = " << T1SitesForFitnessNoMultiplicityCalculation.size() << "\n";
        }
    }
    forceToRecomputeNextTime = false;


    /*
    std::cout << "T1SitesForFitnessMultiplicityCalculation.size() = " << T1SitesForFitnessMultiplicityCalculation.size() << "\n";
    for (auto& elem : T1SitesForFitnessMultiplicityCalculation)
            std::cout << elem.size() << " ";
        std::cout << "\n";
    */
        
        

     /*int SizeBefore = T1SitesForFitnessNoMultiplicityCalculation.size();
    if (SizeBefore >= MinimalSizeForReducingPolymorphicT1Sites) 
    {
        RemoveNoFitnessEffectT1SitesForFitnessNoMultiplicity();

        int SizeAfter = T1SitesForFitnessNoMultiplicityCalculation.size();
        double relativeDiffInSize = (double)(SizeBefore - SizeAfter) / (double) SizeBefore;
        if (relativeDiffInSize < 0.1) // If the size reduction was not considerable enough
        {
            MinimalSizeForReducingPolymorphicT1Sites = SizeAfter * 1.1; // then increase the minimal threshold for doing such size reduction
        } 
    } else if (SizeBefore * 2 < MinimalSizeForReducingPolymorphicT1Sites)
    {
        MinimalSizeForReducingPolymorphicT1Sites /= 1.5;
    }*/
    
}


void SimulationTracker::printT1SitesForFitnessNoMultiplicity()
{
    for (auto& elem : T1SitesForFitnessNoMultiplicityCalculation)
    {
        std::cout << elem.locus << " ";
    }
    std::cout << std::endl;
}

void SimulationTracker::prepareFitnessMapIndexForT1SitesForFitnessMultiplicity()
{
    T1SitesForFitnessMultiplicityCalculation.resize(SSP->NbElementsInFitnessMap);
    T1SitesForFitnessMultiplicityCalculationMutationsForNextGeneration.resize(SSP->NbElementsInFitnessMap);
}


void SimulationTracker::mergeSortRemoveDuplicates(std::vector<T1_locusDescription>& v1, std::vector<T1_locusDescription>& v2)
{
  std::sort(v2.begin(), v2.end());
  v2.erase(
        unique(
            v2.begin(),
            v2.end()
        ),
        v2.end()
    );
  
  v1.insert(
    v1.end(),
    make_move_iterator(v2.begin()),
    make_move_iterator(v2.end())
  );
  v2.clear();

  std::sort(v1.begin(), v1.end());
  v1.erase(
        unique(
            v1.begin(),
            v1.end()
        ),
        v1.end()
    );
}

void SimulationTracker::addMutsSortAndRemoveDuplicateOfT1SitesForFitnessNoMultiplicity()
{
    // sort and remove duplicates of the new mutations
    mergeSortRemoveDuplicates(
        T1SitesForFitnessNoMultiplicityCalculation,
        T1SitesForFitnessNoMultiplicityCalculationMutationsForNextGeneration
    );

}


void SimulationTracker::addMutsSortAndRemoveDuplicateOfT1SitesForFitnessMultiplicity()
{
    assert(T1SitesForFitnessMultiplicityCalculation.size() == SSP->NbElementsInFitnessMap);
    assert(T1SitesForFitnessMultiplicityCalculationMutationsForNextGeneration.size() == SSP->NbElementsInFitnessMap);
    for (int fitnessMapIndex = 0 ; fitnessMapIndex < SSP->NbElementsInFitnessMap; fitnessMapIndex++)
    {
        // add mutations of last generation
        mergeSortRemoveDuplicates(
            T1SitesForFitnessMultiplicityCalculation[fitnessMapIndex],
            T1SitesForFitnessMultiplicityCalculationMutationsForNextGeneration[fitnessMapIndex]
        );
    }
}

void SimulationTracker::setToMinus1LociThatShouldNotBeTracked(Pop& pop, std::vector<T1_locusDescription>& trackedLoci )
{
    int nbElementsFiguredOut = 0;

    // 1) pick one allele for each locus tracked, 2) check if there is any selection and 3) note those who can't be remove anyway because the one allele picked already causes a fitness loss
    std::vector<bool> oneAlleles(trackedLoci.size());
    std::vector<bool> mustBeKeptForSure(trackedLoci.size());
    for (int trackedLocusIndex = 0 ; trackedLocusIndex < trackedLoci.size() ; trackedLocusIndex++)
    {
        auto& T1Site = trackedLoci[trackedLocusIndex];

        // initialize mustBeKeptForSure
        mustBeKeptForSure[trackedLocusIndex] = false;

        // pick a population from which to pick an allele
        int patchForReferenceAllele = 0;
        if (SSP->patchSize[patchForReferenceAllele] == 0)
        {
            for (;patchForReferenceAllele < GP->PatchNumber ; patchForReferenceAllele++)
            {
                if (SSP->patchSize[patchForReferenceAllele] > 0)
                {
                    break;
                }
            }
        }

        bool isUnderSelectionAnywhere = false;
        if (patchForReferenceAllele < GP->PatchNumber)
        {
            // pick an allele
            oneAlleles[trackedLocusIndex] = pop.getPatch(patchForReferenceAllele).getInd(0).getHaplo(0).getT1_Allele(T1Site.byte_index, T1Site.bit_index);

            // check if it is under selection and if the allele found is itself causing fitness loss
            for (int Habitat = 0 ; Habitat <= SSP->MaxHabitat ; Habitat++)
            {
                if (SSP->T1_isMultiplicitySelection)
                {
                    if (SSP->T1_FitnessEffects[Habitat][T1Site.locus] != 1.0)
                    {
                        isUnderSelectionAnywhere = true;
                        
                        if (oneAlleles[trackedLocusIndex] != 0)
                        {
                            mustBeKeptForSure[trackedLocusIndex] = true;
                            //std::cout << "1. trackedLocusIndex = " << trackedLocusIndex << "\n";
                            nbElementsFiguredOut++;
                        }
                    }
                } else
                {
                    if (
                        SSP->T1_FitnessEffects[Habitat][THREE * T1Site.locus + 2] != 1.0 ||
                        SSP->T1_FitnessEffects[Habitat][THREE * T1Site.locus + 1] != 1.0 ||
                        SSP->T1_FitnessEffects[Habitat][THREE * T1Site.locus + 0] != 1.0
                    )
                    {
                        isUnderSelectionAnywhere = true;
                    }
                }
            }        
        }

    
        if (!isUnderSelectionAnywhere)
        {
            //std::cout << "removed " << T1Site.locus << "\n";
            T1Site.locus = -1;
            //std::cout << "2. trackedLocusIndex = " << trackedLocusIndex << "\n";
            nbElementsFiguredOut++;
        }
    }

    // Is it polymorphic
    std::vector<bool> isPolymorphic(trackedLoci.size(), false);
    for (int patch_index = 0 ; patch_index < GP->PatchNumber ; patch_index++)
    {
        Patch& patch = pop.getPatch(patch_index);
        for (int ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ind_index++)
        {
            Individual& Ind = patch.getInd(ind_index);
            for (int haplo_index = 0 ; haplo_index < SSP->ploidy ; haplo_index++)
            {
                if (nbElementsFiguredOut == trackedLoci.size())
                {
                    goto endOfLoops;
                }
                for (int trackedLocusIndex = 0 ; trackedLocusIndex < trackedLoci.size() ; trackedLocusIndex++)
                {
                    auto& T1Site = trackedLoci[trackedLocusIndex];
                     // if the fate has not been figured out yet
                    if (T1Site.locus != -1 && !mustBeKeptForSure[trackedLocusIndex])
                    {
                        if (oneAlleles[trackedLocusIndex] != Ind.getHaplo(haplo_index).getT1_Allele(T1Site.byte_index, T1Site.bit_index))
                        {
                            mustBeKeptForSure[trackedLocusIndex] = true;
                            isPolymorphic[trackedLocusIndex] = true;
                            //std::cout << "3. trackedLocusIndex = " << trackedLocusIndex << "\n";
                            nbElementsFiguredOut++;
                        }
                    }
                }
            }
        }
    }

    endOfLoops:

    // set to -1 those that were not polymorphic
    for (int trackedLocusIndex = 0 ; trackedLocusIndex < trackedLoci.size() ; trackedLocusIndex++)
    {
        if (!mustBeKeptForSure[trackedLocusIndex] && trackedLoci[trackedLocusIndex].locus != -1)
        {
            //std::cout << "removed " << trackedLoci[trackedLocusIndex].locus << "\n";
            trackedLoci[trackedLocusIndex].locus = -1;
            //std::cout << "4. trackedLocusIndex = " << trackedLocusIndex << "\n";
            nbElementsFiguredOut++;
        }
    }

    /*
    std::cout << "nbElementsFiguredOut = " << nbElementsFiguredOut << "\n";
    std::cout << "trackedLoci.size() = " << trackedLoci.size() << "\n";
    std::cout << "loci tracked: ";
    for (int trackedLocusIndex = 0 ; trackedLocusIndex < trackedLoci.size() ; trackedLocusIndex++)
    {
        assert(trackedLoci[trackedLocusIndex].locus == -1 || trackedLoci[trackedLocusIndex].locus == trackedLocusIndex);
        if (trackedLoci[trackedLocusIndex].locus != -1)
        {
            std::cout << trackedLoci[trackedLocusIndex].locus << " ";
        }
    }
    std::cout << "\n";
    */
    assert(nbElementsFiguredOut == trackedLoci.size());

}

void SimulationTracker::setToMinus1LociThatShouldNotBeTracked(Pop& pop, std::vector<T1_locusDescription*>& trackedLoci )
{
    int nbElementsFiguredOut = 0;

    // 1) pick one allele for each locus tracked, 2) check if there is any selection and 3) note those who can't be remove anyway because the one allele picked already causes a fitness loss
    std::vector<bool> oneAlleles(trackedLoci.size());
    std::vector<bool> mustBeKeptForSure(trackedLoci.size());
    for (int trackedLocusIndex = 0 ; trackedLocusIndex < trackedLoci.size() ; trackedLocusIndex++)
    {
        auto& T1Site = trackedLoci[trackedLocusIndex];

        // initialize mustBeKeptForSure
        mustBeKeptForSure[trackedLocusIndex] = false;

        // pick a population from which to pick an allele
        int patchForReferenceAllele = 0;
        if (SSP->patchSize[patchForReferenceAllele] == 0)
        {
            for (;patchForReferenceAllele < GP->PatchNumber ; patchForReferenceAllele++)
            {
                if (SSP->patchSize[patchForReferenceAllele] > 0)
                {
                    break;
                }
            }
        }

        bool isUnderSelectionAnywhere = false;
        if (patchForReferenceAllele < GP->PatchNumber)
        {
            // pick an allele
            oneAlleles[trackedLocusIndex] = pop.getPatch(patchForReferenceAllele).getInd(0).getHaplo(0).getT1_Allele(T1Site->byte_index, T1Site->bit_index);

            // check if it is under selection and if the allele found is itself causing fitness loss
            for (int Habitat = 0 ; Habitat <= SSP->MaxHabitat ; Habitat++)
            {
                if (SSP->T1_isMultiplicitySelection)
                {
                    if (SSP->T1_FitnessEffects[Habitat][T1Site->locus] != 1.0)
                    {
                        isUnderSelectionAnywhere = true;
                        
                        if (oneAlleles[trackedLocusIndex] != 0)
                        {
                            mustBeKeptForSure[trackedLocusIndex] = true;
                            //std::cout << "1. trackedLocusIndex = " << trackedLocusIndex << "\n";
                            nbElementsFiguredOut++;
                        }
                    }
                } else
                {
                    if (
                        SSP->T1_FitnessEffects[Habitat][THREE * T1Site->locus + 2] != 1.0 ||
                        SSP->T1_FitnessEffects[Habitat][THREE * T1Site->locus + 1] != 1.0 ||
                        SSP->T1_FitnessEffects[Habitat][THREE * T1Site->locus + 0] != 1.0
                    )
                    {
                        isUnderSelectionAnywhere = true;
                    }
                }
            }        
        }

    
        if (!isUnderSelectionAnywhere)
        {
            //std::cout << "removed " << T1Site.locus << "\n";
            T1Site->locus = -1;
            //std::cout << "2. trackedLocusIndex = " << trackedLocusIndex << "\n";
            nbElementsFiguredOut++;
        }
    }

    // Is it polymorphic
    std::vector<bool> isPolymorphic(trackedLoci.size(), false);
    for (int patch_index = 0 ; patch_index < GP->PatchNumber ; patch_index++)
    {
        Patch& patch = pop.getPatch(patch_index);
        for (int ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ind_index++)
        {
            Individual& Ind = patch.getInd(ind_index);
            for (int haplo_index = 0 ; haplo_index < SSP->ploidy ; haplo_index++)
            {
                if (nbElementsFiguredOut == trackedLoci.size())
                {
                    goto endOfLoops;
                }
                for (int trackedLocusIndex = 0 ; trackedLocusIndex < trackedLoci.size() ; trackedLocusIndex++)
                {
                    auto& T1Site = trackedLoci[trackedLocusIndex];
                     // if the fate has not been figured out yet
                    if (T1Site->locus != -1 && !mustBeKeptForSure[trackedLocusIndex])
                    {
                        if (oneAlleles[trackedLocusIndex] != Ind.getHaplo(haplo_index).getT1_Allele(T1Site->byte_index, T1Site->bit_index))
                        {
                            mustBeKeptForSure[trackedLocusIndex] = true;
                            isPolymorphic[trackedLocusIndex] = true;
                            //std::cout << "3. trackedLocusIndex = " << trackedLocusIndex << "\n";
                            nbElementsFiguredOut++;
                        }
                    }
                }
            }
        }
    }

    endOfLoops:

    // set to -1 those that were not polymorphic
    for (int trackedLocusIndex = 0 ; trackedLocusIndex < trackedLoci.size() ; trackedLocusIndex++)
    {
        if (!mustBeKeptForSure[trackedLocusIndex] && trackedLoci[trackedLocusIndex]->locus != -1)
        {
            //std::cout << "removed " << trackedLoci[trackedLocusIndex].locus << "\n";
            trackedLoci[trackedLocusIndex]->locus = -1;
            //std::cout << "4. trackedLocusIndex = " << trackedLocusIndex << "\n";
            nbElementsFiguredOut++;
        }
    }

    /*
    std::cout << "nbElementsFiguredOut = " << nbElementsFiguredOut << "\n";
    std::cout << "trackedLoci.size() = " << trackedLoci.size() << "\n";
    std::cout << "loci tracked: ";
    for (int trackedLocusIndex = 0 ; trackedLocusIndex < trackedLoci.size() ; trackedLocusIndex++)
    {
        assert(trackedLoci[trackedLocusIndex].locus == -1 || trackedLoci[trackedLocusIndex].locus == trackedLocusIndex);
        if (trackedLoci[trackedLocusIndex].locus != -1)
        {
            std::cout << trackedLoci[trackedLocusIndex].locus << " ";
        }
    }
    std::cout << "\n";
    */
    assert(nbElementsFiguredOut == trackedLoci.size());
}

void SimulationTracker::RemoveNoFitnessEffectT1SitesForFitnessNoMultiplicity(Pop& pop)
{
    setToMinus1LociThatShouldNotBeTracked(pop, T1SitesForFitnessNoMultiplicityCalculation);

    // remove those set to -1
    this->T1SitesForFitnessNoMultiplicityCalculation.erase(
        std::remove_if(
            this->T1SitesForFitnessNoMultiplicityCalculation.begin(),
            this->T1SitesForFitnessNoMultiplicityCalculation.end(),
            [](T1_locusDescription& x)
            {
                return x.locus == -1;
            }
        ),
        this->T1SitesForFitnessNoMultiplicityCalculation.end()
    );
}

void SimulationTracker::RemoveNoFitnessEffectT1SitesForFitnessMultiplicity(Pop& pop)
{
    assert(T1SitesForFitnessMultiplicityCalculation.size() == SSP->NbElementsInFitnessMap);

    std::vector<T1_locusDescription*> trackedLoci;
    for (int fitnessMapIndex = 0 ; fitnessMapIndex < SSP->NbElementsInFitnessMap; fitnessMapIndex++)
    {
        for (auto& T1Site : T1SitesForFitnessMultiplicityCalculation[fitnessMapIndex])
        {
            trackedLoci.push_back(&T1Site);
        }
    }

    setToMinus1LociThatShouldNotBeTracked(pop, trackedLoci);



    // remove those set to -1
    for (int fitnessMapIndex = 0 ; fitnessMapIndex < SSP->NbElementsInFitnessMap; fitnessMapIndex++)
    {
        this->T1SitesForFitnessMultiplicityCalculation[fitnessMapIndex].erase(
            std::remove_if(
                this->T1SitesForFitnessMultiplicityCalculation[fitnessMapIndex].begin(),
                this->T1SitesForFitnessMultiplicityCalculation[fitnessMapIndex].end(),
                [](T1_locusDescription& x)
                {
                    return x.locus == -1;
                }
            ),
            this->T1SitesForFitnessMultiplicityCalculation[fitnessMapIndex].end()
        );
    }
}

void SimulationTracker::recomputeT1SitesForFitnessNoMultiplicityCalculation(Pop& pop)
{
    #ifdef CALLENTRANCEFUNCTIONS
        std::cout << "Enters in 'SimulationTracker::recomputeT1SitesForFitnessNoMultiplicityCalculation(Pop& pop)'" << std::endl;
    #endif

    //std::cout << "In 'SimulationTracker::recomputeT1SitesForFitnessNoMultiplicityCalculation' T1SitesForFitnessNoMultiplicityCalculation.size() = " << T1SitesForFitnessNoMultiplicityCalculation.size() << "\n";


    // List all sites into T1SitesForFitnessNoMultiplicityCalculation
    T1SitesForFitnessNoMultiplicityCalculation.resize(SSP->T1_nbBits);
    T1SitesForFitnessNoMultiplicityCalculation.shrink_to_fit();
    for (int byte_index = 0 ; byte_index < SSP->T1_nbChars ; byte_index++)
    {
        int bit_index_to;
        if (byte_index == SSP->T1_nbChars - 1)
        {
            bit_index_to = SSP->T1_nbBitsLastByte;
        } else
        {
            bit_index_to = EIGHT;
        }
        for (int bit_index = 0 ; bit_index < bit_index_to ; bit_index++)
        {
            int locus = byte_index * 8 + bit_index;
            T1_locusDescription T1Site(byte_index, bit_index, locus);
            T1SitesForFitnessNoMultiplicityCalculation[locus] = T1Site;
        }
    }

    //std::cout << "In 'SimulationTracker::recomputeT1SitesForFitnessNoMultiplicityCalculation' T1SitesForFitnessNoMultiplicityCalculation.size() = " << T1SitesForFitnessNoMultiplicityCalculation.size() << "\n";


    if (SSP->recomputeLociOverWhichFitnessMustBeComputedEveryHowManyGenerations != -1)
        RemoveNoFitnessEffectT1SitesForFitnessNoMultiplicity(pop);

    //std::cout << "In 'SimulationTracker::recomputeT1SitesForFitnessNoMultiplicityCalculation' T1SitesForFitnessNoMultiplicityCalculation.size() = " << T1SitesForFitnessNoMultiplicityCalculation.size() << "\n";

}



void SimulationTracker::recomputeT1SitesForFitnessMultiplicityCalculation(Pop& pop)
{
    #ifdef CALLENTRANCEFUNCTIONS
        std::cout << "Enters in 'SimulationTracker::recomputeT1SitesForFitnessNoMultiplicityCalculation(Pop& pop)'" << std::endl;
    #endif

    assert(T1SitesForFitnessMultiplicityCalculation.size() == SSP->NbElementsInFitnessMap);

    int T1_locusFrom = 0; // from included
    for (int fitnessMapIndex = 0 ; fitnessMapIndex < SSP->NbElementsInFitnessMap; fitnessMapIndex++)
    {
        T1SitesForFitnessMultiplicityCalculation[fitnessMapIndex].resize(0);
        int T1_locusTo = SSP->FromFitnessMapIndexToTXLocus[fitnessMapIndex].T1; // to excluded
        assert(T1_locusTo >= T1_locusFrom);

        int byte_index = T1_locusFrom / 8;
        int bit_index = T1_locusFrom % 8;
        for (int locus = T1_locusFrom ; locus < T1_locusTo ; locus++)
        {
            T1_locusDescription T1Site(byte_index, bit_index, locus);
            T1SitesForFitnessMultiplicityCalculation[fitnessMapIndex].push_back(T1Site);
            if (bit_index == 7)
            {
                bit_index = 0;
                byte_index++;
            } else
            {
               bit_index++; 
            }
        }

        // set from Locus
        T1_locusFrom = T1_locusTo;
    }

    if (SSP->recomputeLociOverWhichFitnessMustBeComputedEveryHowManyGenerations != -1)
        RemoveNoFitnessEffectT1SitesForFitnessMultiplicity(pop);

}



std::vector<T1_locusDescription> SimulationTracker::listAllPolymorphicT1Sites(Pop& pop)
{
    // Pick one allele for each locus
    std::vector<bool> oneAlleles(SSP->T1_nbBits);
    std::vector<bool> alreadyListed(SSP->T1_nbBits);
    {
        int locus = 0;
        for (int byte_index = 0 ; byte_index < SSP->T1_nbChars ; byte_index++)
        {
            int bit_index_to;
            if (byte_index == SSP->T1_nbChars - 1)
            {
                bit_index_to = SSP->T1_nbBitsLastByte;
            } else
            {
                bit_index_to = EIGHT;
            }
            for (int bit_index = 0 ; bit_index < bit_index_to ; bit_index++)
            {
                int patch_index = SSP->selectNonEmptyPatch(-1, SSP->patchSize, false);
                oneAlleles[locus] = pop.getPatch(patch_index).getInd(0).getHaplo(0).getT1_Allele(byte_index, bit_index);
                alreadyListed[locus] = false;
                locus++;
            }
        }
    }
        
    
    // Find polymorphic loci
    std::vector<T1_locusDescription> polymorphicLoci;
    polymorphicLoci.reserve(SSP->T1_MutationRate.back() * 5 * SSP->TotalpatchCapacity);
    for (int patch_index = 0 ; patch_index < GP->PatchNumber ; patch_index++)
    {
        Patch& patch = pop.getPatch(patch_index);
        for (int ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ind_index++)
        {
            Individual& Ind = patch.getInd(ind_index);
            for (int haplo_index = 0 ; haplo_index < SSP->ploidy ; haplo_index++)
            {
                int locus = 0;
                for (int byte_index = 0 ; byte_index < SSP->T1_nbChars ; byte_index++)
                {
                    int bit_index_to;
                    if (byte_index == SSP->T1_nbChars - 1)
                    {
                        bit_index_to = SSP->T1_nbBitsLastByte;
                    } else
                    {
                        bit_index_to = EIGHT;
                    }
                    for (int bit_index = 0 ; bit_index < bit_index_to ; bit_index++)
                    {                       
                        if (!alreadyListed[locus])
                        {
                            if (oneAlleles[locus] != Ind.getHaplo(haplo_index).getT1_Allele(byte_index, bit_index))
                            {
                                T1_locusDescription T1Site(byte_index, bit_index, locus);
                                //std::cout << "locus " << locus << " Ind " << ind_index << " Haplo " << haplo_index << "\n";
                                polymorphicLoci.push_back(T1Site);
                                alreadyListed[locus] = true;
                            }
                        }
                        locus++;
                    }
                }
            }
        }
    }

    std::sort(polymorphicLoci.begin(),polymorphicLoci.end());
    polymorphicLoci.shrink_to_fit();

    return polymorphicLoci;
}



/*void SimulationTracker::CheckIfT1SitesAreStillPolymorphic(Pop& pop)
{
    this->SortAndRemoveDuplicateOfPolymorphicT1Sites();

    // Then find if the sites that are left are still polymorphic
    // It will set the attribute 'locus' to '-1' if found to not be polymorphic
    for (auto& TrackedMutation : this->PolymorphicT1Sites )
    {
        int byte_index = TrackedMutation.byte_index;
        int bit_index  = TrackedMutation.bit_index;
        bool FoundToBeStillPolymorphic = false;
        
        bool OneAllele = pop.getPatch(...).getInd(0).getHaplo(0).getT1_Allele(byte_index,bit_index);
        for (int patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
        {
            Patch& Patch = pop.getPatch(patch_index);
            for (int ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ++ind_index)
            {
                Individual& Ind = Patch.getInd(ind_index);
                for (int haplo_index = 0 ; haplo_index < SSP->ploidy ; ++haplo_index)
                {
                    if (OneAllele != Ind.getHaplo(haplo_index).getT1_Allele(byte_index,bit_index))
                    {
                        FoundToBeStillPolymorphic = true;
                        break;
                    } 
                }
                if (FoundToBeStillPolymorphic) break;
            }
            if (FoundToBeStillPolymorphic) break;
        }
        if (!FoundToBeStillPolymorphic) // If did not find polymorphism
        {
            TrackedMutation.locus = -1; // In 'RemovePolymorphicT1SitesWhichLocusIsSetToNegative', all sites which locus is set to -1 will be removed
        }
    }
    this->RemovePolymorphicT1SitesWhichLocusIsSetToNegative();
}*/
