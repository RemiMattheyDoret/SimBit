
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


int ResetGeneticsEvent::generalAssertEvent()
{
    assert(SSP != nullptr);
    assert(GP != nullptr);
    assert(GP->__GenerationChange.size() >= 0);
    assert(GP->__PatchNumber.size() == GP->__GenerationChange.size());

    if (generation < 0 || generation > GP->nbGenerations)
    {
        std::cout << "For option --resetGenetics (error found in the constructor of ResetGeneticsEvent) the generation received is either inferior to zero or greater to the number of generations set at option --nbGenerations (" << GP->nbGenerations << "). Generation received is " << generation << ".\n";
        abort();
    }

    int generation_index = std::upper_bound(GP->__GenerationChange.begin(), GP->__GenerationChange.end(), generation) - GP->__GenerationChange.begin() - 1;
    assert(generation_index >= 0 && generation_index < GP->__GenerationChange.size());
    return generation_index;
}

ResetGeneticsEvent_B::ResetGeneticsEvent_B( 
    int g,
    std::vector<std::string> indTypeNames,
    std::vector<unsigned> h,
    int pa
)
:individualTypeNames(indTypeNames), howManyIndividualOfEachType(h), patch_index(pa)
{
    eventType = 'B';
    generation = g;
    (void) generalAssertEvent();
    assert(individualTypeNames.size() == howManyIndividualOfEachType.size());
}


ResetGeneticsEvent_A::ResetGeneticsEvent_A( 
    int g,
    char MT,
    std::vector<int> T1l,
    std::vector<int> T2l,
    std::vector<int> T3l,
    std::vector<int> T5l,
    std::vector<int> pa,
    std::vector<int> hap,
    std::vector<std::vector<int>> inds
    )
    :mutationType(MT),
    T1loci(T1l),
    T2loci(T2l),
    T3loci(T3l),
    T5loci(T5l),
    patches(pa),
    haplotypes(hap),
    individuals(inds)
{
    eventType = 'A';
    generation = g;

    int generation_index = generalAssertEvent();



    if (T1loci.size() + T2loci.size() + T3loci.size() + T5loci.size() == 0)
    {
        std::cout << "For option --resetGenetics (error found in the constructor of ResetGeneticsEvent) no loci (nor T1, T2, T3 or T5) have been indicated for an event happening at generation "<< this->generation<<".\n";
        abort();
    }

    

    for (auto& T3locus : T3loci)
    {
        if (T3locus < 0)
        {
            std::cout << "For option --resetGenetics (error found in the constructor of ResetGeneticsEvent) received T3locus index of "<<T3locus<<" for an event happening at generation "<<generation<<".\n";
            abort();
        }
        if (T3locus >= SSP->T3_nbLoci)
        {
            std::cout << "For option --resetGenetics (error found in the constructor of ResetGeneticsEvent) for an event happening at generation "<<generation<<", received T3 locus index of "<<T3locus<<" while there are only "<<SSP->T3_nbLoci<<" T3 loci in total. As a reminder, the first locus has index 0.\n";
            abort();
        }
    }

    if (haplotypes.size() == 0 || haplotypes.size() > 2)
    {
        std::cout << "For option --resetGenetics (error found in the constructor of ResetGeneticsEvent) "<<haplotypes.size()<<" haplotypes have been indicated for an event happening at generation "<<generation<<". That might be an internal error but you might want to check your input anyway.\n";
        abort();
    }

    int sumOfHaplotypes = std::accumulate(haplotypes.begin(), haplotypes.end(), 0);
    if (
        !(
            ((sumOfHaplotypes == 0 || sumOfHaplotypes == 1) && haplotypes.size() == 1) || 
            (sumOfHaplotypes == 1 && haplotypes.size() == 2)
        )
    )
    {
        if (SSP->ploidy != 2)
        {
            std::cout << "Internal Error. Ooops the option --resetGenetics is not seem to robust to change of ploidy () it assumes ploidy of 2. Please report that to Remi.\n";
            abort();
        }
        std::cout << "For option --resetGenetics (error found in the constructor of ResetGeneticsEvent) for an event happening at generation "<<generation<<", "<<haplotypes.size()<<" haplotypes have been indicated but their sum (sum of indices) is "<<sumOfHaplotypes<<". This indicates a problem!.\n";
        abort();
    }

    if (patches.size() == 0)
    {
        std::cout << "For option --resetGenetics (error found in the constructor of ResetGeneticsEvent) no patch have been indicated for an event happening at generation "<<generation<<".\n";
        abort();
    }

    int minPatchIndex = INT_MAX;
    int maxPatchIndex = -1;

    for (int& patch_index : patches)
    {
        minPatchIndex = std::min(minPatchIndex, patch_index);
        maxPatchIndex = std::max(maxPatchIndex, patch_index);
    }

    if (minPatchIndex < 0)
    {
        std::cout << "For option --resetGenetics (error found in the constructor of ResetGeneticsEvent) the patch index "<<minPatchIndex<<" has been received. The lowest patch index must be positive.\n";
        abort();   
    }
    if (maxPatchIndex >= GP->__PatchNumber[generation_index])
    {
        std::cout << "For option --resetGenetics (error found in the constructor of ResetGeneticsEvent), the patch index "<<maxPatchIndex<<" has been received for an event at generation "<<generation<<". At generation "<<generation<<", there are only "<<GP->__PatchNumber[generation_index]<<" patches. Note that the first patch is designated with 0 and the last patch at this generation with index "<<GP->__PatchNumber[generation_index] - 1<<"\n";
        abort();
    }

    if (individuals.size() != patches.size())
    {
        std::cout << "For option --resetGenetics (error found in the constructor of ResetGeneticsEvent), for an event at generation "<<generation<<". Received " << patches.size() << " patch indices but  description for individuals in "<<individuals.size()<<" patches. That might be an internal error but you might want to check your input parameters.\n";
        abort();
    }

    int totalNbModifiedInds = 0;
    for (int patch_index_index = 0 ; patch_index_index < patches.size() ; patch_index_index++)
    {
        totalNbModifiedInds += individuals[patch_index_index].size();
    }

    if (totalNbModifiedInds == 0)
    {
        std::cout << "For option --resetGenetics (error found in the constructor of ResetGeneticsEvent), for an event at generation "<<generation<<". The total number of individuals affected by the events is 0! That might be an internal error but you might want to check your input parameters.\n";
        abort();
    }

    if (mutationType >= 3 && T1loci.size() > 0)
    {
        std::cout << "Internal Error! For option --resetGenetics (error found in the constructor of ResetGeneticsEvent), for an event at generation "<<generation<<". The mutationType received is stored in SimBit as "<<mutationType<<" which has no meaning. Note that if you inputted a value that cannot be stored in a 8 bits object (outside the range -127 to 127), it is possible that your input has been recasted to "<<mutationType<<" although it is not the actual value you wrote! The value you wrote was different from 0 (meaning set allele to 0), 1 (meaning set allele to 1) or 2 (meaning toggle allele) anyway!\n";
        abort();
    }
}



/*
######################################################################
######################################################################
######################################################################
*/

void ResetGenetics::addEvent(ResetGeneticsEvent_A event)
{
    eventsA.push_back(event);

    std::stable_sort(
        eventsA.begin(), 
        eventsA.end(),
        [](const ResetGeneticsEvent_A& left, const ResetGeneticsEvent_A& right)
        {
          return left.generation < right.generation;
        }
    );
}

void ResetGenetics::addEvent(ResetGeneticsEvent_B event)
{
    eventsB.push_back(event);
    
    std::stable_sort(
        eventsB.begin(), 
        eventsB.end(),
        [](const ResetGeneticsEvent_B& left, const ResetGeneticsEvent_B& right)
        {
          return left.generation < right.generation;
        }
    );   
}

bool ResetGenetics::isGeneration()
{
    return isEventAGeneration() || isEventBGeneration();
}

bool ResetGenetics::isEventAGeneration()
{
    if (eventsA.size() > 0)
    {
        return eventsA[0].generation == GP->CurrentGeneration;
    }
    return false;
}


bool ResetGenetics::isEventBGeneration()
{
    if (eventsB.size() > 0)
    {
        return eventsB[0].generation == GP->CurrentGeneration;
    }
    return false;
}


void ResetGenetics::makeEventAHappen(ResetGeneticsEvent_A& event, Pop& pop)
{
    for (int& patch_index : event.patches)
    {
        assert(patch_index < GP->PatchNumber);
        Patch& patch = pop.getPatch(patch_index);
        assert(event.individuals.size() > patch_index);
        assert(SSP->Habitats.size() > patch_index);
        for (int& ind_index : event.individuals[patch_index])
        {
            assert(ind_index < SSP->patchCapacity[patch_index]);
            if (ind_index < SSP->patchSize[patch_index])
            {
                for (int& haplo_index : event.haplotypes)
                {
                    Haplotype& haplo = patch.getInd(ind_index).getHaplo(haplo_index);

                    // T1
                    for (int& T1locus : event.T1loci)
                    {
                        if (event.mutationType == 2)
                        {
                            haplo.mutateT1_Allele(T1locus, SSP->Habitats[patch_index]);
                        } else if (event.mutationType == 0)
                        {
                            if (haplo.getT1_Allele(T1locus))
                            {
                                haplo.mutateT1_Allele(T1locus, SSP->Habitats[patch_index]);
                            }
                        } else if (event.mutationType == 1)
                        {
                            if (!haplo.getT1_Allele(T1locus))
                            {
                                haplo.mutateT1_Allele(T1locus, SSP->Habitats[patch_index]);
                            }
                        }

                        /*if (SSP->T1mutsDirectional)
                        {
                            haplo.mutateT1_Allele(T1locus, SSP->Habitats[patch_index]);
                        } else
                        {
                            if (event.mutationType == 2)
                            {
                                haplo.mutateT1_Allele(T1locus, SSP->Habitats[patch_index]);
                            } else if (event.mutationType == 0)
                            {
                                if (haplo.getT1_Allele(T1locus))
                                {
                                    haplo.mutateT1_Allele(T1locus, SSP->Habitats[patch_index]);
                                }
                            } else if (event.mutationType == 1)
                            {
                                if (!haplo.getT1_Allele(T1locus))
                                {
                                    haplo.mutateT1_Allele(T1locus, SSP->Habitats[patch_index]);
                                }
                            }
                        }*/
                    }

                    // T2
                    for (int& T2locus : event.T2loci)
                    {
                        haplo.setT2_Allele(T2locus, 0);
                    }

                    // T3
                    for (int& T3locus : event.T3loci)
                    {
                        haplo.setT3_Allele(T3locus, 0);
                    }

                    // T5
                    for (int& T5locus : event.T5loci)
                    {
                        auto& T5genderlocus = SSP->FromT56LocusToT56genderLocus[T5locus];
                        int T5locusInGender = (int) T5genderlocus.second;

                        // Security
                        if (T5genderlocus.first)
                        {
                            assert(T5locusInGender < SSP->T5ntrl_nbLoci);
                        } else
                        {
                            assert(T5locusInGender < SSP->T5sel_nbLoci);
                        }


                        if (SSP->T56_isMultiplicitySelection)
                        {
                            haplo.setW_T56(-1.0, SSP->FromLocusToFitnessMapIndex[SSP->FromT56selLocusToLocus[T5locus]]);
                        }


                        if (event.mutationType == 2)
                        {
                            if (T5genderlocus.first)
                            {
                                haplo.toggleT5ntrl_Allele(T5locusInGender);
                            } else
                            {
                                haplo.toggleT5sel_Allele(T5locusInGender);
                            }
                                
                        } else if (event.mutationType == 0)
                        {
                            if (T5genderlocus.first)
                            {
                                // ntrl
                                haplo.setT5ntrl_AlleleToZero(T5locusInGender); // setT5ntrl_AlleleToZero should be able to deal with flip meaning system thingy
                            } else
                            {   
                                // sel
                                haplo.setT5sel_AlleleToZero(T5locusInGender);
                                if (SSP->T56_isMultiplicitySelection)
                                {
                                    haplo.setW_T56(-1.0, SSP->FromLocusToFitnessMapIndex[SSP->FromT56selLocusToLocus[T5locus]]);
                                }
                            }
                        } else if (event.mutationType == 1)
                        {
                            if (T5genderlocus.first)
                            {
                                // ntrl
                                haplo.setT5ntrl_AlleleToOne(T5locusInGender); // setT5ntrl_AlleleToOne should be able to deal with flip meaning system thingy
                            } else
                            {   
                                // sel
                                haplo.setT5sel_AlleleToOne(T5locusInGender);
                                if (SSP->T56_isMultiplicitySelection)
                                {
                                    haplo.setW_T56(-1.0, SSP->FromLocusToFitnessMapIndex[SSP->FromT56selLocusToLocus[T5locus]]);
                                }
                            }
                        }
                        

                        /*if (SSP->T5mutsDirectional) // Directional
                        {
                            if (T5genderlocus.first)
                            {
                                haplo.mutateT5ntrl_Allele(T5locusInGender);
                            } else
                            {
                                haplo.mutateT5sel_Allele(T5locusInGender, SSP->Habitats[patch_index]);
                            }
                                
                        } else // not directional
                        {
                            if (SSP->T56_isMultiplicitySelection)
                            {
                                haplo.setW_T56(-1.0, SSP->FromLocusToFitnessMapIndex[SSP->FromT56selLocusToLocus[T5locus]]);
                            }


                            if (event.mutationType == 2)
                            {
                                if (T5genderlocus.first)
                                {
                                    haplo.toggleT5ntrl_Allele(T5locusInGender);
                                } else
                                {
                                    haplo.toggleT5sel_Allele(T5locusInGender);
                                }
                                    
                            } else if (event.mutationType == 0)
                            {
                                if (T5genderlocus.first)
                                {
                                    // ntrl
                                    haplo.setT5ntrl_AlleleToZero(T5locusInGender); // setT5ntrl_AlleleToZero should be able to deal with flip meaning system thingy
                                } else
                                {   
                                    // sel
                                    haplo.setT5sel_AlleleToZero(T5locusInGender);
                                    if (SSP->T56_isMultiplicitySelection)
                                    {
                                        haplo.setW_T56(-1.0, SSP->FromLocusToFitnessMapIndex[SSP->FromT56selLocusToLocus[T5locus]]);
                                    }
                                }
                            } else if (event.mutationType == 1)
                            {
                                if (T5genderlocus.first)
                                {
                                    // ntrl
                                    haplo.setT5ntrl_AlleleToOne(T5locusInGender); // setT5ntrl_AlleleToOne should be able to deal with flip meaning system thingy
                                } else
                                {   
                                    // sel
                                    haplo.setT5sel_AlleleToOne(T5locusInGender);
                                    if (SSP->T56_isMultiplicitySelection)
                                    {
                                        haplo.setW_T56(-1.0, SSP->FromLocusToFitnessMapIndex[SSP->FromT56selLocusToLocus[T5locus]]);
                                    }
                                }
                            }
                        }*/
                    }
                }
            }
        }
    }
}

void ResetGenetics::makeEventBHappen(ResetGeneticsEvent_B& event, Pop& pop)
{   

    // Assertions
    assert(event.individualTypeNames.size() > 0);
    assert(event.individualTypeNames.size() == event.howManyIndividualOfEachType.size());


    // Get the patch
    auto patch_index = event.patch_index;
    Patch& ThePatch = pop.getPatch(patch_index);
    assert(ThePatch.getpatchCapacity() == SSP->patchCapacity[event.patch_index]);

    // Where should I start add individuals from?
    unsigned global_ind_index = 0;
    if (SSP->fecundityForFitnessOfOne != -1)
    {
        if (SSP->patchSize[patch_index] == SSP->patchCapacity[patch_index])
        {
            global_ind_index = 0;
        } else
        {
            global_ind_index = SSP->patchSize[patch_index];
        }

    } else
    {
        assert(SSP->patchCapacity[patch_index] == SSP->patchSize[patch_index]);
        global_ind_index = 0;
    }
    unsigned startAddIndividualsAt = global_ind_index; // Just for security check



    // Loop through each indType to add
    assert(event.individualTypeNames.size() == event.howManyIndividualOfEachType.size());
    for(unsigned individualType_index = 0 ; individualType_index < event.individualTypeNames.size() ; ++individualType_index)
    {
        // Gather info
        std::string indTypeName = event.individualTypeNames[individualType_index];
        auto nbIndsToAdd = event.howManyIndividualOfEachType[individualType_index];

        // Modify patch size info
        if (SSP->fecundityForFitnessOfOne != -1)
        {
            SSP->patchSize[patch_index] += nbIndsToAdd;
            if (SSP->patchSize[patch_index] > SSP->patchCapacity[patch_index])
            {
                SSP->patchSize[patch_index] = SSP->patchCapacity[patch_index];   
            }
        }
            
        
        // Get individual Type
        assert(SSP->individualTypes.find(indTypeName) != SSP->individualTypes.end());
        Individual& TheIndType = SSP->individualTypes[indTypeName];
        (void) TheIndType.CalculateFitness(patch_index); // So that it only computes once.
        assert(nbIndsToAdd <= SSP->patchCapacity[patch_index]);

        

        // Prepare Setting individuals
        auto from = global_ind_index;
        auto to = global_ind_index + nbIndsToAdd;
        
        if (to >= SSP->patchCapacity[patch_index]) // 'to' should never be equal to SSP->patchCapacity[patch_index]
        {
            to -= SSP->patchCapacity[patch_index];
            assert(to < SSP->patchCapacity[patch_index]);
            assert(to <= from);
        }
        if (to == from) {assert(nbIndsToAdd == SSP->patchCapacity[patch_index]);}



        // Set individuals
        unsigned ind_index = from;
        do
        {
            // Do not use push_back as inds is actually always intiialized with size at carrying capacity (which is probably stupid)
            assert(ind_index < SSP->patchCapacity[patch_index]);
            ThePatch.setInd(TheIndType, ind_index);

            // incrementation
            ++ind_index;
            if (ind_index == SSP->patchCapacity[patch_index]) // 'to' should never be equal to SSP->patchCapacity[patch_index]
            {
                ind_index = 0;
            }
        } while (ind_index != to);
        global_ind_index = to;


        // security
        if (to > from) // Before went back to zero
        {
            assert(global_ind_index > startAddIndividualsAt); // Make sure it did less or exaclty one full circle
        } else if (to <= from) // After it went back to zero
        {
            assert(global_ind_index <= startAddIndividualsAt);
        }
    }
}

void ResetGenetics::resetPopIfNeeded(Pop& pop)
{
    bool isA = isEventAGeneration();
    bool isB = isEventBGeneration();
    while (isA || isB)
    {
        // Make event
        // always the first even because we destroy them as we read them.

        while (isA)
        {
            makeEventAHappen(eventsA[0], pop);
            eventsA.erase(eventsA.begin());
            isA = isEventAGeneration();
        }

        if (isB)
        {
            makeEventBHappen(eventsB[0], pop);
            eventsB.erase(eventsB.begin());
            isB = isEventBGeneration();
        }
        
        //pop.hasCrazyResettingHappened = true;
    }
}


































/*
void ProgrammedT1Mutation::MakeProgrammedT1Mutation(Pop& pop)
{
    if (SSP->ProgrammedT1Mutations.size() != 0 && SSP->ProgrammedT1MutationsIndexToDo < SSP->ProgrammedT1Mutations.size())
    {
        while (SSP->ProgrammedT1MutationsIndexToDo < SSP->ProgrammedT1Mutations.size())
        {
            if (SSP->ProgrammedT1Mutations[SSP->ProgrammedT1MutationsIndexToDo].generation == GP->CurrentGeneration)
            {
                ProgrammedT1Mutation& mut = SSP->ProgrammedT1Mutations[SSP->ProgrammedT1MutationsIndexToDo];
                int randomChromosome = GP->random_0or1(GP->mt);
                std::uniform_int_distribution<int> runiform_int_0andNbInds(0,SSP->patchSize[mut.patch_index]);
                int randomIndividual = runiform_int_0andNbInds(GP->mt);
                pop.getPatch(mut.patch_index).getInd(randomIndividual).getHaplo(randomChromosome).toggleT1_Allele(mut.locus, SSP->Habitats[mut.patch_index]);

                SSP->ProgrammedT1MutationsIndexToDo++;
            } else
            {
                break;
            }
        }
    }
}

*/


