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



int SpeciesSpecificParameters::selectNonEmptyPatch(int firstPatchToLookAt, std::vector<int>& PSs, bool increasingOrder)
{
    int r = -1;
    if (firstPatchToLookAt != -1 && PSs[firstPatchToLookAt] > 0)
    {
        r = firstPatchToLookAt;
    } else
    {
        if (increasingOrder)
        {
            for (int patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
            {
                if (PSs[patch_index] > 0)
                {
                    r = patch_index;
                    break;
                }
            }
        } else
        {
            for (int patch_index = GP->PatchNumber - 1 ; patch_index >= 0 ; --patch_index)
            {
                if (PSs[patch_index] > 0)
                {
                    r = patch_index;
                    break;
                }
            }
        }
            
    }
    return r;
}

void SpeciesSpecificParameters::setRandomDistributions()
{
    SSP = this;
    geneticSampler.set_nbRecombinations(TotalRecombinationRate);
    geneticSampler.set_T1_nbMuts(T1_Total_Mutation_rate);
    geneticSampler.set_T2_nbMuts(T2_Total_Mutation_rate);
    geneticSampler.set_T3_nbMuts(T3_Total_Mutation_rate);
    geneticSampler.set_T4_nbMuts(T4_Total_Mutation_rate);
    geneticSampler.set_T56_nbMuts(T56_Total_Mutation_rate);
    geneticSampler.set_T8_mutationStuff(T8_Total_Mutation_rate, T8_MutationRate, T8_map);

    geneticSampler.set_recombinationPosition(RecombinationRate);

    geneticSampler.set_T1_mutationPosition(T1_MutationRate);
    geneticSampler.set_T2_mutationPosition(T2_MutationRate);
    geneticSampler.set_T3_mutationPosition(T3_MutationRate);
    geneticSampler.set_T4_mutationPosition(T4_MutationRate);
    geneticSampler.set_T56_mutationPosition(T56_MutationRate);
    
    SSP = nullptr;

#ifdef DEBUG
    std::cout << "Exiting 'setRandomDistributions' " << std::endl;
#endif
}

void SpeciesSpecificParameters::readSegDiversityFile_includeMainColor(InputReader& input)
{
    if (input.PeakNextElementString() == "default")
    {
        input.skipElement();
        SegDiversityFile_includeMainColor = false;
    } else
    {
        SegDiversityFile_includeMainColor = input.GetNextElementBool();
    }
}

void SpeciesSpecificParameters::readSwapInLifeCycle(InputReader& input)
{
    if (input.PeakNextElementString() == "default")
    {
        input.skipElement();
        if (selectionOn == 0 && Gmap.TotalNbLoci > 1000 && TotalRecombinationRate < 10.0)
        {
            SwapInLifeCycle = true;
        } else
        {
            SwapInLifeCycle = false;
        }
    } else
    {
        SwapInLifeCycle = input.GetNextElementBool();
        if (SwapInLifeCycle && selectionOn != 0)
        {
            std::cout << "You asked for selection to be not only on fecundity and for SwapInLifeCycle. Swapping is only possible if selection is only on fecundity.\n";
            abort();
        }
    }
}

void SpeciesSpecificParameters::readDispMat(InputReader& input)
{
    dispersalData.readDispMat(input);
}

void SpeciesSpecificParameters::readGameteDispersal(InputReader& input)
{
    std::string s = input.GetNextElementString();
    if (s == "no" || s == "n" || s == "N" || s == "No" || s == "NO" || s == "0")
    {
        gameteDispersal = false;
    } else if (s == "yes" || s == "y" || s == "Y" || s == "Yes" || s == "YES" || s == "1")
    {
        gameteDispersal = true;
    } else
    {
        std::cout << "In option --gameteDispersal, expected either 'yes' (or 'y','Y', '1' and a few others) or 'no' (or 'n', 'N', '0' and a few others) but instead received " << s << "\n";
        abort();
    }
    
}

void SpeciesSpecificParameters::readOutputSFSbinSize(InputReader& input)
{
    while (input.IsThereMoreToRead())
    {
        if (input.PeakNextElementString() == "default")
        {
            input.skipElement();
            outputSFSbinSizes.push_back(0.01);
        } else
        {
            double x = input.GetNextElementDouble();            
            if (x <= 0 || x > 1)
            {
                std::cout << "While reading option --outputSFSbinSize, a bin size of " << x << " has been received. A bin size must be great than zero and must be lower (or equal) to one.\n";
                abort();
            }
            outputSFSbinSizes.push_back(x);
        }
    }
}


SpeciesSpecificParameters::SpeciesSpecificParameters(std::string sN, int sI)
: speciesName(sN), speciesIndex(sI)
{
    //std::cout << "Constructor of SpeciesSpecificParameters\n";
}

SpeciesSpecificParameters::SpeciesSpecificParameters()
{
    std::cout << "Default cnstructor of SpeciesSpecificParameters called\n";
}

SpeciesSpecificParameters::~SpeciesSpecificParameters()
{
    //std::cout << "Destructor of SpeciesSpecificParameters\n";
}

/*
template <typename T>
T log(T x, double base)
{
    return std::log(x) / std::log(base);
}
*/

bool SpeciesSpecificParameters::setFromLocusToFitnessMapIndex_DecidingFunction(double sumOfProb, int nbLoci)
{
    assert(this->FitnessMapCoefficient == -9.0 || this->FitnessMapProbOfEvent == -9.0);
    assert(nbLoci > 0);
    
    if (this->FitnessMapCoefficient > 0.0)
    {
        return (nbLoci > 80 && (sumOfProb * nbLoci / 500) > this->FitnessMapCoefficient);
    } else if (this->FitnessMapProbOfEvent > 0.0 && nbLoci > this->FitnessMapMinimNbLoci)
    {
        return sumOfProb > this->FitnessMapProbOfEvent;
    }
    return false;
}

void SpeciesSpecificParameters::setFromLocusToFitnessMapIndex()
{
    assert(this->FromLocusToFitnessMapIndex.size() == 0);
    if (this->FitnessMapProbOfEvent < 0.0)
    {
        std::cout << "In 'setFromLocusToFitnessMapIndex', it appears that the inputted parameter for 'FitnessMapProbOfEvent' is negative. Value received is " << this->FitnessMapProbOfEvent << "\n";
        abort();
    }

    /*
    std::cout << "T56_isMultiplicitySelection = " << T56_isMultiplicitySelection << "\n";
    std::cout << "T1_isSelection = " << T1_isSelection << "\n";
    std::cout << "T1_isMultiplicitySelection = " << T1_isMultiplicitySelection << "\n";
    */
    bool ShouldThereBeSeveralFitnessBlocks = (this->T1_isSelection && this->T1_isMultiplicitySelection) || (this->T56_isSelection && this->T56_isMultiplicitySelection) || this->T2_isSelection;

    double sumOfProb = 0.0; // only used if not whole description
    int FitnessMapIndex = 0;
    int nbLociInPiece = 0;
    if (ShouldThereBeSeveralFitnessBlocks)
    {
        for (int interlocus = 0 ; interlocus < (this->Gmap.TotalNbLoci - 1) ; interlocus++)
        {
            nbLociInPiece++; // only used if not whole description

            // Get Locus info
            int T1_locus = this->Gmap.FromLocusToNextT1Locus(interlocus);
            int T2_locus = this->Gmap.FromLocusToNextT2Locus(interlocus);
            int T3_locus = this->Gmap.FromLocusToNextT3Locus(interlocus);
            int T4_locus = this->Gmap.FromLocusToNextT4Locus(interlocus);
            int T56ntrl_locus = this->Gmap.FromLocusToNextT56ntrlLocus(interlocus);
            int T56sel_locus = this->Gmap.FromLocusToNextT56selLocus(interlocus);
            int locusType;

            
            {
                //std::cout << "interlocus " << interlocus << " / " << this->Gmap.TotalNbLoci - 1 << "\n";
                // Get locusType and long security
                if (interlocus == 0)
                {
                    assert(T1_locus + T2_locus + T3_locus + T4_locus + T56ntrl_locus + T56sel_locus == 1);
                    if (T1_locus == 1)
                    {
                        locusType = 1;
                    } else if (T2_locus == 1)
                    {
                        locusType = 2;
                    } else if (T3_locus == 1)
                    {
                        locusType = 3;
                    } else if (T4_locus == 1)
                    {
                        locusType = 4;
                    } else if (T56ntrl_locus == 1)
                    {
                        locusType = 50;
                    } else if (T56sel_locus == 1)
                    {
                        locusType = 51;
                    } else
                    {
                        std::cout << "Internal error in void SpeciesSpecificParameters::setFromLocusToFitnessMapIndex(): unknown locusType received.\n";
                        abort();
                    }
                } else
                {
                    assert(T1_locus + T2_locus + T3_locus + T4_locus + T56ntrl_locus + T56sel_locus - this->Gmap.FromLocusToNextT1Locus(interlocus - 1) - this->Gmap.FromLocusToNextT2Locus(interlocus - 1) - this->Gmap.FromLocusToNextT3Locus(interlocus - 1) - this->Gmap.FromLocusToNextT4Locus(interlocus - 1) - this->Gmap.FromLocusToNextT56ntrlLocus(interlocus - 1) - this->Gmap.FromLocusToNextT56selLocus(interlocus - 1) == 1);
                    if (T1_locus - this->Gmap.FromLocusToNextT1Locus(interlocus - 1) == 1)
                    {
                        locusType = 1;
                    } else if (T2_locus - this->Gmap.FromLocusToNextT2Locus(interlocus - 1) == 1)
                    {
                        locusType = 2;
                    } else if (T3_locus - this->Gmap.FromLocusToNextT3Locus(interlocus - 1) == 1)
                    {
                        locusType = 3;
                    } else if (T4_locus - this->Gmap.FromLocusToNextT4Locus(interlocus - 1) == 1)
                    {
                        locusType = 4;
                    } else if (T56ntrl_locus - this->Gmap.FromLocusToNextT56ntrlLocus(interlocus - 1) == 1)
                    {
                        locusType = 50;
                    } else if (T56sel_locus - this->Gmap.FromLocusToNextT56selLocus(interlocus - 1) == 1)
                    {
                        locusType = 51;
                    } else
                    {
                        std::cout << "Internal error in 'void SpeciesSpecificParameters::setFromLocusToFitnessMapIndex()'.\n";
                        abort();
                    }
                }
                //assert(locusType == this->FromLocusToTXLocus[interlocus].TType);

                // Has there been a whole description of the map
                if (this->FitnessMapInfo_wholeDescription.size())
                {
                    //std::cout << "this->FitnessMapInfo_wholeDescription.size() = " << this->FitnessMapInfo_wholeDescription.size() << "\n";
                    //std::cout << "FitnessMapIndex = " << FitnessMapIndex << "\n";
                    assert(this->FitnessMapInfo_wholeDescription.size() > FitnessMapIndex);
                    //std::cout << "nbLociInPiece = " << nbLociInPiece << "\n";
                    //std::cout << "this->FitnessMapInfo_wholeDescription["<<FitnessMapIndex<<"] = " << this->FitnessMapInfo_wholeDescription[FitnessMapIndex] << "\n";
                    if (nbLociInPiece == this->FitnessMapInfo_wholeDescription[FitnessMapIndex] && FitnessMapIndex != this->FitnessMapInfo_wholeDescription.size() - 1) // Don't add the last element in FromFitnessMapIndexToTXLocus as it will be done afterward.
                    {
                        FitnessMapIndex++;
                        nbLociInPiece = 0;
                        //std::cout << T1_locus << "\n";
                        FromLocusToTXLocusElement E(T1_locus, T2_locus, T3_locus, T4_locus, T56ntrl_locus, T56sel_locus, locusType);
                        FromFitnessMapIndexToTXLocus.push_back(E); // last locus included of each type
                    }
                } else
                {
                    // Get Recombination Rate info
                    double r;
                    if (locusType == 1 || locusType == 2 || locusType == 50  || locusType == 51)
                    {
                        if (this->RecombinationRate.size() == 1) // which is the case if the rate is constant between any two loci.
                        {
                            r = this->RecombinationRate[0];
                        } else
                        {
                            if (interlocus == 0)
                            {
                                r = this->RecombinationRate[0];
                            } else
                            {
                                assert(this->RecombinationRate.size() > interlocus);
                                r = this->RecombinationRate[interlocus] - this->RecombinationRate[interlocus - 1];
                            }
                        }
                    } else
                    {
                        r = 0.0;
                    }
                    if (locusType == 50 || locusType == 51)
                    {
                        r *= this->FitnessMapT5WeightProbOfEvent;
                    }
                    

                    // Sum probabilities of an event
                    sumOfProb += r;

                    // test if needs to make a new boundary
                    bool avoidTooMuchByteSplitting_mustRunDecidingFun = true;
                    if (locusType == 1 && (T1_locus%8) != 7)
                    {
                        avoidTooMuchByteSplitting_mustRunDecidingFun = false;
                    }
                    if (avoidTooMuchByteSplitting_mustRunDecidingFun && setFromLocusToFitnessMapIndex_DecidingFunction(sumOfProb, nbLociInPiece))
                    {
                        FitnessMapIndex++;
                        sumOfProb = 0.0;
                        FromLocusToTXLocusElement E(T1_locus, T2_locus, T3_locus, T4_locus, T56ntrl_locus, T56sel_locus, locusType);
                        FromFitnessMapIndexToTXLocus.push_back(E); // last locus included of each type
                        nbLociInPiece = 0;
                    }
                }
            } 
        }
    }

    ++FitnessMapIndex;
    NbElementsInFitnessMap = FitnessMapIndex;
    

    if (this->FitnessMapInfo_wholeDescription.size())
    {
        //std::cout << "this->FitnessMapInfo_wholeDescription.size() = " << this->FitnessMapInfo_wholeDescription.size() << "\n";
        //std::cout << "NbElementsInFitnessMap = " << NbElementsInFitnessMap << "\n";
        if (ShouldThereBeSeveralFitnessBlocks) assert(this->FitnessMapInfo_wholeDescription.size() == NbElementsInFitnessMap);
        std::vector<int>().swap(this->FitnessMapInfo_wholeDescription);
    }
        

    


    // Add last element to boundaries
    FromLocusToTXLocusElement E(this->Gmap.T1_nbLoci, this->Gmap.T2_nbLoci, this->Gmap.T3_nbLoci, this->Gmap.T4_nbLoci, this->Gmap.T56ntrl_nbLoci, this->Gmap.T56sel_nbLoci, this->Gmap.getLocusType(this->Gmap.TotalNbLoci-1));
    FromFitnessMapIndexToTXLocus.push_back(E);
    assert(FromFitnessMapIndexToTXLocus.size() == NbElementsInFitnessMap);

    // Create FromLocusToFitnessMapIndex from FromFitnessMapIndexToTXLocus
    if ( ShouldThereBeSeveralFitnessBlocks )
    {
        int previous_b_interlocus = 0;
        for (int i = 0 ; i < FromFitnessMapIndexToTXLocus.size() ; i++)
        {
            auto& b = FromFitnessMapIndexToTXLocus[i];
            assert(b.T1 + b.T2 + b.T3 + b.T4 + b.T56ntrl + b.T56sel <= this->Gmap.TotalNbLoci);
            int b_interlocus = b.T1 + b.T2 + b.T3 + b.T4 + b.T56ntrl + b.T56sel;

            for (int j = previous_b_interlocus  ; j < b_interlocus ; j++)
            {
                FromLocusToFitnessMapIndex.push_back(i);        
            }
            previous_b_interlocus = b_interlocus;
        }
        //std::cout << "FromLocusToFitnessMapIndex.size() = " << FromLocusToFitnessMapIndex.size() << " this->Gmap.TotalNbLoci = " << this->Gmap.TotalNbLoci << "\n";
        assert(FromLocusToFitnessMapIndex.size() == this->Gmap.TotalNbLoci);
        //std::cout << "NbElementsInFitnessMap = " << NbElementsInFitnessMap << "\n";
        //std::cout << "FromLocusToFitnessMapIndex.back() = " << FromLocusToFitnessMapIndex.back() << "\n";
        assert(NbElementsInFitnessMap == FromLocusToFitnessMapIndex.back()+1);
    }
    
    // print to console
    if ( ShouldThereBeSeveralFitnessBlocks )
    {
        assert(NbElementsInFitnessMap >= 1);
        std::cout << "\t\tThe fitnessMap of species '"<<this->speciesName<<"' contains " << NbElementsInFitnessMap << " element\n\n";
    }
    

    /*
    std::cout << "FromFitnessMapIndexToTXLocus";
    std::cout << "("<<FromFitnessMapIndexToTXLocus.size()<<"):\n";
    for (auto& e : FromFitnessMapIndexToTXLocus)
        std::cout << e.T1 + e.T2 + e.T3 << " ";
    std::cout << "\n";

    std::cout << "FromLocusToFitnessMapIndex" << std::flush;
    std::cout << "("<<FromLocusToFitnessMapIndex.size()<<"):\n" << std::flush;
    for (int i = 0 ; i < FromLocusToFitnessMapIndex.size() ; i++)
        std::cout << "FromLocusToFitnessMapIndex["<<i<<"] = " << FromLocusToFitnessMapIndex[i] << "\n" << std::flush;
    */
    

}

/*
void SpeciesSpecificParameters::readT4_printTree(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option '--T4_printTree', the std::string that is read is: " << input.print() << std::endl;
#endif
    OutputFile file(input.GetNextElementString(), T4_printTree);
    T4Tree.indicateOutputFile(&file);
}
*/

void SpeciesSpecificParameters::readSNPfreqCalculationAssumption(InputReader& input)
{
    if (input.PeakNextElementString() == "default")
    {
        // Awefully approximative expectation
        input.skipElement();

        size_t maxEverTotalpatchCapacity = 0;
        for (auto& elem : maxEverpatchCapacity) maxEverTotalpatchCapacity += elem;


        if (this->Gmap.T4_nbLoci)
        {            
            T4SNPfreq_assumeMuchDiversity = Gmap.T4_nbLoci < 1e4 ? true : 4 * maxEverTotalpatchCapacity * T4_Total_Mutation_rate / Gmap.T4_nbLoci > 0.001;
        } else
        {
            T4SNPfreq_assumeMuchDiversity = false; // whatever
        }
            
        TxSNPfreq_assumeMuchDiversity = this->Gmap.TotalNbLoci - this->Gmap.T4_nbLoci < 1e4 ? true : maxEverTotalpatchCapacity > 5e3 && ((this->Gmap.TotalNbLoci - this->Gmap.T4_nbLoci) / this->Gmap.TotalNbLoci) < 0.5;
        
    } else
    {
        T4SNPfreq_assumeMuchDiversity = input.GetNextElementBool();
        TxSNPfreq_assumeMuchDiversity = input.GetNextElementBool();
    }
}



void SpeciesSpecificParameters::readT4_nbRunsToPlaceMutations(InputReader& input)
{
    if (input.PeakNextElementString() == "default")
    {
        input.skipElement();
        T4_nbRunsToPlaceMutations = 200;
    } else
    {
        auto x = input.GetNextElementInt();
        
        if (x <= 0)
        {
            std::cout << "For option '--T4_nbRunsToPlaceMutations', only positive values are accepted. Received " << x << "\n";
            abort();
        }

        T4_nbRunsToPlaceMutations = x;
    }
}

void SpeciesSpecificParameters::readGeneticSampling_withWalker(InputReader& input)
{
    geneticSampling_withWalker = input.GetNextElementBool();
}

void SpeciesSpecificParameters::readIndividualSampling_withWalker(InputReader& input)
{
    if (input.PeakNextElementString() == "default")
    {
        input.skipElement();
        double meanPatchCapacity = 0.0;
        for (uint32_t patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
        {
            meanPatchCapacity += SSP->patchCapacity[patch_index];
        }
        meanPatchCapacity /= GP->PatchNumber;
        /*
        The walker seem to be buggy. So until, I fix it, I won't use it by default
        if (meanPatchCapacity > 1e4)
        {
            individualSampling_withWalker = true;       
        } else
        {
            individualSampling_withWalker = false;
        }*/
        individualSampling_withWalker = false;
    } else
    {
        individualSampling_withWalker = input.GetNextElementBool();
    }
}

void SpeciesSpecificParameters::readSampleSeq_info(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option '--sampleSeq_file', the std::string that is read for for a given species is: " << input.print() << std::endl;
#endif 

    sampleSequenceDataContainer.readInput(input);    
}

void SpeciesSpecificParameters::readForcedMigration(InputReader& input)
{
    this->forcedMigration.readInput(input);
}

void SpeciesSpecificParameters::readStochasticMigration(InputReader& input)
{
    if (input.PeakNextElementString() == "default")
    {
        input.skipElement();
        this->isStochasticMigration = true;
    } else
    {
        this->isStochasticMigration = input.GetNextElementBool();
    }
}




void SpeciesSpecificParameters::readLoci(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option '--Loci', the std::string that is read is: " << input.print() << std::endl;
#endif

    if (input.PeakNextElementString() == "NA" || input.PeakNextElementString() == "0")
    {
        std::cout << "In option --L (--Loci), received 'NA' or '0' at the beginning of input. I suppose you meant that you don't want any genetics. Do you really want that? SimBit is a population genetics software but I can understand that you might be interested in demographic processes only. It is probably easy to allow SimBit to run simulations without genetic architecture but I did not made sure it would not cause bugs. So, if you don't want any genetics, I suggest you just use a single locus without selection on it (do '--L T5 1 --T5_mu A 0 --T5_fit multfitA 1' for example) and you just don't ask for outputs concerning the genetics.\n";
        abort();
    }

    Gmap.readLoci(input);


    if (Gmap.T4_nbLoci > 0)
    {
        if (this->nbSubGenerationsPerGeneration != 1)
        {
            std::cout << "You have asked for loci of type T4 via option --Loci (--L). You have also asked for having a number of subGeneration per generation different than 1 via option --nbSubGenerations (--nbSubGens). Sorry, the coalescent tree used to track T4 loci assumes on subgeneration per generation. It would be easy for Remi to get rid of this assumption. Please let me know!\n";
            abort();
        }

        /*int nbDemes = patchCapacity.size();
        for (auto& pc : __patchCapacity)
        {
            if (pc.size() != nbDemes)
            {
                std::cout << "You have asked for loci of type T4 via option --Loci (--L). You have also asked for having change in the number of demes over time (via --PN). Sorry, the coalescent tree used to track T4 loci assumes no change in the carrying capacity or in the number of demes over time. This assumption is simply caused by Remi being too lazy to make it more flexible!\n";
                abort();
            }

            for (uint32_t patch_index = 0 ; patch_index < nbDemes ; ++patch_index)
            {
                if (pc[patch_index] != patchCapacity[patch_index])
                {
                    std::cout << "You have asked for loci of type T4 via option --Loci (--L). You have also asked for having change in the carrying patchCapacity over time.(via --N). Sorry, the coalescent tree used to track T4 loci assumes no change in the carrying capacity or in the number of demes over time. This assumption is simply caused by Remi being too lazy to make it more flexible!\n";
                    abort();
                }
            }
        }*/
    }
    
    

    ///////////////////////////
    ///// T5 and T6 stuff /////
    ///////////////////////////
    // T6 will use T56_FitnessEffects. It is a bit misleading but it is useless to rename it
    assert(T56_FitnessEffects.size() == MaxEverHabitat+1);
    
    // Remove from T5_fit everything that is not under selection
    for (uint32_t habitat = 0 ; habitat <= this->MaxEverHabitat ; habitat++)
    {
        auto& fits = T56_FitnessEffects[habitat]; 

        if (T56_isMultiplicitySelection)
        {
            //std::cout << "fits.size() = " << fits.size() << "\n";
            //std::cout << "T5_nbLoci = " << T5_nbLoci << "\n";
            //std::cout << "T6_nbLoci = " << T6_nbLoci << "\n";
            assert(fits.size() == Gmap.T56_nbLoci);

            //would be faster with erase-remove idiom
            fits.erase(
                std::remove_if(
                    fits.begin(),
                    fits.end(),
                    [&fits, this](fitnesstype& elem)
                    {
                        uint32_t locus = &elem - fits.data();
                        return Gmap.isT56neutral(locus);
                    }   
                ),
                fits.end()
            );

            //std::cout << "fits.size() = " << fits.size() << "\n";
            //std::cout << "T5sel_nbLoci = " << T5sel_nbLoci << "\n";
            assert(fits.size() == Gmap.T56sel_nbLoci);
        } else
        {
            //std::cout << "fits.size() = " << fits.size() << "\n";
            //std::cout << "T5_nbLoci = " << T5_nbLoci << "\n";
            assert(fits.size() == Gmap.T56_nbLoci*2);

            //would be faster with erase-remove idiom
            fits.erase(
                std::remove_if(
                    fits.begin(),
                    fits.end(),
                    [&fits, this](fitnesstype& elem)
                    {
                        uint32_t locus = (&elem - fits.data())/2;
                        return Gmap.isT56neutral(locus);
                    }   
                ),
                fits.end()
            );

            //std::cout << "fits.size() = " << fits.size() << "\n";
            //std::cout << "T5sel_nbLoci = " << T5sel_nbLoci << "\n";
            assert(fits.size() == Gmap.T56sel_nbLoci * 2);
        }
    }
}

void SpeciesSpecificParameters::readT56_compress(InputReader& input)
{
    Gmap.readT56Compression(input);
}


void SpeciesSpecificParameters::readT56_approximationForNtrl(InputReader& input)
{
    if (input.PeakNextElementString() == "default")
    {
        input.skipElement();
        T56_approximationForNtrl = 0.0;
    } else
    {
        T56_approximationForNtrl = input.GetNextElementDouble();
    }


        
    if (T56_approximationForNtrl < 0.0 || T56_approximationForNtrl > 1.0)
    {
        std::cout << "In --T56_approximationForNtrl, the value received is " << T56_approximationForNtrl << ". Sorry only values greater or equal to 0.0 and smaller or equal to 1.0 make sense. Default is 0.0, which means 'no approximation'. The approximation is computed as following. For a fitness value 'fit' greater than 1.0, SimBit will consider the locus to be neutral if 'fit' is greater than '1+T56_approximationForNtrl'. For a fitness value 'fit' lower than 1.0, SimBit will consider the locus to be neutral if 'fit' is lower than '1 / (1+T56_approximationForNtrl)'.\n";
        abort();
    }
}

void SpeciesSpecificParameters::resetT56_freqThresholToDefault()
{
    T56ntrl_frequencyThresholdForFlippingMeaning = 1.0;
    T56sel_frequencyThresholdForFlippingMeaning = 1.0;
/*
    if (TotalpatchCapacity > 50)
    {
        T56ntrl_frequencyThresholdForFlippingMeaning = 0.95;
        T56sel_frequencyThresholdForFlippingMeaning = 0.95;   
    }
    
    if (TotalpatchCapacity > 500)
    {
        T56ntrl_frequencyThresholdForFlippingMeaning = 0.9;
        T56sel_frequencyThresholdForFlippingMeaning = 0.9;
    }
    
    if (TotalpatchCapacity > 1000)
    {
        T56ntrl_frequencyThresholdForFlippingMeaning = 0.8;
        T56sel_frequencyThresholdForFlippingMeaning = 0.8;   
    }
    
    if (TotalpatchCapacity > 3000)
    {
        T56ntrl_frequencyThresholdForFlippingMeaning = 0.7;
        T56sel_frequencyThresholdForFlippingMeaning = 0.7;
    }

    if (TotalpatchCapacity > 10000)
    {
        T56ntrl_frequencyThresholdForFlippingMeaning = 0.6;
        T56sel_frequencyThresholdForFlippingMeaning = 0.6;
    }*/
    //std::cout << "T56ntrl_frequencyThresholdForFlippingMeaning = " << T56ntrl_frequencyThresholdForFlippingMeaning << "\n";
    //std::cout << "T56sel_frequencyThresholdForFlippingMeaning = " << T56sel_frequencyThresholdForFlippingMeaning << "\n";
}

void SpeciesSpecificParameters::readT56_freqThreshold(InputReader& input)
{
    if (input.PeakNextElementString() == "default")
    {
        input.skipElement();
        T5_freqThresholdWasSetToDefault = true;
        this->resetT56_freqThresholToDefault();
    } else
    {
        T5_freqThresholdWasSetToDefault = false;
        T56ntrl_frequencyThresholdForFlippingMeaning = input.GetNextElementDouble();
        T56sel_frequencyThresholdForFlippingMeaning = input.GetNextElementDouble();
    }
        

    // Securities
    if (T56ntrl_frequencyThresholdForFlippingMeaning <= 0.5 )
    {
        std::cout << "In --T5_freqThreshold, received the value " << T56ntrl_frequencyThresholdForFlippingMeaning << " for T5ntrl. Sorry a value lower or equal to 0.5 makes no sense.\n";
        abort();
    }
    if (T56sel_frequencyThresholdForFlippingMeaning <= 0.5 )
    {
        std::cout << "In --T5_freqThreshold, received the value " << T56sel_frequencyThresholdForFlippingMeaning << " for T5sel. Sorry a value lower or equal to 0.5 makes no sense.\n";
        abort();
    }
}

/*void SpeciesSpecificParameters::readreverseFixedT5selMuts(InputReader& input)
{
    std::string x = input.GetNextElementString();
    if (x == "false" || x == "0" || x == "f")
    {
        reverseFixedT5selMuts = false;
    } else if (x == "true" || x == "1" || x == "t")
    {
        reverseFixedT5selMuts = true;
    } else
    {
        std::cout << "In option --reverseFixedT5selMuts, was expecting either 'true' (or '1' or 't') or 'false' (or '0' or 'f') but instead received " << x << "\n";
        abort();
    }
}*/


/*void SpeciesSpecificParameters::readT1mutsDirectional(InputReader& input)
{
    std::string x = input.GetNextElementString();
    if (x == "false" || x == "0" || x == "f")
    {
        T1mutsDirectional = false;
    } else if (x == "true" || x == "1" || x == "t")
    {
        T1mutsDirectional = true;
    } else
    {
        std::cout << "In option --T1mutsDirectional, was expecting either 'true' (or '1' or 't') or 'false' (or '0' or 'f') but instead received " << x << "\n";
        abort();
    }
}*/


/*void SpeciesSpecificParameters::readT5mutsDirectional(InputReader& input)
{
    std::string x = input.GetNextElementString();
    if (x == "false" || x == "0" || x == "f")
    {
        T5mutsDirectional = false;
    } else if (x == "true" || x == "1" || x == "t")
    {
        T5mutsDirectional = true;
        T5sel_knownMutFixed.resize(T5sel_nbLoci, false);

    } else
    {
        std::cout << "In option --T5mutsDirectional, was expecting either 'true' (or '1' or 't') or 'false' (or '0' or 'f') but instead received " << x << "\n";
        abort();
    }
    assert(MaxEverHabitat >= 0);
    T5sel_fitnessEffectOfknownMutFixed.resize(MaxEverHabitat+1,1.0);
}*/


void SpeciesSpecificParameters::readCloningRate(InputReader& input)
{
    this->cloningRate = input.GetNextElementDouble();
    if (this->cloningRate > 1.0 || this->cloningRate < 0.0)
    {
        std::cout << "In --cloningRate, received a cloning rate for species '"<<this->speciesName<<"' that is either lower than zero or greater than one. cloningRate receveid is " << this->cloningRate << "\n";
        abort();
    }
}

void SpeciesSpecificParameters::readSelfingRate(InputReader& input)
{
   this->selfingRate = input.GetNextElementDouble();
 
    if ((this->selfingRate > 1.0 || this->selfingRate < 0.0) && this->selfingRate != -1.0)
    {
        std::cout << "In --selfingRate, received a selfing rate that is either lower than zero or greater than one. selfing rate receveid is " << this->selfingRate << ". Note that a selfingRate of -1 (default) refers to a simple Wright-Fisher model where the selfing rate is at 1/2N (a little bit lower with migration)\n";
        abort();
    }
}

void SpeciesSpecificParameters::readShrinkT56EveryNGeneration(InputReader& input)
{
    if (input.PeakNextElementString() == "default")
    {
        T56_memManager.set_attempt_shrinking_every_N_generation(1000);
        input.skipElement();
    } else
    {
        T56_memManager.set_attempt_shrinking_every_N_generation(input.GetNextElementInt());
    }
}

void SpeciesSpecificParameters::readMatingSystem(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option '--matingSystem', the std::string that is read is: " << input.print() << std::endl;
#endif
    std::string s = input.GetNextElementString();
    if (s == "H" || s == "h")
    {
        this->malesAndFemales = false;
        this->sexRatio = 0.0;
    } else if (s == "MF" || s == "FM" || s == "fm"  || s == "mf")
    {
        this->malesAndFemales = true;
        this->sexRatio = input.GetNextElementDouble();

        if (this->sexRatio >= 1.0 || this->sexRatio <= 0.0)
        {
            std::cout << "For option '--matingSystem', the mating system chosen is 'fm' (males and females). The sex ratio received is " << this->sexRatio << ". A sex ratio must be between 0 and 1 (non-bounded). If set at zero, all individuals are female and no reproduction is possible. If set at one, all individuals are males and no reproduction is possible. If you want hermaphrodites, please indicate 'h' (or 'H') and no sex-ratio\n";
            abort();
        }


        assert(SSP->__patchCapacity.size() == GP->__GenerationChange.size());
        assert(GP->__PatchNumber.size() == GP->__GenerationChange.size());

        if (SSP->fecundityForFitnessOfOne == -1.0)
        {
            for (int generation_index = 0 ; generation_index < GP->__GenerationChange.size() ; generation_index++)
            {
                assert(SSP->__patchCapacity[generation_index].size() == GP->__PatchNumber[generation_index]);
                for (double N : SSP->__patchCapacity[generation_index])
                {
                    int nbMales = (int) (sexRatio * N + 0.5);
                    int nbFemales = N - nbMales;
                    //std::cout << "nbMales = "<<nbMales<<" nbFemales = "<<nbFemales<<"\n";
                    if (nbMales <= 0 || nbFemales <= 0)
                    {
                        std::cout << "For option '--matingSystem', the mating system chosen is 'fm' (males and females). The sex ratio received for species "<< this->speciesName<< " is " << this->sexRatio << ". At the " << generation_index << "th generation index (after generation " << GP->__GenerationChange[generation_index] << "), the patchCapacity multiplied by the sex ratio lead to either the number of males or the number of females to be zero (nbMales = "<<nbMales<<" nbFemales = "<<nbFemales<<"). This error message would not appear if the fecundity (from option --fec) was set to -1.0 (fecundityForFitnessOfOne = "<<SSP->fecundityForFitnessOfOne<<"). If fecundity was not set to -1.0, then the patch size would just reach 0. SimBit enforces that there is at least one female and one male per patch when fecundity is different from -1.0. (this error message was written at a time when sex ratio is deterministic. If at some point SimBit evovles to have stochasticity in the sex ratio, then I should remove this error message and deal with cases where the sex ratio is 0 or 1).\n";
                        abort(); 
                    }
                }
            }
        }
    } else
    {
        std::cout << "For option '--matingSystem', the mating system received is " << s << ". Only types 'fm' (or 'mf', 'FM','MF') and 'h' (or 'H') are recognized. Sorry.\n";
        abort();
    }
}

void SpeciesSpecificParameters::readDispWeightByFitness(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option '--DispWeightByFitness', the std::string that is read is: " << input.print() << std::endl;
#endif

    // Note that --fec has already been received but not --m

    if (input.PeakNextElementString() == "default")
    {
        input.skipElement();
        if (this->fecundityForFitnessOfOne == -1)
        {
            this->DispWeightByFitness = false;
        } else
        {
            this->DispWeightByFitness = true;
        }

    } else
    {
        this->DispWeightByFitness = input.GetNextElementBool();
    }
    
        
    
    if (fecundityForFitnessOfOne != -1 && !DispWeightByFitness)
    {
        std::cout << "For option '--DispWeightByFitness', received 'false' and fecundityForFitnessOfOne is not -1 (is " << fecundityForFitnessOfOne << "). When the fecundity is not -1 (which means infinite), then dispersal must necessarily be weigthed by fitness (or I don't really know what it would mean otherwise).\n";
        abort();
    }
}

void SpeciesSpecificParameters::readPloidy(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option '--ploidy', the std::string that is read is: " << input.print() << std::endl;
#endif    
    this->ploidy = input.GetNextElementInt();
    if (this->ploidy!=2)
    {
        std::cout << "Sorry, for the moment ploidy can only be equal to 2" << std::endl;
        abort();
    }
}

void SpeciesSpecificParameters::readT4_paintedHaplo_ignorePatchSizeSecurityChecks(InputReader& input)
{
    if (input.PeakNextElementString() == "default")
    {
        input.skipElement();
        T4_paintedHaplo_shouldIgnorePatchSizeSecurityChecks = false;
    } else
    {
        T4_paintedHaplo_shouldIgnorePatchSizeSecurityChecks = input.GetNextElementBool();
    }
}


ResetGeneticsEvent_A SpeciesSpecificParameters::readResetGenetics_readEventA(InputReader& input, int& eventIndex)
{
    // generation
    int generation = input.GetNextElementInt();
    if (generation < 0)
    {
        std::cout << "For option '--resetGenetics', the generation must be 0 or larger. Received generation = " << generation << ".\n";
        abort();
    }
    int generation_index = std::upper_bound(GP->__GenerationChange.begin(), GP->__GenerationChange.end(), generation) - GP->__GenerationChange.begin() - 1;

    // loci type
    std::string locusType = input.GetNextElementString();
    if (locusType != "T1" && locusType != "T2" && locusType != "T3"  && locusType != "T5")
    {
        std::cout << "For option '--resetGenetics', " << eventIndex << "th event, for the "<<eventIndex<<"th event (generation "<<generation<<"; species "<<SSP->speciesName<<"), the locus type is  '" << locusType <<"'. Sorry only mutation types 'T1', 'T2', 'T3' and 'T5' (but not T4) are accepted.\n";
        abort();
    }

    // mutation type
    char mutationType;
    if (locusType == "T1" || locusType == "T5")
    {
        std::string string_mutationType = input.GetNextElementString();
        if (string_mutationType == "setTo0")
        {
            mutationType = 0;
            /*if (SSP->T5mutsDirectional)
            {
                std::cout << "For option --resetGenetics, you indicated a mutation type 'setTo0'. It is impossible to ask to set a T5 allele to 0 when the model chosen is direction. The option --T5mutsDirectional was indeed set to 'true' (either by the user or by default). Note that it would be relatively easy to allow that. Just ask Remi to make a few modifications. Otherwise, you could set --T1mutsDirectional (and/or --T5mutsDirectional) to false (which might makes things slower).\n";
                abort();
            }*/
        } else if (string_mutationType == "setTo1")
        {
            mutationType = 1;
        } else if (string_mutationType == "toggle")
        {
            mutationType = 2;
            /*if (SSP->T5mutsDirectional)
            {
                std::cout << "For option --resetGenetics, you indicated a mutation type 'toggle'. This means that it could happen that resetGenetics would attempt to convert a 1 into a 0. It is impossible to ask to set a T5 or T1 allele to 0 when the model chosen is direction. The option --T5mutsDirectional was indeed set to 'true' (either by the user or by default). Note that it would be relatively easy to allow that. Just ask Remi to make a few modifications. Otherwise, you could set --T1mutsDirectional (and/or --T5mutsDirectional) to false (which might makes things slower).\n";
                abort();
            }*/
        } else
        {
            std::cout << "For option '--resetGenetics', " << eventIndex << "th event, for the "<<eventIndex<<"th event (generation "<<generation<<"; species "<<SSP->speciesName<<") for loci of type "<<locusType<<", the mutation type received is '" << string_mutationType <<"'. Sorry only mutation types 'setTo0', 'setTo1' and 'toggle' are accepted.\n";
            abort();
        }
    } else
    {
        mutationType = -1;
    }

    // Loci indices
    std::string lociKeyword = input.GetNextElementString();
    if (lociKeyword != "lociList" && lociKeyword != "allLoci")
    {
        std::cout << "For option '--resetGenetics', " << eventIndex << "th event, for the "<<eventIndex<<"th event (generation "<<generation<<"; species "<<SSP->speciesName<<") expected the keyword 'lociList' or 'allLoci' but instead received "<<lociKeyword<<".\n";
            abort();
    }

    std::vector<int> T1loci;
    std::vector<int> T2loci;
    std::vector<int> T3loci;
    std::vector<int> T5loci;
    if (lociKeyword == "lociList")
    {
        while (input.PeakNextElementString() != "haplo")
        {
            if (locusType == "T1")
            {
                int T1locus = input.GetNextElementInt();
                if (T1locus < 0)
                {
                    std::cout << "For option --resetGenetics received T1locus index of "<<T1locus<<" for an event happening at generation "<<generation<<".\n";
                    abort();
                }
                if (T1locus >= SSP->Gmap.T1_nbLoci)
                {
                    std::cout << "For option --resetGenetics for an event happening at generation "<<generation<<", received T1locus index of "<<T1locus<<" while there are only "<<SSP->Gmap.T1_nbLoci<<" T1 loci (in bits) in total. As a reminder, the first locus has index 0.\n";
                    abort();
                }
                T1loci.push_back(T1locus);
            } else if (locusType == "T2")
            {
                int T2locus = input.GetNextElementInt();
                if (T2locus < 0)
                {
                    std::cout << "For option --resetGenetics received T2 locus index of "<<T2locus<<" for an event happening at generation "<<generation<<".\n";
                    abort();
                }
                if (T2locus >= SSP->Gmap.T2_nbLoci)
                {
                    std::cout << "For option --resetGenetics for an event happening at generation "<<generation<<", received T2 locus index of "<<T2locus<<" while there are only "<<SSP->Gmap.T2_nbLoci<<" T2 loci in total. As a reminder, the first locus has index 0.\n";
                    abort();
                }
                T2loci.push_back(T2locus);
            } else if (locusType == "T3")
            {
                int T3locus = input.GetNextElementInt();
                if (T3locus < 0)
                {
                    std::cout << "For option --resetGenetics received T3locus index of "<<T3locus<<" for an event happening at generation "<<generation<<".\n";
                    abort();
                }
                if (T3locus >= SSP->Gmap.T3_nbLoci)
                {
                    std::cout << "For option --resetGenetics for an event happening at generation "<<generation<<", received T3 locus index of "<<T3locus<<" while there are only "<<SSP->Gmap.T3_nbLoci<<" T3 loci in total. As a reminder, the first locus has index 0.\n";
                    abort();
                }
                T3loci.push_back(T3locus);
            } else if (locusType == "T5")
            {
                int T5locus = input.GetNextElementInt();
                if (T5locus < 0)
                {
                    std::cout << "For option --resetGenetics received T5locus index of "<<T5locus<<" for an event happening at generation "<<generation<<".\n";
                    abort();
                }
                if (T5locus >= SSP->Gmap.T56_nbLoci)
                {
                    std::cout << "For option --resetGenetics for an event happening at generation "<<generation<<", received T5locus index of "<<T5locus<<" while there are only "<<SSP->Gmap.T56_nbLoci<<" T5 loci in total. As a reminder, the first locus has index 0.\n";
                    abort();
                }
                T5loci.push_back(T5locus);
            }
        }
    } else // lociKeyword == "allLoci"
    {
        if (locusType == "T1")
        {
            for (int T1locus = 0 ; T1locus < SSP->Gmap.T1_nbLoci; T1locus++)
            {
                T1loci.push_back(T1locus);
            }
        } else if (locusType == "T2")
        {
            for (int T2locus = 0 ; T2locus < SSP->Gmap.T2_nbLoci; T2locus++)
            {
                T2loci.push_back(T2locus);
            }
        } else if (locusType == "T3")
        {
            for (int T3locus = 0 ; T3locus < SSP->Gmap.T3_nbLoci; T3locus++)
            {
                T3loci.push_back(T3locus);
            }
        } else if (locusType == "T5")
        {
            for (int T5locus = 0 ; T5locus < SSP->Gmap.T56_nbLoci; T5locus++)
            {
                T5loci.push_back(T5locus);
            }
        }
    }

    // haplotypes concerned
    std::string haploString = input.GetNextElementString();
    if (haploString != "haplo")
    {
        std::cout << "For option '--resetGenetics', " << eventIndex << "th event, for the "<<eventIndex<<"th event (generation "<<generation<<"; species "<<SSP->speciesName<<"), expected the keyword 'haplo' but instead received '" << haploString <<"'.\n";
        abort();
    }
    std::string haploDescription = input.GetNextElementString();
    std::vector<int> haplotypes;
    if (haploDescription == "0" || haploDescription == "f" || haploDescription == "F")
    {
        haplotypes.push_back(0);
    } else if (haploDescription == "1" || haploDescription == "m" || haploDescription == "M")
    {
        haplotypes.push_back(1);
    } else if (haploDescription == "both")
    {
        haplotypes.push_back(0);
        haplotypes.push_back(1);
    }



    // patches and individuals
    std::vector<std::vector<int>> individuals;
    std::vector<int> patches;
    while (input.IsThereMoreToRead() && input.PeakNextElementString().substr(0,5) != "event")
    {
        // patch
        std::string patchString = input.GetNextElementString();
        if (patchString != "patch")
        {
            std::cout << "For option '--resetGenetics', " << eventIndex << "th event, for the "<<eventIndex<<"th event (generation "<<generation<<"; species "<<SSP->speciesName<<"), expected the keyword 'patch' but instead received '" << patchString <<"'.\n";
            abort();
        }

        

        int patch_index = input.GetNextElementInt();
        if (patch_index < 0 || patch_index > GP->__PatchNumber[generation_index])
        {
            std::cout << "For option '--resetGenetics', " << eventIndex << "th event, for the "<<eventIndex<<"th event (generation "<<generation<<"; species "<<SSP->speciesName<<"), after the keyword 'patch' received the patch index '" << patch_index <<"'. At the generation "<<generation<<" there are "<<GP->__PatchNumber[generation_index]<< "patches. The patch index must be lower than the number of patches (because it is a zero based counting) and positive.\n";
            abort();
        }
        patches.push_back(patch_index);


        // individuals
        std::string modeOfEntryOfIndividuals = input.GetNextElementString();

        if (modeOfEntryOfIndividuals == "allInds")
        {

            assert(SSP->__patchCapacity[generation_index].size() == GP->__PatchNumber[generation_index]);

            std::vector<int> onePatchIndviduals;
            onePatchIndviduals.reserve(SSP->__patchCapacity[generation_index][patch_index]);
            for (int ind_index = 0 ; ind_index < SSP->__patchCapacity[generation_index][patch_index]; ind_index++)
            {
                onePatchIndviduals.push_back(ind_index);
            }
            individuals.push_back(onePatchIndviduals);

        } else if (modeOfEntryOfIndividuals == "indsList")
        {
            std::vector<int> onePatchIndviduals;
            while (input.IsThereMoreToRead() && input.PeakNextElementString().substr(0,5) != "event" && input.PeakNextElementString() != "patch")
            {
                onePatchIndviduals.push_back(input.GetNextElementInt());
            }
            if (onePatchIndviduals.size() == 0)
            {
                std::cout << "For option '--resetGenetics', " << eventIndex << "th event, for the "<<eventIndex<<"th event (generation "<<generation<<"; species "<<SSP->speciesName<<"), for patch index "<< patch_index <<" there seems to have no individuals indicated.\n";
                abort();
            }
            individuals.push_back(onePatchIndviduals);


        } else
        {
            std::cout << "For option '--resetGenetics', " << eventIndex << "th event, for the "<<eventIndex<<"th event (generation "<<generation<<"; species "<<SSP->speciesName<<"), for patch index "<< patch_index <<" expected mode of entry of individuals (either 'indsList' or 'allInds') but instead received '" << modeOfEntryOfIndividuals <<"'.\n";
            abort();
        }
    }
    if (patches.size() == 0)
    {
        std::cout << "For option '--resetGenetics', " << eventIndex << "th event, for the "<<eventIndex<<"th event (generation "<<generation<<"; species "<<SSP->speciesName<<"), there seems to have no patch indicated.\n";
        abort();
    }


    return ResetGeneticsEvent_A(
        generation,
        mutationType,
        T1loci,
        T2loci,
        T3loci,
        T5loci,
        patches,
        haplotypes,
        individuals
    );
}

ResetGeneticsEvent_B SpeciesSpecificParameters::readResetGenetics_readEventB(InputReader& input, int& eventIndex)
{
    // generation
    int generation = input.GetNextElementInt();
    if (generation < 0)
    {
        std::cout << "For option '--resetGenetics', the generation must be 0 or larger. Received generation = " << generation << ".\n";
        abort();
    }
    int generation_index = std::upper_bound(GP->__GenerationChange.begin(), GP->__GenerationChange.end(), generation) - GP->__GenerationChange.begin() - 1;

    // patch index
    int patch_index = input.GetNextElementInt();
    
    if (patch_index < 0 || patch_index >= GP->__PatchNumber[generation_index])
    {
        std::cout << "For option '--resetGenetics', " << eventIndex << "th event, for the "<<eventIndex<<"th event (generation "<<generation<<"; species "<<SSP->speciesName<<"), after the keyword 'patch' received the patch index '" << patch_index <<"'. At the generation "<<generation<<" there are "<<GP->__PatchNumber[generation_index]<< "patches. The patch index must be lower than the number of patches (because it is a zero based counting) and positive.\n";
        abort();
    }

    // Type
    std::vector<std::string> typeNames;
    std::vector<unsigned> howManys;
    while (input.IsThereMoreToRead() && input.PeakNextElementString().substr(0,5) != "event")
    {
        std::string typeName = input.GetNextElementString();
        if (this->individualTypes.find(typeName) == this->individualTypes.end())
        {
            std::cout << "For option '--resetGenetics', " << eventIndex << "th event, for the "<<eventIndex<<"th event (generation "<<generation<<"; species "<<SSP->speciesName<<"), indicated individualType "<< typeName << " but no individual type with this name has been defined in option --indTypes\n";
            abort();
        }
        typeNames.push_back(typeName);
        howManys.push_back((unsigned)input.GetNextElementInt());
    }

    // Security
    auto sumNbInds = std::accumulate(howManys.begin(), howManys.end(),0);
    if (sumNbInds > this->__patchCapacity[generation_index][patch_index])
    {
        std::cout << "For option '--resetGenetics', " << eventIndex << "th event, for the "<<eventIndex<<"th event (generation "<<generation<<"; species "<<SSP->speciesName<<"). Indicated a total of "<<sumNbInds<<" individuals but the carrying capacity will only be " << this->__patchCapacity[generation_index][patch_index] << " at generation "<< generation << " (generation index: "<< generation_index <<")\n";
        abort();
    }

    // return
    return ResetGeneticsEvent_B(
        generation,
        typeNames,
        howManys,
        patch_index
    );


}

void SpeciesSpecificParameters::readResetGenetics(InputReader& input)
{
    #ifdef DEBUG
    std::cout << "For option '--resetGenetics', the std::string that is read is: " << input.print() << std::endl;
#endif 

    if (input.PeakNextElementString() == "default")
    {
        input.skipElement();
        return;
    }


    // input example with eventA's --resetGenetics @S0 event 100 T1 setTo0 lociList 2 4 6 8 haplo both patch 0 allInds patch 2 indsList 1 2 3 4 5 
    //                                          ^generation

    int eventIndex = 0;
    while (input.IsThereMoreToRead())
    {
        std::string eventString = input.GetNextElementString();
        if (eventString == "eventA")
        {
            resetGenetics.addEvent(readResetGenetics_readEventA(input, eventIndex));
        } else if (eventString == "eventB")
        {
            resetGenetics.addEvent(readResetGenetics_readEventB(input, eventIndex));
        } else
        {
            std::cout << "For option '--resetGenetics', for the "<<eventIndex<<"th event (species "<<SSP->speciesName<<"), expected the keyword 'eventA' or 'eventB' but got "<<eventString<<"\n";
            abort();
        }
        eventIndex++;
    }
}

void SpeciesSpecificParameters::readfecundityDependentOfFitness(InputReader& input)
{
    if (input.PeakNextElementString() == "default")
    {
        input.skipElement();
        fecundityDependentOfFitness = true;
    } else
    {
        fecundityDependentOfFitness = input.GetNextElementBool();
    }
}

void SpeciesSpecificParameters::readfecundityForFitnessOfOne(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option '--fec', the std::string that is read is: " << input.print() << std::endl;
#endif
    this->fecundityForFitnessOfOne = input.GetNextElementDouble();
    //std::cout << "fecundity set to " << this->fecundityForFitnessOfOne << "\n";

    if (this->fecundityForFitnessOfOne > 10000.0)
    {
        std::cout << "this->fecundityForFitnessOfOne received is greater than " << 10000.0 << " (is " << this->fecundityForFitnessOfOne << "). If you want the patch to always be at carrying capacity, just use fecundityForFitnessOfOne = -1. It will make everything faster. Change the code if you really want to use such high fecundity\n";
        abort();
    } else if (this->fecundityForFitnessOfOne < 0.0)
    {
        if (this->fecundityForFitnessOfOne != -1.0)
        {
            std::cout << "this->fecundityForFitnessOfOne received is negative but different from -1.0 (which would mean patch size is always at carrying capacity).\n";
            abort();
        }
    }
}

/*void SpeciesSpecificParameters::readResetTrackedT1Muts(InputReader& input)
{
    if (input.PeakNextElementString() == "default")
    {
        input.skipElement();
        if (T1_isSelection)
        {
            if (T1_isMultiplicitySelection)
            {
                recomputeLociOverWhichFitnessMustBeComputedEveryHowManyGenerations = -1;
            } else
            {
                
                double averageMu = T1_MutationRate.back() / T1_nbLoci;
                
                double a = 0.0;
                for (int i = 1 ; i < TotalpatchCapacity ; i++)
                {
                    a += 1/i;
                }
                //std::cout << "TotalpatchCapacity = " << TotalpatchCapacity << " averageMu = " << averageMu << " a = " <<  a << "\n";
                //std::cout << "TotalpatchCapacity * averageMu * 4 * a = " << TotalpatchCapacity * averageMu * 4 * a << "\n";
                if (TotalpatchCapacity * averageMu * 4 * a < 1.0)
                {
                    recomputeLociOverWhichFitnessMustBeComputedEveryHowManyGenerations = 20;
                } else
                {
                    recomputeLociOverWhichFitnessMustBeComputedEveryHowManyGenerations = -1;
                }
            }
        } else
        {
            recomputeLociOverWhichFitnessMustBeComputedEveryHowManyGenerations = -1;
        }  
    } else
    {
        recomputeLociOverWhichFitnessMustBeComputedEveryHowManyGenerations = input.GetNextElementInt();
        if (recomputeLociOverWhichFitnessMustBeComputedEveryHowManyGenerations == 0 || recomputeLociOverWhichFitnessMustBeComputedEveryHowManyGenerations < -1)
        {
            std::cout << "In --T1_resetTrackedT1Muts, the number of generations received is " << recomputeLociOverWhichFitnessMustBeComputedEveryHowManyGenerations << ". Please use a number that is strictly positive or -1 (-1 means, to loop through all the sites systematically).\n";
            abort();
        }
    }
        

    std::string s = input.GetNextElementString();
    if (s == "n" || s == "N" || s == "no" || s == "0")
    {
        allowToCorrectRecomputeLociOverWhichFitnessMustBeComputedEveryHowManyGenerations = false;
    } else if (s == "y" || s == "Y" || s == "yes" || s == "1")
    {
        std::cout << "In --T1_resetTrackedT1Muts, the string received for the second argument is " << s << ". This string is valid. However, the current version of SimBit is not able to modify the generations at which the listing of tracked mutations is done. Sorry! The only possible entry for the second argument is therefore 'no' (or 'n' or '0' or a few other equivalents).\n";
        abort();
        allowToCorrectRecomputeLociOverWhichFitnessMustBeComputedEveryHowManyGenerations = true;
    } else
    {
        std::cout << "In --T1_resetTrackedT1Muts, expected yes or no (or a few equivalents) as second argument but instead received " << s << ".\n";
        abort();
    }

    
    //std::cout << "recomputeLociOverWhichFitnessMustBeComputedEveryHowManyGenerations = "<< recomputeLociOverWhichFitnessMustBeComputedEveryHowManyGenerations << "\n";
}*/



template<typename INT>
double SpeciesSpecificParameters::getRecombinationRatePositionOfLocus(INT locus)
{
    assert(locus < Gmap.TotalNbLoci-1);


    if (locus == 0) return 0.0;

    assert(RecombinationRate.size());
    
    if (RecombinationRate.size() == 1)
    {
        return RecombinationRate.front() * locus;
    } else
    {
        assert(RecombinationRate.size() == Gmap.TotalNbLoci-1);
        return RecombinationRate[locus];
    }
}


void SpeciesSpecificParameters::readT8_mapInfo(InputReader& input)
{
    std::string mode = input.GetNextElementString();

    if (mode == "default" || mode == "prob")
    {
        double prob;
        if (mode == "default")
            prob = 0.05;
        else
            prob = input.GetNextElementDouble();

        if (prob <= 0.0)
        {
            std::cout << "For option --T8_mapInfo, mode 'prob' (or 'default', but 'default' should not cause any trouble), received value '" << prob << "'. The value received must be greater than 0.\n";
            abort();
        }

        if (Gmap.T8_nbLoci)
        {
            auto totalProb = (double)getRecombinationRatePositionOfLocus(Gmap.FromT8LocusToLocus(Gmap.T8_nbLoci-2));

            double lastPosition = 0.0;
            for (long T8locus = 0 ; T8locus < Gmap.T8_nbLoci - 1 ; ++T8locus)
            {
                //std::cout << "Gmap.FromT8LocusToLocus(T8locus) = " << Gmap.FromT8LocusToLocus(T8locus) << "\n";
                auto position = getRecombinationRatePositionOfLocus(Gmap.FromT8LocusToLocus(T8locus));
                //std::cout << "position = " << position << "\n";

                if (position > totalProb - prob / 2) break;

                auto sizeOfBlock = position - lastPosition;
                if (sizeOfBlock > prob)
                {
                    T8_map.push_back(T8locus);
                    lastPosition = position;
                }
            }
            if (T8_map.size() == 0 ) T8_map.push_back(Gmap.T8_nbLoci-1); 
            else if (T8_map.back() != Gmap.T8_nbLoci-1) T8_map.push_back(Gmap.T8_nbLoci-1);

            //T8_map = {20000000, 40000000, 60000000, 80000000, 99999999};

            std::cout << "nbBlocks = " << T8_map.size() << "\n";
            std::cout << "prob = " << prob << "\n";
        }

    } else if (mode == "descr")
    {
        uint32_t sum = 0;
        while (input.IsThereMoreToRead())
        {
            auto x = input.GetNextElementInt();
            if (x <= 0)
            {
                std::cout << "For option --T8_mapInfo, mode 'descr', received valule '" << x << "'. Values must be greater than 0 (no blocks of T8 loci can contain less than a singe T8 locus).\n";
                abort();
            }
            sum += x;
            if (sum > Gmap.T8_nbLoci)
            {
                std::cout << "For option --T8_mapInfo, mode 'descr', the sum of values is greater than the total number of T8 loci. There are " << Gmap.T8_nbLoci << " T8 loci. The last value read is  '" << x << "' and it caused the sum of values read so far to be equal " << sum << ".\n";
                abort();
            }
            T8_map.push_back(sum-1);
        }
            
    } else
    {
        std::cout << "For option --T8_mapInfo, received mode '" << mode << "'. Only modes 'default', 'prob' and 'descr' are accepted.\n";
        abort();
    }

    // Test
    if (Gmap.T8_nbLoci)
    {
        if (T8_map.back() != Gmap.T8_nbLoci-1)
        {
            std::cout << "Problem when creating the T8 map (--T8_mapInfo). The last value is " << T8_map.back() << " and it does not equal the total number of T8 loci minus 1 (which is "<< Gmap.T8_nbLoci-1<< ").\n";
            abort();
        }

        for (size_t i = 1 ; i < T8_map.size() ; ++i)
        {
            if (T8_map[i] <= T8_map[i-1])
            {
                std::cout << "Problem when creating the T8 map (--T8_mapInfo). One of the value (the "<< i << "th value which equals "<<T8_map[i]<<") is not strictly greater than the previous (the "<< i-1 << "th value which equals "<<T8_map[i-1]<<"). Note by the way that there are "<< Gmap.T8_nbLoci<< " T8 loci in total.\n";
                abort();
            }
        }

    }
    
    std::cout << "T8_map: ";
    for (size_t i = 0 ; i < T8_map.size() ; ++i) std::cout << T8_map[i] << " ";
    std::cout <<  "\n";
    
}


void SpeciesSpecificParameters::readFitnessMapInfo(InputReader& input)
{
    this->FitnessMapMinimNbLoci = 0; // This is set to a different value (in allParameters.cpp) only if used default
    this->FitnessMapProbOfEvent = 0.0;
    std::string mode;
    if (!this->T1_isMultiplicitySelection && !this->T56_isMultiplicitySelection &&  this->Gmap.T2_nbLoci == 0)
    {
        mode = "prob";
        this->FitnessMapProbOfEvent = DBL_MAX;
        while (input.IsThereMoreToRead())
        {
            input.skipElement();
        }
        this->FitnessMapT5WeightProbOfEvent = 1.0; // whatever
    } else
    {
        mode = input.GetNextElementString();
        if (mode == "prob")
        {
            this->FitnessMapProbOfEvent = input.GetNextElementDouble();
            if (this->FitnessMapProbOfEvent <= 0.0)
            {
                std::cout << "For species " << this->speciesName << " when reading --FitnessMapInfo ( mode = " << mode << ") received not strictly positive value for 'FitnessMapProbOfEvent' (received "<<this->FitnessMapProbOfEvent<<").\n";
                abort();        
            }

            this->FitnessMapT5WeightProbOfEvent = input.GetNextElementDouble();
            if (this->FitnessMapT5WeightProbOfEvent <= 0.0)
            {
                std::cout << "For species " << this->speciesName << " when reading --FitnessMapInfo ( mode = " << mode << ") received not strictly positive value for 'FitnessMapT5WeightProbOfEvent'(received "<<this->FitnessMapT5WeightProbOfEvent<<").\n";
                abort();        
            }

            this->FitnessMapCoefficient = -9.0;
        } else if (mode == "coef")
        {
            std::cout << "Sorry for --FitnessMapInfo, mode coef is not available on this version.\n";
            abort();
            this->FitnessMapProbOfEvent = -9.0;
            this->FitnessMapCoefficient = input.GetNextElementDouble();
            if (this->FitnessMapCoefficient <= 0.0)
            {
                std::cout << "For species " << this->speciesName << " when reading --FitnessMapInfo ( mode = " << mode << ") received a not strictly positive value (received "<<this->FitnessMapCoefficient<<").\n";
                abort();        
            }

            this->FitnessMapT5WeightProbOfEvent = input.GetNextElementDouble();
            if (this->FitnessMapT5WeightProbOfEvent <= 0.0)
            {
                std::cout << "For species " << this->speciesName << " when reading --FitnessMapInfo ( mode = " << mode << ") received not strictly positive value for 'FitnessMapT5WeightProbOfEvent'(received "<<this->FitnessMapT5WeightProbOfEvent<<").\n";
                abort();        
            }
        } else if (mode == "descr")
        {
            assert(this->FitnessMapInfo_wholeDescription.size() == 0);
            int sum = 0;
            while (input.IsThereMoreToRead())
            {
                auto x = input.GetNextElementInt();
                if (x < 1)
                {
                    std::cout << "For species " << this->speciesName << " when reading --FitnessMapInfo ( mode = " << mode << ") received a number of loci that is lower than 1 (received "<<x<<").\n";
                    abort();
                }
                sum += x;
                this->FitnessMapInfo_wholeDescription.push_back(x);
            }
            if (sum != this->Gmap.TotalNbLoci)
            {
                std::cout << "For species " << this->speciesName << " when reading --FitnessMapInfo ( mode = " << mode << ") received a total of "<<sum<<" loci but there is a total of "<< this->Gmap.TotalNbLoci << " loci asked for option --L. Please note that all loci must be allocated to a fitness map index whether or not the fitness of that locus will be saved for use in next generation. In other words, even for a T3 locus who cannot do the assumption of multiplciity, you have to tell which index map this locus will fall on. It is a little confusing maybe but that's how SimBit works and I've got my reasons :)\n";
                abort();
            }
        }
        else
        {
            std::cout << "For species " << this->speciesName << " when reading --FitnessMapInfo, received mode = " << mode << ". Sorry only modes 'prob', 'descr' and 'coef' are allowed\n";
            std::cout << "Hum.... actually even mode 'coef' is not allowed on this version due to potential misfunctioning and slow down. Sorry! Please use 'prob' or 'descr'.\n";
            abort(); 
        }
    }
    //std::cout << "\n\n\n\n\n\nin SSP : FitnessMapT5WeightProbOfEvent = " << FitnessMapT5WeightProbOfEvent << "\n\n\n\n\n\n\n";
}

void SpeciesSpecificParameters::readInitialpatchSize(InputReader& input)
{
    int TotalNbIndsAtStart = 0;
    assert(this->patchSize.size() == 0);

    std::string Mode = input.GetNextElementString();

    if (Mode == "A")
    {
        for (int patch_index = 0 ; patch_index < GP->__PatchNumber[0] ; patch_index++)
        {
            int x = input.GetNextElementInt();
            int ps;

            if (x < 0)
            {
                ps = this->__patchCapacity[0][patch_index];
            } else
            {
                ps = x;
            }

            if (ps > this->__patchCapacity[0][patch_index])
            {
                std::cout << "For patch " << patch_index << " (zero-based counting), the capacity is " << this->__patchCapacity[0][patch_index] << " but the initial patch size is " << ps << "\n";
                abort();
            }
            if (this->fecundityForFitnessOfOne == -1 && ps != this->__patchCapacity[0][patch_index])
            {
                std::cout << "For patch " << patch_index << " (zero-based counting), the capacity is " << this->__patchCapacity[0][patch_index] << " and the initial patch size is " << ps << ". When the fecundityForFitnessOfOne does not differ from -1 (-1 is the default and means infinite fecundity so that patch sizes are always at carrying capacity), the initial patch size must be at carrying capacity. This error message should only appear if you set both --InitialpatchSize and --fec yourself and made a mismatch between the initial patch size and the carrying capacity at the zeroth generation (indicated by option --N)\n";
                abort();
            }
            this->patchSize.push_back(ps);
            TotalNbIndsAtStart+=ps;
        }
    } else if (Mode == "unif")
    {
        int x = input.GetNextElementInt();
        
        for (int patch_index = 0 ; patch_index < GP->__PatchNumber[0] ; patch_index++)
        {
            int ps; // stands for patchSize
            if (x < 0)
            {
                ps = this->__patchCapacity[0][patch_index];
            } else
            {
                ps = x;
            }

            if (ps > this->__patchCapacity[0][patch_index])
            {
                std::cout << "For patch " << patch_index << " (zero-based counting), the capacity is " << this->__patchCapacity[0][patch_index] << " but the initial patch size is " << ps << "\n";
            }
            if (this->fecundityForFitnessOfOne == -1 && ps != this->__patchCapacity[0][patch_index])
            {
                std::cout << "For patch " << patch_index << " (zero-based counting), the capacity is " << this->__patchCapacity[0][patch_index] << " and the initial patch size is " << ps << ". When the fecundityForFitnessOfOne does not differs from -1 (-1 is the default and means infinite fecundity so that the patch size is always at carrying capacity), the initial patch size must be at carrying capacity. This error message should only appear if you set both --InitialpatchSize and --fec yourself and made a mismatch between the initial patch size and the carrying capacity at the zeroth generation (indicated by option --N)\n";
                abort();
            }
            this->patchSize.push_back(ps);
            TotalNbIndsAtStart+=ps;
        }

    } else
    {
        std::cout << "For option --InitialpatchSize, mode received is " << Mode << ". Sorry only modes 'A' and 'unif' are recognized.\n";
        abort();
    }

    assert(this->patchSize.size() == GP->__PatchNumber[0]);

    if (TotalNbIndsAtStart == 0)
    {
        std::cout << "You asked for a total of exactly 0 individual (over all patches) at initialization (via option '--InitialpatchSize')! Nothing can be simulated for this species!\n" << std::endl;
        abort();
    }
    assert(TotalNbIndsAtStart > 0);

    // Just an extra security
    TotalpatchSize = 0;
    for (uint32_t patch_index = 0 ; patch_index < GP->PatchNumber ; patch_index++)
    {
        TotalpatchSize += patchSize[patch_index];
    }
    assert(TotalpatchSize <= TotalpatchCapacity);        
}


Haplotype SpeciesSpecificParameters::getHaplotypeForIndividualType(InputReader& input, bool haploIndex, std::string& IndividualTypeName)
{
    // IsThereSelection() is called in resetGenetics that comes first.
    std::string beginKeyword;
    std::string endKeyword;
    if (!haploIndex) // if haplo0
    {
        beginKeyword = "haplo0";
        endKeyword = "haplo1";
    } else
    {
        beginKeyword = "haplo1";
        endKeyword = "ind";
    }

    if (!input.IsThereMoreToRead())
    {
        std::cout << "In --indIni (--IndividualInitialization), for individualType '" << IndividualTypeName << "'. expected keyword "<<beginKeyword<< " but the input has been completely read already.\n";
        abort();
    }
    
    if (input.PeakNextElementString() != beginKeyword && input.PeakNextElementString() != "bothHaplo")
    {
        std::cout << "In --indIni (--IndividualInitialization), for individualType '" << IndividualTypeName << "'. expected keyword '"<<beginKeyword<< "'' (or 'bothHaplo') but got " << input.PeakNextElementString() << " instead.\n";
        abort();
    }
    input.skipElement();


    if (input.PeakNextElementString() == "empty")
    {
        input.skipElement();
        std::vector<unsigned char> T1_info(Gmap.T1_nbChars,0);
        std::vector<unsigned char> T2_info(Gmap.T2_nbLoci, 0);
        std::vector<T3type> T3_info(Gmap.T3_nbLoci, 0);
        uint32_t T4ID = std::numeric_limits<uint32_t>::max();
        std::vector<uint32_t> T56_info;
        return Haplotype(T1_info, T2_info, T3_info, T4ID, T56_info);
    }


    // Gather info
    bool receivedT1_info = false;
    bool receivedT2_info = false;
    bool receivedT3_info = false;
    bool receivedT4_info = false;
    bool receivedT56_info = false;
    std::vector<unsigned char> T1_info;
    std::vector<unsigned char> T2_info;
    std::vector<T3type> T3_info;
    uint32_t T4ID = std::numeric_limits<uint32_t>::max();
    std::vector<uint32_t> T56_info;

    while (input.IsThereMoreToRead() && input.PeakNextElementString().substr(0,5) != "haplo" && input.PeakNextElementString() != "ind")
    {
        auto TT = input.GetNextElementString();
        if (TT == "T1" || TT == "t1")
        {
            // Security
            if (receivedT1_info)
            {
                std::cout << "In --indTypes, for individualType '" << IndividualTypeName << "'. Received info about "<<TT<<" twice!\n";
                abort(); 
            }
            receivedT1_info = true;

            // read info
            T1_info.resize(this->Gmap.T1_nbChars,0); // put all zeros for every bit
            for (unsigned locus = 0 ; locus < this->Gmap.T1_nbLoci ; ++locus)
            {
                auto char_index = locus / 8;
                auto bit_index = locus % 8;
                assert(char_index * 8 + bit_index == locus);
                assert(char_index < T1_info.size());
                bool value = input.GetNextElementBool();

                // Set value (code just copy pasted from haplotype.cpp; not a great design to copy paste!)
                //std::cout << "set locus " << locus << " to " << value << "\n";
                //std::cout << "char_index " << char_index << " bit_index " << bit_index << "\n";

                if (value)
                {
                    T1_info[char_index] |= 1 << bit_index;
                } else
                {
                    T1_info[char_index] &= ~(1 << bit_index);
                }
            }
        } else if (TT == "T2" || TT == "t2")
        {
            // Security
            if (receivedT2_info)
            {
                std::cout << "In --indTypes, for individualType '" << IndividualTypeName << "'. Received info about "<<TT<<" twice!\n";
                abort(); 
            }
            receivedT2_info = true;

            // read info
            T2_info.reserve(this->Gmap.T2_nbLoci);
            for (unsigned T2locus = 0 ; T2locus < this->Gmap.T2_nbLoci ; ++T2locus)
            {
                T2_info.push_back(input.GetNextElementInt()); // implicit cast
            }
        } else if (TT == "T3" || TT == "t3")
        {
            // Security
            if (receivedT3_info)
            {
                std::cout << "In --indTypes, for individualType '" << IndividualTypeName << "'. Received info about "<<TT<<" twice!\n";
                abort(); 
            }
            receivedT3_info = true;

            // read info
            T3_info.reserve(this->Gmap.T3_nbLoci);
            for (unsigned T3locus = 0 ; T3locus < this->Gmap.T3_nbLoci ; ++T3locus)
            {
                T3_info.push_back(input.GetNextElementInt()); 
            }

        } else if (TT == "T4" || TT == "t4")
        {
            // Security
            if (receivedT4_info)
            {
                std::cout << "In --indTypes, for individualType '" << IndividualTypeName << "'. Received info about "<<TT<<" twice!\n";
                abort(); 
            }
            receivedT4_info = true;

            // read info
            T4ID = input.GetNextElementInt();
        } else if (TT == "T5" || TT == "t5" || TT == "T6" || TT == "t6"  || TT == "T56" || TT == "t56")
        {
            // Security
            if (receivedT56_info)
            {
                std::cout << "In --indTypes, for individualType '" << IndividualTypeName << "'. Received info about "<<TT<<" twice!\n";
                abort(); 
            }
            receivedT56_info = true;

            // read info
            // Watch out T56 get index of mutations
            while (input.IsThereMoreToRead() && input.PeakNextElementString().substr(0,5) != "patch" && input.PeakNextElementString() != endKeyword && input.PeakNextElementString().at(0) != 'T' && input.PeakNextElementString().at(0) != 't')
            {
                auto newMutPosition = input.GetNextElementInt();
                if (T56_info.size() != 0 && T56_info.back() >= newMutPosition)
                {
                    std::cout << "In --indTypes, for individualType '" << IndividualTypeName << "'. When reading info about T5 loci, got the mutation position " << newMutPosition << " after mutPosition " << T56_info.back() << ". Please make sure that the info value are strictly increased. Note that for T5 loci, you are expected to indicate mutation positions (unlike for T1 loci that expects T1_nbLoci boolean values).\n";
                    abort();
                }
                T56_info.push_back(newMutPosition);
            }
        } else if (TT == "T7" || TT == "t7") 
        {
            std::cout << "In --indTypes, received info about type of locus T7. Sorry, --indIni is currently not able to initialize the T7 loci.\n";
            abort();
        } else
        {
            std::cout << "In --indTypes, received unknown type of locus (" << TT << "). Please if you want to indicate locus of type 1, write either T1 or t1. Same logic apply to other type of loci.\nMaybe you tried to write 'haplo1' but wrote 'Haplo1' or something like that!\nMaybe, did inputted the correct number of values for a given type of loci (give the number of loci for this type at option --Loci (--L))\n";
            abort();
        }
    }

    // Security
    if (Gmap.T1_nbLoci && !receivedT1_info)
    {
        std::cout << "In --indTypes, for individualType "<< IndividualTypeName<< " expected info for T1 loci but did not receive it!\n";
        abort();
    }
    if (Gmap.T2_nbLoci && !receivedT2_info)
    {
        std::cout << "In --indTypes, for individualType "<< IndividualTypeName<< " expected info for T2 loci but did not receive it!\n";
        abort();
    }
    if (Gmap.T3_nbLoci && !receivedT3_info)
    {
        std::cout << "In --indTypes, for individualType "<< IndividualTypeName<< " expected info for T3 loci but did not receive it!\n";
        abort();
    }
    if (Gmap.T4_nbLoci && !receivedT4_info)
    {
        std::cout << "In --indTypes, for individualType "<< IndividualTypeName<< " expected info for T4 loci but did not receive it!\n";
        std::cout << "Note btw that you cannot specify what mutations the haplotype is carrying for the T4 loci as the mutations are placed on the tree. You can only tell what ID the haplotype has (that is its position in the tree).\n";
        abort();
    }
    if (Gmap.T56_nbLoci && !receivedT56_info)
    {
        std::cout << "In --indTypes, for individualType "<< IndividualTypeName<< " expected info for T56 (aka T5; aka T6) loci but did not receive it!\n";
        abort();
    }
    assert(T1_info.size() == Gmap.T1_nbChars);
    assert(T2_info.size() == Gmap.T2_nbLoci);

    return Haplotype(T1_info, T2_info, T3_info, T4ID, T56_info);
}


void SpeciesSpecificParameters::readRedefIndTypes(InputReader& input)
{
    if (input.PeakNextElementString() == "default")
    {
        input.skipElement();
        return;
    }


    while( input.IsThereMoreToRead() )
    {
        // Expecting 'redef'
        {
            auto x = input.GetNextElementString();
            if (x != "redef")
            {
                std::cout << "For option --redefIndTypes, expected keyword 'redef' but received "<<x<<" instead.\n";
                abort();
            }
        }

        // Name of individual type to modify
        {
            auto x = input.GetNextElementString();
            bool foundIt = false;
            for (auto it = individualTypes.begin() ; it != individualTypes.end() ;++it)
            {
                if (it->first == x)
                {
                    foundIt = true;
                    break;
                }
            }
            if (!foundIt)
            {
                std::cout << "For option --redefIndTypes, received a individual type named "<<x<<" but SimBit could not find any previously defined individual types with the same name given in option --indTypes.\n";
                abort();
            }
            
            redefIndTypes_types.push_back(x);
        }


        // Get haplo to modify
        std::string haploToModify = input.GetNextElementString();
        if (haploToModify != "haplo0" && haploToModify != "haplo1"  && haploToModify != "bothHaplo")
        {
            std::cout << "For option --redefIndTypes, expected keyword 'haplo0', 'haplo1' or 'bothHaplo' (indicating what haplotype to redefine) but received "<<haploToModify<<" instead.\n";
                abort();
        }


        // Expecting 'at'
        {
            auto x = input.GetNextElementString();
            if (x != "at")
            {
                std::cout << "For option --redefIndTypes, expected keyword 'at' but received "<<x<<" instead.\n";
                abort();
            }
        }            

        // Generation
        {
            auto x = input.GetNextElementInt();
            if (x <= 0 || x > GP->nbGenerations)
            {
                std::cout << "For option --redefIndTypes, received a time value that is either not strictly positive or greater than the total number of generations. Received " << x << ".\n";
                abort();
            }
            redefIndTypes_times.push_back(x);
        }


        // Expecting 'basedOn'
        {
            auto x = input.GetNextElementString();
            if (x != "basedOn")
            {
                std::cout << "For option --redefIndTypes, expected keyword 'basedOn' but received "<<x<<" instead.\n";
                abort();
            }
        }


        // definition of what haplotpes to copy
        int patch;
        int ind;
        std::string haplo; 
        {
            int generation_index = std::upper_bound(GP->__GenerationChange.begin(), GP->__GenerationChange.end(), redefIndTypes_times.back()) - GP->__GenerationChange.begin() - 1;

            patch = input.GetNextElementInt();
            auto currentPatchNumber = GP->__PatchNumber[generation_index];
            if (patch < 0 || patch >= currentPatchNumber)
            {
                std::cout << "For option --redefIndTypes, got patch index " << patch << " at generation " << redefIndTypes_times.back() << ". It appears that at this generation, there are only " << currentPatchNumber << " patches in the world. As a reminder, the patch index is zero-based counting and can therefore not equal the number of patches.\n";
                abort();
            }

            ind   = input.GetNextElementInt();
            auto capacity = SSP->__patchCapacity[generation_index][patch];
            if (patch < 0 || patch >= capacity)
            {
                std::cout << "For option --redefIndTypes, got individual index " << ind << " in patch index " << patch << " at generation " << redefIndTypes_times.back() << ". It appears that at this generation, this patch has a carrying capacity of " << capacity << " only.\n";
                abort();
            }

            haplo = input.GetNextElementString();
            if (haplo != "haplo0" && haplo != "haplo1" && haplo != "bothHaplo")
            {
                std::cout << "For option --redefIndTypes, expected keywords 'haplo0', 'haplo1' or 'bothHaplo' for the haplotypes the redefinition is based on (after keyword 'basedOn').\n";
                abort();
            }

            if (haplo == "bothHaplo" && haploToModify != "bothHaplo")
            {
                std::cout << "For option --redefIndTypes, redefining indType " << redefIndTypes_types.back() << "  at generation " << redefIndTypes_times.back() << ". You indicated that the redefinition must based based on both haplotypes of individual " << ind << " in patch " << patch << " but to set only one haplotype (you used keyword '" << haploToModify << "') in the individual type. If you base your redefinition on two haplotypes, you need to modify two haplotypes (use 'bothHaplo' for both the indType to redefine and the haplotype you base redefinition on.\n";
                abort();
            }
        }


        // Set how to redefine things
        if (haplo == "bothHaplo")
        {
            assert(haploToModify == "bothHaplo");
            redefIndTypes_whereFromInfo.push_back({{patch, ind, 0},{patch, ind, 1}});
        } else
        {
            std::array<int, 3> haploArray;
            if (haplo == "haplo0")
            {
                haploArray = {patch, ind, 0};
            } else
            {
                assert(haplo == "haplo1");
                haploArray = {patch, ind, 1};
            }

            if (haploToModify == "haplo0")
            {
                redefIndTypes_whereFromInfo.push_back({haploArray,{-1,-1,-1}});
            } else if (haploToModify == "haplo1")
            {
                redefIndTypes_whereFromInfo.push_back({{-1,-1,-1},haploArray});
            } else
            {
                assert(haploToModify == "bothHaplo");
                redefIndTypes_whereFromInfo.push_back({haploArray,haploArray});
            }
        }
    } // Finished reading input


    // Reorder stuff to go in decreasing generation order
    assert(redefIndTypes_times.size() == redefIndTypes_types.size());
    assert(redefIndTypes_times.size() == redefIndTypes_whereFromInfo.size());
    
    if (redefIndTypes_times.size() >= 2)
    {
        auto order = reverse_sort_indexes(redefIndTypes_times);
        reorderNoAssertions(redefIndTypes_times, order);
        reorderNoAssertions(redefIndTypes_types, order);
        reorderNoAssertions(redefIndTypes_whereFromInfo, order);        

        assert(redefIndTypes_times[0] >= redefIndTypes_times[1]);
    }    
}


void SpeciesSpecificParameters::redefineIndividualTypesIfNeeded(Pop& pop)
{
    while (redefIndTypes_times.size() && redefIndTypes_times.back() == GP->CurrentGeneration)
    {
        auto indTypeName = redefIndTypes_types.back();
        auto info = redefIndTypes_whereFromInfo.back();
        redefIndTypes_times.pop_back();
        redefIndTypes_types.pop_back();
        redefIndTypes_whereFromInfo.pop_back();

        auto indTypeIt = SSP->individualTypes.find(indTypeName);
        assert(indTypeIt != SSP->individualTypes.end());
        Individual& indType = indTypeIt->second;

        // Zeroth (first) haplotype
        if (info.first[0] == -1)
        {
            assert(info.first[1] == -1);
            assert(info.first[2] == -1);
        } else
        {
            auto patchi = info.first[0];
            auto indi = info.first[1];
            auto haploi = info.first[2];
            assert(patchi >= 0 && patchi < GP->PatchNumber);
            assert(indi >= 0 && indi < patchCapacity[patchi]);
            assert(haploi == 0 || haploi == 1);

            if (indi < patchSize[patchi])
            {
                indType.getHaplo(0) = pop.getPatch(patchi).getInd(indi).getHaplo(haploi);
            }
        }

        // first (second) haplotype
        if (info.second[0] == -1)
        {
            assert(info.second[1] == -1);
            assert(info.second[2] == -1);
        } else
        {
            auto patchi = info.second[0];
            auto indi = info.second[1];
            auto haploi = info.second[2];
            assert(patchi >= 0 && patchi < GP->PatchNumber);
            assert(indi >= 0 && indi < patchCapacity[patchi]);
            assert(haploi == 0 || haploi == 1);

            if (indi < patchSize[patchi])
            {
                indType.getHaplo(1) = pop.getPatch(patchi).getInd(indi).getHaplo(haploi);
            }
        }

    }

    assert(GP->CurrentGeneration >= 0);
    if (redefIndTypes_times.size()) assert(redefIndTypes_times.back() > GP->CurrentGeneration);
}



void SpeciesSpecificParameters::readIndividualTypes(InputReader& input)
{
    IsThereSelection();              // Need that to initialize individual types
    setFromLocusToFitnessMapIndex(); // Need that to initialize individual types

    // Default 
    if (input.PeakNextElementString() == "default")
    {
        input.skipElement();
        return;
    }
    

    // Read and build individual types
    while (input.IsThereMoreToRead())
    {
        if (input.PeakNextElementString() != "ind")
        {
            std::cout << "In --indIni, expected keyword 'ind' but instead received '" << input.PeakNextElementString() << "'\n";
            abort();
        }
        input.skipElement(); // skipping "ind"

        // Get name and make sure it is a new name
        auto IndividualTypeName = input.GetNextElementString();
        if (individualTypes.find(IndividualTypeName) != individualTypes.end())
        {
            std::cout << "In --indIni, you start a new individualType (with keyword 'ind') and you name it " << IndividualTypeName << ". This name has already been used by a previous IndiivdualType. Please use unique names! Just call them with number if you are lazy to give them names!\n";
            abort();
        }

        // Get Haplotypes
        Haplotype haplo0;
        Haplotype haplo1;
        if (input.PeakNextElementString() == "bothHaplo")
        {
            haplo0 = getHaplotypeForIndividualType(input, 1, IndividualTypeName); // second argument is 1 so that it stops at keyword 'ind'
            haplo1 = haplo0;
        } else
        {
            haplo0 = getHaplotypeForIndividualType(input, 0, IndividualTypeName);
            haplo1 = getHaplotypeForIndividualType(input, 1, IndividualTypeName);
        }
            
        // Add IndividualType
        individualTypes[IndividualTypeName] = Individual(haplo0, haplo1);
    }
}






void SpeciesSpecificParameters::readIndividualInitialization(InputReader& input)
{
    // Set isIndividualInitialization to inform other initialization option that this one has been used
    if (input.PeakNextElementString() == "default")
    {
        input.skipElement();
        this->isIndividualInitialization = false;
        return;
    } else
    {
        this->isIndividualInitialization = true;
    }

    // security
    if (individualTypes.size() == 0)
    {
        std::cout << "You explicitely used the option '--indIni (--IndividualInitialization)', but no individual types have been defined (with option --indTypes).\n";
        abort();
    }
        


    ////////////////////////////
    // Assign individual type //
    ////////////////////////////
    
    individualTypesForInitialization.resize(GP->PatchNumber);
    for (unsigned patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
    {
        std::string patchKeyword = "patch" + std::to_string(patch_index);
        if (input.PeakNextElementString() != patchKeyword)
        {
            std::cout << "In --indIni (--IndividualInitialization), expected keyword "<<patchKeyword<< " but got " << input.PeakNextElementString() << " instead\n";
            abort();
        }
        input.skipElement(); // skipping patch<number>

        unsigned ind_index = 0;
        while (ind_index < this->patchSize[patch_index])
        {
            auto IndividualTypeName = input.GetNextElementString();
            if (individualTypes.find(IndividualTypeName) == individualTypes.end())
            {
                std::cout << "In --indIni (--IndividualInitialization), asking for individual type " << IndividualTypeName << " in patch " << patch_index << " but individual type " << IndividualTypeName << " has not been defined previously. It is probably caused by a typo when naming individual types.\n";
                abort();
            }
            unsigned howMany = input.GetNextElementInt();

            // security
            ind_index += howMany;
            if (ind_index > this->patchSize[patch_index])
            {
                std::cout << "In --indIni (--IndividualInitialization), for patch " << patch_index << " received more individuals than was set as initial patch size (initial patch size for this patch is " << this->patchSize[patch_index] << "). Note that if you have not specified the initial patch size yourself, then it is set to carrying capacity of this patch at generation 0.\n";
                abort();
            }

            for (unsigned howMany_index = 0 ; howMany_index < howMany ; ++howMany_index )
            {
                individualTypesForInitialization[patch_index].push_back(IndividualTypeName);
            }
        }
        assert(ind_index == this->patchSize[patch_index]);
        assert(individualTypesForInitialization[patch_index].size() == this->patchSize[patch_index]);
    }
}


void SpeciesSpecificParameters::readT1_Initial_AlleleFreqs(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option --T1_Initial_AlleleFreqs, the std::string that is read is: " << input.print() << std::endl;
#endif
    if (input.PeakNextElementString() == "NA")
    {
        if (this->Gmap.T1_nbLoci != 0)
        {
            std::cout << "For option --T1_Initial_AlleleFreqs, received 'NA' however there are T1 loci as indiciated by --L (--Loci) (this->Gmap.T1_nbLoci = "<<this->Gmap.T1_nbLoci<<")" << "\n";
            abort();
        }
        input.skipElement();
    } else
    {
        std::string Mode = input.GetNextElementString();

        if (Mode != "default" && isIndividualInitialization)
        {
            std::cout << "You cannot use both --individualInitialization (--indIni) and --T1_Initial_AlleleFreqs\n";
            abort();
        }


        this->T1_Initial_AlleleFreqs_AllZeros = false; // Set to true laster if received "AllZeros"
        this->T1_Initial_AlleleFreqs_AllOnes = false;  // Set to true laster if received "AllOnes"
        

        if (Mode.compare("freqs") == 0)
        {
            for (int patch_index = 0 ; patch_index < GP->__PatchNumber[0] ; ++patch_index)
            {
                int nbRepeatsLeft = 0;
                double freq;
                std::vector<double> OneLine;
                for (int T1Locus = 0 ; T1Locus < this->Gmap.T1_nbLoci ; T1Locus++)
                {
                    if (nbRepeatsLeft == 0 && input.PeakNextElementString() == "R")
                    {
                        input.skipElement();
                        freq = input.GetNextElementDouble();
                        nbRepeatsLeft = input.GetNextElementInt();
                        if (nbRepeatsLeft < 0)
                        {
                            std::cout << "In '--T1_Initial_AlleleFreqs', received a number of repeats lower than zero. Number of repeats received is " << nbRepeatsLeft << "\n";
                            abort();
                        }
                    }
                    if (nbRepeatsLeft == 0)
                    {
                        freq = input.GetNextElementDouble();
                    } else
                    {
                        nbRepeatsLeft--;
                    }
                    
                    if (freq < 0.0 || freq > 1.0)
                    {
                        std::cout << "In '--T1_Initial_AlleleFreqs', received an impossible allele frequency. Received '" << freq << "'\n";
                        abort();
                    }

                    //std::cout << "T1Locus = " << T1Locus << " x = " << x << " nbRepeatsLeft = "<<nbRepeatsLeft<<"\n";
                    OneLine.push_back(freq);
                }
                assert(nbRepeatsLeft >= 0);
                if (nbRepeatsLeft != 0)
                {
                    std::cout << "In '--T1_Initial_AlleleFreqs', too many values received!'\n";
                    abort();
                }
                assert(OneLine.size() == this->Gmap.T1_nbLoci);
                this->T1_Initial_AlleleFreqs.push_back(OneLine);
            }
        } else if (Mode.compare("Shift")==0)
        {
            int shiftPosition = input.GetNextElementInt();
            if (shiftPosition < 0 || shiftPosition > GP->__PatchNumber[0])
            {
                std::cout << "In '--T1_Initial_AlleleFreqs', received an impossible shiftPosition. Received '" << shiftPosition << ". For information, the simulation is starting with " << GP->__PatchNumber[0] << " as indicated in option '--PN (--PatchNumber)'\n";
                abort();
            }
            
            for (int patch_index = 0 ; patch_index < GP->__PatchNumber[0] ; ++patch_index)
            {
                std::vector<double> line;
                for (int T1Locus = 0 ; T1Locus < this->Gmap.T1_nbLoci ; T1Locus++)
                {
                    if (patch_index < shiftPosition)
                    {
                        line.push_back(0.0);
                    } else
                    {
                        line.push_back(1.0);
                    }
                }
                assert(line.size() == this->Gmap.T1_nbLoci);
                this->T1_Initial_AlleleFreqs.push_back(line);
            }
        } else if (Mode.compare("AllZeros") == 0 || Mode.compare("default") == 0)
        {
            this->T1_Initial_AlleleFreqs_AllZeros = true;
        } else if (Mode.compare("AllOnes") == 0)
        {
            this->T1_Initial_AlleleFreqs_AllOnes = true;
        } else
        {
            std::cout << "Sorry, for '--T1_Initial_AlleleFreqs', the Mode " << Mode << " has not been implemented yet. Only Modes 'freqs' (note that 'freqs' was previously called 'A'), 'AllZeros', 'AllOnes' and 'Shift' are accepted for the moment.";
            abort();
        }
    }

    // Security check
    if (!this->T1_Initial_AlleleFreqs_AllZeros && !this->T1_Initial_AlleleFreqs_AllOnes)
    {
        assert(this->T1_Initial_AlleleFreqs.size() == GP->__PatchNumber[0]);
        for (int patch_index = 0 ; patch_index < GP->__PatchNumber[0] ; ++patch_index)
        {
            assert(this->T1_Initial_AlleleFreqs[patch_index].size() == this->Gmap.T1_nbLoci);
        }
    }




    if (!T1_Initial_AlleleFreqs_AllZeros && !T1_Initial_AlleleFreqs_AllOnes)
    {
        funkyMathForQuasiRandomT1AllFreqInitialization.resize(GP->PatchNumber);
        for (int patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
        {
            funkyMathForQuasiRandomT1AllFreqInitialization[patch_index].reserve(Gmap.T1_nbLoci);
            for (int locus = 0 ; locus < Gmap.T1_nbLoci ; ++locus)
            {
                funkyMathForQuasiRandomT1AllFreqInitialization[patch_index].push_back(
                    (int)(((long)locus*103) % (long)(2 * SSP->TotalpatchCapacity))
                );    
            }    
        }
    }
}

void SpeciesSpecificParameters::readT56_toggleMutsEveryNGeneration(InputReader& input)
{
    T56_toggleMutsEveryNGeneration = input.GetNextElementInt();
    if (T56_toggleMutsEveryNGeneration <= 0)
    {
        T56_toggleMutsEveryNGeneration_nextGeneration = -1;
    } else
    {
        T56_toggleMutsEveryNGeneration_nextGeneration = T56_toggleMutsEveryNGeneration;
    }
}


void SpeciesSpecificParameters::readT56_Initial_AlleleFreqs(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option --T5_Initial_AlleleFreqs, the std::string that is read is: " << input.print() << std::endl;
#endif

    if (input.PeakNextElementString() == "default")
    {
        input.skipElement();
        this->T56_Initial_AlleleFreqs_AllZeros = true;
        this->T56_Initial_AlleleFreqs_AllOnes = false;
    } else if (input.PeakNextElementString() == "NA")
    {
        if (this->Gmap.T56_nbLoci != 0)
        {
            std::cout << "For option --T5_Initial_AlleleFreqs, received 'NA' however there are T5 loci as indiciated by --L (--Loci) (this->Gmap.T5_nbLoci = "<<this->Gmap.T56_nbLoci<<")" << "\n";
            abort();
        }
        input.skipElement();
    } else
    {
        std::string Mode;

        Mode = input.GetNextElementString();

        if (Mode != "default" && isIndividualInitialization)
        {
            std::cout << "You cannot use both --individualInitialization (--indIni) and --T56_Initial_AlleleFreqs\n";
            abort();
        }

        this->T56_Initial_AlleleFreqs_AllZeros = false; // Set to true later if received "AllZeros"
        this->T56_Initial_AlleleFreqs_AllOnes = false;  // Set to true later if received "AllOnes"
        
        if (Mode.compare("A") == 0)
        {
            for (int patch_index = 0 ; patch_index < GP->__PatchNumber[0] ; ++patch_index)
            {
                int nbRepeatsLeft = 0;
                double freq;
                std::vector<double> OneLine;
                for (int T5Locus = 0 ; T5Locus < this->Gmap.T56_nbLoci ; T5Locus++)
                {
                    if (nbRepeatsLeft == 0 && input.PeakNextElementString() == "R")
                    {
                        input.skipElement();
                        freq = input.GetNextElementDouble();
                        nbRepeatsLeft = input.GetNextElementInt();
                        if (nbRepeatsLeft < 0)
                        {
                            std::cout << "In '--T56_Initial_AlleleFreqs', received a number of repeats lower than zero. Number of repeats received is " << nbRepeatsLeft << "\n";
                            abort();
                        }
                    }
                    if (nbRepeatsLeft == 0)
                    {
                        freq = input.GetNextElementDouble();
                    } else
                    {
                        nbRepeatsLeft--;
                    }
                    
                    if (freq < 0.0 || freq > 1.0)
                    {
                        std::cout << "In '--T56_Initial_AlleleFreqs', received an impossible allele frequency. Received '" << freq << "'\n";
                        abort();
                    }

                    //std::cout << "T1Locus = " << T1Locus << " x = " << x << " nbRepeatsLeft = "<<nbRepeatsLeft<<"\n";
                    OneLine.push_back(freq);
                }
                assert(nbRepeatsLeft >= 0);
                if (nbRepeatsLeft != 0)
                {
                    std::cout << "In '--T1_Initial_AlleleFreqs', too many values received!'\n";
                    abort();
                }
                assert(OneLine.size() == this->Gmap.T56_nbLoci);
                this->T56_Initial_AlleleFreqs.push_back(OneLine);
            }
        } else if (Mode.compare("Shift")==0)
        {
            int shiftPosition = input.GetNextElementInt();
            if (shiftPosition < 0 || shiftPosition > GP->__PatchNumber[0])
            {
                std::cout << "In '--T56_Initial_AlleleFreqs', received an impossible shiftPosition. Received '" << shiftPosition << ". For information, the simulation is starting with " << GP->__PatchNumber[0] << " as indicated in option '--PN (--PatchNumber)'\n";
                abort();
            }
            
            for (int patch_index = 0 ; patch_index < GP->__PatchNumber[0] ; ++patch_index)
            {
                std::vector<double> line;
                for (int T5Locus = 0 ; T5Locus < this->Gmap.T56_nbLoci ; T5Locus++)
                {
                    if (patch_index < shiftPosition)
                    {
                        line.push_back(0.0);
                    } else
                    {
                        line.push_back(1.0);
                    }
                }
                assert(line.size() == this->Gmap.T56_nbLoci);
                this->T56_Initial_AlleleFreqs.push_back(line);
            }
        } else if (Mode.compare("AllZeros") == 0 || Mode.compare("default") == 0)
        {
            this->T56_Initial_AlleleFreqs_AllZeros = true;
        } else if (Mode.compare("AllOnes") == 0)
        {
            std::cout << "Sorry the option AllOnes is currently not available for T56_Initial_AlleleFreqs.\n";
            abort();
            this->T56_Initial_AlleleFreqs_AllOnes = true;
        } else
        {
            std::cout << "Sorry, for '--T56_Initial_AlleleFreqs', the Mode " << Mode << " has not been implemented yet. Only Modes 'A', 'AllZeros', 'AllOnes' and 'Shift' are accepted for the moment.";
            abort();
        }
    }

    // Security check
    if (!this->T56_Initial_AlleleFreqs_AllZeros && !this->T56_Initial_AlleleFreqs_AllOnes)
    {
        assert(this->T56_Initial_AlleleFreqs.size() == GP->__PatchNumber[0]);
        for (int patch_index = 0 ; patch_index < GP->__PatchNumber[0] ; ++patch_index)
        {
            assert(this->T56_Initial_AlleleFreqs[patch_index].size() == this->Gmap.T56_nbLoci);
        }
    }

    // Initialize flipped if T6
    if (SSP->Gmap.T6ntrl_nbLoci)
    {
        T6ntrl_flipped = CompressedSortedDeque(SSP->Gmap.T6ntrl_nbLoci);
    }

    /*if (T56sel_compress && SSP->Gmap.T6sel_nbLoci)
    {
        T6sel_flipped = CompressedSortedDeque(SSP->Gmap.T6sel_nbLoci);
    }*/

    // If all fixed
    if (T56_Initial_AlleleFreqs_AllOnes)
    {
        if (Gmap.isT56ntrlCompress)
        {
            for (auto locus = 0 ; locus < this->Gmap.T6ntrl_nbLoci ; ++locus)
                T6ntrl_flipped.push_back(locus);
        } else
        {
            for (auto locus = 0 ; locus < this->Gmap.T5ntrl_nbLoci ; ++locus)
                T5ntrl_flipped.push_back(locus);
        }

        /*if (T56sel_compress)
        {
            for (auto locus = 0 ; locus < this->Gmap.T6sel_nbLoci ; ++locus)
                T6sel_flipped.push_back(locus);
        } else
        {
            for (auto locus = 0 ; locus < this->Gmap.T5sel_nbLoci ; ++locus)
                T5sel_flipped.push_back(locus);
        }*/  
    } 

    // The following sounds awefully sily but it is just to ensure that begin and end are not nullptr (so as to not give nullptr to ZipIterator<std::vector<uint32_t>, std::vector<uint32_t>::iterator>). It probably serves no purpose though!
    if (T5ntrl_flipped.size() == 0)
    {
        T5ntrl_flipped.push_back(0);
        T5ntrl_flipped.resize(0);
    }
    /*if (T5sel_flipped.size() == 0)
    {
        T5sel_flipped.push_back(0);
        T5sel_flipped.resize(0);
    }*/


    if (!T56_Initial_AlleleFreqs_AllZeros && !T56_Initial_AlleleFreqs_AllOnes)
    {
        funkyMathForQuasiRandomT56AllFreqInitialization.resize(GP->PatchNumber);
        for (int patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
        {
            funkyMathForQuasiRandomT56AllFreqInitialization[patch_index].reserve(Gmap.T56_nbLoci);
            for (int locus = 0 ; locus < Gmap.T56_nbLoci ; ++locus)
            {
                funkyMathForQuasiRandomT56AllFreqInitialization[patch_index].push_back(
                    (int)(((long)locus*107) % (long)(2 * SSP->TotalpatchCapacity))
                );    
            }    
        }
    }
    
    
}

void SpeciesSpecificParameters::readSelectionOn(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option '--selectionOn', the std::string that is read is: " << input.print() << std::endl;
#endif    

    std::string s = input.GetNextElementString();

    if (s == "fertility")
    {
        selectionOn = 0;
    } else if (s == "viability")
    {
        selectionOn = 1;
        if (fecundityForFitnessOfOne != -1)
        {
            std::cout << "You asked for fecundityForFitnessOfOne different from -1 (" << fecundityForFitnessOfOne << ") and for selection on viability only. Sorry, variation of patch size must be computed from selection on fertility and therefore both can not be simulated simultaneously.\n";
            abort();
        }
    } else if (s == "both" || s == "fertilityAndViability")
    {
        selectionOn = 2;
    } else
    {
        std::cout << "For option '--selectionOn', expected either 'fertility' (which should be the default), 'viability' or 'fertilityAndViability' (aka. 'both'). Instead it received '" << s << "'.\n";
        abort();
    }
}

void SpeciesSpecificParameters::readPopGrowthModel(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option '--popGrowthModel', the std::string that is read is: " << input.print() << std::endl;
#endif


    bool isAnyLogisticModel = false;// To test user input


    this->__growthK.resize(GP->__GenerationChange.size());
    for ( int generation_index = 0; generation_index < GP->__GenerationChange.size() ; generation_index++)
    {
        (void) input.GetNextGenerationMarker(generation_index);
        int currentPatchNumber = GP->__PatchNumber[generation_index];
        
        this->__growthK[generation_index].resize(currentPatchNumber);

        std::string Mode = input.GetNextElementString();

        if (Mode == "A")
        {
            for (int patch_index = 0 ; patch_index < currentPatchNumber ; patch_index++)
            {
                int x;
                std::string s = input.PeakNextElementString();
                if (s == "logistic")
                {
                    input.skipElement();
                    x = -2;
                        
                } else if (s == "exponential")
                {
                    input.skipElement(); 
                    x = -1;
                } else
                {
                    x = input.GetNextElementInt();
                    if (x < -2)
                    {
                        std::cout << "In --popGrowthModel, received a negative K value different from -1 and -2 (-1 is aka 'exponential'; -2 (aka. 'logistic') means logistic with the carrying capacity as the growthK; k value received = "<<x << ")\n";
                        abort();
                    }
                }
                if (x != -1) isAnyLogisticModel = true;
                this->__growthK[generation_index][patch_index] = x;
            }
        } else if (Mode == "unif")
        {
            std::string s = input.PeakNextElementString();
            int x;
            if (s == "logistic")
            {
                input.skipElement();                
                x = -2;
            } else if (s == "exponential")
            {
                input.skipElement(); 
                x = -1;
            } else
            {
                x = input.GetNextElementInt();
                if (x < -2)
                {
                    std::cout << "In --popGrowthModel, received a negative K value different from -1 and -2 (-1 is aka 'exponential'; -2 (aka. 'logistic') means logistic with the carrying capacity as the growthK; k value received = "<<x << ")\n";
                    abort();
                }
            }

            if (x != -1) isAnyLogisticModel = true;
            for (int patch_index = 0 ; patch_index < currentPatchNumber ; patch_index++)
            {
                assert(this->__patchCapacity[generation_index][patch_index] >= 0);
                this->__growthK[generation_index][patch_index] = x;
            }
        } else
        {
            std::cout << "In --popGrowthModel, received " << Mode << " for mode. Only modes 'A' and 'unif' are recognized. Sorry\n";
            abort(); 
        }
        
            
    }
    this->growthK = this->__growthK[0];


    // Test user input
    if (!isAnyLogisticModel)
    {
        // Then make sure user did not specify competitiono matrix
        for (int speciesIndexA = 0 ; speciesIndexA < GP->nbSpecies ; ++speciesIndexA)
        {
            for (int speciesIndexB = 0 ; speciesIndexB < GP->nbSpecies ; ++speciesIndexB)
            {
                if (speciesIndexA == speciesIndexB)
                {
                    assert(GP->speciesCompetition[speciesIndexA][speciesIndexB] == 1.0);
                } else
                {
                    if (GP->speciesCompetition[speciesIndexA][speciesIndexB] != 0.0)
                    {
                        std::cout << "When received --popGrowthModel, some values differ from -1 (-1 is aka 'exponential' in the input). So far, so good! The problem is that at the option --eco you specified a competition matrix that does have some non-zero values on the off-diagonal. Competition (unlike 'interaction') is computed through the logistic model and therefore such non-zeros off the diagonal cannot be simulated if all of the growthK are set to -1 (aka. exponential). I hope that makes sense to you! So either, remove the non-zeros values from the off-diagonal of the competition matrix or explicitely specify a logistic growth model.\n";
                        abort();
                    }
                }
                
            }
        }
    }
}

void SpeciesSpecificParameters::readHabitats(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option '--H (--Habitats)', the std::string that is read is: " << input.print() << std::endl;
#endif

    assert(GP->__GenerationChange.size() > 0);
    for ( int generation_index = 0; generation_index < GP->__GenerationChange.size() ; generation_index++)
    {
        int Generation = input.GetNextGenerationMarker(generation_index);
        //std::cout << "Generation = " << Generation << "\n";
        std::string Mode = input.GetNextElementString();
        //std::cout << "Mode = " << Mode << "\n";

        std::vector<int> ForASingleGeneration;

        if (Mode.compare("unif")==0)
        {
            int UniqueHabitat = input.GetNextElementInt();
            assert(GP->__PatchNumber.size() > generation_index);
            for (int patch_index = 0 ; patch_index < GP->__PatchNumber[generation_index] ; ++patch_index)
            {
                ForASingleGeneration.push_back(UniqueHabitat);
            }
            this->__MaxHabitat.push_back(UniqueHabitat);
        } else if (Mode.compare("A")==0)
        {
            assert(GP->__PatchNumber[generation_index] > 0);
            int max = 0;
            assert(GP->__PatchNumber.size() > generation_index);
            for (int patch_index = 0 ; patch_index < GP->__PatchNumber[generation_index] ; ++patch_index)
            {
                int h = input.GetNextElementInt();
                if (max < h) {max = h;}                
                ForASingleGeneration.push_back(h);
            }
            this->__MaxHabitat.push_back(max);
        } else
        {
            std::cout << "In '--H (--Habitats)', at generation " << Generation << ", Mode received is " << Mode << ". Sorry, only Modes 'unif' and 'A' are implemented for the moment.\n";
            abort();
        }
        assert(ForASingleGeneration.size() == GP->__PatchNumber[generation_index]);
        
        // Assign to __Habitats
        this->__Habitats.push_back(ForASingleGeneration);
    }
    assert(this->__Habitats.size() == GP->__GenerationChange.size());
    assert(this->__Habitats.size() == this->__MaxHabitat.size());
    

    // flatten, sort and remove duplicates from the matrix to 1) get the maximum number and 2) make sure no habitat index is missing.
    std::vector<int> flat;
    for (int generation_index = 0 ; generation_index < GP->__GenerationChange.size() ; generation_index++)
    {
        auto& ForASingleGeneration = this->__Habitats[generation_index];
        assert(ForASingleGeneration.size() == GP->__PatchNumber[generation_index]);
        for (int patch_index = 0 ; patch_index < GP->__PatchNumber[generation_index] ; patch_index++)
        {
            flat.push_back(ForASingleGeneration[patch_index]);
        }
    }
    std::sort(flat.begin(), flat.end());
    flat.erase(std::unique(flat.begin(), flat.end()), flat.end());
    assert(this->__MaxHabitat.size() > 0);
    assert(this->__Habitats.size() > 0);
    this->MaxHabitat = this->__MaxHabitat[0];
    this->Habitats = this->__Habitats[0];

    if (flat[0] != 0)
    {
        std::cout << "In option '--H (--Habitats)', the lowest element received is " << flat[0] << ". The lowest element has to be 0 (zero based counting). Please also note that one cannot input an habitat index of n if all [0,n-1] have have not been included.\n";
        abort();
    }

    // 2) make sure no habitat index is missing.
    int previous = flat[0] - 1;
    for (int current : flat)
    {
        if (previous + 1 != current)
        {
            std::cout << "In option '--H (--Habitats)', you included the index " << current << " while the previous index has not been included. The last index included was " << previous << ".  Please also note that you used. Please note --Habitats use zero based counting.\n";
            abort();
        }
        previous = current;
    }

    // set MaxEverHabitat
    MaxEverHabitat = 0;
    assert(MaxHabitat >= 0);
    assert(__MaxHabitat.size() > 0);
    assert(__MaxHabitat[0] == MaxHabitat);
    for (auto& MaxHabitatInSpecificGenerationRange : __MaxHabitat)
    {
        assert(MaxHabitatInSpecificGenerationRange >= 0);
        if (MaxHabitatInSpecificGenerationRange > MaxEverHabitat)
            MaxEverHabitat = MaxHabitatInSpecificGenerationRange;
    }
    assert(MaxEverHabitat >= 0);
    
}

void SpeciesSpecificParameters::readpatchCapacity(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option '--N (--patchCapacity)', the std::string that is read is: " << input.print() << std::endl;
#endif

    this->TotalpatchCapacity = 0;
    maxEverpatchCapacity.resize(GP->maxEverPatchNumber,-1);

    for (int generation_index = 0 ; generation_index < GP->__GenerationChange.size() ; generation_index++)
    {
        (void) input.GetNextGenerationMarker(generation_index);
        int CurrentPatchNumber = GP->__PatchNumber[generation_index];
        assert(CurrentPatchNumber <= GP->maxEverPatchNumber);

        std::string Mode = input.GetNextElementString();
        std::vector<int> patchCapacity_line; // Will contain a single entry (for one interval of time) that will be pushed to '__patchCapacity'
        if (Mode.compare("A")==0)
        {
            for (int patch_index = 0 ; patch_index < CurrentPatchNumber; patch_index++)
            {
                int PC = input.GetNextElementInt();
                if (PC < 0)
                {
                    std::cout << "For option 'patchCapacity', received a patch capacity that is negative for the patch index "<<patch_index<<". Value received is "<<PC<<"."<< std::endl;
                    abort();
                }
                if (generation_index == 0) this->TotalpatchCapacity += PC;
                patchCapacity_line.push_back(PC);
            }
        } else if (Mode.compare("unif")==0)
        {
            int uniquepatchCapacity = input.GetNextElementInt();
            if (uniquepatchCapacity < 0)
            {
                std::cout << "For option 'patchCapacity', received a unique patch capacity that is negative. Value received is "<<uniquepatchCapacity<<"."<< std::endl;
                abort();
            }
            for (int patch_index=0 ; patch_index < CurrentPatchNumber ; ++patch_index)
            {
                if (generation_index == 0) this->TotalpatchCapacity += uniquepatchCapacity; 
                patchCapacity_line.push_back(uniquepatchCapacity);
            }
            
        } else
        {
            std::cout << "Sorry, for option 'patchCapacity', only Modes 'A' and 'unif' are implemented for the moment. Mode received is " << Mode << std::endl;
            abort();
        }

        for (int patch_index=0 ; patch_index < CurrentPatchNumber ; ++patch_index)
        {
            maxEverpatchCapacity[patch_index] = std::max(maxEverpatchCapacity[patch_index], patchCapacity_line[patch_index]);
        }

        this->__patchCapacity.push_back(patchCapacity_line);
            

        if (generation_index == 0) 
        {
            if (this->TotalpatchCapacity <= 0)
            {
                std::cout << "The total patch capacity over all patches for species " << this->speciesName << " is " << this->TotalpatchCapacity << ". Nothing can be simulated.\n";
                abort();
            }
        }
    }
    
    assert(this->__patchCapacity.size() == GP->__GenerationChange.size());
    this->patchCapacity = this->__patchCapacity[0];
    

    for (int patch_index=0 ; patch_index < GP->maxEverPatchNumber ; ++patch_index)
    {
        if (maxEverpatchCapacity[patch_index] == 0)
        {
            std::cout << "patch_index " << patch_index << " never gets a carrying capacity greater than " << 0 << ". This is a weird unseless patch!\n";
            abort();
        }
    }
}

void SpeciesSpecificParameters::readnbSubGenerations(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option '--nbSubGens (--nbSubGenerations)', the std::string that is read is: " << input.print() << std::endl;
#endif

    this->nbSubGenerationsPerGeneration = input.GetNextElementInt();
    if (this->nbSubGenerationsPerGeneration < 1)
    {
        std::cout << "In option '--nbSubGens (--nbSubGenerations)', for species " << this->speciesName << " received " << this->nbSubGenerationsPerGeneration << ". The number of subGeneration must be equal or greater than 1.\n";
        abort();
    }
    
    
}


void SpeciesSpecificParameters::readT1_EpistaticFitnessEffects(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option '--T1_epistasis (--T1_EpistaticFitnessEffects)', the std::string that is read is: " << input.print() << std::endl;
#endif 

    if (input.PeakNextElementString() == "NA")
    {
        if (this->Gmap.T1_nbLoci != 0)
        {
            std::cout << "In option '--T1_epistasis (--T1_EpistaticFitnessEffects)', for species "<<this->speciesName<<", received 'NA' but there are T1 loci as indiciated in --L (--Loci) (this->Gmap.T1_nbLoci = "<<this->Gmap.T1_nbLoci<<").\n";
            abort();   
        }
        input.skipElement();
    } else
    {
        for (int habitat = 0; habitat <= this->MaxEverHabitat ; habitat++)
        {
            input.GetNextHabitatMarker(habitat);

            // Get info for this habitat
            std::vector<std::vector<T1_locusDescription>> LociIndices_ForASingleHabitat;
            std::vector<std::vector<fitnesstype>>              FitnessEffects_ForASingleHabitat;


            while (input.IsThereMoreToRead() && input.PeakNextElementString().at(0) != '@')
            {
                // Initialize vector to fill
                std::vector<T1_locusDescription> LociIndices_ForASingleGroupOfLoci;
                std::vector<fitnesstype>         FitnessEffects_ForASingleGroupOfLoci;


                // 'loci' keyword
                if (input.PeakNextElementString() != "loci")
                {
                    std::cout << "In option '--T1_epistasis (--T1_EpistaticFitnessEffects)', for species "<<this->speciesName<<", expected string 'loci' but received "<<input.PeakNextElementString()<<" instead.\n";
                }
                input.skipElement(); // skip expectedString (set<number>)

                // Read interacting loci
                while (input.PeakNextElementString() != "fit")
                {
                    int LocusPosition = input.GetNextElementInt();
                    if ( LocusPosition < 0 || LocusPosition >= this->Gmap.T1_nbLoci )
                    {
                        std::cout << "In option '--T1_epistasis (--T1_EpistaticFitnessEffects)', received locus position "<<LocusPosition<<" (indicated for habitat " << habitat << " that is either lower than zero or greater or equal to the number of T1 sites (" << this->Gmap.T1_nbLoci <<"). As a reminder the first locus is the zeroth locus (zero based counting).\n";
                        abort();  
                    }
                    //std::cout << "LocusPosition = " << LocusPosition << "\n";
                    T1_locusDescription T1_locus(LocusPosition / EIGHT, LocusPosition % EIGHT, LocusPosition);
                    LociIndices_ForASingleGroupOfLoci.push_back(T1_locus);
                }
                auto nbLociUnderEpistasis = LociIndices_ForASingleGroupOfLoci.size();
                input.skipElement(); // skip "fit"
                if (nbLociUnderEpistasis <= 1)
                {
                    std::cout << "In option '--T1_epistasis (--T1_EpistaticFitnessEffects)', received "<<nbLociUnderEpistasis<<" loci but it takes at least 2 loci to have an epistatic interaction.\n";
                    abort();  
                }



                // Read fitness effects
                /*
                SSP->T1_Epistasis_LociIndices[Habitat][0].
                locus 1:    00 00 00 00 00 00 00 00 00 01 01 01 01 01 01 01 01 01 11 11 11 11 11 11 11 11 11
                locus 2:    00 00 00 01 01 01 11 11 11 00 00 00 01 01 01 11 11 11 00 00 00 01 01 01 11 11 11
                locus 3:    00 01 11 00 01 11 00 01 11 00 01 11 00 01 11 00 01 11 00 01 11 00 01 11 00 01 11
                Index:      0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 
                */

                int TotalNbFitnessesToIndicate = pow(3,nbLociUnderEpistasis);
                for (int i = 0 ; i < TotalNbFitnessesToIndicate; i++)
                {
                    double fit = input.GetNextElementDouble();
                    if ( fit < 0.0)
                    {
                      std::cout << "In option '--T1_epistasis (--T1_EpistaticFitnessEffects)', the " << i << "th fitness value indicated after habitat " << habitat << "is either lower than zero. The fitness value indicated is  " << fit << ".\n";
                        abort();  
                    }
                    FitnessEffects_ForASingleGroupOfLoci.push_back(fit);               
                }

                // Push to "for a single habitat"
                LociIndices_ForASingleHabitat.push_back(LociIndices_ForASingleGroupOfLoci);
                FitnessEffects_ForASingleHabitat.push_back(FitnessEffects_ForASingleGroupOfLoci); 
            }

            this->T1_Epistasis_LociIndices.push_back(LociIndices_ForASingleHabitat);
            this->T1_Epistasis_FitnessEffects.push_back(FitnessEffects_ForASingleHabitat); 
        } // end of for habitat
        
        assert(this->T1_Epistasis_LociIndices.size() == this->MaxEverHabitat + 1);
        assert(this->T1_Epistasis_FitnessEffects.size() == this->MaxEverHabitat + 1);

        //get fitness value with this->FitnessEffects_ForASingleHabitat[habitat][groupOfLoci][fitnessValueIndex]
    } // end of if "NA"  
}

void SpeciesSpecificParameters::readT56_FitnessEffects(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option '--T5_FitnessEffects', the std::string that is read is: " << input.print() << std::endl;
#endif
    // WATCH OUT: This option is read before --L. Only a quick screen of L has been done to set quickScreenAtL_T56_nbLoci

    assert(SSP != nullptr);
    this->T56_isMultiplicitySelection = false;
    char isMultfit = -1; // -1 is undefeined. 0 is false, 1 is true
    if (input.PeakNextElementString() == "NA")
    {
        if (this->quickScreenAtL_T56_nbLoci != 0)
        {
            std::cout << "In option '--T5_FitnessEffects', for species " << this->speciesName << ", received 'NA' but there are T5 loci as indiciated in --L (--Loci) (this->Gmap.T5_nbLoci = " << quickScreenAtL_T56_nbLoci << ").\n";
            abort();   
        }
        isMultfit = 0; // whatever but not -1
        input.skipElement();
    } else
    {
        //std::cout << "this->MaxEverHabitat = "<<this->MaxEverHabitat << std::endl ;
        assert(this->MaxEverHabitat >= 0);
        std::string firstMode;
        for (int habitat = 0; habitat <= this->MaxEverHabitat ; habitat++)
        {
            (void) input.GetNextHabitatMarker(habitat);
            std::string Mode = input.GetNextElementString();

            if ((Mode == "MultiplicityA" || Mode == "MultiplicityGamma" || Mode == "MultiplicityUnif" || Mode == "multfitA" || Mode == "multfitGamma" || Mode == "multfitUnif"))
            {
                isMultfit = 1;
            } else
            {
                isMultfit = 0;
            }
                

            if (habitat == 0)
            {
                firstMode = Mode;
            } else
            {
                bool wasFirstMultFit = !(Mode == "MultiplicityA" || Mode == "MultiplicityGamma" || Mode == "MultiplicityUnif" || Mode == "multfitA" || Mode == "multfitGamma" || Mode == "multfitUnif");
                bool isCurrentMultFit = !(firstMode == "MultiplicityA" || firstMode == "MultiplicityGamma" || firstMode == "MultiplicityUnif" || firstMode == "multfitA" || firstMode == "multfitGamma" || firstMode == "multfitUnif");
                if  (wasFirstMultFit != isCurrentMultFit)
                {
                    std::cout << "For --T5_fit, the Mode received for the zeroth habitat is '" << firstMode << "'. The Mode received for the " << habitat << "th habitat is '" << Mode << "'. In the current version, if one patch make the assumption of 'multfit', then all patches must make this assumption.\n";
                    abort();
                }
            }
       
            std::vector<fitnesstype> ForASingleHabitat;
                   
            if (Mode.compare("A")==0)
            {
                this->T56_isMultiplicitySelection = false;
                ForASingleHabitat.reserve(2 * this->quickScreenAtL_T56_nbLoci);
                for (int entry_index = 0 ; entry_index < (2 * this->quickScreenAtL_T56_nbLoci) ; entry_index++)    
                {
                    double fit = input.GetNextElementDouble();
                    if (fit < 0.0)
                    {
                        std::cout << "In option '--T5_FitnessEffects', habitat " << habitat << ", Mode 'A', one 'fit' value is lower than 0 ( is " << fit << " ).\n";
                        abort();
                    }
                    if (isLocusConsideredNeutral(fit, T56_approximationForNtrl)) fit = 1.0;
                    
                    ForASingleHabitat.push_back(fit);
                }
                if (ForASingleHabitat.size() != 2 * this->quickScreenAtL_T56_nbLoci)
                {
                    std::cout << "In option '--T5_FitnessEffects', habitat "  << habitat << ", Mode 'A' " << ForASingleHabitat.size() << " have been received while " << 2 * (this->quickScreenAtL_T56_nbLoci) << " were expected\n";
                    abort();
                }
            } else if (Mode.compare("cstHgamma") == 0)
            {
                this->T56_isMultiplicitySelection = false;
                double h = input.GetNextElementDouble();
                bool isHomoFollows;
                {
                    std::string s = input.GetNextElementString();
                    if (s == "homo" || s == "Homo" || s == "HOMO")
                    {
                        isHomoFollows = true;
                    } else if (s == "hetero" || s == "Hetero" || s == "HETERO")
                    {
                        if (h == 0.0)
                        {
                            std::cout << "In option '--T5_FitnessEffects', habitat " << habitat << ", Mode 'cstHgamma', receive an h value of 0.0, which means that heterozygotes have a fitness of 1.0. But the values are being specified for the heterozygotes (as indicated by keyword "<<s<<").\n";
                            abort();
                        }
                        isHomoFollows = false;
                    } else
                    {
                        std::cout << "In option '--T5_FitnessEffects', habitat " << habitat << ", Mode 'cstHgamma', after the 'h' value, expected either keywor 'homo' or keyword 'hetero' but got '" << s << "' instead.\n";
                        abort();
                    }
                }

                double alpha = input.GetNextElementDouble();
                double beta = input.GetNextElementDouble();
                if (alpha <= 0 || beta <= 0)
                {
                    std::cout << "Either alpha ("<<alpha<<") or beta ("<<beta<<") have non-positive values.\n";
                    abort();
                }
                std::gamma_distribution<double> dist(alpha, beta);

                for (int entry_index = 0 ; entry_index < this->quickScreenAtL_T56_nbLoci ; entry_index++)    
                {
                    double fitHomo;
                    double fitHetero;
                    if (isHomoFollows)
                    {
                        fitHomo = dist(GP->rngw.getRNG());
                        fitHetero = 1 - h * (1 - fitHomo);
                    } else 
                    {
                        fitHetero = dist(GP->rngw.getRNG());
                        fitHomo   = 1 - (1-fitHetero)/h;
                    }
                    //if (fitHomo > 1.0) fitHomo = 1.0;
                    if (fitHomo < 0.0) fitHomo = 0.0;
                    //if (fitHetero > 1.0) fitHetero = 1.0;
                    if (fitHetero < 0.0) fitHetero = 0.0;

                    if (isLocusConsideredNeutral(fitHetero, T56_approximationForNtrl)) fitHetero = 1.0;
                    if (isLocusConsideredNeutral(fitHomo, T56_approximationForNtrl)) fitHomo = 1.0;

                    ForASingleHabitat.push_back(fitHetero);
                    ForASingleHabitat.push_back(fitHomo);
                }
                assert(ForASingleHabitat.size() == 2 * this->quickScreenAtL_T56_nbLoci);


            } else if (Mode.compare("cstH")==0)
            {
                this->T56_isMultiplicitySelection = false;
                double h = input.GetNextElementDouble();
                bool isHomoFollows;
                {
                    std::string s = input.GetNextElementString();
                    if (s == "homo" || s == "Homo" || s == "HOMO")
                    {
                        isHomoFollows = true;
                    } else if (s == "hetero" || s == "Hetero" || s == "HETERO")
                    {
                        if (h == 0.0)
                        {
                            std::cout << "In option '--T5_FitnessEffects', habitat " << habitat << ", Mode 'cstH', receive an h value of 0.0, which means that heterozygotes have a fitness of 1.0. But the values are being specified for the heterozygotes (as indicated by keyword "<<s<<").\n";
                            abort();
                        }
                        isHomoFollows = false;
                    } else
                    {
                        std::cout << "In option '--T5_FitnessEffects', habitat " << habitat << ", Mode 'cstH', after the 'h' value, expected either keywor 'homo' or keyword 'hetero' but got '" << s << "' instead.\n";
                        abort();
                    }
                }
                    
                for (int entry_index = 0 ; entry_index < this->quickScreenAtL_T56_nbLoci ; entry_index++)    
                {
                    double fitHomo;
                    double fitHetero;
                    if (isHomoFollows)
                    {
                        fitHomo = input.GetNextElementDouble();
                        fitHetero = 1 - h * (1 - fitHomo);
                    } else 
                    {
                        fitHetero = input.GetNextElementDouble();
                        fitHomo   = 1 - (1-fitHetero)/h;
                    }

                    if (fitHomo < 0 )
                    {
                        std::cout << "In option '--T5_FitnessEffects', habitat " << habitat << ", Mode 'cstH', one 'fitHomo' value is lower than 0 ( is " << fitHomo << " ).\n";
                        abort();
                    }
                    if (fitHetero < 0 )
                    {
                        std::cout << "In option '--T5_FitnessEffects', habitat " << habitat << ", Mode 'cstH', one 'fitHetero' value is lower than 0 ( is " << fitHetero << ", fitHomo = "<<fitHomo<<", h = "<<h<<" ).\n";
                        abort();
                    }

                    if (isLocusConsideredNeutral(fitHetero, T56_approximationForNtrl)) fitHetero = 1.0;
                    if (isLocusConsideredNeutral(fitHomo, T56_approximationForNtrl)) fitHomo = 1.0;

                    
                    ForASingleHabitat.push_back(fitHetero);
                    ForASingleHabitat.push_back(fitHomo);
                }
                if (ForASingleHabitat.size() != 2 * this->quickScreenAtL_T56_nbLoci)
                {
                    std::cout << "In option '--T5_FitnessEffects', habitat "  << habitat << ", Mode 'cstH' " << ForASingleHabitat.size() << " have been received while " << 2 * (this->quickScreenAtL_T56_nbLoci) << " were expected\n";
                    abort();
                }
            } else if (Mode.compare("MultiplicityA")==0 || Mode.compare("multfitA")==0)
            {
                this->T56_isMultiplicitySelection = true;
                ForASingleHabitat.reserve(this->quickScreenAtL_T56_nbLoci);
                for (int entry_index = 0 ; entry_index < this->quickScreenAtL_T56_nbLoci ; )
                {
                    double fit01 = input.GetNextElementDouble();
                    if (fit01 < 0)
                    {
                        std::cout << "In option '--T5_FitnessEffects', habitat " << habitat << ", Mode 'multfitA', one 'fit01' value is lower than 0 ( is " << fit01 << " ).\n";
                        abort();
                    }

                    if (isLocusConsideredNeutral(fit01, T56_approximationForNtrl)) fit01 = 1.0;

                    ForASingleHabitat.push_back(fit01);
                    entry_index++;
                }
                if (ForASingleHabitat.size() != this->quickScreenAtL_T56_nbLoci)
                {
                    std::cout << "In option '--T5_FitnessEffects', habitat " << habitat << ", Mode 'multfitA' " << ForASingleHabitat.size() << " have been received while " << this->quickScreenAtL_T56_nbLoci << "  were expected\n";
                    abort();
                }
            } else if (Mode.compare("unif")==0)
            {
                ForASingleHabitat.reserve(2 * this->quickScreenAtL_T56_nbLoci);
                this->T56_isMultiplicitySelection = false;
                
                double fitHet = input.GetNextElementDouble();
                double fitHom = input.GetNextElementDouble();
                
                if (fitHet<0 || fitHom<0)
                {
                    std::cout << "In option '--T5_FitnessEffects', habitat " << habitat << ", Mode 'unif', value for fitness are negative (they are " << fitHet << " " << fitHom << " )\n";
                    abort();
                }
                

                if (isLocusConsideredNeutral(fitHet, T56_approximationForNtrl)) fitHet = 1.0;
                if (isLocusConsideredNeutral(fitHom, T56_approximationForNtrl)) fitHom = 1.0;

                for (int locus = 0 ; locus < this->quickScreenAtL_T56_nbLoci ; ++locus)
                {
                    ForASingleHabitat.push_back(fitHet);
                    ForASingleHabitat.push_back(fitHom);
                }
                if (ForASingleHabitat.size() != 2 * this->quickScreenAtL_T56_nbLoci)
                {
                    std::cout << "In option '--T5_FitnessEffects', habitat " << habitat << ", Mode 'unif' " << ForASingleHabitat.size() << " have been received while " << 2 * (this->quickScreenAtL_T56_nbLoci) << " were expected\n";
                    abort();
                }
            } else if (Mode.compare("MultiplicityUnif")==0 || Mode.compare("multfitUnif")==0)
            {
                this->T56_isMultiplicitySelection = true;
                
                ForASingleHabitat.reserve(this->quickScreenAtL_T56_nbLoci);
                double fit01 = input.GetNextElementDouble();
                if (fit01<0)
                {
                    std::cout << "In option '--T5_FitnessEffects' Mode 'multfitUnif', habitat " << habitat << ", value for 'fit01' is negative (is " << fit01 << ")\n";
                    abort();
                }
                if (isLocusConsideredNeutral(fit01, T56_approximationForNtrl)) fit01 = 1.0;
                
                for (int locus = 0 ; locus < this->quickScreenAtL_T56_nbLoci ; locus++)
                {
                    ForASingleHabitat.push_back(fit01);
                }
                if (ForASingleHabitat.size() != this->quickScreenAtL_T56_nbLoci)
                {
                    std::cout << "In option '--T5_FitnessEffects', habitat " << habitat << ", Mode 'multfitUnif' " << ForASingleHabitat.size() << " have been received while " << this->quickScreenAtL_T56_nbLoci << "  were expected (likely an internal error)\n";
                    abort();
                }
            } else if (Mode.compare("multfitGamma")==0 || Mode.compare("MultiplicityGamma")==0)
            {
                this->T56_isMultiplicitySelection = true;
                double alpha = input.GetNextElementDouble();
                double beta = input.GetNextElementDouble();
                if (alpha <= 0 || beta <= 0)
                {
                    std::cout << "Either alpha ("<<alpha<<") or beta ("<<beta<<") have non-positive values.\n";
                    abort();
                }
                std::gamma_distribution<double> dist(alpha, beta);
                ForASingleHabitat.reserve(this->quickScreenAtL_T56_nbLoci);
                
                for (int locus = 0 ; locus < this->quickScreenAtL_T56_nbLoci ; locus++)
                {
                    double fit = 1 - dist(GP->rngw.getRNG());
                    if (fit < 0.0) fit = 0.0;
                    //if (fit > 1.0) fit = 1.0;
                    if (isLocusConsideredNeutral(fit, T56_approximationForNtrl)) fit = 1.0;
                    ForASingleHabitat.push_back(fit);
                }
                assert(ForASingleHabitat.size() == this->quickScreenAtL_T56_nbLoci);
            } else
            {
                std::cout << "Sorry, for option '--T5_fit (--T56_FitnessEffects)', only Modes 'A', 'unif', 'cstHgamma', 'multfitA' (aka. 'MultiplicityA'), 'multfitUnif' (aka. 'MultiplicityUnif') and 'multfitGamma' (aka. 'MultiplicityGamma') are implemented for the moment (received Mode " << Mode << ")" << std::endl;
                abort();
            } // end of ifelse Mode
            this->T56_FitnessEffects.push_back(ForASingleHabitat);
        } // end of for habitat
    
        assert(this->T56_FitnessEffects.size() == this->MaxEverHabitat + 1);
    }

    assert(isMultfit == 1 || isMultfit == 0);
    
    Gmap.setT56GenderLoci(T56_FitnessEffects, isMultfit == 1);
    

    // Set default values for T56_compress that will be changed later if the user wants
    uint32_t localCount_T56ntrl_nbLoci = 0;
    uint32_t localCount_T56sel_nbLoci = 0;
    for (uint32_t locus = 0 ; locus < this->quickScreenAtL_T56_nbLoci ; locus++)
    {
        //std::cout << "line 3188 locus " << locus << "\n";
        if (!Gmap.isT56neutral(locus))
        {
            ++localCount_T56sel_nbLoci;
        } else
        {
            ++localCount_T56ntrl_nbLoci;
        }
    }

    if (localCount_T56ntrl_nbLoci <= 65534 && localCount_T56ntrl_nbLoci > 5e3) // 65535 is the largest unsigned short (65535 = 2^16-1) and we are actually wasting a bit at every operation in order to gain .
    {
        Gmap.isT56ntrlCompress = true;
    } else
    {
        Gmap.isT56ntrlCompress = false;
    }

    Gmap.isT56selCompress = false;
    /*if (localCount_T56sel_nbLoci <= 65534) // 65535 is the largest unsigned short (65535 = 2^16-1) and we are actually wasting a bit at every operation in order to gain .
    {
        T56sel_compress = true;
    } else
    {
        T56sel_compress = false;
    }*/
        
}

/*bool SpeciesSpecificParameters::isT56LocusUnderSelection(uint32_t locus)
{
    assert(this->MaxEverHabitat >= 0);
    assert(this->T56_FitnessEffects.size() == this->MaxEverHabitat + 1);
    for (uint32_t habitat = 0 ; habitat <= this->MaxEverHabitat ; habitat++)
    {
        if (this->T56_isMultiplicitySelection)
        {
            assert(this->T56_FitnessEffects[habitat].size() > locus);
            if (this->T56_FitnessEffects[habitat][locus] < T56_approximationForNtrl)
            {
                return true;
            }
        } else
        {
            assert(this->T56_FitnessEffects[habitat].size() > 2*(locus) + 1);
            if (this->T56_FitnessEffects[habitat][2*(locus) + 1] < T56_approximationForNtrl || this->T56_FitnessEffects[habitat][2*(locus)] < T56_approximationForNtrl)
            {
                return true;
            }
        }
    }
    return false;
}*/


void SpeciesSpecificParameters::readT8_FitnessEffects(InputReader& input)
{
   #ifdef DEBUG
    std::cout << "For option '--T8_FitnessEffects', the std::string that is read is: " << input.print() << std::endl;
    #endif
    assert(SSP != nullptr);
    if (input.PeakNextElementString() == "NA")
    {
        if (this->Gmap.T8_nbLoci != 0)
        {
            std::cout << "In option '--T8_FitnessEffects', for species "<<this->speciesName<<", received 'NA' but there are T8 loci as indiciated in --L (--Loci) (this->Gmap.T8_nbLoci = "<<this->Gmap.T8_nbLoci<<").\n";
            abort();   
        }
        input.skipElement();
    } else
    {       
        std::string Mode = input.GetNextElementString();
           
        if (Mode.compare("MultiplicityA")==0 || Mode.compare("multfitA")==0 || Mode.compare("A")==0)
        {
            for (int entry_index = 0 ; entry_index < this->Gmap.T8_nbLoci ; )    
            {
                if (input.isNextRKeyword())
                {
                    double fit01 = input.GetNextElementDouble();
                    int nbRepeats = input.GetNextElementInt();
                    if (fit01 < 0 )
                    {
                        std::cout << "In option '--T8_FitnessEffects', Mode 'multfitA', one 'fit01' value (coming after keyword 'R') is lower than 0 ( is " << fit01 << " ).\n";
                        abort();
                    }
                    for (int i = 0 ; i < nbRepeats; i++)
                    {
                        this->T8_FitnessEffects.push_back(fit01);
                        entry_index++;
                        if (entry_index > this->Gmap.T8_nbLoci)
                        {
                            std::cout << "In 'T8_FitnessEffects', mode 'multfitA', when reading after keyword 'R', you asked to repeat the value " << fit01 << " for " << nbRepeats << " times. This (on top with previous entries) results of a number of entry greater than " << this->Gmap.T8_nbLoci << " as requested given the mode here (multfitA) and the arguments given to the option '--Loci'\n";
                            abort();
                        }
                    }
                } else
                {
                    double fit01 = input.GetNextElementDouble();
                    if (fit01 < 0 )
                    {
                        std::cout << "In option '--T8_FitnessEffects', Mode 'MultiplicityA', one 'fit01' value is lower than 0 ( is " << fit01 << " ).\n";
                        abort();
                    }
                    this->T8_FitnessEffects.push_back(fit01);
                    entry_index++;
                }
            }
            if (this->T8_FitnessEffects.size() != this->Gmap.T8_nbLoci)
            {
                std::cout << "In option '--T8_FitnessEffects', Mode 'multfitA' " << this->T8_FitnessEffects.size() << " have been received while " << this->Gmap.T8_nbLoci << " were expected\n";
                abort();
            }
        } else if (Mode.compare("MultiplicityUnif")==0 || Mode.compare("multfitUnif")==0 || Mode.compare("unif")==0)
        {
            double fit01 = input.GetNextElementDouble();
            if (fit01<0)
            {
                std::cout << "In option '--T8_FitnessEffects' Mode 'multfitUnif', value for 'fit01' is negative (is " << fit01 << ")\n";
                abort();
            }
            for (int locus = 0 ; locus < this->Gmap.T8_nbLoci ; locus++)
            {
                this->T8_FitnessEffects.push_back(fit01);
            }
            if (this->T8_FitnessEffects.size() != this->Gmap.T8_nbLoci)
            {
                std::cout << "In option '--T8_FitnessEffects', Mode 'multfitUnif' " << this->T8_FitnessEffects.size() << " have been received while " << this->Gmap.T8_nbLoci << " were expected\n";
                abort();
            }
        } else if (Mode.compare("MultiplicityGamma")==0 || Mode.compare("multfitGamma")==0 || Mode.compare("gamma")==0)
        {
            double alpha = input.GetNextElementDouble();
            double beta = input.GetNextElementDouble();
            if (alpha <= 0 || beta <= 0)
            {
                std::cout << "Either alpha ("<<alpha<<") or beta ("<<beta<<") have non-positive values.\n";
                abort();
            }
            std::gamma_distribution<double> dist(alpha, beta);
            for (int locus = 0 ; locus < this->Gmap.T8_nbLoci ; locus++)
            {
                double fit = 1 - dist(GP->rngw.getRNG());
                if (fit < 0.0) fit = 0.0;
                //if (fit > 1.0) fit = 1.0;
                this->T8_FitnessEffects.push_back(fit);
            }
        } else
        {
            std::cout << "Sorry, for option '--T8_fit (--T8_fitnessEffects)', only Modes 'A' (or 'multfitA' or 'MultiplicityA'), 'unif' (or 'multfitUnif' or 'MultiplicityUnif'). Mode received is '" << Mode << "'\n";
            abort();
        } // end of ifelse Mode
    }
}



void SpeciesSpecificParameters::readT1_FitnessEffects(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option '--T1_FitnessEffects', the std::string that is read is: " << input.print() << std::endl;
#endif
    assert(SSP != nullptr);
    if (input.PeakNextElementString() == "NA")
    {
        if (this->Gmap.T1_nbLoci != 0)
        {
            std::cout << "In option '--T1_FitnessEffects', for species "<<this->speciesName<<", received 'NA' but there are T1 loci as indiciated in --L (--Loci) (this->Gmap.T1_nbLoci = "<<this->Gmap.T1_nbLoci<<").\n";
            abort();   
        }
        input.skipElement();
    } else
    {
        std::string firstMode;
        for (int habitat = 0; habitat <= this->MaxEverHabitat ; habitat++)
        {
            (void) input.GetNextHabitatMarker(habitat);
            std::string Mode = input.GetNextElementString();
            if (habitat == 0)
            {
                firstMode = Mode;
            } else
            {
                bool wasFirstMultFit = !(Mode == "MultiplicityA" || Mode == "MultiplicityGamma" || Mode == "MultiplicityUnif" || Mode == "multfitA" || Mode == "multfitGamma" || Mode == "multfitUnif");
                bool isCurrentMultFit = !(firstMode == "MultiplicityA" || firstMode == "MultiplicityGamma" || firstMode == "MultiplicityUnif" || firstMode == "multfitA" || firstMode == "multfitGamma" || firstMode == "multfitUnif");
                if  (wasFirstMultFit != isCurrentMultFit)
                {
                    std::cout << "For --T1_fit, the Mode received for the zeroth habitat is '" << firstMode << "'. The Mode received for the " << habitat << "th habitat is '" << Mode << "'. In the current version, if one patch make the assumption of 'multfit', then all patches must make this assumption.\n";
                    abort();
                }
            }
       
            std::vector<fitnesstype> ForASingleHabitat;
            if (Mode.compare("DomA")==0 || Mode.compare("domA")==0)
            {
            
                this->T1_isMultiplicitySelection = false;
                std::string cst_or_fun = input.GetNextElementString();
                if (cst_or_fun.compare("cst") != 0 && cst_or_fun.compare("fun") != 0)
                {
                    std::cout << "In option '--T1_FitnessEffects', habitat " << habitat << ", Mode 'domA', you had to indicate 'cst' or 'fun' for how dominance coefficient is defined but instead, SimBit received " << cst_or_fun << ".\n";
                    abort();
                }
            

                double dominance = input.GetNextElementDouble(); // either cst or mean dominance
                if (dominance < 0.0 || dominance > 1.0)
                {
                    std::cout << "In option '--T1_FitnessEffects', habitat " << habitat << ", Mode 'domA', the dominance value is either smaller than 0.0 or larger than 1.0. dominance recevied = " << dominance << std::endl;
                    abort();
                }
                double mean_s = 0; // used only if cst_or_fun is 'fun'
            
                for (int entry_index = 0 ; entry_index < this->Gmap.T1_nbLoci ; )    
                {
                
                    if (input.isNextRKeyword())
                    {
                        double fit11 = input.GetNextElementDouble();
                        if (fit11 < 0)
                        {
                            std::cout << "In option '--T1_FitnessEffects', habitat " << habitat << ", Mode 'domA', one 'fit' value (directly after keyword 'R') is lower than 0 ( is " << fit11 << " ).\n";
                            abort();
                        }
                        double fit01;
                        if (cst_or_fun.compare("cst") == 0)
                        {
                            fit01 = 1 - ((1-fit11) * dominance);
                        } else if (cst_or_fun.compare("fun") == 0)
                        {
                            fit01 = -1.0; // this is temporary because we need to read all s, before being able to compute the locus specific dominance
                            mean_s += 1 - fit11;
                        } else
                        {
                            std::cout << "Internal error in Parameters::readT1_FitnessEffects.\n";
                            abort();
                        }
                            
                        int nbRepeats = input.GetNextElementInt();
                        for (int i = 0 ; i < nbRepeats ; i++)
                        {
                            ForASingleHabitat.push_back(1.0);
                            ForASingleHabitat.push_back(fit01);
                            ForASingleHabitat.push_back(fit11);
                            entry_index++;
                            if (entry_index > this->Gmap.T1_nbLoci)
                            {
                                std::cout << "In 'T1_FitnessEffects', mode 'DomA', when reading after keyword 'R', you asked to repeat the value " << fit11 << " for " << nbRepeats << " times. This (on top with previous entries) results of a number of entry greater than " << this->Gmap.T1_nbLoci << " as requested given the mode here ('CstDomA') and the arguments given to the option '--Loci'\n";
                                abort();
                            }
                        }
                    } else
                    {
                    
                        double fit11 = input.GetNextElementDouble();
                        if (fit11 < 0 || fit11 > 1 )
                        {
                            std::cout << "In option '--T1_FitnessEffects', habitat " << habitat << ", Mode 'A', one 'fit' value is either greater than 1 or lower than 0 ( is " << fit11 << " ).\n";
                            abort();
                        }
                    
                        double fit01;
                        if (cst_or_fun.compare("cst") == 0)
                        {
                            fit01 = 1 - ((1-fit11) * dominance);
                        } else if (cst_or_fun.compare("fun") == 0)
                        {
                            fit01 = -1.0; // this is temporary because we need to read all s, before being able to compute the locus specific dominance
                            mean_s += 1 - fit11;
                        } else
                        {
                            std::cout << "Internal error in Parameters::readT1_FitnessEffects.\n";
                            abort();
                        }

                        ForASingleHabitat.push_back(1.0);
                        ForASingleHabitat.push_back(fit01);
                        ForASingleHabitat.push_back(fit11);
                        entry_index++;
                    
                    }
                
                }  
                if (ForASingleHabitat.size() != THREE * this->Gmap.T1_nbLoci)
                {
                    std::cout << "In option '--T1_FitnessEffects', habitat "  << habitat << ", Mode 'DomA' " << ForASingleHabitat.size() / 3 << " have been received while " << this->Gmap.T1_nbLoci << " (= 8 * ( " << this->Gmap.T1_nbChars << "))  were expected\n";
                    abort();
                }

                // Now if used the function, compute locus-specific dominance and insert them in place of the current -1.0
                if (cst_or_fun.compare("fun") == 0)
                {
                    mean_s /= this->Gmap.T1_nbLoci;
                    assert(mean_s >= 0.0 && mean_s <= 1.0);
                    if (mean_s == 0.0) {std::cout << "All T1 are neutral!\n";}
                    double k = -1.0;
                    if (mean_s > 0.0)
                    {
                        k = log(2 * dominance) / mean_s; // See Nemo manual for definition
                    }
                    int securityCount = 0;
                    for (int i = 2 ; i < ForASingleHabitat.size() ; i += 3)
                    {
                        double fit11 = ForASingleHabitat[i];
                        assert(fit11 >= 0.0 && fit11 <= 1.0);
                        double s11 = 1 - fit11;
                        assert(s11 >= 0.0 && s11 <= 1.0);
                        double dom = -1.0;
                        double fit01 = 1.0;
                        if (mean_s > 0.0)
                        {
                            dom = exp(k * s11) / 2;   
                            assert(dom >= 0.0 && dom <= 1.0);
                            fit01 = 1 - (s11 * dom);
                        }
                        assert(fit01 >= 0.0 && fit01 <= 1.0);
                        assert(ForASingleHabitat[i-2] == 1.0);  // homozygote wildtype
                        assert(ForASingleHabitat[i-1] == -1.0); // heterozygote should currently be set to -1.0
                        ForASingleHabitat[i-1] = fit01;
                        securityCount++;
                    }
                    assert(securityCount == this->Gmap.T1_nbLoci);
                }
                    
            } else if (Mode.compare("A")==0)
            {
                this->T1_isMultiplicitySelection = false;
                for (int entry_index = 0 ; entry_index < THREE * this->Gmap.T1_nbLoci ; )    
                {
                    if (input.isNextRKeyword())
                    {
                        double fit00 = input.GetNextElementDouble();
                        double fit01 = input.GetNextElementDouble();
                        double fit11 = input.GetNextElementDouble();
                        int nbRepeats = input.GetNextElementInt();
                        for (int i = 0 ; i < nbRepeats ; i++)
                        {
                            ForASingleHabitat.push_back(fit00);
                            ForASingleHabitat.push_back(fit01);
                            ForASingleHabitat.push_back(fit11);
                            entry_index += 3;
                            if (entry_index > (THREE * this->Gmap.T1_nbLoci))
                            {
                                std::cout << "In 'T1_FitnessEffects', mode 'A', when reading after keyword 'R', you asked to repeat the values " << fit00 << " " << fit01 << " " << fit11 << " for " << nbRepeats << " times. This (on top with previous entries) results of a number of entry greater than " << THREE * this->Gmap.T1_nbLoci << " as requested given the mode here ('A') and the arguments given to the optinn '--Loci'\n";
                                abort();
                            }
                        }
                    } else
                    {
                        double fit = input.GetNextElementDouble();
                        if (fit < 0)
                        {
                            std::cout << "In option '--T1_FitnessEffects', habitat " << habitat << ", Mode 'A', one 'fit' value is lower than 0 ( is " << fit << " ).\n";
                            abort();
                        }
                        ForASingleHabitat.push_back(fit);
                        entry_index++;
                    }
                }
                if (ForASingleHabitat.size() != THREE * this->Gmap.T1_nbLoci)
                {
                    std::cout << "In option '--T1_FitnessEffects', habitat "  << habitat << ", Mode 'A' " << ForASingleHabitat.size() << " have been received while " << THREE * (this->Gmap.T1_nbLoci) << " (= 3 * 8 * ( " << this->Gmap.T1_nbChars << "))  were expected\n";
                    abort();
                }
            } else if (Mode.compare("MultiplicityA")==0 || Mode.compare("multfitA")==0)
            {
                this->T1_isMultiplicitySelection = true;
                for (int entry_index = 0 ; entry_index < this->Gmap.T1_nbLoci ; )    
                {
                    if (input.isNextRKeyword())
                    {
                        double fit01 = input.GetNextElementDouble();
                        int nbRepeats = input.GetNextElementInt();
                        if (fit01 < 0 )
                        {
                            std::cout << "In option '--T1_FitnessEffects', habitat " << habitat << ", Mode 'multfitA', one 'fit01' value (coming after keyword 'R') is lower than 0 ( is " << fit01 << " ).\n";
                            abort();
                        }
                        for (int i = 0 ; i < nbRepeats; i++)
                        {
                            ForASingleHabitat.push_back(fit01);
                            entry_index++;
                            if (entry_index > this->Gmap.T1_nbLoci)
                            {
                                std::cout << "In 'T1_FitnessEffects', mode 'multfitA', when reading after keyword 'R', you asked to repeat the value " << fit01 << " for " << nbRepeats << " times. This (on top with previous entries) results of a number of entry greater than " << this->Gmap.T1_nbLoci << " as requested given the mode here (multfitA) and the arguments given to the option '--Loci'\n";
                                abort();
                            }
                        }
                    } else
                    {
                        double fit01 = input.GetNextElementDouble();
                        if (fit01 < 0 )
                        {
                            std::cout << "In option '--T1_FitnessEffects', habitat " << habitat << ", Mode 'MultiplicityA', one 'fit01' value is lower than 0 ( is " << fit01 << " ).\n";
                            abort();
                        }
                        ForASingleHabitat.push_back(fit01);
                        entry_index++;
                    }
                }
                if (ForASingleHabitat.size() != this->Gmap.T1_nbLoci)
                {
                    std::cout << "In option '--T1_FitnessEffects', habitat " << habitat << ", Mode 'multfitA' " << ForASingleHabitat.size() << " have been received while " << this->Gmap.T1_nbLoci << " (= 8 * ( " << this->Gmap.T1_nbChars << "))  were expected\n";
                    abort();
                }
            } else if (Mode.compare("unif")==0)
            {
                this->T1_isMultiplicitySelection = false;
                
                double fit00 = input.GetNextElementDouble();
                double fit01 = input.GetNextElementDouble();
                double fit11 = input.GetNextElementDouble();
                if (fit00<0 || fit01<0 || fit11<0)
                {
                    std::cout << "In option '--T1_FitnessEffects', habitat " << habitat << ", Mode 'unif', value for fitness are negative (they are " << fit00 << " " << fit01 << " " << fit11 << " )\n";
                    abort();
                }
                
                for (int locus = 0 ; locus < this->Gmap.T1_nbLoci ; ++locus)
                {
                    ForASingleHabitat.push_back(fit00);
                    ForASingleHabitat.push_back(fit01);
                    ForASingleHabitat.push_back(fit11);
                }
                if (ForASingleHabitat.size() != THREE * this->Gmap.T1_nbLoci)
                {
                    std::cout << "In option '--T1_FitnessEffects', habitat " << habitat << ", Mode 'unif' " << ForASingleHabitat.size() << " have been received while " << THREE * (this->Gmap.T1_nbLoci) << " (= 3 * 8 * ( " << this->Gmap.T1_nbChars << "))  were expected\n";
                    abort();
                }
            } else if (Mode.compare("MultiplicityUnif")==0 || Mode.compare("multfitUnif")==0)
            {
                this->T1_isMultiplicitySelection = true;
                
                double fit01 = input.GetNextElementDouble();
                if (fit01<0)
                {
                    std::cout << "In option '--T1_FitnessEffects' Mode 'multfitUnif', habitat " << habitat << ", value for 'fit01' is negative (is " << fit01 << ")\n";
                    abort();
                }
                for (int locus = 0 ; locus < this->Gmap.T1_nbLoci ; locus++)
                {
                    ForASingleHabitat.push_back(fit01);
                }
                if (ForASingleHabitat.size() != this->Gmap.T1_nbLoci)
                {
                    std::cout << "In option '--T1_FitnessEffects', habitat " << habitat << ", Mode 'multfitUnif' " << ForASingleHabitat.size() << " have been received while " << this->Gmap.T1_nbLoci << " (= 8 * ( " << this->Gmap.T1_nbChars << "))  were expected\n";
                    abort();
                }
            } else if (Mode.compare("MultiplicityGamma")==0 || Mode.compare("multfitGamma")==0)
            {
                double alpha = input.GetNextElementDouble();
                double beta = input.GetNextElementDouble();
                if (alpha <= 0 || beta <= 0)
                {
                    std::cout << "Either alpha ("<<alpha<<") or beta ("<<beta<<") have non-positive values.\n";
                    abort();
                }
                std::gamma_distribution<double> dist(alpha, beta);
                for (int locus = 0 ; locus < this->Gmap.T1_nbLoci ; locus++)
                {
                    double fit = 1 - dist(GP->rngw.getRNG());
                    if (fit < 0.0) fit = 0.0;
                    //if (fit > 1.0) fit = 1.0;
                    ForASingleHabitat.push_back(fit);
                }
            } else if (Mode.compare("cstHgamma") == 0)
            {
                this->T1_isMultiplicitySelection = false;
                double h = input.GetNextElementDouble();
                bool isHomoFollows;
                {
                    std::string s = input.GetNextElementString();
                    if (s == "homo" || s == "Homo" || s == "HOMO")
                    {
                        isHomoFollows = true;
                    } else if (s == "hetero" || s == "Hetero" || s == "HETERO")
                    {
                        if (h == 0.0)
                        {
                            std::cout << "In option '--T1_FitnessEffects', habitat " << habitat << ", Mode 'cstHgamma', receive an h value of 0.0, which means that heterozygotes have a fitness of 1.0. But the values are being specified for the heterozygotes (as indicated by keyword "<<s<<").\n";
                            abort();
                        }
                        isHomoFollows = false;
                    } else
                    {
                        std::cout << "In option '--T1_FitnessEffects', habitat " << habitat << ", Mode 'cstHgamma', after the 'h' value, expected either keywor 'homo' or keyword 'hetero' but got '" << s << "' instead.\n";
                        abort();
                    }
                }

                double alpha = input.GetNextElementDouble();
                double beta = input.GetNextElementDouble();
                if (alpha <= 0 || beta <= 0)
                {
                    std::cout << "Either alpha ("<<alpha<<") or beta ("<<beta<<") have non-positive values.\n";
                    abort();
                }
                std::gamma_distribution<double> dist(alpha, beta);

                for (int entry_index = 0 ; entry_index < this->Gmap.T1_nbLoci ; entry_index++)    
                {
                    double fitHomo;
                    double fitHetero;
                    if (isHomoFollows)
                    {
                        fitHomo = dist(GP->rngw.getRNG());
                        fitHetero = 1 - h * (1 - fitHomo);
                    } else 
                    {
                        fitHetero = dist(GP->rngw.getRNG());
                        fitHomo   = 1 - (1-fitHetero)/h;
                    }
                    //if (fitHomo > 1.0) fitHomo = 1.0;
                    if (fitHomo < 0.0) fitHomo = 0.0;
                    //if (fitHetero > 1.0) fitHetero = 1.0;
                    if (fitHetero < 0.0) fitHetero = 0.0;

                    
                    ForASingleHabitat.push_back(1.0);
                    ForASingleHabitat.push_back(fitHetero);
                    ForASingleHabitat.push_back(fitHomo);
                }
                assert(ForASingleHabitat.size() == 3 * this->Gmap.T1_nbLoci);


            } else
            {
                std::cout << "Sorry, for option '--T1_fit (--T1_fitnessEffects)', only Modes 'A', 'unif', 'domA' (or 'DomA'), 'multfitA' (aka. 'MultiplicityA'), 'multfitUnif' (aka. 'MultiplicityUnif'), cstHgamma and 'multfitGamma' (aka. 'MultiplicityGamma') are implemented for the moment (received Mode " << Mode << ")" << std::endl;
                abort();
            } // end of ifelse Mode
            this->T1_FitnessEffects.push_back(ForASingleHabitat);
        } // end of for habitat
    
        assert(this->T1_FitnessEffects.size() == this->MaxEverHabitat + 1);
    }
        
    
}


void SpeciesSpecificParameters::readT2_FitnessEffects(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option '--T2_FitnessEffects', the std::string that is read is: " << input.print() << std::endl;
#endif

    if (input.PeakNextElementString() == "NA")
    {
        if (this->Gmap.T2_nbLoci != 0)
        {
            std::cout << "In option '--T2_FitnessEffects', for species "<<this->speciesName<<", received 'NA' but there are T2 loci as indiciated in --L (--Loci) (this->Gmap.T2_nbLoci = "<<this->Gmap.T2_nbLoci<<").\n";
            abort();   
        }
        input.skipElement();
    } else
    {
        for (int habitat = 0 ; habitat <= this->MaxEverHabitat ; habitat++)
        {
            (void) input.GetNextHabitatMarker(habitat);

            std::vector<fitnesstype> ForASingleHabitat;
            
            std::string Mode = input.GetNextElementString();
            if (Mode.compare("unif") == 0)
            {
                double fit = input.GetNextElementDouble();
                if (fit < 0.0)
                {
                    std::cout <<  "In option '--T2_FitnessEffects', habitat " << habitat << ", Mode 'unif', 'fit' value received is " << fit << " which is lower than zero.\n";
                    abort();
                }
                for (int char_index = 0 ; char_index < this->Gmap.T2_nbLoci ; char_index++)
                {
                    ForASingleHabitat.push_back(fit);
                }
            } else if (Mode.compare("A") == 0)
            {
                for (int entry_index = 0 ; entry_index < this->Gmap.T2_nbLoci ; ++entry_index)
                {
                    double fit = input.GetNextElementDouble();
                    if (fit < 0.0)
                    {
                        std::cout <<  "In option '--T2_FitnessEffects' Mode 'A', 'fit' value received is " << fit << " which is lower than zero.\n";
                        abort();
                    }
                    ForASingleHabitat.push_back(fit);
                }
                if (ForASingleHabitat.size() != this->Gmap.T2_nbLoci)
                {
                    std::cout << "In option '--T2_FitnessEffects', habitat " << habitat << ", Mode 'A', " << ForASingleHabitat.size() << " values received while " << this->Gmap.T2_nbLoci << " were expected." << std::endl;
                    abort();
                }
            } else if (Mode.compare("Gamma") == 0)
            {
                double alpha = input.GetNextElementDouble();
                double beta = input.GetNextElementDouble();
                if (alpha <= 0 || beta <= 0)
                {
                    std::cout << "Either alpha ("<<alpha<<") or beta ("<<beta<<") have non-positive values.\n";
                    abort();
                }
                std::gamma_distribution<double> dist(alpha, beta);
                for (int locus = 0 ; locus < this->Gmap.T2_nbLoci ; locus++)
                {
                    double fit = 1 - dist(GP->rngw.getRNG());
                    if (fit < 0.0) fit = 0.0;
                    //if (fit > 1.0) fit = 1.0;
                    ForASingleHabitat.push_back(fit);
                }
            } else
            {
                std::cout << "In option '--T2_FitnessEffects', habitat " << habitat << ", Mode received is " << Mode << ". Sorry only Modes 'A' and 'unif' are recognized for the moment!\n";
                abort();
            } // end of ifelse Mode
            assert(this->T2_FitnessEffects.size() == habitat);
            this->T2_FitnessEffects.push_back(ForASingleHabitat);
        } // end of Habitat loop

        assert(this->T2_FitnessEffects.size() == this->MaxEverHabitat + 1);
    }
}

void SpeciesSpecificParameters::readT56_MutationRate(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option '--T5_MutationRate', the std::string that is read is: " << input.print() << std::endl;
#endif
    //std::cout << "T5_nbLoci = " << T5_nbLoci << "\n";

    this->T56_Total_Mutation_rate = 0.0; // just for safety but it is not necessary

    if (input.PeakNextElementString() == "NA" || input.PeakNextElementString() == "default")
    {
        this->T56_Total_Mutation_rate = 0.0;
        if (this->Gmap.T56_nbLoci != 0)
        {
            std::cout << "In option '--T5_MutationRate', for species "<<this->speciesName<<", received 'default' (default input). The default input is only valid when there are no T5 loci indicated in option --L (--Loci). In short, if you want T5 loci, please indicate a mutation rate for these loci (with option --T5_mu). Note SimBit received "<<this->Gmap.T56_nbLoci<<" T5 loci from option --L (--Loci)\n";
            abort();   
        }
        input.skipElement();
    } else
    {
        std::string Mode = input.GetNextElementString();
        if (Mode.compare("A")==0)
        {
            double currentSum = 0.0;
            while (input.IsThereMoreToRead())
            {
                if (input.isNextRKeyword())
                {
                    double inputValue = input.GetNextElementDouble();
                    int nbRepeats     = input.GetNextElementInt();

                    if (inputValue < 0.0)
                    {
                        std::cout << "For option 'T5_MutationRate', received a negative rate (received " << inputValue << ").\n";
                        abort();
                    }

                    if (inputValue > 0.5)
                    {
                        std::cout << "For option 'T5_MutationRate', received a value greater than 0.5 (received " << inputValue << ").\n";
                        abort();
                    }

                    for (uint32_t repeat = 0 ; repeat < nbRepeats ; ++repeat)
                    {
                        currentSum += inputValue;
                        T56_MutationRate.push_back(currentSum);
                    }
                } else
                {
                    double inputValue = input.GetNextElementDouble();
                    if (inputValue < 0.0)
                    {
                        std::cout << "For option 'T5_MutationRate', received a negative rate (received " << inputValue << ").\n";
                        abort();
                    }
                    currentSum += inputValue;
                    T56_MutationRate.push_back(currentSum);
                }
            }

            if (T56_MutationRate.size() != Gmap.T56_nbLoci)
            {
                std::cout << "For option 'T5_MutationRate', received " << T56_MutationRate.size() << " values but it was expected "<< Gmap.T56_nbLoci <<" values.\n";
                        abort();
            } 

            T56_Total_Mutation_rate = T56_MutationRate.back();

        } else if (Mode.compare("unif")==0)
        {
            double perLocusRate = input.GetNextElementDouble();

            if (perLocusRate < 0)
            {
                std::cout << "T56_MutationRate of Mode 'unif' is lower than zero" << std::endl;
                abort();
            }

            if (perLocusRate > 0.5)
            {
                std::cout << "For option 'T5_MutationRate', received a value greater than 0.5 (received " << perLocusRate << ").\n";
                abort();
            }

            if (Gmap.T56_nbLoci > 0)
            {
                this->T56_MutationRate.push_back(perLocusRate * this->Gmap.T56_nbLoci);
            }

            T56_Total_Mutation_rate = T56_MutationRate.back();
                
        } else
        {
            std::cout << "Sorry, for 'T56_MutationRate' only Mode 'unif' and 'A' are recognized so far. Mode received (" << Mode << ") is unknown!" << std::endl;
            abort();
        }
        
        assert(this->T56_MutationRate.size() == this->Gmap.T56_nbLoci || this->T56_MutationRate.size() == 1);
    }    

    //std::cout << "T5ntrl_Total_Mutation_rate = "<< T5ntrl_Total_Mutation_rate << "\n";
    //std::cout << "T5sel_Total_Mutation_rate = "<< T5sel_Total_Mutation_rate << "\n";
}

void SpeciesSpecificParameters::readT1_mutDirection(InputReader& input)
{
    if (input.PeakNextElementString() == "default")
    {
        T1_mutDirection = 2;
        input.skipElement();
    } else
    {
        std::string s = input.GetNextElementString();
        if (s == "0to1")
        {
            T1_mutDirection = 1;
        } else if (s == "1to0")
        {
            T1_mutDirection = 0;
        } else if (s == "both")
        {
            T1_mutDirection = 2;
        } else
        {
            std::cout << "For option --T1_mutDirection, only values '1to0', '0to1' and 'both' (or 'default') are accepted. Received '" << s << "'\n";
            abort();
        }
    }
    assert(T1_mutDirection >= 0 && T1_mutDirection <= 2);
}

void SpeciesSpecificParameters::readT4_nbMutationPlacingsPerOutput(InputReader& input)
{
    if (input.PeakNextElementString() == "default")
    {
        input.skipElement();
        T4_nbMutationPlacingsPerOutput = 1;
    } else
    {
        auto x = input.GetNextElementInt();
        if (x < 1)
        {
            std::cout << "For option --T4_nbMutationPlacingsPerOutput, you asked for " << x << " mutations placings. Sorry, the number of mutation placings must be great or equal to 1.\n";
            abort();
        }

        T4_nbMutationPlacingsPerOutput = x; 
    }
}

void SpeciesSpecificParameters::readT4_respectPreviousMutationsOutputs(InputReader& input)
{
    if (input.PeakNextElementString() == "default")
    {
        input.skipElement();
        if (T4_nbMutationPlacingsPerOutput == 1 && !outputWriter.isFile(T4_paintedHaplo_file) && !outputWriter.isFile(T4_paintedHaploSegmentsDiversity_file) && !outputWriter.isFile(T4_paintedHaploSegmentsOrigin_file))
        {
            T4_respectPreviousMutationsOutputs = true;
        }
        else
        {
            T4_respectPreviousMutationsOutputs = false;
        }
    } else
    {
        T4_respectPreviousMutationsOutputs = input.GetNextElementBool();
        if (T4_respectPreviousMutationsOutputs)
        {
            if (T4_nbMutationPlacingsPerOutput > 1)
            {
                std::cout << "For options --T4_respectPreviousMutationsOutputs, received 'true' (or equivalent).  T4_respectPreviousMutationsOutputs can only be set to true if you only ask for a single mutation placing per T4 output (as set with option --T4_nbMutationPlacingsPerOutput; because it would not make any sense otherwise). However, you appear to have asked for " << T4_nbMutationPlacingsPerOutput << " mutation placings per T4 outputs.\n";
                abort();
            }

            if (outputWriter.isFile(T4_paintedHaplo_file) || outputWriter.isFile(T4_paintedHaploSegmentsDiversity_file) || outputWriter.isFile(T4_paintedHaploSegmentsOrigin_file))
            {
                std::cout << "For options --T4_respectPreviousMutationsOutputs, received 'true' (or equivalent).  T4_respectPreviousMutationsOutputs can only be set to true if you do not use haplotype painted (but you seem to use haploytpe painting). Sorry, that's a little bit of a silly limitation of SimBit.\n";
                abort();
            }
        }
    }
}




void SpeciesSpecificParameters::readT4_mutDirection(InputReader& input)
{
    if (input.PeakNextElementString() == "default")
    {
        T4_mutDirection = 2;
        input.skipElement();
    } else
    {
        std::string s = input.GetNextElementString();
        if (s == "0to1")
        {
            T4_mutDirection = 1;
        } else if (s == "1to0")
        {
            T4_mutDirection = 0;
        } else if (s == "both")
        {
            T4_mutDirection = 2;
        } else
        {
            std::cout << "For option --T4_mutDirection, only values '1to0', '0to1' and 'both' (or 'default') are accepted. Redeived " << s << "\n";
            abort();
        }
    }
    assert(T4_mutDirection >= 0 && T4_mutDirection <= 2);
}


void SpeciesSpecificParameters::readT8_mutDirection(InputReader& input)
{
    if (input.PeakNextElementString() == "default")
    {
        T8_mutDirection = 2;
        input.skipElement();
    } else
    {
        std::string s = input.GetNextElementString();
        if (s == "0to1")
        {
            T8_mutDirection = 1;
        } else if (s == "1to0")
        {
            T8_mutDirection = 0;
        } else if (s == "both")
        {
            T8_mutDirection = 2;
        } else
        {
            std::cout << "For option --T8_mutDirection, only values '1to0', '0to1' and 'both' (or 'default') are accepted. Redeived " << s << "\n";
            abort();
        }
    }
    assert(T8_mutDirection >= 0 && T8_mutDirection <= 2);
}


void SpeciesSpecificParameters::readT56_mutDirection(InputReader& input)
{
    if (input.PeakNextElementString() == "default")
    {
        T56_mutDirection = 2;
        input.skipElement();
    } else
    {
        std::string s = input.GetNextElementString();
        if (s == "0to1")
        {
            T56_mutDirection = 1;
        } else if (s == "1to0")
        {
            T56_mutDirection = 0;
        } else if (s == "both")
        {
            T56_mutDirection = 2;
        } else
        {
            std::cout << "For option --T56_mutDirection, only values '1to0', '0to1' and 'both' (or 'default') are accepted. Redeived " << s << "\n";
            abort();
        }
    }
    assert(T56_mutDirection >= 0 && T56_mutDirection <= 2);
}

void SpeciesSpecificParameters::readT1_MutationRate(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option '--T1_MutationRate', the std::string that is read is: " << input.print() << std::endl;
#endif

    if (input.PeakNextElementString() == "NA")
    {
        if (this->Gmap.T1_nbLoci != 0)
        {
            std::cout << "In option '--T1_MutationRate', for species "<<this->speciesName<<", received 'NA' but there are T1 loci as indiciated in --L (--Loci) (this->Gmap.T1_nbLoci = "<<this->Gmap.T1_nbLoci<<").\n";
            abort();   
        }
        input.skipElement();
    } else
    {
        std::string Mode = input.GetNextElementString();
        if (Mode.compare("A")==0)
        {
            double currentSum = 0;
            int T1_char_index = 0;
            int bit_index = 0;
            while (input.IsThereMoreToRead())
            {
                if (input.isNextRKeyword())
                {
                    double inputValue = input.GetNextElementDouble();
                    int NbRepeats  = input.GetNextElementInt();

                    for (int repeat = 0 ; repeat < NbRepeats ; ++repeat)
                    {
                        currentSum += inputValue;
                        this->T1_MutationRate.push_back(currentSum); 

                        // Adjust the bit_index and T1_char_index
                        if (bit_index == 7)
                        {
                            bit_index = 0;
                            ++T1_char_index;
                        } else
                        {
                            ++bit_index;
                        }   
                    }

                    // security
                    if (inputValue < 0)
                    {
                        std::cout << "In 'T1_MutationRate' of Mode 'A' the value for the " << T1_char_index << "th T1_char_index is negative (is " << inputValue << ")" << std::endl;
                        abort();
                    }
                    if (inputValue >= 0.1)
                    {
                        std::cout << "In 'T1_MutationRate' of Mode 'A' the value for the " << T1_char_index << "th T1_char_index is greater or equal to 0.1 (is " << inputValue << "). Is it really what you want to do?" << std::endl;
                        abort();
                    }

                } else
                {
                    double inputValue = input.GetNextElementDouble();
                    currentSum += inputValue;
                    this->T1_MutationRate.push_back(currentSum);

                    // security
                    if (inputValue < 0)
                    {
                        std::cout << "In 'T1_MutationRate' of Mode 'A' the value for the " << T1_char_index << "th T1_char_index is negative (is " << inputValue << ")" << std::endl;
                        abort();
                    }
                    if (inputValue >= 0.1)
                    {
                        std::cout << "In 'T1_MutationRate' of Mode 'A' the value for the " << T1_char_index << "th T1_char_index is greater than 0.1 (is " << inputValue << "). Is it really what you want to do?" << std::endl;
                        abort();
                    }

                    // Adjust the bit_index and T1_char_index
                    if (bit_index == 7)
                    {
                        bit_index = 0;
                        ++T1_char_index;
                    } else
                    {
                        ++bit_index;
                    }             
                }
            }
            // Did I get the info I expected?
            if ((T1_char_index * EIGHT + bit_index) != this->Gmap.T1_nbLoci)
            {
                std::cout << "In '--T1_MutationRate', expected " << this->Gmap.T1_nbLoci << " elements. But " << T1_char_index * EIGHT + bit_index << " elements were received.\n";
                abort();
            }

        } else if (Mode.compare("unif")==0)
        {
            this->T1_MutationRate.push_back(input.GetNextElementDouble() * this->Gmap.T1_nbLoci);
            if (this->T1_MutationRate.back() < 0)
            {
                std::cout << "T1_MutationRate of Mode 'unif' is lower than zero" << std::endl;
                abort();
            }
        } else
        {
            std::cout << "Sorry, for 'T1_MutationRate' only Mode 'unif' and 'A' are recognized so far. Mode received (" << Mode << ") is unknown!" << std::endl;
            abort();
        }
        
        assert(this->T1_MutationRate.size() == this->Gmap.T1_nbLoci || this->T1_MutationRate.size() == 1);
        this->T1_Total_Mutation_rate = this->T1_MutationRate.back();
    }

    
}

void SpeciesSpecificParameters::readT2_MutationRate(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option '--T2_MutationRate', the std::string that is read is: " << input.print() << std::endl;
#endif

    if (input.PeakNextElementString() == "NA")
    {
        if (this->Gmap.T2_nbLoci != 0)
        {
            std::cout << "In option '--T2_MutationRate', for species "<<this->speciesName<<", received 'NA' but there are T2 loci as indiciated in --L (--Loci) (this->Gmap.T2_nbLoci = "<<this->Gmap.T2_nbLoci<<").\n";
            abort();   
        }
        input.skipElement();
    } else
    {
        std::string Mode = input.GetNextElementString();
        if (Mode.compare("A")==0)
        {
            double currentSum = 0;
            int T2_char_index = 0;
            while (input.IsThereMoreToRead())
            {
                if (input.isNextRKeyword())
                {
                    double inputValue = input.GetNextElementDouble();
                    int nbRepeats     = input.GetNextElementInt();
                    T2_char_index += nbRepeats;
                    for (int repeat = 0 ; repeat < nbRepeats ; ++repeat)
                    {
                        currentSum += inputValue;
                        this->T2_MutationRate.push_back(currentSum);
                    }

                    // security
                    if (inputValue < 0)
                    {
                        std::cout << "In 'T2_MutationRate' of Mode 'A' the value for the " << T2_char_index << "th T2_char_index is negative (is " << inputValue << ")" << std::endl;
                        abort();
                    }
                    if (inputValue > 0.5)
                    {
                        std::cout << "In 'T2_MutationRate' of Mode 'A' the value for the " << T2_char_index << "th T2_char_index is greater than 0.5 (is " << inputValue << "). Is it really what you want to do?" << std::endl;
                        abort();
                    }
                } else
                {
                    double inputValue = input.GetNextElementDouble();
                    currentSum += inputValue;
                    this->T2_MutationRate.push_back(currentSum);
                    ++T2_char_index;

                    // security
                    if (inputValue < 0)
                    {
                        std::cout << "In 'T2_MutationRate' of Mode 'A' the value for the " << T2_char_index << "th T2_char_index is negative (is " << this->T2_MutationRate.back() << ")" << std::endl;
                        abort();
                    }
                    if (inputValue > 0.5)
                    {
                        std::cout << "In 'T2_MutationRate' of Mode 'A' the value for the " << T2_char_index << "th T2_char_index is greater than 0.5 (is " << this->T2_MutationRate.back() << "). Is it really what you want to do?" << std::endl;
                        abort();
                    }
                }
            }
        } else if (Mode.compare("unif")==0)
        {
            this->T2_MutationRate.push_back(input.GetNextElementDouble() * this->Gmap.T2_nbLoci);
            if (this->T2_MutationRate.back() < 0)
            {
                std::cout << "T2_MutationRate of Mode 'unif' is lower than zero" << std::endl;
                abort();
            }
        } else
        {
            std::cout << "Sorry, for 'T2_MutationRate' only Mode 'unif' and 'A' are recognized so far. Mode received (" << Mode << ") is unknown!" << std::endl;
            abort();
        }
        this->T2_Total_Mutation_rate = this->T2_MutationRate.back();
        assert(this->T2_MutationRate.size() == this->Gmap.T2_nbLoci || this->T2_MutationRate.size() == 1);
        
    }
        
    
}


void SpeciesSpecificParameters::readT3_mutationalEffect(InputReader& input)
{
    #ifdef DEBUG
        std::cout << "For option 'T3_mutationalEffect', the std::string that is read is: " << input.print() << std::endl;
    #endif

    auto type = input.GetNextElementString();
    if (type == "default")
    {
        T3_mutationType = '1';
        T3_mutationType_effectSizeInfo = 1.0;
    }
    else if (type == "cst")
    {
        T3_mutationType = '1';
        T3_mutationType_effectSizeInfo = input.GetNextElementDouble();
    }
    else if (type == "gauss1")
    {
        T3_mutationType = 'g';
        T3_mutationType_effectSizeInfo = input.GetNextElementDouble();
    }
    else if (type == "gauss2") 
    {
        T3_mutationType = 'G';
        T3_mutationType_effectSizeInfo = input.GetNextElementDouble();
    }
    else
    {
        std::cout << "For option 'T3_mutationType', received type '"<< type << "'. Sorry, only types 'default', 'cst', 'gauss1' and 'gauss2' are accepted. Note that 'default' and 'cst' are the same and mean that at each mutation a constant value is either added or removed (with equal probability). For the types 'gauss1' and 'gauss2', at each mutation, SimBit draws an mutation effect from a normal (Gaussian) distribution with specified variance and mean of zero. With 'gauss1' this mutational effect is added to the current allelic value while with 'gauss2', the mutational effect replaces the current allelic value." << std::endl;
        abort();
    }

    if (T3_mutationType_effectSizeInfo < 0.0)
    {
        std::cout << "For option 'T3_mutationType', received type '"<< type << "' followed by the value "<<T3_mutationType_effectSizeInfo<<". This value cannot be negative." << std::endl;
        abort();
    }
}


void SpeciesSpecificParameters::readT3_MutationRate(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option 'T3_MutationRate', the std::string that is read is: " << input.print() << std::endl;
#endif

    if (input.PeakNextElementString() == "NA")
    {
        if (this->Gmap.T3_nbLoci != 0)
        {
            std::cout << "In option '--T3_MutationRate', for species "<<this->speciesName<<", received 'NA' but there are T3 loci as indiciated in --L (--Loci) (this->Gmap.T3_nbLoci = "<<this->Gmap.T3_nbLoci<<").\n";
            abort();   
        }
        input.skipElement();
    } else
    {
        std::string Mode = input.GetNextElementString();
        
        if (Mode.compare("A")==0)
        {
            double currentSum = 0;
            int T3_char_index = 0;
            while (input.IsThereMoreToRead())
            {
                if (input.isNextRKeyword())
                {
                    double inputValue = input.GetNextElementDouble();
                    int nbRepeats     = input.GetNextElementInt();
                    T3_char_index += nbRepeats;
                    for (int repeat = 0 ; repeat < nbRepeats ; ++repeat)
                    {
                        currentSum += inputValue;
                        this->T3_MutationRate.push_back(currentSum);
                    }

                    // security
                    if (inputValue < 0)
                    {
                        std::cout << "In 'T3_MutationRate' of Mode 'A' the value for the " << T3_char_index << "th T3_char_index is negative (is " << inputValue << ")" << std::endl;
                        abort();
                    }
                    if (inputValue > 0.5)
                    {
                        std::cout << "In 'T3_MutationRate' of Mode 'A' the value for the " << T3_char_index << "th T3_char_index is greater than 0.5 (is " << inputValue << "). Is it really what you want to do?" << std::endl;
                        abort();
                    }
                } else
                {
                    double inputValue = input.GetNextElementDouble();
                    currentSum += inputValue;
                    this->T3_MutationRate.push_back(currentSum);
                    ++T3_char_index;

                    // security
                    if (inputValue < 0)
                    {
                        std::cout << "In 'T3_MutationRate' of Mode 'A' the value for the " << T3_char_index << "th T3_char_index is negative (is " << this->T3_MutationRate.back() << ")" << std::endl;
                        abort();
                    }
                    if (inputValue > 0.5)
                    {
                        std::cout << "In 'T3_MutationRate' of Mode 'A' the value for the " << T3_char_index << "th T3_char_index is greater than 0.5 (is " << this->T3_MutationRate.back() << "). Is it really what you want to do?" << std::endl;
                        abort();
                    }
                }
            }
        } else if (Mode.compare("unif")==0)
        {
            this->T3_MutationRate.push_back(input.GetNextElementDouble() * this->Gmap.T3_nbLoci);
            if (this->T3_MutationRate.back() < 0)
            {
                std::cout << "T3_MutationRate of Mode 'unif' is lower than zero" << std::endl;
                abort();
            }
        } else
        {
            std::cout << "Sorry, for 'T3_MutationRate' only Mode 'unif' and 'A' are recognized so far. Mode received (" << Mode << ") is unknown!" << std::endl;
            abort();
        }
        this->T3_Total_Mutation_rate = this->T3_MutationRate.back();
        assert(this->T3_MutationRate.size() == this->Gmap.T3_nbLoci || this->T3_MutationRate.size() == 1);
        
    }
    
    
}


void SpeciesSpecificParameters::readT3_PhenotypicEffects(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option 'T3_PhenotypicEffects', the std::string that is read is: " << input.print() << std::endl;
#endif

    if (input.PeakNextElementString() == "NA")
    {
        if (this->Gmap.T3_nbLoci != 0)
        {
            std::cout << "In option 'T3_PhenotypicEffects', for species "<<this->speciesName<<", received 'NA' but there are T3 loci as indiciated in --L (--Loci) (this->Gmap.T3_nbLoci = "<<this->Gmap.T3_nbLoci<<").\n";
            abort();   
        }
        input.skipElement();
    } else
    {
        this->T3_PhenoNbDimensions = input.GetNextElementInt();
        if (this->T3_PhenoNbDimensions < 1)
        {
            std::cout << "In option 'T3_PhenotypicEffects', the first entry should be the number of dimensions of the fitness space. The minimal possible value is 1. Received " << this->T3_PhenoNbDimensions  << "\n";
            abort();
        }
        if (this->T3_PhenoNbDimensions > this->Gmap.T3_nbLoci)
        {
            std::cout << "In option 'T3_PhenotypicEffects', the first entry should be the number of dimensions of the fitness space. The number of dimensions received ("<<this->T3_PhenoNbDimensions <<") is higher than the number of T3 loci ("<<this->Gmap.T3_nbLoci<<"). Is it really what you want? If yes, turn this security check off.\n";
            abort();
        }
        
        for (int i = 0 ; i < this->T3_PhenoNbDimensions; i++)
        {
            Individual::T3_IndPhenotype.push_back(0.0); // This is a static variable that is used by all inds to calculate fitness on T3 trait
        }

        input = InputReader(input, 1); // remove first element to avoid messing up when calling 'GetNextHabitatMarker'

        for (int habitat = 0; habitat <= this->MaxEverHabitat ; habitat++)
        {
            (void) input.GetNextHabitatMarker(habitat);

            std::vector<double> T3_PhenotypicEffects_OneHabitat;
            
            std::string Mode = input.GetNextElementString();

            if (Mode.compare("A") == 0)
            {
                for (int T3_locus = 0 ; T3_locus < this->Gmap.T3_nbLoci ; ++T3_locus)
                {
                    for (int dim = 0 ; dim < this->T3_PhenoNbDimensions ; dim++)
                    {
                        T3_PhenotypicEffects_OneHabitat.push_back(input.GetNextElementDouble());
                    }
                }  
            } else if (Mode.compare("unif") == 0)
            {
                std::vector<double> v;
                for (int dim = 0 ; dim < this->T3_PhenoNbDimensions ; dim++)
                {
                    v.push_back(input.GetNextElementDouble());
                }
                assert(v.size() == this->T3_PhenoNbDimensions );
                for (int T3_locus = 0 ; T3_locus < this->Gmap.T3_nbLoci ; ++T3_locus)
                {
                    for (int dim = 0 ; dim < this->T3_PhenoNbDimensions ; dim++)
                    {
                        T3_PhenotypicEffects_OneHabitat.push_back(v[dim]);
                    }
                }  
            } else
            {
                std::cout << "In option 'T3_PhenotypicEffects', received Mode "<< Mode << ". Sorry only modes 'A' and 'unif' are recognized for the moment. \n";
                abort();
            }
            this->T3_PhenotypicEffects.push_back(T3_PhenotypicEffects_OneHabitat);
        } // end of habitat for loop
    }
 
   
    
}

void SpeciesSpecificParameters::readCentralT1LocusForExtraGeneticInfo(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option '--centralT1LocusForExtraGeneticInfo', the std::string that is read is: " << input.print() << std::endl;
#endif
    this->centralT1LocusForExtraGeneticInfo = input.GetNextElementInt();
      
}


void SpeciesSpecificParameters::readT3_FitnessLandscape(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option '--T3_fit (--T3_fitnessLandscape)', the std::string that is read is: " << input.print() << std::endl;
#endif

    if (input.PeakNextElementString() == "NA")
    {
        if (this->Gmap.T3_nbLoci != 0)
        {
            std::cout << "In option 'T3_FitnessLandscape', for species "<<this->speciesName<<", received 'NA' (default value if nothing is specified) but there are T3 loci as indiciated in --L (--Loci) (this->Gmap.T3_nbLoci = "<<this->Gmap.T3_nbLoci<<").\n";
            abort();   
        }
        input.skipElement();
    } else
    {
        std::string SelectionMode = input.GetNextElementString();
        input = InputReader(input, 1); // remove first element to avoid messing up with GetNextHabitatMarker


        for (int habitat = 0; habitat <= this->MaxEverHabitat ; habitat++)
        {
            (void) input.GetNextHabitatMarker(habitat);
        
            std::string EntryMode = input.GetNextElementString();

            std::vector<float> T3_fitnessLandscapeOptimum_OneHabitat;
            std::vector<float> T3_fitnessLandscapeLinearGradient_OneHabitat;
            std::vector<float> T3_fitnessLandscapeGaussStrength_OneHabitat;

            if (SelectionMode.compare("linear")==0 || SelectionMode.compare("simple")==0)
            {
                this->T3_fitnessLandscapeType='L';
                if (EntryMode.compare("A")==0)
                {
                    for (int dim = 0 ; dim < this->T3_PhenoNbDimensions  ; dim++)
                    {
                        double mean = input.GetNextElementDouble();
                        double gradient = input.GetNextElementDouble();

                        if (gradient < 0.0)
                        {
                            std::cout << "In '--T3_FitnessLandscape', SelectionMode 'linear', EntryMode 'A', for dimension "<<dim<<" (counting zero-based) received a gradient of "<<gradient<<". Gradient must be non-negative.\n";
                            abort();
                        }

                        T3_fitnessLandscapeOptimum_OneHabitat.push_back(mean);
                        T3_fitnessLandscapeLinearGradient_OneHabitat.push_back(gradient);

                    }
                } else if (EntryMode.compare("unif")==0) 
                {
                    double mean = input.GetNextElementDouble();
                    double gradient = input.GetNextElementDouble();
                    
                    if (gradient < 0.0)
                    {
                        std::cout << "In '--T3_FitnessLandscape', SelectionMode 'linear', EntryMode 'unif', received a gradient of "<<gradient<<". Gradient must be non-negative.\n";
                        abort();
                    }

                    for (int dim = 0 ; dim < this->T3_PhenoNbDimensions  ; dim++)
                    {
                        T3_fitnessLandscapeOptimum_OneHabitat.push_back(mean);
                        T3_fitnessLandscapeLinearGradient_OneHabitat.push_back(gradient);
                    }
                } else
                {
                    std::cout << "In option '--T3_fit (--T3_fitnessLandscape)', received EntryMode " << EntryMode << ". The only EntryMode accepted are 'A' and 'unif'. Sorry.\n";
                    abort();       
                }
            } else if (SelectionMode.compare("gaussian")==0 || SelectionMode.compare("gauss")==0)
            {
                this->T3_fitnessLandscapeType='G';
                if (EntryMode.compare("A")==0)
                {
                    for (int dim = 0 ; dim < this->T3_PhenoNbDimensions  ; dim++)
                    {
                        double mean = input.GetNextElementDouble();
                        double strength = input.GetNextElementDouble();

                        if (strength <= 0.0)
                        {
                            std::cout << "In '--T3_FitnessLandscape', SelectionMode 'gauss', EntryMode 'A', for dimension "<<dim<<" (counting zero-based) received a strength of "<<strength<<". strength must be non-negative.\n";
                            abort();
                        }

                        T3_fitnessLandscapeOptimum_OneHabitat.push_back(mean);
                        T3_fitnessLandscapeGaussStrength_OneHabitat.push_back(strength);

                    }
                } else if (EntryMode.compare("unif")==0) 
                {
                    double mean = input.GetNextElementDouble();
                    double strength = input.GetNextElementDouble();
                    
                    if (strength <= 0.0)
                    {
                        std::cout << "In '--T3_FitnessLandscape', SelectionMode 'gauss', EntryMode 'unif', received a strength of "<<strength<<". strength must be non-negative.\n";
                        abort();
                    }

                    assert(this->T3_PhenoNbDimensions > 0);
                    for (int dim = 0 ; dim < this->T3_PhenoNbDimensions ; dim++)
                    {
                        T3_fitnessLandscapeOptimum_OneHabitat.push_back(mean);
                        T3_fitnessLandscapeGaussStrength_OneHabitat.push_back(strength);
                    }
                } else
                {
                    std::cout << "In option '--T3_fit (--T3_fitnessLandscape)', received EntryMode " << EntryMode << ". The only EntryMode accepted are 'A' and 'unif'.\n";
                    abort();       
                }
            } else
            {
                std::cout << "In option '--T3_fit (--T3_fitnessLandscape)', the only SelectionMode accepted are 'linear' and 'gaussian' (or 'gauss'). Mode received = '"<<SelectionMode<< "'.\n";
                abort();
            }
            this->T3_fitnessLandscapeOptimum.push_back(T3_fitnessLandscapeOptimum_OneHabitat);

            if (SelectionMode.compare("gaussian")==0 || SelectionMode.compare("gauss")==0)
            {
                this->T3_fitnessLandscapeGaussStrength.push_back(T3_fitnessLandscapeGaussStrength_OneHabitat);
            } else if (SelectionMode.compare("linear")==0)
            {
                this->T3_fitnessLandscapeLinearGradient.push_back(T3_fitnessLandscapeLinearGradient_OneHabitat);    
            } else if (SelectionMode != "NA")
            {
                std::cout << "Unknown SelectionMode in T3_FitnessLandscape. The error should have been detected earlier in the process. There is therefore an internal error but you might want to check your input parameters anyway.\n";
                abort();
            }
            
        } // end of for habitat
        assert(this->T3_fitnessLandscapeOptimum.size() == this->MaxEverHabitat+1);
        if (SelectionMode.compare("gaussian")==0 || SelectionMode.compare("gauss")==0)
        {
            assert(this->T3_fitnessLandscapeGaussStrength.size() == this->MaxEverHabitat+1);
        } else if (SelectionMode.compare("linear")==0)
        {
            assert(this->T3_fitnessLandscapeLinearGradient.size() == this->MaxEverHabitat+1);
        }       
    }
}



void SpeciesSpecificParameters::readT3_DevelopmentalNoise(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option '--T3_DN (--T3_DevelopmentalNoise)', the std::string that is read is: " << input.print() << std::endl;
#endif

    if (this->Gmap.T3_nbLoci)
    {
        for (int habitat = 0; habitat <= this->MaxEverHabitat ; habitat++)
        {
            (void) input.GetNextHabitatMarker(habitat);
        
            std::string EntryMode = input.GetNextElementString();

            std::vector<double> DN_OneHabitat; // each value is for one dimension of the phenotypic space
            

            if (EntryMode.compare("A")==0)
            {
                for (int dim = 0 ; dim < this->T3_PhenoNbDimensions  ; dim++)
                {
                    double sd = input.GetNextElementDouble();

                    if (sd < 0)
                    {
                        std::cout << "In '--T3_DevelopmentalNoise', mode = A. Received a negative standard deviation. sd = "<< sd << "\n";
                        abort();
                    }

                    DN_OneHabitat.push_back(sd);
                }
            } else if (EntryMode.compare("unif")==0) 
            {
                double sd = input.GetNextElementDouble();
                if (sd < 0)
                {
                    std::cout << "In '--T3_DevelopmentalNoise', mode = unif. Received a negative standard deviation. sd = "<< sd << "\n";
                    abort();
                }
                for (int dim = 0 ; dim < this->T3_PhenoNbDimensions  ; dim++)
                {
                    DN_OneHabitat.push_back(sd);
                }
            } else
            {
                std::cout << "In option '--T3_fit (--T3_fitnessLandscape)', received EntryMode " << EntryMode << ". The only EntryMode accepted are 'A' and 'unif'. Sorry.\n";
                abort();       
            }

            this->T3_DevelopmentalNoiseStandardDeviation.push_back(DN_OneHabitat);
            
        } // end of for habitat

        assert(this->T3_DevelopmentalNoiseStandardDeviation.size() == this->MaxEverHabitat+1);
        
    } else
    {
        this->T3_PhenoNbDimensions = 0;
        input.consideredFullyRead();
    }
}


void SpeciesSpecificParameters::readT4_MutationRate(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option 'T4_MutationRate', the std::string that is read is: " << input.print() << std::endl;
#endif
    std::string mode = input.GetNextElementString();
    if (mode == "unif")
    {
        double constantRate = input.GetNextElementDouble();
        T4_MutationRate.push_back(constantRate); // yes, just one value. It will be taken care of. Not very clean though as a design
        T4_Total_Mutation_rate = constantRate * Gmap.T4_nbLoci;
    } else if (mode == "A")
    {
        double currentSum = 0.0;
        while (input.IsThereMoreToRead())
        {
            if (input.isNextRKeyword())
            {
                double inputValue = input.GetNextElementDouble();
                int nbRepeats     = input.GetNextElementInt();

                if (inputValue < 0.0)
                {
                    std::cout << "For option 'T4_MutationRate', received a negative rate (received " << inputValue << ").\n";
                    abort();
                }

                if (inputValue > 0.5)
                {
                    std::cout << "For option 'T4_MutationRate', received a value greater than 0.5 (received " << inputValue << ").\n";
                    abort();
                }

                for (uint32_t repeat = 0 ; repeat < nbRepeats ; ++repeat)
                {
                    currentSum += inputValue;
                    T4_MutationRate.push_back(currentSum);
                }
            } else
            {
                double inputValue = input.GetNextElementDouble();
                if (inputValue < 0.0)
                {
                    std::cout << "For option 'T4_MutationRate', received a negative rate (received " << inputValue << ").\n";
                    abort();
                }
                currentSum += inputValue;
                T4_MutationRate.push_back(currentSum);
            }
        }

        if (T4_MutationRate.size() != Gmap.T4_nbLoci)
        {
            std::cout << "For option 'T4_MutationRate', received " << T4_MutationRate.size() << " values but it was expected "<< Gmap.T4_nbLoci <<" values.\n";
                    abort();
        } 

        T4_Total_Mutation_rate = T4_MutationRate.back();
    } else
    {
        std::cout << "For option 'T4_MutationRate', received the mode of entry '"<<mode<<"'. Sorry, only modes 'unif' and 'A' are accepted.\n";
        abort();
    }

    //std::cout << "T4_Total_Mutation_rate = " << T4_Total_Mutation_rate << "\n";
}


void SpeciesSpecificParameters::readT8_MutationRate(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option 'T8_MutationRate', the std::string that is read is: " << input.print() << std::endl;
#endif
    std::string mode = input.GetNextElementString();
    if (mode == "unif")
    {
        double constantRate = input.GetNextElementDouble();
        T8_MutationRate.push_back(constantRate); // yes, just one value. It will be taken care of. Not very clean though as a design
        T8_Total_Mutation_rate = constantRate * Gmap.T8_nbLoci;
    } else if (mode == "A")
    {
        double currentSum = 0.0;
        while (input.IsThereMoreToRead())
        {
            if (input.isNextRKeyword())
            {
                double inputValue = input.GetNextElementDouble();
                int nbRepeats     = input.GetNextElementInt();

                if (inputValue < 0.0)
                {
                    std::cout << "For option --T8_mu, received a negative mutation rate (" << inputValue << ").\n";
                    abort();
                }

                if (nbRepeats > Gmap.T8_nbLoci)
                {
                    std::cout << "For option --T8_mu, received more valules than expected. Expected " << Gmap.T8_nbLoci << " values but received at least " << nbRepeats << " values\n";
                    abort();
                }

                for (uint32_t repeat = 0 ; repeat < nbRepeats ; ++repeat)
                {
                    currentSum += inputValue;
                    T8_MutationRate.push_back(currentSum);
                }
            } else
            {
                double inputValue = input.GetNextElementDouble();
                if (inputValue < 0.0)
                {
                    std::cout << "For option --T8_mu, received a negative mutation rate (" << inputValue << ").\n";
                    abort();
                }
                currentSum += inputValue;
                T8_MutationRate.push_back(currentSum);
            }
        }

        if (T8_MutationRate.size() != Gmap.T8_nbLoci)
        {
            std::cout << "For option --T8_mu, received " << T8_MutationRate.size() << " values while it was expecting " << Gmap.T8_nbLoci << " values.\n";
            abort();
        }

        T8_Total_Mutation_rate = T8_MutationRate.back();
    } else
    {
        std::cout << "For option 'T8_MutationRate', received the mode of entry '"<<mode<<"'. Sorry, only modes 'unif' and 'A' are accepted.\n";
        abort();
    }


    //std::cout << "T8_Total_Mutation_rate = " << T8_Total_Mutation_rate << "\n";
}


void SpeciesSpecificParameters::readT7fitnessParameters(InputReader& input)
{
    /*
    T7phenpars.nbDimensions = input.GetNextElementInt();
    if (T7phenpars.nbDimensions < 1)
    {
        std::cout << "For option --T7_fit, received a negative number of dimensions. Received " << T7phenpars.nbDimensions << "\n";
        abort();
    }

    auto type = input.GetNextElementString();
    if (type != "linear" && type != "gaussian")
    {
        std::cout << "For option --T7_fit, expected a type of fitness lansdcape. Types accepted are 'linear' and 'gaussian'. Type received is " << type << "\n";
        abort();
    }
    T7phenpars.fitnessLandscapeType = type == "linear" ?  'L' : 'G';

    T7phenpars.fitnessLandscapeOptimum = input.GetNextElementDouble();
    if (T7phenpars.fitnessLandscapeOptimum < 0.0)
    {
        std::cout << "For option --T7_fit, expected a negative fitnessLandscapeOptimum\n";
        abort();
    }

    T7phenpars.fitnessLandscapeLinearGradient = input.GetNextElementDouble();
    if (T7phenpars.fitnessLandscapeLinearGradient < 0.0)
    {
        std::cout << "For option --T7_fit, expected a negative fitnessLandscapeLinearGradients\n";
        abort();
    }

    T7phenpars.fitnessLandscapeGaussStrength = input.GetNextElementDouble();
    if (T7phenpars.fitnessLandscapeLinearGradient < 0.0)
    {
        std::cout << "For option --T7_fit, expected a negative fitnessLandscapeGaussStrength\n";
        abort();
    }



    while (input.IsThereMoreToRead())
    {
        auto x = input.GetNextElementInt();
        if (x <= 0 || x > T7devpars.maxAge)
        {
            std::cout << "Invalid age in T7_fit.\n";
            abort();
        }
        T7phenpars.agesAtwhichPhenotypeIsSampled.push_back(x);
    }
    */
}
    



void SpeciesSpecificParameters::readT7developmentParameters(InputReader& input)
{
    T7devpars.maxAge = input.GetNextElementInt();
    if (T7devpars.maxAge <= 0)
    {
        std::cout << "For option '--T7developmentParameters', first argument (maxAge) is not great than zero. Value received ("<<T7devpars.maxAge<<")\n";
        abort();
    }

    T7devpars.maxDeltaT = input.GetNextElementInt();
    if (T7devpars.maxDeltaT <= 0)
    {
        std::cout << "For option '--T7developmentParameters', second argument (maxDeltaT) is not great than zero. Value received ("<<T7devpars.maxDeltaT<<")\n";
        abort();
    }


    T7devpars.basal_transcription_rate = input.GetNextElementDouble();
    if (T7devpars.basal_transcription_rate < 0)
    {
        std::cout << "For option '--T7developmentParameters', third argument (basal_transcription_rate) is not great than zero. Value received ("<<T7devpars.basal_transcription_rate<<")\n";
        abort();
    }

    T7devpars.EPSILON = input.GetNextElementDouble();
    if (T7devpars.EPSILON < 0)
    {
        std::cout << "For option '--T7developmentParameters', fourth argument (EPSILON) is not great than zero. Value received ("<<T7devpars.EPSILON<<")\n";
        abort();
    }

    T7devpars.stochasticDevelopment = input.GetNextElementBool();
    T7devpars.fitnessOverTime = input.GetNextElementBool();


    T7devpars.nbCisSites = input.GetNextElementInt();
    if (T7devpars.nbCisSites < 0)
    {
        std::cout << "For option '--T7developmentParameters', fifth argument (T7devpars.nbCisSites) is not great than zero. Value received ("<<T7devpars.nbCisSites<<")\n";
        abort();
    }


    T7devpars.basic_signal_trans_effect = input.GetNextElementDouble();
    if (T7devpars.basic_signal_trans_effect < 0)
    {
        std::cout << "For option '--T7developmentParameters', sixth argument (T7devpars.basic_signal_trans_effect) is not great than zero. Value received ("<<T7devpars.basic_signal_trans_effect<<")\n";
        abort();
    }

    T7devpars.basic_signal_conc = input.GetNextElementInt();
    if (T7devpars.basic_signal_conc < 0)
    {
        std::cout << "For option '--T7developmentParameters', seventh argument (T7devpars.basic_signal_conc) is not great than zero. Value received ("<<T7devpars.basic_signal_conc<<")\n";
        abort();
    }

    T7devpars.protein_decayRate = input.GetNextElementDouble();
    if (T7devpars.basic_signal_conc < 0)
    {
        std::cout << "For option '--T7developmentParameters', eighth argument (T7devpars.protein_decayRate) is not great than zero. Value received ("<<T7devpars.protein_decayRate<<")\n";
        abort();
    }

    T7devpars.mRNA_decayRate = input.GetNextElementDouble();
    if (T7devpars.mRNA_decayRate < 0)
    {
        std::cout << "For option '--T7developmentParameters', nineth argument (T7devpars.mRNA_decayRate) is not great than zero. Value received ("<<T7devpars.mRNA_decayRate<<")\n";
        abort();
    }


    T7devpars.translationRate = input.GetNextElementDouble();
    if (T7devpars.translationRate < 0)
    {
        std::cout << "For option '--T7developmentParameters', nineth argument (T7devpars.translationRate) is not great than zero. Value received ("<<T7devpars.translationRate<<")\n";
        abort();
    }


}


void SpeciesSpecificParameters::readT7_MutationRate(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option 'T7_MutationRate', the std::string that is read is: " << input.print() << std::endl;
#endif
   
    T7mutpars.duplication = input.GetNextElementDouble();
    if (T7mutpars.duplication < 0)
    {
        std::cout << "For option 'T7_MutationRate', received a negative duplication mutation rate (received " << T7mutpars.duplication << "\n";
        abort();
    }
    T7mutpars.deletion = input.GetNextElementDouble();
    if (T7mutpars.deletion < 0)
    {
        std::cout << "For option 'T7_MutationRate', received a negative deletion mutation rate (received " << T7mutpars.deletion << "\n";
        abort();
    }

    T7mutpars.cisEffectMu = input.GetNextElementDouble();
    if (T7mutpars.cisEffectMu < 0)
    {
        std::cout << "For option 'T7_MutationRate', received a negative value for the cis effect mutation rate (received " << T7mutpars.cisEffectMu << "\n";
        abort();
    }

    T7mutpars.cisEffectSD = input.GetNextElementDouble();
    if (T7mutpars.cisEffectSD < 0)
    {
        std::cout << "For option 'T7_MutationRate', received a negative value for 'cisEffectSD' (received " << T7mutpars.cisEffectSD << "\n";
        abort();
    }

    T7mutpars.transEffectMu = input.GetNextElementDouble();
    if (T7mutpars.transEffectSD < 0)
    {
        std::cout << "For option 'T7_MutationRate', received a negative value for the trans effect mutation rate (received " << T7mutpars.transEffectMu << "\n";
        abort();
    }

    T7mutpars.transEffectSD = input.GetNextElementDouble();
    if (T7mutpars.transEffectSD < 0)
    {
        std::cout << "For option 'T7_MutationRate', received a negative transEffectSD mutation rate (received " << T7mutpars.transEffectSD << "\n";
        abort();
    }

     T7mutpars.changeTargetMu = input.GetNextElementDouble();
    if (T7mutpars.changeTargetMu < 0)
    {
        std::cout << "For option 'T7_MutationRate', received a negative value for thee mutation to change the target (aka. cisSite). (received " << T7mutpars.changeTargetMu << "\n";
        abort();
    }

   
    
    T7mutpars.totalMutationRatePerGene = T7mutpars.duplication + T7mutpars.deletion + T7mutpars.cisEffectMu + T7mutpars.transEffectMu;
    assert(T7mutpars.totalMutationRatePerGene >= 0);
}




void SpeciesSpecificParameters::readT4_simplifyEveryNGenerations(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option 'T4_simplifyEveryNGenerations', the std::string that is read is: " << input.print() << std::endl;
#endif 
    if (input.PeakNextElementString() == "default")
    {
        input.skipElement();
        //std::cout << "SSP->TotalpatchCapacity = " << SSP->TotalpatchCapacity << "\n";
        //std::cout << "SSP->TotalRecombinationRate = " << SSP->TotalRecombinationRate << "\n";

        size_t maxEverTotalpatchCapacity = 0;
        for (auto& elem : maxEverpatchCapacity) maxEverTotalpatchCapacity += elem;
        auto nbNewEdgesPerGenerations = maxEverTotalpatchCapacity * (1 + this->TotalRecombinationRate);
        auto nbBytesPerGeneration = nbNewEdgesPerGenerations * 4 * 5;
        T4_simplifyEveryNGenerations = 1e9 / nbBytesPerGeneration; // This is a very arbitrary choice and will strongly affect the RAM vs CPU time trade-off. has not been thoroughly tested and might be a poor default choice. I could use something like (SSP->TotalpatchCapacity + SSP->TotalpatchCapacity * SSP->TotalRecombinationRate * 100000) * 100000; maybe
        if (T4_simplifyEveryNGenerations < 30) T4_simplifyEveryNGenerations = 30;

    } else
    {
        T4_simplifyEveryNGenerations = input.GetNextElementInt();
    }
    if (T4_simplifyEveryNGenerations <= 0)
    {
        std::cout << "T4_simplifyEveryNGenerations was set to " << T4_simplifyEveryNGenerations << ".\nIt makes no sense to have a value lower than 1.";
        abort();
    }
}

void SpeciesSpecificParameters::readT8_simplifyEveryNGenerations(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option 'T4_simplifyEveryNGenerations', the std::string that is read is: " << input.print() << std::endl;
#endif 
    if (input.PeakNextElementString() == "default")
    {
        input.skipElement();
        //std::cout << "SSP->TotalpatchCapacity = " << SSP->TotalpatchCapacity << "\n";
        //std::cout << "SSP->TotalRecombinationRate = " << SSP->TotalRecombinationRate << "\n";

        size_t maxEverTotalpatchCapacity = 0;
        for (auto& elem : maxEverpatchCapacity) maxEverTotalpatchCapacity += elem;
        
        T8_simplifyEveryNGenerations = 2*maxEverTotalpatchCapacity;

        if (T8_simplifyEveryNGenerations < 30) T8_simplifyEveryNGenerations = 30;

    } else
    {
        T8_simplifyEveryNGenerations = input.GetNextElementInt();
    }
    if (T8_simplifyEveryNGenerations <= 0)
    {
        std::cout << "T8_simplifyEveryNGenerations was set to " << T8_simplifyEveryNGenerations << ".\nIt makes no sense to have a value lower than 1.";
        abort();
    }
}


void SpeciesSpecificParameters::readRecRateOnMismatch(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option 'recRateOnMismatch', the std::string that is read is: " << input.print() << std::endl;
#endif


    if (input.PeakNextElementString() == "NA")
    {
        recRateOnMismatch_bool = false;
        input.skipElement();
    } else
    {
        std::cout << "Sorry, the option 'readRecRateOnMismatch' has been disabled.\n";
        if (this->Gmap.T1_nbLoci != this->Gmap.TotalNbLoci)
        {
            std::cout << "Sorry, recRateOnMismatch can only be used if you only use T1 loci.\n";
            abort();
        }

        recRateOnMismatch_bool = true;
        recRateOnMismatch_halfWindow = input.GetNextElementInt() / 2;
        recRateOnMismatch_factor = input.GetNextElementDouble();
        
        if (recRateOnMismatch_halfWindow*2 > this->Gmap.T1_nbLoci )
        {
            std::cout << "recRateOnMismatch_Window must be smaller or equal to the total number ot T1 loci. Received recRateOnMismatch_Window = " << 2 * recRateOnMismatch_halfWindow << " while there are " << this->Gmap.T1_nbLoci << " T1 loci." << std::endl;
            abort();   
        }
        if (recRateOnMismatch_factor >= 0.0)
        {
            std::cout << "recRateOnMismatch_factor must be negative. Received recRateOnMismatch_factor = " << recRateOnMismatch_factor << std::endl;
            abort();
        }
    }


    
}

void SpeciesSpecificParameters::readRecombinationRate(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option 'RecombinationRate', the std::string that is read is: " << input.print() << std::endl;
#endif

    std::string unit = input.GetNextElementString();
    if (unit != "M" && unit != "cM" && unit != "rate")
    {
        std::cout << "Error found while reading option --r (--RecombinationRate): First entry for a given species is the unit. There are only three units possible, 'M', 'cM' and 'rate'. Unit received is " << unit << "\n";
        abort();
    }


    std::string Mode = input.GetNextElementString();
    
    if (Mode.compare("unif")==0)
    {
        double InputValue = input.GetNextElementDouble();
        //std::cout << "InputValue = " << InputValue << "\n";
        double AddToParam;
        double cumsum = 0.0;
        
        if (InputValue == -1.0)
        {
            AddToParam = -1.0;
        } else
        {
            if (unit.compare("M")==0)
            {
                AddToParam = InputValue;
            } else if (unit.compare("cM")==0)
            {
                AddToParam = InputValue / 100;
            } else if (unit.compare("rate")==0)
            {
                if (InputValue > 0.5)
                {
                    std::cout << "Error found while reading option --r (--RecombinationRate): Received unit 'rate' and value '"<< InputValue <<"'. When you indicate a recombination rate in 'rate' (as opposed to 'cM' or 'M'), you cannot specify a value greater than 0.5 because it makes no sense. Independent chromosome have a recombination rate of 0.5. If you meant to express a recombination distance in centimorgans or morgans, please use 'cM' or 'M' as a unit instead.\n";
                    abort();
                } else if (InputValue == 0.5)
                {
                    AddToParam = -1;
                } else
                    AddToParam = - log(1 - 2 * InputValue)/2;
            } else
            {
                std::cout << "First element of option 'RecombinationRate' should be the unit (either 'M' (Morgan), 'cM' (centiMorgan) or 'c'(crossOver rate)) should be inputed. Value received is " << unit << ".\n";
                abort();
            }
        }

        if (AddToParam > 10)
        {
            std::cout << "In option '--RecombinationRate' received a value which equivalent 'M' is higher than 10 (received " << AddToParam << ". Is it really what you wanted? It is not going to be very efficient. Note you can have independent chromosomes with '-1' if this is what you want.\n";
            abort();
        }

        if (AddToParam == -1.0)
        {
            for (int pos=0 ; pos < (this->Gmap.TotalNbLoci - 1); pos++)
            {
                this->ChromosomeBoundaries.push_back(pos);
            }
            cumsum += 0.5;
            this->RecombinationRate.push_back(cumsum);
        } else
        {
            assert(AddToParam >= 0.0);
            cumsum += AddToParam;
            this->RecombinationRate.push_back(cumsum);
        }
    } else if (Mode.compare("A")==0)
    {
        double cumsum = 0.0;
        int pos=0;
        
        while (input.IsThereMoreToRead())
        {
            int NbRepeats;
            double InputValue;
            double AddToParam;
            if (input.isNextRKeyword())
            {
                InputValue = input.GetNextElementDouble();
                NbRepeats = input.GetNextElementDouble();
                if (NbRepeats < 0)
                {
                    std::cout << "In option '--RecombinationRate', NbRepeats (not directly after indication 'R') is negative\n";
                    abort();
                }
            } else
            {
                InputValue = input.GetNextElementDouble();
                NbRepeats = 1;
            }
            if (unit.compare("M")==0)
            {
                if (InputValue==-1)
                {
                    AddToParam = -1;
                } else if (InputValue >= 0)
                {
                    AddToParam = InputValue;
                } else
                {
                    std::cout << "In option '--RecombinationRate', Unit 'M' received a value which is negative but different from -1. Value received is " << InputValue << ".\n";
                    abort();
                }
            } else if (unit.compare("cM")==0)
            {
                if (InputValue==-1)
                {
                    AddToParam = -1;
                } else if (InputValue >= 0)
                {
                    AddToParam = InputValue / 100;
                } else
                {
                    std::cout << "In option '--RecombinationRate', Unit 'cM' received a value which is negative but different from -1. Value received is " << InputValue << ".\n";
                    abort();
                }
            } else if (unit.compare("rate")==0)
            {
                if (InputValue==-1 || InputValue==0.5)
                {
                    AddToParam = -1;
                } else if (InputValue >= 0 && InputValue < 0.5) 
                {
                    AddToParam = - log(1 - 2 * InputValue)/2;
                } else
                {
                    std::cout << "In option '--RecombinationRate', Unit 'rate' received a value which is either negative but different from -1 or greater than 0.5. Value received is " << InputValue << ".\n";
                    abort();
                }
            } else
            {
                std::cout << "First element of option 'RecombinationRate' should be the unit (either 'M' (Morgan), 'cM' (centiMorgan) or 'c'(crossOver rate)) should be inputed. Value received is " << unit << ".\n";
                abort();
            }
            if (AddToParam > 10)
            {
                std::cout << "In option '--RecombinationRate' received a value which equivalent 'M' is higher than 10 (received " << AddToParam << ". Is it really what you wanted? It is not going to be very efficient!  Note you can have independent chromosomes with '-1' if this is what you want.\n";
                abort();
            }
            for (int rep = 0 ; rep < NbRepeats ; rep++)
            {
                if (AddToParam == -1)
                {
                    cumsum += 0.5;
                    this->ChromosomeBoundaries.push_back(pos);
                    this->RecombinationRate.push_back(cumsum);
                } else if (AddToParam >= 0)
                {
                    cumsum += AddToParam;
                    this->RecombinationRate.push_back(cumsum);
                } else
                {
                    std::cout << "Internal Error: Should not get here in function 'readRecombinationRate'\n";
                    abort();
                }
                pos++;
            }
        }
    } else
    {
        std::cout << "Sorry, for 'RecombinationRate' only Mode 'unif' and 'A' are recognized so far. Mode received (" << Mode << " and unit " << unit << ") is unknown!" << std::endl;
        abort();
    }
    if (this->RecombinationRate.size() == 1)
    {
        this->TotalRecombinationRate = this->RecombinationRate[0] * (this->Gmap.TotalNbLoci - 1);
    } else if (this->RecombinationRate.size() > 1)
    {
        this->TotalRecombinationRate = this->RecombinationRate.back();
    } else
    {
        abort();
    }
        
    assert(this->TotalRecombinationRate >= 0);
    
    if (
        (
            this->RecombinationRate.size() != this->Gmap.TotalNbLoci - 1
            &&
            this->RecombinationRate.size() != 1
        )
        ||
        (this->Gmap.TotalNbLoci != 1 && this->RecombinationRate.size() == 1 && Mode.compare("A")==0)
    )
    {
        std::cout << "In '--RecombinationRate' did not receive the expected number of elements!\n";

        std::cout << "Nb values expected was " << this->Gmap.TotalNbLoci-1 << " (that is the total number of loci minus 1)\n\n";

        std::cout << "\nthis->RecombinationRate.size() = " << this->RecombinationRate.size() << "    this->Gmap.T1_nbLoci = " << this->Gmap.T1_nbLoci << "    this->Gmap.T2_nbLoci = " << this->Gmap.T2_nbLoci  << "    this->Gmap.T3_nbLoci = " << this->Gmap.T3_nbLoci << "    this->Gmap.T4_nbLoci = " << this->Gmap.T4_nbLoci << "    this->Gmap.T56_nbLoci = " << this->Gmap.T56_nbLoci << "    this->Gmap.TotalNbLoci = " << this->Gmap.TotalNbLoci << "\n\n";

        abort();
    }
    assert(this->ChromosomeBoundaries.size() < this->Gmap.TotalNbLoci);

    if (this->Gmap.TotalNbLoci < 2)
    {
        this->TotalRecombinationRate = 0;
    }
    
    
}


/*void SpeciesSpecificParameters::readT1_vcfOutput_sequence(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option 'T1_vcfOutput_sequence', the std::string that is read is: " << input.print() << std::endl;
#endif

    std::string Mode = input.GetNextElementString();
    if (Mode.compare("range") == 0)
    {
        int from = input.GetNextElementInt();
        int to = input.GetNextElementInt();
        if (from > to)
        {
            std::cout << "In option 'T1_vcfOutput_sequence', 'from' is greater than 'to' from = " << from << ", to = " << to << ".\n";
            abort();
        }
        if (from < 0)
        {
            std::cout << "In option 'T1_vcfOutput_sequence', 'from' is lower than zero from = " << from << ", to = " << to << ".\n";
            abort();   
        }
        if (to >= this->Gmap.T1_nbLoci)
        {
            std::cout << "In option 'T1_vcfOutput_sequence', 'to' is greater or equal to the number of T1 loci from = " << from << ", to = " << to << ", number of T1 loci = " <<  this->Gmap.T1_nbLoci << ".\n";
            abort();   
        }
        for (int locus = from ; locus <= to ; locus++)
        {

            T1_locusDescription T1_locus(locus / EIGHT, locus % EIGHT, locus);
            T1_vcfOutput_sequenceList.push_back(T1_locus);
        }

        T1_vcfOutput_sequenceIsRange = true;
    } else if (Mode.compare("A") == 0)
    {
        std::string sub;
        T1_vcfOutput_sequenceIsRange = true; // will be set to false if loci dont follow in a range
        while (input.IsThereMoreToRead())
        {
            int locus = input.GetNextElementInt();
            if (locus < 0)
            {
                std::cout << "In option 'T1_vcfOutput_sequence', 'locus' is lower than zero. locus = " << locus << ".\n";
                abort();   
            }
            if (locus >= this->Gmap.T1_nbLoci)
            {
                std::cout << "In option 'T1_vcfOutput_sequence', 'locus' is equal or greater than the number the T1 loci. locus = " << locus << ". Number of T1 loci = " << this->Gmap.T1_nbLoci << ". Please note that the loci indices here must be zero-based counting. The first locus has index '0'. The last locus has index 'Number of T1 loci - 1'\n";
                abort();   
            }
            T1_locusDescription T1_locus(locus / EIGHT, locus % EIGHT, locus);
            // Test if it is in a range
            if (T1_vcfOutput_sequenceList.size() > 0)
            {
                if (T1_vcfOutput_sequenceList.back().locus !=  T1_locus.locus - 1)
                {
                    T1_vcfOutput_sequenceIsRange = false;
                }
            }
            // Add to T1_vcfOutput_sequenceList
            T1_vcfOutput_sequenceList.push_back(T1_locus);     
        }
    } else
    {
        std::cout << "In option '--T1_vcfOutput_sequence', the only two Modes accepted are 'range' and 'A'. Mode received is '" << Mode << "\n";
        abort();
    }
    std::sort(
        T1_vcfOutput_sequenceList.begin(),
        T1_vcfOutput_sequenceList.end(),
        [](const T1_locusDescription& left, const T1_locusDescription& right)
        {
          return left.locus < right.locus;
        }
    );

    if (T1_vcfOutput_sequenceList.size() > 1)
    {
        assert(T1_vcfOutput_sequenceList[0].locus < T1_vcfOutput_sequenceList[1].locus);
    }

    

}*/

void SpeciesSpecificParameters::readKillOnDemand(InputReader& input)
{
    killOnDemand.readUserInput(input);
}


void SpeciesSpecificParameters::IsThereSelection()
{
#ifdef CALLENTRANCEFUNCTIONS
    std::cout << "In 'IsThereSelection' start of T1_isSelection\n";
#endif 
    T1_FitnessEffects.shrink_to_fit();
    T2_FitnessEffects.shrink_to_fit();
    T56_FitnessEffects.shrink_to_fit();
    T8_FitnessEffects.shrink_to_fit();

    
    if (Gmap.T1_nbLoci)
    {
        T1_isSelection=false;
        assert(this->T1_FitnessEffects.size() == this->MaxEverHabitat + 1);
        for (int Habitat = 0 ; Habitat <= this->MaxEverHabitat ; Habitat++)
        {
            for (int locus = 0 ; locus < this->Gmap.T1_nbLoci ; ++locus)
            {
                if (this->T1_isMultiplicitySelection)
                {
                    
                    assert(this->T1_FitnessEffects[Habitat].size() > locus);
                    assert(this->T1_FitnessEffects[Habitat][locus] >= 0.0);
                    if (this->T1_FitnessEffects[Habitat][locus] != 1.0)
                    {
                        this->T1_isSelection = true;
                        goto T1isSelectionDone;
                    }   
                } else
                {
                    for (int geno = 0;geno < this->ploidy+1;geno++)
                    {
                        assert(this->T1_FitnessEffects[Habitat].size() > (THREE * locus + geno));
                        assert(this->T1_FitnessEffects[Habitat][THREE * locus + geno] >= 0.0);
                        if (this->T1_FitnessEffects[Habitat][THREE * locus + geno] != 1.0)
                        {
                            this->T1_isSelection = true;
                            goto T1isSelectionDone;
                        }
                    }
                } 
            }
        }
        T1isSelectionDone:


        // Local selection
        if (this->T1_isSelection && this->MaxEverHabitat > 1)
        {
            for (int locus = 0 ; locus < this->Gmap.T1_nbLoci ; ++locus)
            {
                for (int Habitat = 0 ; Habitat <= this->MaxEverHabitat ; Habitat++)
                {
                    if (this->T1_isMultiplicitySelection)
                    {
                        if (this->T1_FitnessEffects[Habitat][locus] != this->T1_FitnessEffects[0][locus])
                        {
                            this->T1_isLocalSelection = true;
                            goto T1isLocalSelectionDone;
                        }
                    } else
                    {
                        for (int geno = 0; geno < this->ploidy+1;geno++)
                        {
                            if (this->T1_FitnessEffects[Habitat][THREE * locus + geno] != this->T1_FitnessEffects[0][THREE * locus + geno])
                            {
                                this->T1_isLocalSelection = true;
                                goto T1isLocalSelectionDone;
                            }
                        }
                    }
                }
            }
        }
        
    } else
    {
        this->T1_isMultiplicitySelection = false;
        this->T1_isSelection = false;
        this->T1_isEpistasis = false;
        this->T1_isLocalSelection = false;
    }
    T1isLocalSelectionDone:

    if (!this->T1_isSelection)
        this->T1_isMultiplicitySelection = false;



#ifdef CALLENTRANCEFUNCTIONS
    std::cout << "In 'IsThereSelection' start of T1_isEpistasis\n";
#endif
    if (this->Gmap.T1_nbLoci)
    {
        this->T1_isEpistasis = false;
        //std::cout << "this->Gmap.T1_nbChars = " << this->Gmap.T1_nbChars << "\n";
        int howManyWarningGiven = 0;
        //std::cout << "this->T1_Epistasis_LociIndices.size() = " << this->T1_Epistasis_LociIndices.size() << "\n";
        if (this->T1_Epistasis_LociIndices.size() != 0)
        {    
            //std::cout << "this->MaxEverHabitat = " << this->MaxEverHabitat << "\n";
            //std::cout << "this->T1_Epistasis_FitnessEffects.size() = " << this->T1_Epistasis_FitnessEffects.size() << "\n";
            assert(this->T1_Epistasis_FitnessEffects.size() == this->MaxEverHabitat + 1);
            for (int Habitat = 0 ; Habitat <= this->MaxEverHabitat ; Habitat++)
            {
                //std::cout << "this->T1_Epistasis_FitnessEffects[Habitat].size() = " << this->T1_Epistasis_FitnessEffects[Habitat].size() << "\n";
                for (int groupOfLociIndex = 0 ; groupOfLociIndex < this->T1_Epistasis_FitnessEffects[Habitat].size() ; groupOfLociIndex++)
                {
                    //std::cout << "groupOfLociIndex = " << groupOfLociIndex << "\n";
                    assert(this->T1_Epistasis_LociIndices[Habitat].size() > groupOfLociIndex);

                    // Send warnings if a locus is under both types of selection
                    for (auto& T1_locus : this->T1_Epistasis_LociIndices[Habitat][groupOfLociIndex])
                    {
                        //std::cout << "T1_locus.locus = " << T1_locus.locus << "\n";

                        
                        if (this->T1_FitnessEffects.size() != 0)
                        {
                            if (this->T1_isMultiplicitySelection)
                            {
                                assert(this->T1_FitnessEffects.size() > Habitat);
                                assert(this->T1_FitnessEffects[Habitat].size() > T1_locus.locus);
                                if (this->T1_FitnessEffects[Habitat][T1_locus.locus] != 1.0)
                                {
                                    if (howManyWarningGiven < 50)
                                    {
                                        std::cout << "\tWARNING: The " << T1_locus.locus << "th locus is under both epistatic selection (--T1_epistasis (--T1_EpistaticFitnessEffects)) and regular selection (--T1_fit (--T1_FitnessEffects)) under habitat " << Habitat << ".\n";
                                        howManyWarningGiven++;
                                    } else
                                    {
                                        if (howManyWarningGiven == 50)
                                        {
                                            std::cout << "You should have already received 50 warnings for having a locus under both types of selection. Further warnings will be silenced.\n";
                                            howManyWarningGiven++;
                                        }
                                    }
                                }
                            } else
                            {
                                bool ShouldIGiveWarning = false;
                                for (int geno=0 ; geno <= this->ploidy ; geno++)
                                {
                                    //std::cout << "geno = " << geno << "\n";
                                    assert(this->T1_FitnessEffects.size() > Habitat);
                                    assert(this->T1_FitnessEffects[Habitat].size() > THREE * T1_locus.locus + geno);
                                    //std::cout << "this->T1_FitnessEffects[Habitat][THREE * T1_locus.locus + geno] = " << this->T1_FitnessEffects[Habitat][THREE * T1_locus.locus + geno] << "\n";
                                    //std::cout << "this->T1_FitnessEffects[Habitat][THREE * T1_locus.locus] = " << this->T1_FitnessEffects[Habitat][THREE * T1_locus.locus] << "\n";
                                    if (this->T1_FitnessEffects[Habitat][THREE * T1_locus.locus + geno] != this->T1_FitnessEffects[Habitat][THREE * T1_locus.locus])
                                    {
                                        ShouldIGiveWarning = true;
                                        break;
                                    }
                                }
                                if (ShouldIGiveWarning)
                                {
                                    if (howManyWarningGiven < 50)
                                    {
                                        std::cout << "\tWARNING: The " << T1_locus.locus << "th locus is under both epistatic selection (--T1_epistasis (--T1_EpistaticFitnessEffects)) and regular selection (--T1_fit (--T1_FitnessEffects)) under habitat " << Habitat << ".\n";
                                        howManyWarningGiven++;
                                    } else
                                    {
                                        if (howManyWarningGiven == 50)
                                        {
                                            std::cout << "You should have already received 50 warnings for having a locus under both types of selection. Further warnings will be silenced.\n";
                                            howManyWarningGiven++;
                                        }
                                    }    
                                }
                            }
                        }
                    }

                    // Check if there is any selection
                    assert(this->T1_Epistasis_FitnessEffects.size() > Habitat);
                    assert(this->T1_Epistasis_FitnessEffects[Habitat].size() > groupOfLociIndex);
                    for (auto& fit : this->T1_Epistasis_FitnessEffects[Habitat][groupOfLociIndex])
                    {
                        if (fit != this->T1_Epistasis_FitnessEffects[Habitat][groupOfLociIndex][0])
                        {
                            this->T1_isEpistasis = true;
                        }
                    }
                } // end of groupOfLoci
            }
        } else
        {
            assert(this->T1_Epistasis_FitnessEffects.size() == 0);
        }
    } else
    {
        this->T1_isMultiplicitySelection = false;
        this->T1_isSelection = false;
        this->T1_isEpistasis = false;
        this->T1_isLocalSelection = false;
    }

#ifdef CALLENTRANCEFUNCTIONS
    std::cout << "In 'IsThereSelection' start of T8_isSelection\n";
#endif     

    if (Gmap.T8_nbLoci)
    {
        T8_isSelection=false;
        
        for (int locus = 0 ; locus < this->Gmap.T8_nbLoci ; ++locus)
        {    
            assert(this->T8_FitnessEffects.size() > locus);
            assert(this->T8_FitnessEffects[locus] > 0.0);
            if (this->T8_FitnessEffects[locus] != 1.0)
            {
                this->T8_isSelection = true;
                goto T8isSelectionDone;
            }   
        }  
    } else
    {
        this->T8_isSelection = false;
    }
    T8isSelectionDone:
  

#ifdef CALLENTRANCEFUNCTIONS
    std::cout << "In 'IsThereSelection' start of T2_isSelection\n";
#endif     

    if (this->Gmap.T2_nbLoci)
    {
        this->T2_isSelection = false;
        assert(T2_FitnessEffects.size() == this->MaxEverHabitat + 1);
        for (int Habitat = 0 ; Habitat <= this->MaxEverHabitat ; Habitat++)
        {
            for (int T2_char_index=0;T2_char_index < this->Gmap.T2_nbLoci;T2_char_index++)
            {
                assert(this->T2_FitnessEffects[Habitat].size() > T2_char_index);
                if (this->T2_FitnessEffects[Habitat][T2_char_index] != 1.0)
                {
                    this->T2_isSelection = true;
                    goto T2isSelectionDone;
                }
            }
        }
        T2isSelectionDone:


        // Local selection
        if (this->T2_isSelection && this->MaxEverHabitat > 1)
        {
            for (int locus = 0 ; locus < this->Gmap.T2_nbLoci ; ++locus)
            {
                for (int Habitat = 0 ; Habitat <= this->MaxEverHabitat ; Habitat++)
                {
                    if (this->T2_FitnessEffects[Habitat][locus] != this->T2_FitnessEffects[0][locus])
                    {
                        this->T2_isLocalSelection = true;
                        goto T2isLocalSelectionDone;
                    }
                    
                }
            }
        }
    } else
    {
        this->T2_isSelection = false;
    }
    T2isLocalSelectionDone:


#ifdef CALLENTRANCEFUNCTIONS
    std::cout << "In 'IsThereSelection' start of T3_isSelection\n";
#endif     

    
    
    if (this->Gmap.T3_nbLoci)
    {
        this->T3_isSelection=false;
        if (this->T3_fitnessLandscapeType!='L' && this->T3_fitnessLandscapeType!='G')
        {
            std::cout << "Internal error: function 'IsThereSelection' only knows how to deal with T3_fitnessLandscapeType=='L' and T3_fitnessLandscapeType=='G' but it appears that T3_fitnessLandscapeType=="<<this->T3_fitnessLandscapeType<<"\n";
            abort();
        }
        if (this->T3_PhenoNbDimensions < 1)
        {
            std::cout << "Internal error: In 'IsThereSelection' the number of dimensions for T3 is lower than 1. It is "<<this->T3_PhenoNbDimensions<<".\n";
            abort();   
        }

        if (this->T3_fitnessLandscapeType=='L')
        {
            assert(this->T3_fitnessLandscapeLinearGradient.size() == this->T3_fitnessLandscapeOptimum.size());
            assert(this->T3_fitnessLandscapeLinearGradient.size() == this->MaxEverHabitat + 1);

            for (int Habitat = 0 ; Habitat <= this->MaxEverHabitat ; Habitat++)
            {
                assert(this->T3_fitnessLandscapeLinearGradient[Habitat].size() == this->T3_fitnessLandscapeOptimum[Habitat].size());
                assert(this->T3_fitnessLandscapeLinearGradient[Habitat].size() == this->T3_PhenoNbDimensions);
                for (int dim = 0 ; dim < this->T3_PhenoNbDimensions ; dim++)
                {
                    if (this->T3_fitnessLandscapeLinearGradient[Habitat][dim] != 0.0)
                    {
                        this->T3_isSelection = true;
                        goto T3isSelectionDone;
                    }
                }
            }
        } else if (this->T3_fitnessLandscapeType=='G')
        {
            for (auto& elem1 : this->T3_fitnessLandscapeGaussStrength)
            {
                for (auto& elem2 : elem1)
                {
                    assert(elem2 > 0.0);
                }
            }      

            // There is necessarily some selection unless this->T3_fitnessLandscapeGaussian is infinity which is impossible.
            this->T3_isSelection = true;
        } else
        {
            std::cout << "Internal error in funciton IsThereSelection\n";
            abort();
        }
        T3isSelectionDone:


        if (T3_isSelection) T3_isLocalSelection = true; // This could really be improved. A priori I don't think it will have any consequence for the moment though
    } else
    {
        this->T3_isSelection = false;
    }
    
#ifdef CALLENTRANCEFUNCTIONS
    std::cout << "In 'IsThereSelection' end of T3_isSelection\n";
#endif
    


#ifdef CALLENTRANCEFUNCTIONS
    std::cout << "In 'IsThereSelection' start of T56_isSelection\n";
#endif        
    
    //std::cout << "T5sel_nbLoci = " << T5sel_nbLoci << "\n";
    if (Gmap.T56sel_nbLoci)
    {
        this->T56_isSelection = true;
    } else
    {
        this->T56_isSelection = false;
        this->T56_isMultiplicitySelection = false;
        this->T56_isLocalSelection = false;
    }
        

    //std::cout << "this->T56_isSelection = " << this->T56_isSelection << "\n";

    if (this->T56_isSelection)
    {
        assert(this->T56_FitnessEffects.size() == this->MaxEverHabitat + 1);
        for (int locus = 0 ; locus < this->Gmap.T5sel_nbLoci ; ++locus)
        {
            for (int Habitat = 0 ; Habitat <= this->MaxEverHabitat ; Habitat++)
            {
                if (this->T56_isMultiplicitySelection)
                {
                    
                    assert(this->T56_FitnessEffects[Habitat].size() > locus);
                    assert(this->T56_FitnessEffects[Habitat][locus] >= 0.0);
                    if (this->T56_FitnessEffects[Habitat][locus] != this->T56_FitnessEffects[0][locus])
                    {
                        this->T56_isLocalSelection = true;
                        // no need for goto as I actually want the assert statements to run by security
                    }
                } else
                {
                    for (int geno = 0; geno < 2; geno++) // only values for het and double mutant
                    {
                        assert(this->T56_FitnessEffects[Habitat].size() > (2 * locus + geno));
                        assert(this->T56_FitnessEffects[Habitat][2 * locus + geno] >= 0.0);
                        if (this->T56_FitnessEffects[Habitat][2 * locus + geno] != this->T56_FitnessEffects[0][2 * locus + geno])
                        {
                            this->T56_isLocalSelection = true;
                            // no need for goto as I actually want the assert statements to run by security
                        }
                    }
                } 
            }
        }
    }

    if (!this->T56_isSelection)
    {
        this->T56_isMultiplicitySelection = false;
        this->T56_isLocalSelection = false;
    }

#ifdef CALLENTRANCEFUNCTIONS
    std::cout << "In 'IsThereSelection' end of T56_isSelection\n";
#endif


   

    // if fecundity is different from -1, then dispersal rate necessarily depends on fitness
    if (this->fecundityForFitnessOfOne != -1.0)
    {
        assert(this->DispWeightByFitness);
    }

    // Dispersal probability cannot be a function of fitness if there is no selection
    if (!this->T1_isSelection && !this->T1_isEpistasis && !this->T2_isSelection && !this->T3_isSelection  && !this->T56_isSelection  && !this->T8_isSelection )
    {
        isAnySelection = false;
        this->selectionOn = 0; // This is to avoid some bullshit in lifeCycle.pp

        if (this->fecundityForFitnessOfOne == -1.0)
            this->DispWeightByFitness = false;
    } else
    {
        isAnySelection = true;
    }


    if (this->MaxEverHabitat == 1)
    {
        this->T1_isLocalSelection = false;
        this->T2_isLocalSelection = false;
        this->T3_isLocalSelection = false;
        this->T56_isLocalSelection = false;
    }

    if (additiveEffectAmongLoci && (Gmap.T2_nbLoci || T1_isMultiplicitySelection || T56_isMultiplicitySelection))
    {
        std::cout << "For option --additiveEffectAmongLoci, received 'true' (selective effects are additive among loci, including loci of different types). However, you used the multfit (multiplicative fitness) assumption, either for T1, T2 or T5 loci. This is not compatible.\n";
        abort();
    }
}

/*void SpeciesSpecificParameters::ClearT1_Initial_AlleleFreqs()
{
    for ( int patch_index = 0 ; patch_index < this->T1_Initial_AlleleFreqs.size() ; ++patch_index )
    {
        this->T1_Initial_AlleleFreqs[patch_index].clear();
    }
    this->T1_Initial_AlleleFreqs.clear();
}*/


/*void SpeciesSpecificParameters::setInitialT1_AlleleFreqTo(const int uniqueFreq)
{
    std::vector<double> line(this->Gmap.T1_nbChars);
    for (int T1_char_index=0;T1_char_index<this->Gmap.T1_nbChars;T1_char_index++)
    {
        line[T1_char_index] = uniqueFreq;
    }
    for (int patch_index=0;patch_index<GP->PatchNumber;++patch_index)
    {
        this->T1_Initial_AlleleFreqs.push_back(line);
    }
}*/




std::vector<int> SpeciesSpecificParameters::UpdateParameters(int generation_index)
{
#ifdef DEBUG
    std::cout << "Enters in SpeciesSpecificParameters::UpdateParameters\n";
#endif

    /*
    std::cout << "patchSize when enterin SpeciesSpecificParameters::UpdateParameters: ";
    for (auto& e : patchSize) std::cout << e << " ";
    std::cout << "\n";
    */    

    // Change patchCapacity Do not change patchSize yet as it will be used later
    assert(generation_index < this->__patchCapacity.size());
    this->patchCapacity = this->__patchCapacity[generation_index];
    this->patchCapacity.shrink_to_fit();
    assert(this->patchCapacity.size() == GP->PatchNumber);

    // Change Patch size
    assert(this->patchSize.size()>0);
    this->patchSize.resize(GP->PatchNumber, 0); // complete with zeros
    this->patchSize.shrink_to_fit();
    assert(this->patchSize.size() == GP->PatchNumber);
    std::vector<int> previousPatchSizes = this->patchSize;
    for (int patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
    {
        if (this->patchSize[patch_index] > this->patchCapacity[patch_index]) this->patchSize[patch_index] = this->patchCapacity[patch_index];
    }


    // Set patch size to carrying capacity if it is not allowed to differ from it. Also change TotalpatchSize
    // Also change TotalpatchCapacity
    TotalpatchSize = 0;
    if (this->fecundityForFitnessOfOne == -1.0)
    {
        for (int patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
        {
            this->patchSize[patch_index] = this->patchCapacity[patch_index];
            TotalpatchSize += this->patchCapacity[patch_index];
        }
    } else
    {
        this->TotalpatchCapacity = 0;
        for ( auto& OnepatchCapacity : this->patchCapacity )
        {
            this->TotalpatchCapacity += OnepatchCapacity;
        }
    }

    
    if (T5_freqThresholdWasSetToDefault)
        this->resetT56_freqThresholToDefault();

    // change Habitat
    this->Habitats = this->__Habitats[generation_index];
    this->MaxHabitat = this->__MaxHabitat[generation_index];
    assert(this->Habitats.size() == GP->PatchNumber);

    // Change growth model
    growthK = __growthK[generation_index];
    assert(growthK.size() == GP->PatchNumber);

     // Change DispMat
    this->dispersalData.nextGenerationPatchSizes = SSP->patchCapacity;
    assert(this->dispersalData.__forwardMigration.size() > generation_index);
    assert(this->dispersalData.__forwardMigrationIndex.size() > generation_index);
    this->dispersalData.forwardMigration = this->dispersalData.__forwardMigration[generation_index];
    this->dispersalData.forwardMigrationIndex = this->dispersalData.__forwardMigrationIndex[generation_index];
    this->dispersalData.BackwardMigration.resize(GP->PatchNumber);
    this->dispersalData.BackwardMigrationIndex.resize(GP->PatchNumber);

    std::vector<unsigned> howManySendToMe(GP->PatchNumber, 0); //howManySendToMe[patch_to] returns number of patch sending migrants to "patch_to"
    assert(howManySendToMe.size() == GP->PatchNumber);
    assert(dispersalData.forwardMigrationIndex.size() == GP->PatchNumber);
    assert(dispersalData.forwardMigrationIndex.size() == dispersalData.forwardMigration.size());
    for (int patch_from = 0 ; patch_from < GP->PatchNumber ; ++patch_from)
    {
        for (auto& patch_to : dispersalData.forwardMigrationIndex[patch_from])
        {
            assert(patch_to < GP->PatchNumber);
            howManySendToMe[patch_to]++;
        }
    }
    for (int patch_to = 0 ; patch_to < GP->PatchNumber ; ++patch_to)
    {
        //std::cout << "howManySendToMe["<<patch_to<<"] = " << howManySendToMe[patch_to] << "\n";
        this->dispersalData.BackwardMigration[patch_to].resize(howManySendToMe[patch_to]);
        this->dispersalData.BackwardMigrationIndex[patch_to].resize(howManySendToMe[patch_to]);
    }
    
    (void) this->dispersalData.setOriginalBackwardMigrationIfNeeded();

    assert(this->dispersalData.forwardMigration.size() == GP->PatchNumber);
    /*
    std::cout << "patchSize when exiting SpeciesSpecificParameters::UpdateParameters: ";
    for (auto& e : patchSize) std::cout << e << " ";
    std::cout << "\n";    
    */
    return previousPatchSizes; // This return value will be used to know what individuals to duplicate
}

void SpeciesSpecificParameters::readStochasticGrowth(InputReader& input)
{
    #ifdef DEBUG
    std::cout << "Enters in SpeciesSpecificParameters::readStochasticGrowth\n";
    #endif  

    if (input.PeakNextElementString() == "default")
    {
        input.skipElement();
        stochasticGrowth = true;
    } else
    {
        stochasticGrowth = input.GetNextElementBool();
    }
}

void SpeciesSpecificParameters::readadditiveEffectAmongLoci(InputReader& input)
{
#ifdef DEBUG
    std::cout << "Enters in SpeciesSpecificParameters::readadditiveEffectAmongLoci\n";
#endif

    if (input.PeakNextElementString() == "default")
    {
        input.skipElement();
        additiveEffectAmongLoci = false;
    } else
    {
        additiveEffectAmongLoci = input.GetNextElementBool();
    }

    if (additiveEffectAmongLoci)
    {
        std::cout << "For option --additiveEffectAmongLoci, received 'true' (or equivalent). In the current version effects among loci can only be multiplicative. Sorry!\n";
        abort();
    }
}

void SpeciesSpecificParameters::readReadPopFromBinary(InputReader& input)
{
#ifdef DEBUG
    std::cout << "Enters in SpeciesSpecificParameters::readReadPopFromBinary\n";
#endif    

    std::string yesNo = input.GetNextElementString();
    if (yesNo == "yes" || yesNo == "y" || yesNo == "Y" || yesNo == "YES" || yesNo == "true" || yesNo == "1")
    {
        readPopFromBinary = true;
        
        // read file path
        readPopFromBinaryPath = input.GetNextElementString();

        // test that the file exists
        std::ifstream f(readPopFromBinaryPath);
        if (!f.good())
        {
            std::cout << "For species " << speciesName << ", in function 'SpeciesSpecificParameters::readReadPopFromBinary', the file " << readPopFromBinaryPath << " was not found. Note that paths to read binary populations from are NOT relative to the general path (the general path only applies to outputs).\n";
            abort();
        }        
    } else if (yesNo == "no" || yesNo == "n" || yesNo == "N" || yesNo == "NO" || yesNo == "false" || yesNo == "0")
    {
        readPopFromBinary = false;
        readPopFromBinaryPath = "";
    } else
    {
        std::cout << "For option '--readReadPopFromBinary', the first element received is " << yesNo << ". Expected 'yes' (or 'y', 'Y', 'true', '1', etc...), 'no' (or 'n', 'NO', '0', 'false', etc...) as to whether the population must be initialized based on a binary file.\n";
        abort();
    }


    
}


void SpeciesSpecificParameters::readT8_propagationMethod(InputReader& input)
{
    if (input.PeakNextElementString() == "default")
    {
        input.skipElement();
        T8_propagationMethod = 1;
        if (SSP->RecombinationRate.back() == 0)
            T8_WhenToSortData = '0';
        else
            T8_WhenToSortData = 'P';
    } else
    {
        T8_propagationMethod = input.GetNextElementInt();

        if (T8_propagationMethod != 1 && T8_propagationMethod != 2)
        {
            std::cout << "For option --T8_propagationMethod, two elements are expected the method of propagation and how mutattions should be sorted. Received value '" << T8_propagationMethod << "'. Sorry, only '1' and '2' are valid methods of propagation.\n";
            abort();
        }

        auto s = input.GetNextElementString();
        if (s == "endOfPropagation" || s == "p")
        {
            T8_WhenToSortData = 'p';
        } else if (s == "duringPropagation" || s == "P" )
        {
            T8_WhenToSortData = 'P';
        } else if (s == "notSorted" || s == "0")
        {
            T8_WhenToSortData = '0';
        } else
        {
            std::cout << "For option --T8_propagationMethod, two elements are expected the method of propagation and how mutations should be sorted. Received value '" << s << "'. Sorry, only 'endOfProbagation', 'duringProbagation' and 'notSorted' valid methods of sorting for T8 loci.\n";
            abort();
        }
    }


    /*
    if (T8_propagationMethod == 2 && T8_WhenToSortData == 'P')
    {
        std::cout << "For option --T8_propagationMethod, two elements are expected the method of propagation and whether ohw values must be sorted. Sorry, with method of propagation 2, sorting cannot happen during propagation.\n";
        abort();
    }
    */
        
    std::cout << "T8_propagationMethod = " << static_cast<unsigned>(T8_propagationMethod) << "\n";
    std::cout << "T8_WhenToSortData = " << T8_WhenToSortData << "\n";
}

void SpeciesSpecificParameters::readSubsetLociForfitnessSubsetLoci_file(InputReader& input)
{
    // Make sure it starts by LociSet!
    if (input.PeakNextElementString() != "LociSet")
    {
        std::cout << "For option --fitnessStats_file, in function SpeciesSpecificParameters::readSubsetLociForfitnessSubsetLoci_file(InputReader& input), the very first word expected (after the time information and species specific marker) was the keyword 'LociSet'. Instead, SimBit received received "<<input.PeakNextElementString()<<".\n";
        abort();
    }

    int nbLociSets = 0;
    std::string whichT;
    while (input.IsThereMoreToRead())
    {
        if (input.PeakNextElementString() == "LociSet")
        {
            input.skipElement();
            std::string whichTstring = input.PeakNextElementString();
            if (whichTstring != "T1" && whichTstring != "T2" && whichTstring != "T3" && whichTstring != "T1epistasis")
            {
                std::cout << "For option --fitnessStats_file, in function SpeciesSpecificParameters::readSubsetLociForfitnessSubsetLoci_file(InputReader& input), expected a description of the type of locus (T1, T2, T3 or T1epistasis) after the keyword 'LociSet' but instead it received "<<whichTstring<<".\n";
                abort();
            }
            nbLociSets++;

             // prepare vectors
            subsetT1LociForfitnessSubsetLoci_file.push_back({});
            subsetT2LociForfitnessSubsetLoci_file.push_back({});
            subsetT3LociForfitnessSubsetLoci_file.push_back({});
            subsetT5LociForfitnessSubsetLoci_file.push_back({});
            subsetT1epistasisLociForfitnessSubsetLoci_file.push_back({});
            assert(subsetT1LociForfitnessSubsetLoci_file.size() == nbLociSets);
            assert(subsetT2LociForfitnessSubsetLoci_file.size() == nbLociSets);
            assert(subsetT3LociForfitnessSubsetLoci_file.size() == nbLociSets);
            assert(subsetT5LociForfitnessSubsetLoci_file.size() == nbLociSets);
            assert(subsetT1epistasisLociForfitnessSubsetLoci_file.size() == nbLociSets);

        } else
        {
            std::string eventualTstring = input.PeakNextElementString();
            if (eventualTstring.at(0) == 'T')
            {
                input.skipElement();
                whichT = eventualTstring;
                
                if (whichT != "T1" && whichT != "T2" && whichT != "T3" && whichT != "T1epistasis")
                {
                    std::cout << "For option --fitnessStats_file, in function SpeciesSpecificParameters::readSubsetLociForfitnessSubsetLoci_file(InputReader& input), received a type of locus descriptor starting with the letter 'T', however it was not followed by '1', '2', '3' or '1epistasis'. The type of locus descriptor received is "<<whichT<<".\n";
                    abort();
                }
            }
            

            int locus = input.GetNextElementInt();
            if (whichT == "T1")
            {
                //std::cout << "locus = " << locus << "\n";
                //std::cout << "T1_nbLoci = " << T1_nbLoci << "\n";
                if (locus < 0 || locus >= Gmap.T1_nbLoci)
                {
                    std::cout << "For option --fitnessStats_file, in function SpeciesSpecificParameters::readSubsetLociForfitnessSubsetLoci_file(InputReader& input), received the locus index " << locus << " for trait of type "<<whichT<<". Only positive number lower than the total number of loci is accepted. As a reminder, the first locus has index 0 and the last locus has index nbLociForThisTrait-1 (T1_nbLoci = "<<Gmap.T1_nbLoci<<", T2_nbLoci = "<<Gmap.T2_nbLoci<<", T3_nbLoci = "<<Gmap.T3_nbLoci<<", T56_nbLoci = "<<Gmap.T56_nbLoci<<")\n";
                    abort();
                }
                if (
                    std::find(
                        subsetT1LociForfitnessSubsetLoci_file.back().begin(),
                        subsetT1LociForfitnessSubsetLoci_file.back().end(),
                        locus
                    ) != subsetT1LociForfitnessSubsetLoci_file.back().end()
                )
                {
                    std::cout << "For option --fitnessStats_file, in function SpeciesSpecificParameters::readSubsetLociForfitnessSubsetLoci_file(InputReader& input), received the locus index " << locus << " for trait type "<<whichT<<" several times within a single LociSet.\n";
                    abort();
                }
                subsetT1LociForfitnessSubsetLoci_file.back().push_back(locus);
            } else if (whichT == "T2")
            {
                if (locus < 0 || locus >= Gmap.T2_nbLoci)
                {
                    std::cout << "For option --fitnessStats_file, in function SpeciesSpecificParameters::readSubsetLociForfitnessSubsetLoci_file(InputReader& input), received the locus index " << locus << " for trait of type "<<whichT<<". Only positive number lower than the total number of loci is accepted. As a reminder, the first locus has index 0 and the last locus has index nbLociForThisTrait-1 (T1_nbLoci = "<<Gmap.T1_nbLoci<<", T2_nbLoci = "<<Gmap.T2_nbLoci<<", T3_nbLoci = "<<Gmap.T3_nbLoci<<", T56_nbLoci = "<<Gmap.T56_nbLoci<<")\n";
                    abort();
                }
                if (
                    std::find(
                        subsetT2LociForfitnessSubsetLoci_file.back().begin(),
                        subsetT2LociForfitnessSubsetLoci_file.back().end(),
                        locus
                    ) != subsetT2LociForfitnessSubsetLoci_file.back().end()
                )
                {
                    std::cout << "For option --fitnessStats_file, in function SpeciesSpecificParameters::readSubsetLociForfitnessSubsetLoci_file(InputReader& input), received the locus index " << locus << " for trait type "<<whichT<<" several times within a single LociSet.\n";
                    abort();
                }
                subsetT2LociForfitnessSubsetLoci_file.back().push_back(locus);
            } else if (whichT == "T3")
            {
                if (locus < 0 || locus >= Gmap.T3_nbLoci)
                {
                    std::cout << "For option --fitnessStats_file, in function SpeciesSpecificParameters::readSubsetLociForfitnessSubsetLoci_file(InputReader& input), received the locus index " << locus << " for trait of type "<<whichT<<". Only positive number lower than the total number of loci is accepted. As a reminder, the first locus has index 0 and the last locus has index nbLociForThisTrait-1 (T1_nbLoci = "<<Gmap.T1_nbLoci<<", T2_nbLoci = "<<Gmap.T2_nbLoci<<", T3_nbLoci = "<<Gmap.T3_nbLoci<<", T56_nbLoci = "<<Gmap.T56_nbLoci<<")\n";
                    abort();
                }  
                if (
                    std::find(
                        subsetT3LociForfitnessSubsetLoci_file.back().begin(),
                        subsetT3LociForfitnessSubsetLoci_file.back().end(),
                        locus
                    ) != subsetT3LociForfitnessSubsetLoci_file.back().end()
                )
                {
                    std::cout << "For option --fitnessStats_file, in function SpeciesSpecificParameters::readSubsetLociForfitnessSubsetLoci_file(InputReader& input), received the locus index " << locus << " for trait type "<<whichT<<" several times within a single LociSet.\n";
                    abort();
                }
                subsetT3LociForfitnessSubsetLoci_file.back().push_back(locus); 
            } else if (whichT == "T1epistasis")
            {
                if (locus < 0 || locus >= Gmap.T1_nbLoci)
                {
                    std::cout << "For option --fitnessStats_file, in function SpeciesSpecificParameters::readSubsetLociForfitnessSubsetLoci_file(InputReader& input), received the locus index " << locus << " for trait of type "<<whichT<<". Only positive number lower than the total number of loci is accepted. As a reminder, the first locus has index 0 and the last locus has index nbLociForThisTrait-1 (T1_nbLoci = "<<Gmap.T1_nbLoci<<", T2_nbLoci = "<<Gmap.T2_nbLoci<<", T3_nbLoci = "<<Gmap.T3_nbLoci<<", T56_nbLoci = "<<Gmap.T56_nbLoci<<")\n";
                    abort();
                }
                if (
                    std::find(
                        subsetT1epistasisLociForfitnessSubsetLoci_file.back().begin(),
                        subsetT1epistasisLociForfitnessSubsetLoci_file.back().end(),
                        locus
                    ) != subsetT1epistasisLociForfitnessSubsetLoci_file.back().end()
                )
                {
                    std::cout << "For option --fitnessStats_file, in function SpeciesSpecificParameters::readSubsetLociForfitnessSubsetLoci_file(InputReader& input), received the locus index " << locus << " for trait type "<<whichT<<" several times within a single LociSet.\n";
                    abort();
                }
                subsetT1epistasisLociForfitnessSubsetLoci_file.back().push_back(locus);
            } else if (whichT == "T5")
            {
                //std::cout << "locus = " << locus << "\n";
                //std::cout << "T1_nbLoci = " << T1_nbLoci << "\n";
                if (locus < 0 || locus >= Gmap.T56_nbLoci)
                {
                    std::cout << "For option --fitnessStats_file, in function SpeciesSpecificParameters::readSubsetLociForfitnessSubsetLoci_file(InputReader& input), received the locus index " << locus << " for trait of type "<<whichT<<". Only positive number lower than the total number of loci is accepted. As a reminder, the first locus has index 0 and the last locus has index nbLociForThisTrait-1 (T1_nbLoci = "<<Gmap.T1_nbLoci<<", T2_nbLoci = "<<Gmap.T2_nbLoci<<", T3_nbLoci = "<<Gmap.T3_nbLoci<<", T56_nbLoci = "<<Gmap.T56_nbLoci<<")\n";
                    abort();
                }
                if (
                    std::find(
                        subsetT5LociForfitnessSubsetLoci_file.back().begin(),
                        subsetT5LociForfitnessSubsetLoci_file.back().end(),
                        locus
                    ) != subsetT5LociForfitnessSubsetLoci_file.back().end()
                )
                {
                    std::cout << "For option --fitnessStats_file, in function SpeciesSpecificParameters::readSubsetLociForfitnessSubsetLoci_file(InputReader& input), received the locus index " << locus << " for trait type "<<whichT<<" several times within a single LociSet.\n";
                    abort();
                }
                subsetT5LociForfitnessSubsetLoci_file.back().push_back(locus);
            } else
            {
                std::cout << "For option --fitnessStats_file, in function SpeciesSpecificParameters::readSubsetLociForfitnessSubsetLoci_file(InputReader& input), received a type of locus descriptor starting with the letter 'T', however it was not followed by '1', '2', '3', '5' or '1epistasis'. The type of locus descriptor received is "<<whichT<<". Note that this error message comes at a bit of an unexpected time and might therefore be due to an internal error. You might want to check your input nevertheless.\n";
                abort();
            }
        }
    }

    // Make sure at least one set has been received
    assert(subsetT1LociForfitnessSubsetLoci_file.size() == nbLociSets);
    assert(subsetT2LociForfitnessSubsetLoci_file.size() == nbLociSets);
    assert(subsetT3LociForfitnessSubsetLoci_file.size() == nbLociSets);
    assert(subsetT5LociForfitnessSubsetLoci_file.size() == nbLociSets);
    assert(subsetT1epistasisLociForfitnessSubsetLoci_file.size() == nbLociSets);
    if (nbLociSets == 0)
    {
        std::cout << "For option --fitnessStats_file, in function SpeciesSpecificParameters::readSubsetLociForfitnessSubsetLoci_file(InputReader& input), received 0 LociSet.\n";
        abort();
    }

    // Make sure there is not empty set and also just sort them for fun!
    for (int lociSetIndex = 0 ; lociSetIndex < nbLociSets; lociSetIndex++)
    {
        if (
            subsetT1LociForfitnessSubsetLoci_file[lociSetIndex].size() == 0 &&
            subsetT2LociForfitnessSubsetLoci_file[lociSetIndex].size() == 0 &&
            subsetT3LociForfitnessSubsetLoci_file[lociSetIndex].size() == 0 &&
            subsetT5LociForfitnessSubsetLoci_file[lociSetIndex].size() == 0 &&
            subsetT1epistasisLociForfitnessSubsetLoci_file[lociSetIndex].size() == 0
        )
        {
            std::cout << "For option --fitnessStats_file, in function SpeciesSpecificParameters::readSubsetLociForfitnessSubsetLoci_file(InputReader& input), the "<<lociSetIndex<<"th LociSet seem to contain no locus of any type! Maybe the keyword 'LociSet' has been repeated twice (or more times) without indicating loci indices in between. Maybe the option argument ends with the keyword 'LociSet'.\n";
            abort();
        }
        std::sort(
            subsetT1LociForfitnessSubsetLoci_file[lociSetIndex].begin(),
            subsetT1LociForfitnessSubsetLoci_file[lociSetIndex].end()
        );
        std::sort(
            subsetT2LociForfitnessSubsetLoci_file[lociSetIndex].begin(),
            subsetT2LociForfitnessSubsetLoci_file[lociSetIndex].end()
        );
        std::sort(
            subsetT3LociForfitnessSubsetLoci_file[lociSetIndex].begin(),
            subsetT3LociForfitnessSubsetLoci_file[lociSetIndex].end()
        );
        std::sort(
            subsetT5LociForfitnessSubsetLoci_file[lociSetIndex].begin(),
            subsetT5LociForfitnessSubsetLoci_file[lociSetIndex].end()
        );
        std::sort(
            subsetT1epistasisLociForfitnessSubsetLoci_file[lociSetIndex].begin(),
            subsetT1epistasisLociForfitnessSubsetLoci_file[lociSetIndex].end()
        );
    }
}




void SpeciesSpecificParameters::readQuickScreenOfOptionL(InputReader& input)
{
    quickScreenAtL_T56_nbLoci = 0;
    while( input.IsThereMoreToRead() )
    {
        auto tmp = input.GetNextLocusInfo();
        auto type = tmp.first;
        auto nbElements = tmp.second;
        
        if (type == 5)
        {
            quickScreenAtL_T56_nbLoci += nbElements;
        } else
        {
            if (
                type != 1
                &&
                type != 2
                &&
                type != 3
                &&
                type != 4
                &&
                type != 7
                &&
                type != 8
                )
            {
                std::cout << "While making a quick screen through the option --L (--Loci), SimBit has found a unexpected type of trait (type "<< type << ").\n";
                abort();
            }
        }
    }
    //std::cout << "quickScreenAtL_T56_nbLoci = " << quickScreenAtL_T56_nbLoci << "\n";
}



void SpeciesSpecificParameters::readkillIndividuals(InputReader& input)
{
    if (input.PeakNextElementString() == "default")
    {
        input.skipElement();
        return;
    }

    if (input.IsThereMoreToRead() && fecundityForFitnessOfOne == -1)
    {
        std::cout << "Option --killIndividuals can only be used if fecundity is different from -1 (SimBit cannot kill arbitrary individuals if the patch size cannot vary from the carrying capacity).\n";
        abort();
    }

    while (input.IsThereMoreToRead())
    {
        // Expecting 'kill'
        {
            auto s = input.GetNextElementString();
            if (s != "kill")
            {
                std::cout << "For option --killIndividuals, expected keyword 'kill' but got " << s << " instead.\n";
                abort();
            }
        }

        // Check if allBut
        bool isAllBut = false;
        {
            auto s = input.PeakNextElementString();
            if (s == "allBut" || s == "allbut" || s == "AllBut" || s == "Allbut")
            {
                isAllBut = true;
                input.skipElement();
            }
        }

        // Read number of individuals
        int nbInds;
        {
            nbInds = input.GetNextElementInt();
            if (nbInds < 0)
            {
                if (isAllBut)
                {
                    std::cout << "For option --killIndividuals, you asked to kill all but ('allBut') " << nbInds << " individuals. Cannot leave a negative number of indivduals.\n";
                    abort();
                } else
                {
                    std::cout << "For option --killIndividuals, you asked to kill " << nbInds << " individuals. Cannot kill a negative number of indivduals.\n";
                    abort();
                }
            }
        }

        // Expecting 'patch'
        {
            auto s = input.GetNextElementString();
            if (s != "patch")
            {
                std::cout << "For option --killIndividuals, expected keyword 'patch' but got " << s << " instead.\n";
                abort();
            }
        }

        // Read patch index
        int patch;
        {
            patch = input.GetNextElementInt();
            if (patch < 0)
            {
                std::cout << "For option --killIndividuals, received a negative patch index. Patch index received is " << patch << ".\n";
                abort();
            }
        }

        // Expecting 'at'
        {
            auto s = input.GetNextElementString();
            if (s != "at")
            {
                std::cout << "For option --killIndividuals, expected keyword 'at' but got " << s << " instead.\n";
                abort();
            }
        }

        // read generations
        std::vector<size_t> generations;
        while (input.IsThereMoreToRead() && input.PeakNextElementString() != "kill")
        {
            auto g = input.GetNextElementInt();
            if (g <= 0 || g > GP->nbGenerations)
            {
                std::cout << "For option --killIndividuals, received generation " << g << ". This generation is either not greater than zero or is greater than the total number of generations (nbGenerations = "<< GP->nbGenerations <<").\n";
                abort();
            }
            generations.push_back(g);
        }




        for (const auto& generation : generations)
        {
            int generation_index = std::upper_bound(GP->__GenerationChange.begin(), GP->__GenerationChange.end(), generation) - GP->__GenerationChange.begin() - 1;

            {
                auto x = GP->__PatchNumber[generation_index];
                if (patch >= x)
                {
                    std::cout << "For option --killIndividuals, received patch " << patch << " at generation " << generation << " but there are only " << x << " patches in the world at that generation. As a ereminder all indices in SimBit are zero based counting\n";
                    abort();
                }
            }

            {
                auto x = __patchCapacity[generation_index][patch];
                if (nbInds > x)
                {
                    if (isAllBut)
                    {
                        std::cout << "WARNING: For option --killIndividuals, received command to kill all but " << nbInds << " individuals in patch " << patch << " at generation " << generation << ". However, at this generation, in this patch, the carrying capacity is only " << x << ".\n";
                    } else
                    {
                        std::cout << "WARNING: For option --killIndividuals, received command to kill " << nbInds << " individuals in patch " << patch << " at generation " << generation << ". However, at this generation, in this patch, the carrying capacity is only " << x << ".\n";
                    }
                }
            }

            killIndividuals_times.push_back(generation);
            killIndividuals_patch.push_back(patch);
            killIndividuals_isAllBut.push_back(isAllBut);
            killIndividuals_numberInds.push_back(nbInds);
        }
        assert(killIndividuals_times.size() == killIndividuals_patch.size());
    }

    assert(killIndividuals_times.size() == killIndividuals_patch.size());
    assert(killIndividuals_times.size() == killIndividuals_isAllBut.size());
    assert(killIndividuals_times.size() == killIndividuals_numberInds.size());

    if (killIndividuals_times.size() >= 2)
    {
        auto order = reverse_sort_indexes(killIndividuals_times);
        reorderNoAssertions(killIndividuals_times, order);
        reorderNoAssertions(killIndividuals_patch, order);
        reorderNoAssertions(killIndividuals_isAllBut, order);
        reorderNoAssertions(killIndividuals_numberInds, order);
    }

    assert(killIndividuals_times.size() == killIndividuals_patch.size());
    assert(killIndividuals_times.size() == killIndividuals_isAllBut.size());
    assert(killIndividuals_times.size() == killIndividuals_numberInds.size());
    
    /*
    for (auto& t : killIndividuals_times)
        std::cout << t << " ";
    std::cout << "\n";
    */
}


void SpeciesSpecificParameters::killIndividualsIfAskedFor()
{
    //std::cout << "GP->CurrentGeneration = " << GP->CurrentGeneration << "\n";
    //std::cout << "before killIndividuals_times.back() = " << killIndividuals_times.back() << "\n";
    while (killIndividuals_times.size() && killIndividuals_times.back() == GP->CurrentGeneration)
    {
        assert(fecundityForFitnessOfOne != -1);

        auto patch = killIndividuals_patch.back();

        assert(patch >= 0 && patch < GP->PatchNumber);
        assert(patchSize.size() > patch);
        
        int newSize = 
            killIndividuals_isAllBut.back() 
            ?
            killIndividuals_numberInds.back() 
            :
            patchSize[patch] - killIndividuals_numberInds.back() < 0
        ;

        assert(newSize <= patchCapacity[patch]);
        if (newSize >= 0 && newSize < patchSize[patch])
        {
            patchSize[patch] = newSize;
        }
        

        killIndividuals_times.pop_back();
        killIndividuals_patch.pop_back();
        killIndividuals_isAllBut.pop_back();
        killIndividuals_numberInds.pop_back();
    }
    //std::cout << "after killIndividuals_times.back() = " << killIndividuals_times.back() << "\n";

    if (killIndividuals_times.size()) assert(killIndividuals_times.back() > GP->CurrentGeneration);
}
