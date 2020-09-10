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


double DispersalData::computeNbOffspringsProducedInPatch(const unsigned patch_from, const double n_t, const double rn_t, const double r)
{
    
    double nbOffs = 0.0;
    if (SSP->patchSize[patch_from] == 0) {return 0;}
    //std::cout << "SSP->growthK["<<patch_from<<"] = " << SSP->growthK[patch_from] << "\n";
    assert(SSP->growthK.size() > patch_from);

    if (GP->nbSpecies == 1)
    {
        if (SSP->growthK[patch_from] == -1.0) // Simple exponential growth
        {
            nbOffs = rn_t;
        }
        else if (SSP->growthK[patch_from] == -2.0) // -2 means logistic with true carrying capacity
        {
            if (r <= 1) // if decline
            {
                nbOffs = rn_t;
            } else
            {
                nbOffs = n_t + (r-1) * n_t * (1 - n_t / SSP->patchCapacity[patch_from]);
            }
        } else
        {
            if (r <= 1) // if decline
            {
                nbOffs = rn_t;
            } else // if growth
            {
                nbOffs = n_t + (r-1) * n_t * (1 - n_t / SSP->growthK[patch_from]);
            }
        }

    } else // if several species
    {
        // Compute the species competition and interaction
        double sumOfAlphasProd_interaction = 0.0;
        double sumOfAlphasProd_competition = 0.0; // Only used if logistic growth. Will compute itself
        for (int speciesIndex = 0 ; speciesIndex < GP->nbSpecies ; ++speciesIndex)
        {
            //// Competition
            // if there are no competition, then just the self should be computed here (all other will have a speciesComputation of 0)
            sumOfAlphasProd_competition += GP->speciesCompetition[SSP->speciesIndex][speciesIndex] * GP->allSpeciesPatchSizePreviousGeneration[patch_from][speciesIndex];


            // Interaction
            auto interaction = GP->speciesInteraction[SSP->speciesIndex][speciesIndex];
            if (interaction.type == '0')
            {
                // nothing to do
            } else if (interaction.type == 'A')
            {
                sumOfAlphasProd_interaction += interaction.magnitude;
            } else if (interaction.type == 'B')
            {
                auto causalSpPS = GP->allSpeciesPatchSizePreviousGeneration[patch_from][speciesIndex];
                sumOfAlphasProd_interaction += interaction.magnitude * causalSpPS;
            } else if (interaction.type == 'C')
            {
                auto recipientSpPS = GP->allSpeciesPatchSizePreviousGeneration[patch_from][SSP->speciesIndex];
                sumOfAlphasProd_interaction += interaction.magnitude * recipientSpPS;
            } else if (interaction.type == 'D') 
            {
                auto causalSpPS = GP->allSpeciesPatchSizePreviousGeneration[patch_from][speciesIndex];
                auto recipientSpPS = GP->allSpeciesPatchSizePreviousGeneration[patch_from][SSP->speciesIndex];
                sumOfAlphasProd_interaction += interaction.magnitude * causalSpPS * recipientSpPS;
            }
        }


        // Competition
        if (SSP->growthK[patch_from] == -1.0) // exponential growth
        {
            nbOffs = rn_t;   // No competition possible when exponential growth
        } else if (SSP->growthK[patch_from] == -2) // -2 means always at true carrying capacity
        {
            if (r <= 1) // if decline
            {
                nbOffs = rn_t;
            } else // if growth
            {
                nbOffs = n_t + (r-1) * n_t * (1 - sumOfAlphasProd_competition / SSP->patchCapacity[patch_from]);
            }
        } else
        {
            if (r <= 1) // if decline
            {
                nbOffs = rn_t;
            } else // if growth
            {
                nbOffs = n_t + (r-1) * n_t * (1 - sumOfAlphasProd_competition / SSP->growthK[patch_from]);
            }
        }
        
        // Adjust for interaction
        //std::cout << "sumOfAlphasProd_interaction = " << sumOfAlphasProd_interaction << "\n";
        nbOffs += sumOfAlphasProd_interaction;
    }

    // Correct if it went below zero
    if (nbOffs < 0.0) nbOffs = 0.0; // if r is very large, then behaviour is chaotic and it can lead to lower than zero values. Species interaction can also lead to lower than zero values

    if (SSP->stochasticGrowth)
    {
        std::poisson_distribution<> d(nbOffs);
        nbOffs = d(GP->rngw.getRNG());
    }

    /*if (patch_from == 1003)
    {
        std::cout << "\n";
        std::cout << "patch_from = " << patch_from << "\n";
        std::cout << "r = " << r << "\n";
        std::cout << "rn_t = " << rn_t << "\n";
        std::cout << "n_t = " << n_t << "\n";
        std::cout << "nbOffs = " << nbOffs << "\n";
    }*/

    assert(nbOffs >= 0);
    assert(!std::isnan(nbOffs));

    return nbOffs;
}

void DispersalData::getMigrationEventsForEachDestinationPatch(std::vector<std::vector<std::vector<double>>>* CumSumFits_p, std::vector<ListMigrationEvents>& migrationEventsForEachDestinationPatch)
{
#ifdef DEBUG
    std::cout << "Enters in DispersalData::getMigrationEventsForEachDestinationPatch\n";
#endif

    
    //////////////////////////////////////
    /// Offspring production per patch ///
    //////////////////////////////////////

    // compute production of offspring of each patch
    migrationEventsForEachDestinationPatch.resize(GP->PatchNumber);
    for (unsigned patch_from = 0 ; patch_from < GP->PatchNumber ; patch_from++)
    {
        /////////////////////////
        /// If noone in patch ///
        /////////////////////////
        if (SSP->patchSize[patch_from] == 0) // Of course, that can onnly be true if fec != -1
        {
            continue;
        }


        //////////////////////////////////
        /// Get mean and total fitness ///
        //////////////////////////////////
        double totalFitness;
        double meanFitness;
        if (SSP->DispWeightByFitness && SSP->isAnySelection)
        {
            assert(CumSumFits_p != nullptr);
            auto& CumSumFits = *CumSumFits_p;
            if (SSP->malesAndFemales)
            {
                assert(SSP->patchSize[patch_from] == CumSumFits[patch_from][0].size() + CumSumFits[patch_from][1].size());
            } else
            {
                assert(SSP->patchSize[patch_from] == CumSumFits[patch_from][0].size());
            }

            /*if (patch_from == 1003)
            {
                std::cout << "CumSumFits[patch_from][0].back() = " << CumSumFits[patch_from][0].back() << "\n";
                for (unsigned ind = 0 ; ind < CumSumFits[patch_from][0].size() ; ++ind)
                {
                    std::cout << "CumSumFits["<<patch_from<<"]["<<ind<<"].back() = " << CumSumFits[patch_from][0].back() << "\n";    
                }
                
                std::cout << "SSP->patchSize[patch_from] = " << SSP->patchSize[patch_from] << "\n";
            }*/

            totalFitness = CumSumFits[patch_from][0].back();
            if (SSP->patchSize[patch_from] > 0)
                meanFitness = totalFitness / SSP->patchSize[patch_from];
            else
                meanFitness = 0;
        } else
        {
            totalFitness = SSP->patchSize[patch_from];
            meanFitness = 1.0;
        }
            

            

        ////////////////////////////
        /// Compute nbOffsprings ///
        ////////////////////////////
        double nbOffs; // nbOffs will be used even if fec == -1 but it will then just be a total fitness then. Hence the usage of a double.
        if (SSP->fecundityForFitnessOfOne != -1)
        {
            double n_t  = SSP->patchSize[patch_from];
            double r, rn_t;
            if (SSP->fecundityDependentOfFitness)
            {
                rn_t = totalFitness * SSP->fecundityForFitnessOfOne;
                r    = meanFitness * SSP->fecundityForFitnessOfOne;
            } else
            {
                rn_t = n_t * SSP->fecundityForFitnessOfOne;
                r    = SSP->fecundityForFitnessOfOne;
            }
            
            
            /*if (patch_from == 1003)
            {
                std::cout << "meanFitness = " << meanFitness << "\n";
                std::cout << "totalFitness = " << totalFitness << "\n";
            }*/

            nbOffs = computeNbOffspringsProducedInPatch(patch_from, n_t, rn_t, r); // returns a double for performance reasons (among other reasons)
        } else
        {
            nbOffs = totalFitness;
            assert(nbOffs == SSP->patchCapacity[patch_from]);
        }
            

        //if (patch_from == 1003) std::cout << "nbOffs of patch "<< patch_from << " = " << nbOffs << ":\n";

        assert(forwardMigration.size() == forwardMigrationIndex.size());
        assert(forwardMigration.size() > patch_from);
        for (unsigned FMpatchToFakeIndex = 0 ; FMpatchToFakeIndex < forwardMigration[patch_from].size() ; ++FMpatchToFakeIndex)
        {
            auto patch_to = forwardMigrationIndex[patch_from][FMpatchToFakeIndex];
            assert(patch_to >= 0 && patch_to < GP->PatchNumber);
            double forwardMigrationRate = (double) forwardMigration[patch_from][FMpatchToFakeIndex];
            assert(forwardMigrationRate >= 0.0 && forwardMigrationRate <= 1.0);

            assert(!std::isnan(forwardMigrationRate));
            //if (patch_from == 1003) std::cout << "\t" << forwardMigrationRate * nbOffs << " of them go to " << patch_to << "\n";
            migrationEventsForEachDestinationPatch[patch_to].addMigrationEvent(forwardMigrationRate * nbOffs, patch_from); // Will only be considered if greater than 0 migrants
        }
    }

    assert(migrationEventsForEachDestinationPatch.size() == GP->PatchNumber);
#ifdef DEBUG
    std::cout << "Exit in DispersalData::getMigrationEventsForEachDestinationPatch\n";
#endif
}

void DispersalData::setOriginalBackwardMigrationIfNeeded()
{
    // memory
    assert(forwardMigration.size() == forwardMigrationIndex.size());
    assert(forwardMigration.size() == GP->PatchNumber);
    BackwardMigration.resize(GP->PatchNumber);
    BackwardMigrationIndex.resize(GP->PatchNumber);



    if (SSP->DispWeightByFitness)
    {
        return;
    } else
    {
        assert(SSP->fecundityForFitnessOfOne == -1);

        std::vector<ListMigrationEvents> migrationEventsForEachDestinationPatch; // migrationEventsForEachDestinationPatch[patch_to]
        
        getMigrationEventsForEachDestinationPatch(nullptr, migrationEventsForEachDestinationPatch);
        computeBackwardMigrationRates_from_migrationEventsForEachDestinationPatch(migrationEventsForEachDestinationPatch);

        /*for (size_t patch_to = 0 ; patch_to < GP->PatchNumber; ++patch_to)
        {
            std::cout << "From patch " << BackwardMigrationIndex[patch_to][0] << " to patch " << patch_to << " = " << BackwardMigration[patch_to][0] << "\n";
        }*/

    }
}


const std::vector<int>& DispersalData::setBackwardMigrationIfNeededAndGetNextGenerationPatchSizes(std::vector<std::vector<std::vector<double>>>& CumSumFits)
{
#ifdef DEBUG
    std::cout << "Enters in DispersalData::SetBackwardMigrationAndGetNextGenerationPatchSizes\n";
#endif

    // Always think about:
    //    SSP->DispWeightByFitness (must be true if fec != 1)
    //    SSP->fecundityForFitnessOfOne
    //    GP->nbSpecies == 1 // Assume if several species, then user wants none-default ecoology
    //    computingOriginalBackwardMigration
    

    ////////////////////////////////////////////////////
    /// No need to recompute backward migration rate ///
    ////////////////////////////////////////////////////

    if (!SSP->DispWeightByFitness)
    {
        // patch sizes must always be at carrying capacity
        // Then backward migration rate has been previously computed
        assert(SSP->fecundityForFitnessOfOne == -1);
        // original backward migration rate doesn't need to be recomputed
        return SSP->patchCapacity; // The next generation patch size is simply the patchCapacity
    }


    //////////////////
    /// Assertions ///
    //////////////////

    assert(this->nextGenerationPatchSizes.size() == GP->PatchNumber);
    assert(this->forwardMigration.size() == GP->PatchNumber);
    assert(this->forwardMigrationIndex.size() == GP->PatchNumber);
    assert(this->BackwardMigration.size() == GP->PatchNumber);
    assert(this->BackwardMigrationIndex.size() == GP->PatchNumber);
    if (SSP->isAnySelection)
    {
        assert(CumSumFits.size() == GP->PatchNumber);
        // just to test the first patch
        if (SSP->malesAndFemales)
        {
            assert(CumSumFits[0].size() == 2);
        } else
        {
            assert(CumSumFits[0].size() == 1);
        }
    }


    ////////////////////////////////////////////////////
    /// Compute offspring production and destination ///
    ////////////////////////////////////////////////////

    std::vector<ListMigrationEvents> migrationEventsForEachDestinationPatch; // migrationEventsForEachDestinationPatch[patch_to]
    
    getMigrationEventsForEachDestinationPatch(&CumSumFits, migrationEventsForEachDestinationPatch);
    computeBackwardMigrationRates_from_migrationEventsForEachDestinationPatch(migrationEventsForEachDestinationPatch); // will set nextGenerationPatchSizes if fec != -1

    
    if (SSP->fecundityForFitnessOfOne != -1.0)
    {
        /*for (size_t patch_index = 0 ; patch_index < nextGenerationPatchSizes.size(); ++patch_index)
        {
            bool isAnyoneComeFromPatch = false;
            for (size_t i = 0 ; i < BackwardMigrationIndex.size(); ++i)
            {
                for (size_t j = 0 ; j < BackwardMigrationIndex[i].size(); ++j)   
                {
                    if (BackwardMigrationIndex[i][j] == patch_index)
                    {
                        assert(BackwardMigration[i][j] >= 0.0 && BackwardMigration[i][j] <= 1.0);
                        isAnyoneComeFromPatch = true;
                        break;
                    }
                }
            } 
            if (nextGenerationPatchSizes[patch_index] == 0)
            {
                assert(!isAnyoneComeFromPatch);
            } else
            {
                assert(isAnyoneComeFromPatch); // This does not always have to be true but it is just for a test.
            }
        }*/

        return nextGenerationPatchSizes;
    } else
    {
        return SSP->patchCapacity;
    }
}

void DispersalData::computeBackwardMigrationRates_from_migrationEventsForEachDestinationPatch(std::vector<ListMigrationEvents>& migrationEventsForEachDestinationPatch)
{
#ifdef DEBUG
    std::cout << "Enters in DispersalData::computeBackwardMigrationRates_from_migrationEventsForEachDestinationPatch\n";
#endif    
    // Loop through each patch_to
    assert(nextGenerationPatchSizes.size() == GP->PatchNumber);
    for (unsigned patch_to = 0 ; patch_to < GP->PatchNumber ; patch_to++)
    {   
        // sum number of migrants to patch_to
        double sumOfIncomingMigrants = migrationEventsForEachDestinationPatch[patch_to].getSumOfIncomingMigrants();

        /*if (patch_to == 0 || (patch_to > 1e3 && patch_to < 1e3+5))
        {
            std::cout << sumOfIncomingMigrants << " incoming migrants in patch " << patch_to << "\n";   
        }*/

        // Set next generation patch size
        if (SSP->fecundityForFitnessOfOne != -1.0)
        {
            if (sumOfIncomingMigrants < SSP->patchCapacity[patch_to])
            {
                long intPart = (long)sumOfIncomingMigrants;
                double fracPart = sumOfIncomingMigrants - intPart;
                if (fracPart != 0.0 && GP->rngw.uniform_real_distribution(1.0) < fracPart)
                {
                    nextGenerationPatchSizes[patch_to] = intPart+1;
                } else
                {
                    nextGenerationPatchSizes[patch_to] = intPart;
                }
                    
            } else
            {
                nextGenerationPatchSizes[patch_to] = SSP->patchCapacity[patch_to];
            }
        }
            

        // resize
        auto nbEvents = migrationEventsForEachDestinationPatch[patch_to].size();
        BackwardMigration[patch_to].resize(nbEvents);
        BackwardMigrationIndex[patch_to].resize(nbEvents);
        

        // Set backward migration rates
        //if (patch_to == 0)
            //std::cout << sumOfIncomingMigrants << " migrants incoming to " << patch_to << "\n";
        for (unsigned migrantEventIndex = 0 ; migrantEventIndex < nbEvents ; ++migrantEventIndex)
        {
            MigrationEvent& migrationEvent = migrationEventsForEachDestinationPatch[patch_to].getMigrationEvent(migrantEventIndex);

            //std::cout << "from patch " << migrationEvent.comingFrom << " to patch " << patch_to << ", migration event announces a migration rate of " << migrationEvent.nbInds / sumOfIncomingMigrants << "\n";
    
                BackwardMigration[patch_to][migrantEventIndex] = migrationEvent.nbInds / sumOfIncomingMigrants;
                BackwardMigrationIndex[patch_to][migrantEventIndex] = migrationEvent.comingFrom;

            /*{
                auto rate = BackwardMigration[patch_to][migrantEventIndex];
                bool cond = BackwardMigration[patch_to][migrantEventIndex] >= 0.0 && BackwardMigration[patch_to][migrantEventIndex] <= 1.0;
                if (!cond)
                {
                    std::cout << "BackwardMigration["<<patch_to<<"]["<<migrantEventIndex<<"] = " << BackwardMigration[patch_to][migrantEventIndex] << "\n";
                    std::cout << "sumOfIncomingMigrants = " << sumOfIncomingMigrants << "\n";
                    std::cout << "migrationEvent.nbInds = " << migrationEvent.nbInds << "\n";
                }

                assert(cond);
            }*/

            

            assert(BackwardMigration[patch_to][migrantEventIndex] >= 0.0 && BackwardMigration[patch_to][migrantEventIndex] <= 1.0);
            assert(BackwardMigrationIndex[patch_to][migrantEventIndex] >= 0.0 && BackwardMigrationIndex[patch_to][migrantEventIndex] < GP->PatchNumber);
        
        }

        // Reorder
        /*std::cout << "---\n";
        for (size_t fakeFrom = 0 ; fakeFrom < BackwardMigrationIndex[patch_to].size() ; ++fakeFrom)
        {
            std::cout << patch_to << "->" << BackwardMigrationIndex[patch_to][fakeFrom] << " = " << BackwardMigration[patch_to][fakeFrom] << "\n";
        }*/
            
        auto idx = reverse_sort_indexes(BackwardMigration[patch_to]);
        reorder(BackwardMigration[patch_to], idx, 1);
        reorder(BackwardMigrationIndex[patch_to], idx, 2, GP->PatchNumber);
    }

    // Security
    assert(BackwardMigration.size() == GP->PatchNumber);

#ifdef DEBUG
    std::cout << "Exit in DispersalData::computeBackwardMigrationRates_from_migrationEventsForEachDestinationPatch\n";
#endif    
}


std::vector<std::vector<double>> DispersalData::FromProbaLineToFullFormForwardMigration(
                                       std::vector<double>& probs,
                                       int center,
                                       int CurrentPatchNumber
                                       )
{
    std::vector<std::vector<double>> FullFormForwardMigration;



    // security
    if (probs.size() <= center)
    {
        std::cout << "center (first element given to DispMat LSS if LSS was used) is larger (" << center + 1 << ") than the number of probabilities indicated (" << probs.size() << ")" << std::endl;
        abort();
    }
    if (probs.size() > CurrentPatchNumber)
    {
        std::cout << "There are more dispersal probabilities (" << probs.size() << ") than patches (" << GP->PatchNumber << ")" << std::endl;
        abort();
    }
    double sum = 0;
    for (int i=0;i<probs.size();i++)
    {
        sum += probs[i];
    }
    if (std::abs(sum - 1) > 0.0001)
    {
        std::cout << "In option --m (--DispMat), mode LSS, the sum of probabilities is different from one. It sums up to " <<sum<<". Of course, it can sometimes be hard to get a value that is exactly equal to one (typically 0.333333 + 0.333333 + 0.333333 is not equal to 1.0) but SimBit will allow some approximation. (Oops, that was found only at the second security catch)!\n";
        abort();
    }

    if (probs.size() > CurrentPatchNumber)
    {
        std::cout << "There are more dispersal probabilities (" << probs.size() << ") than patches (" << GP->PatchNumber << ") (Oops, that was found only at the second security catch which means there is a bug)!" << std::endl;
        abort();
    }


    // Loop through all patch from
    for (uint32_t patch_from=0 ; patch_from < CurrentPatchNumber ; patch_from++)
    {
        std::vector<double> probsFromPatch(CurrentPatchNumber,-1.0);
        for (uint32_t patch_to = 0 ; patch_to < CurrentPatchNumber ; ++patch_to)
        {
            int positionRelativeToCenter = patch_to - patch_from;
            int indexInProbs = center + positionRelativeToCenter;
            if (indexInProbs >= 0 && indexInProbs < probs.size())
            {
                probsFromPatch[patch_to] = probs[indexInProbs];
            } else
            {
                probsFromPatch[patch_to] = 0.0;
            }
        }

        // Rescale
        double currentSumProbs = 0.0;
        for (uint32_t patch_to = 0 ; patch_to < CurrentPatchNumber ; ++patch_to)
        {
            assert(probsFromPatch[patch_to] >= 0.0 && probsFromPatch[patch_to] <= 1.0);
            currentSumProbs += probsFromPatch[patch_to];
        }
        if (currentSumProbs != 1.0)
        {
            assert(currentSumProbs < 1.0);
            assert(currentSumProbs > 0.0);
            for (uint32_t patch_to = 0 ; patch_to < CurrentPatchNumber ; ++patch_to)
            {
                probsFromPatch[patch_to] /= currentSumProbs;
            }
        }   

        // assert
        double sumOfProbs = 0.0;
        for (uint32_t patch_to = 0 ; patch_to < CurrentPatchNumber ; ++patch_to)
        {
            sumOfProbs += probsFromPatch[patch_to];
        }
        assert(fabs(sumOfProbs - 1.0) < 0.00000001);

        // build object
        FullFormForwardMigration.push_back(probsFromPatch);
    }
    
    // security
    assert(FullFormForwardMigration.size() == CurrentPatchNumber);
    for (int patch_to=0 ; patch_to < CurrentPatchNumber ; patch_to++)
    {
        assert(FullFormForwardMigration[patch_to].size() == CurrentPatchNumber);
    }

    return FullFormForwardMigration;
}

void DispersalData::readDispMat(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option '--m (--DispMat)', the std::string that is read is: " << input.print() << std::endl;
#endif

    for (int generation_index = 0 ; generation_index < GP->__GenerationChange.size() ; generation_index++)
    {
        int generation = input.GetNextGenerationMarker(generation_index);

        int CurrentPatchNumber = GP->__PatchNumber[generation_index];
        assert(CurrentPatchNumber > 0);

        std::string Mode = input.GetNextElementString();

        std::vector<std::vector<double>> FFFM; // Stands for FullFormForwardMigration
        FFFM.reserve(CurrentPatchNumber);
        
        if (Mode.compare("isolate") == 0)
        {
            for (size_t patch_from = 0 ; patch_from < CurrentPatchNumber ; ++patch_from)
            {
                std::vector<double> FFMi;
                FFMi.reserve(CurrentPatchNumber);
                for (size_t patch_to = 0 ; patch_to < CurrentPatchNumber ; ++patch_to)
                {
                    if (patch_to == patch_from)
                        FFMi.push_back(1.0);
                    else
                        FFMi.push_back(0.0);
                }
                FFFM.push_back(FFMi);
            }
        } else if (Mode.compare("2DSS") == 0)
        {
            size_t nbRows = input.GetNextElementInt();
            size_t nbCols = input.GetNextElementInt();
            double migRate = input.GetNextElementDouble();
            
            
            // Diagonal migration is currently not allowed
            /*
            auto diagonalString = input.GetNextElementString();
            bool isDiagonal;
            if (diagonalString == "diag")
            {
                isDiagonal = true;
            } else
            {
                if (diagonalString != "nodiag")
                {
                    std::cout << "In option --m (--DispMat), expected information about whether migration should be allowed in diagonal. Expected either 'diag' or 'nodiag' but received '" << diagonalString << "' instead\n";
                    abort();
                }
                isDiagonal = false;
            }
            size_t maxNbNeighbours = isDiagonal ? 8 : 4;
            */
                

            size_t maxNbNeighbours = 4;
            double notMigRate = 1 - maxNbNeighbours * migRate;

            if (nbRows * nbCols != CurrentPatchNumber)
            {
                std::cout << "In option --m (--DispMat), mode 2DSS, received " << nbRows << " rows and " << nbCols << " columns with a migration rate to an adjacent patch of " << migRate << ". Note btw that patch are numbered by column.\nThe number of patches asked with option PN at this generation is " << CurrentPatchNumber << " and is not equal to "  << nbRows << " X " << nbCols << " = " << nbRows * nbCols << "\n";
                abort();
            }

            if (!(notMigRate >= 0.0 && notMigRate <= 1.0))
            {
                std::cout << "In option --m (--DispMat), mode 2DSS, received " << nbRows << " rows and " << nbCols << " columns with a migration rate to an adjacent patch of " << migRate << ". Note btw that patch are numbered by column.\n";
                /*if (isDiagonal)
                {
                    std::cout << "Because you chose to allow diagonal migration, the migration rate should be bounded between 0 and 1/8 = 0.125\n";
                } else
                {
                    std::cout << "Because you chose to not allow diagonal migration, the migration rate should be bounded between 0 and 1/4 = 0.25\n";
                }*/
                std::cout << "The migration rate should be bounded between 0 (included) and 1/4 = 0.25 (included)\n";
                    
                abort();
            }


            /*
                3 X 10
                    0  3  6  9  12 15 18 21 24 27
                    1  4  7  10 13 16 19 22 25 28
                    2  5  8  11 14 17 20 23 26 29
            */


            
            for (size_t patch_from = 0 ; patch_from < CurrentPatchNumber ; ++patch_from)
            {
                std::vector<double> FFFMi;
                FFFMi.reserve(CurrentPatchNumber);
                for (size_t patch_to = 0 ; patch_to < CurrentPatchNumber ; ++patch_to)
                {
                    double m = 0;

                    size_t patch_from_col = patch_from / nbRows;
                    size_t patch_from_row = patch_from % nbRows;
                    size_t patch_to_col = patch_to / nbRows;
                    size_t patch_to_row = patch_to % nbRows;

                    if (patch_from == patch_to)
                    {
                        size_t nbBorders = 0;
                        if (patch_from_col == 0) ++nbBorders;
                        if (patch_from_row == 0) ++nbBorders;
                        if (patch_from_col == nbCols - 1) ++nbBorders;
                        if (patch_from_row == nbRows - 1) ++nbBorders;
                        m = notMigRate + nbBorders * migRate;
                    } else if (
                            (
                                myAbs(patch_from_col, patch_to_col) == 1
                                &&
                                myAbs(patch_from_row, patch_to_row) == 0
                            )
                            ||
                            (
                                myAbs(patch_from_col, patch_to_col) == 0
                                &&
                                myAbs(patch_from_row, patch_to_row) == 1
                            )
                        )
                    {
                        m = migRate;
                    }
                    FFFMi.push_back(m);
                }
                FFFM.push_back(FFFMi); 
            }
        } else if (Mode.compare("LSS") == 0) // Linear Stepping Stone (but with as many stones as we want and potential assymetry). First value is the number of probabilities that will follow. Then are the probabilities which must sum to one, the last value is an integer which indicate which of the probabilities (0 based counting) correspond to the probability of not migrating
        {
            // Gather values
            int center = input.GetNextElementInt();
            std::vector<double> probs;
            double sum = 0.0;

            while (input.IsThereMoreToRead() && input.PeakNextElementString().at(0) != '@')
            {
                double x = input.GetNextElementDouble();
                if (x < 0.0)
                {
                    std::cout << "In option --m (--DispMat), mode LSS, received a negative probability (received " <<x<<")\n";
                    abort();
                }
                sum+=x;
                probs.push_back(x);
            }
            
            if (std::abs(sum - 1.0) > 0.0001 )
            {
                std::cout << "In option --m (--DispMat), mode LSS, the sum of probabilities is different from one. It is " <<sum<<". Of course, it can sometimes be hard to get a value that is exactly equal to one (typically 0.333333 + 0.333333 + 0.333333 is not equal to 1.0) but SimBit will allow some approximation.\n";
                abort();
            }

            if (center < 0 || center >= probs.size())
            {
                std::cout << "In option '--m (--DispMat)', Mode 'LSS', " << probs.size() << " probabilities were received. Center is " << center << " which is either lower than zero or equal or bigger to " << probs.size() << ". The center is indicated by zero based counting.\n";
                abort();
            }
            
            
            FFFM = this->FromProbaLineToFullFormForwardMigration(
                    probs,
                    center,
                    CurrentPatchNumber
                );
            
        } else if (Mode == "OnePatch" || Mode == "onePatch" || Mode == "onepatch" || Mode == "Onepatch" || Mode == "ONEPATCH" || Mode == "NA") // A single population, no dispersal
        {
            if (CurrentPatchNumber != 1)
            {
                std::cout << "In option '--DispMat', from generation " << GP->__GenerationChange[generation_index] << " received 'OnePatch' (or euivalently 'NA') (either by default of via command line input) but there is more than one patch from this generation ( there are " << CurrentPatchNumber << "\n";
                    abort();
            }
            
            std::vector<double> v;
            v.push_back(1.0);
            FFFM.push_back(v);
        } else if (Mode.compare("Island") == 0 || Mode.compare("island") == 0) // Simple island model. Take one value, the probability of migrating anywhere
        {

            double totalMigration = input.GetNextElementDouble();
            if (totalMigration < 0.0 || totalMigration > 1.0)
            {
                std::cout << "In '--DispMat', mode 'island' the received probability of dispersing from one patch to any other patch is " << totalMigration << ". Sorry, a probability must be bourned between 0 and 1.\n";
                abort();   
            }

            double noMigration = 1 - totalMigration;
            double Migration = totalMigration / (CurrentPatchNumber - 1);


            if (CurrentPatchNumber==1)
            {
                if (noMigration != 1.0)
                {
                    std::cout << "In '--DispMat', mode 'island' the received probability of dispersing from one patch to any other patch is " << totalMigration << " but there is only one patch in the simulation from generation " << generation << ". The probability of dispersing can only be 0. If there is only one patch, instead of using --m island 0, better use --m OnePatch though.\n";
                    abort();
                }
            } else
            {
                if (fabs((noMigration + Migration * (CurrentPatchNumber-1)) - 1.0) > 0.00000001)
                {
                    std::cout << "In '--DispMat', mode 'island' the received probability of dispersing from one patch to any other patch is " << totalMigration << ". From that SimBit computed, from generation "<< generation <<", the probability of migrating from any patch to any other specific patch as " << Migration << " because there are " << CurrentPatchNumber << " patches. But SimBit computed a resulting probability of not migrating of "<<noMigration<<" leading to 'noMigration + Migration * (CurrentPatchNumber-1)' to not equal one. Sounds like an internal error but you might want to check your parameters anyway.\n";
                    abort();
                }
            }

            if (noMigration < 0.0 || noMigration > 1.0)
            {
                std::cout << "In '--DispMat', mode 'island' the received probability of dispersing from one patch to any other specific patch is " << Migration << ". Given that there are "<< CurrentPatchNumber <<" from generation " << generation << ", this means that the probability of not dispersing is " << noMigration << " ( 'noMigration = 1 - (Migration * (CurrentPatchNumber - 1))' ). Sorry, a probability must be bourned between 0 and 1.\n";
                std::cout << "That might an internal error because this mistake should have been caught earlier. But please check your input parameters.\n";
                abort();
            }


            

            for ( int patch_from = 0 ; patch_from < CurrentPatchNumber ; ++patch_from )
            {
                std::vector<double> v;
                for ( int patch_to = 0 ; patch_to < CurrentPatchNumber ; ++patch_to )
                {
                    if (patch_from == patch_to)
                    {
                        v.push_back(noMigration);
                    } else
                    {
                        v.push_back(Migration);
                    }
                }
                assert(v.size() == CurrentPatchNumber);
                FFFM.push_back(v);
            }
        } else if (Mode.compare("LinearNormal")==0) // normal distribution kernel. First value is 'sd', second value is how many sd are considered before truncating the distribution
        {
            double sd(input.GetNextElementDouble());
            double HowMany_sd(input.GetNextElementDouble());
            
            std::vector<double> OneSide;
            
            int patch_index = 0;
            double up;
            double down;
            double proba;
            
            while(1)
            {
                down = patch_index - 0.5;
                up = down + 1;
                
                proba =  0.5 * (  erf(up/(sd*sqrt(2))) - erf(down/(sd*sqrt(2)))  );
                if (patch_index <= HowMany_sd*sd)
                {
                    OneSide.push_back(  proba  );
                } else
                {
                    break;
                }
                ++patch_index;
                if (patch_index*2 > CurrentPatchNumber) // to ensure the list of probabilities is not longer than the number of patches
                {
                    break;
                }
            }
            std::vector<double> OtherSide = OneSide;  //Probabilities for migrating left
            reverse(OtherSide.begin(),OtherSide.end());
            OtherSide.pop_back(); // remove the probability for remaining in the same patch.
            std::vector<double> probas;
            // Merge both vectors
            assert(OneSide.size() == OtherSide.size() + 1);
            probas.reserve(OneSide.size() + OtherSide.size());
            probas.insert( probas.end(), OtherSide.begin(), OtherSide.end() );
            probas.insert( probas.end(), OneSide.begin(), OneSide.end() );
            assert(probas.size() == OneSide.size() + OtherSide.size());
            assert(probas.size()%2 != 0);
            
            double sum =0;
            for (int i = 0; i < probas.size() ; i++)
            {
                sum += probas[i];
            }
            
            assert(sum < 1);
            double newsum =0;
            for (int i = 0 ; i < probas.size() ; i++)
            {
                probas[i] /= sum;
                newsum += probas[i];
            }
            
            if (std::abs(newsum - 1) > 0.000001)
            {
                std::cout << "Internal error! newsum = " << newsum << " | sum =" << sum << " | The probas are" << std::endl;
                for (int i = 0 ; i < probas.size() ; i++)
                {
                    std::cout << probas[i] << std::endl;
                }
                abort();
            }
            assert(probas.size() == probas.capacity());
            int center = (probas.size()-1)/2;
            FFFM = this->FromProbaLineToFullFormForwardMigration(probas, center, CurrentPatchNumber);
        } else if (Mode.compare("A")==0)
        {
            for (int patch_from = 0 ; patch_from < CurrentPatchNumber; patch_from++)
            {
                std::vector<double> v;
                for (int patch_to = 0 ; patch_to < CurrentPatchNumber; patch_to++)
                {
                    v.push_back(input.GetNextElementDouble());
                }
                FFFM.push_back(v);
            }
        }
        else
        {
            std::cout << "Sorry, for DispMat Mode only Mode 'A', 'LSS', '2DSS', 'isolate', 'OnePatch' (or 'NA'), 'LinearNormal' and 'Island' are recognized so far. Mode received (" << Mode << ") is unknown" << std::endl;
            abort();
        }

        // Security
        assert(FFFM.size() == CurrentPatchNumber);
        for (int patch_index = 0 ; patch_index < CurrentPatchNumber ; patch_index++)
        {
            assert(FFFM[patch_index].size() == CurrentPatchNumber);
        }

        // Set this->__ForwardMigration
        this->pushBack__ForwardMigrationRate(FFFM, CurrentPatchNumber);
    }
#ifdef DEBUG
    std::cout << "readDispMat is finished" << std::endl;
#endif
}


void DispersalData::pushBack__ForwardMigrationRate(const std::vector<std::vector<double>>& FFFM, const int& CurrentPatchNumber)
{
    assert(FFFM.size() == CurrentPatchNumber);
    forwardMigration.resize(CurrentPatchNumber);
    forwardMigrationIndex.resize(CurrentPatchNumber);
    for (unsigned patch_from = 0 ; patch_from < CurrentPatchNumber ; ++patch_from)
    {
        assert(FFFM[patch_from].size() == CurrentPatchNumber);
        forwardMigration[patch_from].resize(0);
        forwardMigrationIndex[patch_from].resize(0);
        double totalRateFrom = 0.0;
        for (unsigned patch_to = 0 ; patch_to < CurrentPatchNumber ; ++patch_to)
        {
            auto& rate = FFFM[patch_from][patch_to];
            //std::cout << "from " << patch_from << " to " << patch_to << ": " << rate << "\n";
            totalRateFrom += rate;
            assert(rate >= 0.0 && rate <= 1.0);
            if (rate > 0.0)
            {
                forwardMigration[patch_from].push_back(rate); // No need to order from highest to lowest probability unlike backward migration rate
                forwardMigrationIndex[patch_from].push_back(patch_to);
            }
        }


        // Assert
        assert(fabs(totalRateFrom - 1.0) < 0.00000001);


        assert(forwardMigration[patch_from].size() > 0);
        assert(forwardMigration[patch_from].size() == forwardMigrationIndex[patch_from].size());
    }


    __forwardMigration.push_back(forwardMigration);
    __forwardMigrationIndex.push_back(forwardMigrationIndex);

    assert(__forwardMigration.size() == __forwardMigrationIndex.size());

    // below what is commented out is not needed as SSP update should do it before the first generation
    //if (__forwardMigration.size() > 1)
    //{
    //    forwardMigration = __forwardMigration[0];
    //    forwardMigrationIndex = __forwardMigrationIndex[0];
    //}
}


void DispersalData::print()
{

    /*
    std::cout << "\n\nConstant forward migration states:\n";
    {
        // initialize F
        std::vector<std::vector<double>> F(GP->PatchNumber);
        for (size_t patch_from = 0 ; patch_from < GP->PatchNumber ; ++patch_from)
        {
            F[patch_from].resize(GP->PatchNumber,0.0);
        }

        // compute F
        for (size_t patch_from = 0 ; patch_from < GP->PatchNumber ; ++patch_from)
        {
            for (size_t fakePatchTo = 0 ; fakePatchTo < forwardMigration[patch_from].size() ; ++fakePatchTo)   
            {
                auto patch_to = forwardMigrationIndex[patch_from][fakePatchTo];
                assert(F[patch_from][patch_to] == 0.0);
                F[patch_from][patch_to] = forwardMigration[patch_from][patch_to];
            }
        }


        // print F
        for (size_t patch_from = 0 ; patch_from < GP->PatchNumber ; ++patch_from)
        {
            for (size_t patch_to = 0 ; patch_to < GP->PatchNumber ; ++patch_to)
            {
                auto& m = F[patch_from][patch_to];
                if (m > 0.0)
                    std::cout << patch_from << " -> " << patch_to << " = " << m << "\n";
                else
                    assert(m==0.0);
            }
        } 
    }
    */
        






    std::cout << "\n\nBackward migration states:\n";
    {
        // initialize F
        std::vector<std::vector<double>> F(GP->PatchNumber);
        for (size_t patch_from = 0 ; patch_from < GP->PatchNumber ; ++patch_from)
        {
            F[patch_from].resize(GP->PatchNumber,0.0);
        }

        // compute F
        for (size_t patch_to = 0 ; patch_to < GP->PatchNumber ; ++patch_to)
        {
            for (size_t fakePatchFrom = 0 ; fakePatchFrom < BackwardMigration[patch_to].size() ; ++fakePatchFrom)   
            {
                auto patch_from = BackwardMigrationIndex[patch_to][fakePatchFrom];
                assert(F[patch_from][patch_to] == 0.0);
                F[patch_from][patch_to] = BackwardMigration[patch_to][fakePatchFrom];
            }
        }

        // print F
        for (size_t patch_to = 0 ; patch_to < GP->PatchNumber ; ++patch_to)
        {
            for (size_t patch_from = 0 ; patch_from < GP->PatchNumber ; ++patch_from)
            {
                auto& m = F[patch_from][patch_to];
                if (m > 0.0)
                    std::cout << patch_from << " -> " << patch_to << " = " << m << "\n";
                else
                    assert(m==0.0);
            }
        }   
    }



    /*
    std::cout << "\n\nBackward migration structure:\n";
    {
        for (size_t patch_to = 0 ; patch_to < GP->PatchNumber ; ++patch_to)
        {
            for (size_t fakePatchFrom = 0 ; fakePatchFrom < BackwardMigration[patch_to].size() ; ++fakePatchFrom)   
            {
                auto patch_from = BackwardMigrationIndex[patch_to][fakePatchFrom];
                std::cout << patch_from << " -> " << patch_to << ": " << BackwardMigration[patch_to][fakePatchFrom] << "\n";
            }
            std::cout << "---\n";
        }
    }
    */
}
