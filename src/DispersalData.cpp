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
    
    //std::cout << "SSP->growthK["<<patch_index<<"] = " << SSP->growthK[patch_index] << "\n";

    if (GP->nbSpecies == 1)
    {
        if (SSP->growthK[patch_from] == -1.0) // Simple exponential growth
        {
            nbOffs = rn_t;
        }
        else if (SSP->growthK[patch_from] == -2.0) // -2 means always at true carrying capacity
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
        nbOffs = d(GP->mt);
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

    return nbOffs;
}

void DispersalData::getMigrationEventsForEachDestinationPatch(std::vector<std::vector<std::vector<double>>>& CumSumFits, std::vector<ListMigrationEvents>& migrationEventsForEachDestinationPatch)
{
    
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
        if (SSP->isAnySelection)
        {
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
            meanFitness = totalFitness / SSP->patchSize[patch_from];
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
            double rn_t = totalFitness * SSP->fecundityForFitnessOfOne;
            double r    = meanFitness * SSP->fecundityForFitnessOfOne;
            /*if (patch_from == 1003)
            {
                std::cout << "meanFitness = " << meanFitness << "\n";
                std::cout << "totalFitness = " << totalFitness << "\n";
            }*/

            nbOffs = computeNbOffspringsProducedInPatch(patch_from, n_t, rn_t, r); // returns a double for performance reasons (among other reasons)
        } else
        {
            nbOffs = totalFitness;
        }
            
        
        

        //if (patch_from == 1003) std::cout << "nbOffs of patch "<< patch_from << " = " << nbOffs << ":\n";

        assert(forwardMigration.size() == forwardMigrationIndex.size());
        assert(forwardMigration.size() > patch_from);
        for (unsigned FMpatchToFakeIndex = 0 ; FMpatchToFakeIndex < forwardMigration[patch_from].size() ; ++FMpatchToFakeIndex)
        {
            auto patch_to = forwardMigrationIndex[patch_from][FMpatchToFakeIndex];
            assert(patch_to >= 0 && patch_to < GP->PatchNumber);
            double forwardMigrationRate = (double) forwardMigration[patch_from][FMpatchToFakeIndex];
            assert(forwardMigrationRate >= 0.0);
            //if (patch_from == 1003) std::cout << "\t" << forwardMigrationRate * nbOffs << " of them go to " << patch_to << "\n";
            migrationEventsForEachDestinationPatch[patch_to].addMigrationEvent(forwardMigrationRate * nbOffs, patch_from); // Will only be considered if greater than 0 migrants
        }
    }

    assert(migrationEventsForEachDestinationPatch.size() == GP->PatchNumber);
}

void DispersalData::setOriginalBackwardMigrationIfNeeded()
{
    // memory
    assert(forwardMigration.size() == GP->PatchNumber);
    BackwardMigration.resize(GP->PatchNumber);
    BackwardMigrationIndex.resize(GP->PatchNumber);



    if (SSP->DispWeightByFitness)
    {
        return;
    } else
    {
        assert(SSP->fecundityForFitnessOfOne == -1);

        std::vector<ListMigrationEvents> migrationEventsForEachDestinationPatch(GP->PatchNumber);


        // Set backward migration rates
        for (unsigned patch_from = 0 ; patch_from < GP->PatchNumber ; patch_from++)
        {
            assert(forwardMigrationIndex.size() > patch_from);

            for (unsigned FMpatchToFakeIndex = 0 ; FMpatchToFakeIndex < forwardMigration[patch_from].size() ; ++FMpatchToFakeIndex)
            {
                auto patch_to = forwardMigrationIndex[patch_from][FMpatchToFakeIndex];
                assert(patch_to >= 0 && patch_to < GP->PatchNumber);
                double forwardMigrationRate = (double) forwardMigration[patch_from][patch_to];
                migrationEventsForEachDestinationPatch[patch_to].addMigrationEvent(forwardMigrationRate, patch_from); // Will only be considered if rate is greater than 0
            }
        }

        computeBackwardMigrationRates_from_migrationEventsForEachDestinationPatch(migrationEventsForEachDestinationPatch);
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
    getMigrationEventsForEachDestinationPatch(CumSumFits, migrationEventsForEachDestinationPatch);
    computeBackwardMigrationRates_from_migrationEventsForEachDestinationPatch(migrationEventsForEachDestinationPatch); // will set nextGenerationPatchSizes if fec != -1

    if (SSP->fecundityForFitnessOfOne != -1.0)
    {
        return nextGenerationPatchSizes;
    } else
    {
        return SSP->patchCapacity;
    }
}

void DispersalData::computeBackwardMigrationRates_from_migrationEventsForEachDestinationPatch(std::vector<ListMigrationEvents>& migrationEventsForEachDestinationPatch)
{
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
                if (GP->random_0and1(GP->mt) < fracPart)
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
            

        // Reset to zero
        BackwardMigration[patch_to].resize(migrationEventsForEachDestinationPatch[patch_to].size());
        BackwardMigrationIndex[patch_to].resize(migrationEventsForEachDestinationPatch[patch_to].size());

        // Set backward migration rates
        //if (patch_to == 0)
            //std::cout << sumOfIncomingMigrants << " migrants incoming to " << patch_to << "\n";
        for (unsigned migrantEventIndex = 0 ; migrantEventIndex < migrationEventsForEachDestinationPatch[patch_to].size() ; ++migrantEventIndex)
        {
            MigrationEvent& migrationEvent = migrationEventsForEachDestinationPatch[patch_to].getMigrationEvent(migrantEventIndex);

            //if (patch_to == 0)
            //    std::cout << "\t" << migrationEvent.nbInds << " migrants incoming from " << migrationEvent.comingFrom << " ("<<(double)migrationEvent.nbInds / (double)sumOfIncomingMigrants<<")\n";

            BackwardMigration[patch_to][migrantEventIndex] =  migrationEvent.nbInds / sumOfIncomingMigrants;
            BackwardMigrationIndex[patch_to][migrantEventIndex] = migrationEvent.comingFrom;
        }
            

        // Reorder
        auto idx = reverse_sort_indexes(BackwardMigration[patch_to]);
        reorder(BackwardMigration[patch_to], idx);
        reorder(BackwardMigrationIndex[patch_to], idx);
    }

    // Security
    assert(BackwardMigration.size() == GP->PatchNumber);
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
        std::cout << "center (last element given to DispMat LSS if LSS was used) is larger (" << center + 1 << ") than the number of probabilities indicated (" << probs.size() << ")" << std::endl;
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
        std::cout << "There are more dispersal probabilities (" << probs.size() << ") than patches (" << GP->PatchNumber << ") (Oops, that was found only at the second security catch)!" << std::endl;
        abort();
    }

    int howManyElementsBeforeCenter = center;
    int howManyElementsAfterCenter = probs.size() - center - 1;
    assert(howManyElementsAfterCenter >= 0);

    // Loop through all destination patch
    for (int patch_from=0 ; patch_from < CurrentPatchNumber ; patch_from++)
    {
        // probsTmp wil be a truncated version of probs in case we hit a border of the world
        auto probsTmp = probs;

        
        // truncate beginning of probsTmp
        int nbElementsToRemoveBefore = howManyElementsBeforeCenter - patch_from;
        assert(nbElementsToRemoveBefore < (int) probsTmp.size());
        int newCenter;
        if (nbElementsToRemoveBefore > 0)
        {
            newCenter = center - nbElementsToRemoveBefore;
            assert(newCenter >= 0);
            probsTmp.erase(probsTmp.begin(), probsTmp.begin() + nbElementsToRemoveBefore);
        } else
        {
            newCenter = center;
        }

        // truncate end of probsTmp
        int nbElementsToRemoveAfter = patch_from + howManyElementsAfterCenter - CurrentPatchNumber + 1;
        if (nbElementsToRemoveAfter >= (int) probsTmp.size())
        {
            std::cout << "Internal error in DispersalData::FromProbaLineToFullFormForwardMigration. nbElementsToRemoveBefore = " << nbElementsToRemoveBefore << " nbElementsToRemoveAfter = " << nbElementsToRemoveAfter << " patch_from = "<< patch_from << " CurrentPatchNumber = " << CurrentPatchNumber << " howManyElementsAfterCenter = " << howManyElementsAfterCenter << " howManyElementsBeforeCenter = "<< howManyElementsBeforeCenter << " probsTmp.size() = " << probsTmp.size() << "\n";
            abort();
        }
        if (nbElementsToRemoveAfter > 0)
        {
            probsTmp.resize(probsTmp.size() - nbElementsToRemoveAfter);
        }
        
        // If vTmp diffes from 'v', then we need to scale probabilities up
        if (probs.size() != probsTmp.size())
        {

            double sum=0;
            int i;
            for ( i = 0 ; i < probsTmp.size() ; i++ )
            {
                sum += probsTmp[i];
            }
            assert(sum <= 1 && sum >= 0);
            if (sum!=1)
            {
                for ( i = 0 ; i < probsTmp.size() ; i++ )
                {
                    probsTmp[i] /= sum;
                }
            }
        }
        assert(probsTmp.size());


        // insert in completeProbs
        std::vector<double> completeProbs(CurrentPatchNumber,0.0);
        int probsIndex = 0;
        for (int patch_to = patch_from - newCenter ; patch_to < patch_from - newCenter + probsTmp.size() ; patch_to++)
        {
            assert(patch_to >= 0);
            assert(patch_to < CurrentPatchNumber);
            completeProbs[patch_to] = probsTmp[probsIndex];
            probsIndex++;
        }
        assert(probsIndex == probsTmp.size());
            

        // build object
        FullFormForwardMigration.push_back(completeProbs);
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

        std::vector<std::vector<double>> FFFM; // Stands fir FullFormForwardMigration
        
        if (Mode.compare("LSS") == 0) // Linear Stepping Stone (but with as many stones as we want and potential assymetry). First value is the number of probabilities that will follow. Then are the probabilities which must sum to one, the last value is an integer which indicate which of the probabilities (0 based counting) correspond to the probability of not migrating
        {
            // Gather values
            int nbValues = input.GetNextElementInt();
            std::vector<double> probs;
            double sum = 0.0;
            for (int i = 0 ; i < nbValues ; i++)
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
            int center = input.GetNextElementInt();

            if (center < 0 || center >= nbValues)
            {
                std::cout << "In option '--m (--DispMat)', Mode 'LSS', " << nbValues << " probabilities were expected. Center is " << center << " which is either lower than zero or equal or bigger to " << nbValues << ". The center is indicated by zero based counting.\n";
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
            std::cout << "Sorry, for DispMat Mode only Mode 'A', 'LSS', 'OnePatch' (or 'NA'), 'LinearNormal' and 'Island' are recognized so far. Mode received (" << Mode << ") is unknown" << std::endl;
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
        for (unsigned patch_to = 0 ; patch_to < CurrentPatchNumber ; ++patch_to)
        {
            auto& rate = FFFM[patch_from][patch_to];
            assert(rate >= 0.0 && rate <= 1.0);
            if (rate > 0.0)
            {
                forwardMigration[patch_from].push_back(rate); // No need to order from highest to lowest probability unlike backward migration rate
                forwardMigrationIndex[patch_from].push_back(patch_to);
            }
        }
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


