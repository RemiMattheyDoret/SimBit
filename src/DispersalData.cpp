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


std::vector<int> DispersalData::getNbOffspringProducedPerPatch(std::vector<std::vector<std::vector<double>>>& CumSumFits)
{
    // This function is only called when setting the backward migration rate if fec != -1.0

    assert(SSP->fecundityForFitnessOfOne != -1.0);

    
    //////////////////////////////////////
    /// Offspring production per patch ///
    //////////////////////////////////////

    // compute production of offspring of each patch
    std::vector<int> nbOffspringProduced(GP->PatchNumber);
    for (int patch_index = 0 ; patch_index < GP->PatchNumber ; patch_index++)
    {
        //////////////////////////////////
        /// Get mean and total fitness ///
        //////////////////////////////////
        double totalFitness;
        double meanFitness;
        if (SSP->isAnySelection)
        {
            if (SSP->malesAndFemales)
            {
                assert(SSP->patchSize[patch_index] == CumSumFits[patch_index][0].size() + CumSumFits[patch_index][1].size());
            } else
            {
                assert(SSP->patchSize[patch_index] == CumSumFits[patch_index][0].size());
            }

            totalFitness = CumSumFits[patch_index][0].back();
            meanFitness = totalFitness / SSP->patchSize[patch_index];
        } else
        {
            totalFitness = SSP->patchSize[patch_index];
            meanFitness = 1.0;
        }

            

        ////////////////////////////
        /// Compute nbOffsprings ///
        ////////////////////////////
        double n_t  = SSP->patchSize[patch_index];
        double rn_t = totalFitness * SSP->fecundityForFitnessOfOne;
        double r    = meanFitness * SSP->fecundityForFitnessOfOne;
        
        /*
        std::cout << "\n";
        std::cout << "r = " << r << "\n";
        std::cout << "rn_t = " << rn_t << "\n";
        std::cout << "n_t = " << n_t << "\n";
        */
        

        
        double nbOffs = 0.0;
        

        if (GP->nbSpecies == 1)
        {
            if (SSP->growthK[patch_index] == -1.0) // Simple exponential growth
            {
                nbOffs = rn_t;
            }
            else if (SSP->growthK[patch_index] == -2.0) // -2 means always at true carrying capacity
            {
                if (r <= 1) // if decline
                {
                    nbOffs = rn_t;
                } else
                {
                    nbOffs = n_t + (r-1) * n_t * (1 - n_t / SSP->patchCapacity[patch_index]);
                }
            } else
            {
                if (r <= 1) // if decline
                {
                    nbOffs = rn_t;
                } else // if growth
                {
                    nbOffs = n_t + (r-1) * n_t * (1 - n_t / SSP->growthK[patch_index]);
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
                sumOfAlphasProd_competition += GP->speciesCompetition[SSP->speciesIndex][speciesIndex] * GP->allSpeciesPatchSizePreviousGeneration[patch_index][speciesIndex];


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
                    auto causalSpPS = GP->allSpeciesPatchSizePreviousGeneration[patch_index][speciesIndex];
                    sumOfAlphasProd_interaction += interaction.magnitude * causalSpPS;
                } else if (interaction.type == 'C')
                {
                    auto recipientSpPS = GP->allSpeciesPatchSizePreviousGeneration[patch_index][SSP->speciesIndex];
                    sumOfAlphasProd_interaction += interaction.magnitude * recipientSpPS;
                } else if (interaction.type == 'D') 
                {
                    auto causalSpPS = GP->allSpeciesPatchSizePreviousGeneration[patch_index][speciesIndex];
                    auto recipientSpPS = GP->allSpeciesPatchSizePreviousGeneration[patch_index][SSP->speciesIndex];
                    sumOfAlphasProd_interaction += interaction.magnitude * causalSpPS * recipientSpPS;
                }
            }


            // Competition
            if (SSP->growthK[patch_index] == -1.0) // exponential growth
            {
                nbOffs = rn_t;   // No competition possible when exponential growth
            } else if (SSP->growthK[patch_index] == -2) // -2 means always at true carrying capacity
            {
                if (r <= 1) // if decline
                {
                    nbOffs = rn_t;
                } else // if growth
                {
                    nbOffs = n_t + (r-1) * n_t * (1 - sumOfAlphasProd_competition / SSP->patchCapacity[patch_index]);
                }
            } else
            {
                if (r <= 1) // if decline
                {
                    nbOffs = rn_t;
                } else // if growth
                {
                    nbOffs = n_t + (r-1) * n_t * (1 - sumOfAlphasProd_competition / SSP->growthK[patch_index]);
                }
            }
            
            // Adjust for interaction
            //std::cout << "sumOfAlphasProd_interaction = " << sumOfAlphasProd_interaction << "\n";
            nbOffs += sumOfAlphasProd_interaction;

            // Correct if it went below zero
            if (nbOffs < 0.0) nbOffs = 0.0; // if r is very large, then behaviour is chaotic and it can lead to lower than zero values. Species interaction can also lead to lower than zero values
        }

        assert(nbOffs >= 0.0);

        if (SSP->stochasticGrowth)
        {
            std::poisson_distribution<int> d(nbOffs);
            nbOffspringProduced[patch_index] = d(GP->mt);
        } else
        {
            nbOffspringProduced[patch_index] = nbOffs;
        }
    }

    assert(nbOffspringProduced.size() == GP->PatchNumber);
    return nbOffspringProduced;
}

void DispersalData::setOriginalBackwardMigrationIfNeeded()
{
    // memory
    assert(FullFormForwardMigration.size() == GP->PatchNumber);
    BackwardMigration.resize(GP->PatchNumber);
    BackwardMigrationIndex.resize(GP->PatchNumber);

    for (int patch_to = 0 ; patch_to < GP->PatchNumber ; patch_to++)
    {
        BackwardMigration[patch_to].resize(GP->PatchNumber);
        BackwardMigrationIndex[patch_to].resize(GP->PatchNumber);
    }




    if (SSP->DispWeightByFitness)
    {
        return;
    } else
    {
        assert(SSP->fecundityForFitnessOfOne == -1);


        // Set backward migration rates
        for (int patch_to = 0 ; patch_to < GP->PatchNumber ; patch_to++)
        {   
            // assert 
            assert(this->FullFormForwardMigration[patch_to].size() == GP->PatchNumber); // Here patch_to is used as patch_from but just for an assertion

                    
            // Construct sumOfRatesOfMigrant
            double sumOfImmigrationRates = 0.0;
            for (int patch_from = 0 ; patch_from < GP->PatchNumber ; patch_from++)
            {   
                sumOfImmigrationRates += this->FullFormForwardMigration[patch_from][patch_to];
            }


            // Compute fraction of migrants from patch_from and backward migration rate
            for (int patch_from = 0 ; patch_from < GP->PatchNumber ; patch_from++)
            {
                BackwardMigration[patch_to][patch_from] = this->FullFormForwardMigration[patch_from][patch_to] / sumOfImmigrationRates;
                BackwardMigrationIndex[patch_to][patch_from] = patch_from;
            }


            // Reorder
            auto idx = reverse_sort_indexes(BackwardMigration[patch_to]);
            reorder(BackwardMigration[patch_to], idx);
            reorder(BackwardMigrationIndex[patch_to], idx);
        }
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
    assert(this->FullFormForwardMigration.size() == GP->PatchNumber);
    assert(BackwardMigration.size() == GP->PatchNumber);
    assert(BackwardMigrationIndex.size() == GP->PatchNumber);
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


    ////////////////////////////////////
    /// Compute offspring production ///
    ////////////////////////////////////

    std::vector<int> nbOffspringProduced;
    if (SSP->fecundityForFitnessOfOne != -1.0)
    {
        nbOffspringProduced = getNbOffspringProducedPerPatch(CumSumFits);
    } else
    {
        nbOffspringProduced.resize(GP->PatchNumber, 1.0); // They must all be set to any non-zero value so that they don't affect what comes below
    }


    ////////////////////////////////////////////////
    /// Figure out where all these offsprings go ///
    ////////////////////////////////////////////////
    
    // Loop through each patch_to
    for (int patch_to = 0 ; patch_to < GP->PatchNumber ; patch_to++)
    {   
        // assert 
        assert(this->FullFormForwardMigration[patch_to].size() == GP->PatchNumber); // Here patch_to is used as patch_from but just for an assertion

        
        // Construct numberOfMigrantsConstr and Sum of migrants
        // numberOfMigrantsConstr uses rate when fecunditry is infinite
        double sumOfNumberOfMigrants = 0.0;
        std::vector<double> numberOfMigrantsConstr(GP->PatchNumber); // Here individual coming from the same patch are also called "migrants" which is confusing
        for (int patch_from = 0 ; patch_from < GP->PatchNumber ; patch_from++)
        {   
            //std::cout << "nbOffspringProduced["<<patch_from<<"] = " << nbOffspringProduced[patch_from] << "\n";
            numberOfMigrantsConstr[patch_from] = this->FullFormForwardMigration[patch_from][patch_to] * nbOffspringProduced[patch_from];
            sumOfNumberOfMigrants += numberOfMigrantsConstr[patch_from];
        }


        // Compute fraction of migrants from patch_from and backward migration rate
        nextGenerationPatchSizes[patch_to] = 0;
        for (int patch_from = 0 ; patch_from < GP->PatchNumber ; patch_from++)
        {
            if (SSP->fecundityForFitnessOfOne != -1)
            {
                nextGenerationPatchSizes[patch_to] += (int)(numberOfMigrantsConstr[patch_from] + 0.5);
                if (nextGenerationPatchSizes[patch_to] > SSP->patchCapacity[patch_to])
                {
                    nextGenerationPatchSizes[patch_to] = SSP->patchCapacity[patch_to];
                }
            }
                
            BackwardMigration[patch_to][patch_from] = numberOfMigrantsConstr[patch_from] / sumOfNumberOfMigrants;
            BackwardMigrationIndex[patch_to][patch_from] = patch_from;
        }
            

        // Reorder
        auto idx = reverse_sort_indexes(BackwardMigration[patch_to]);
        reorder(BackwardMigration[patch_to], idx);
        reorder(BackwardMigrationIndex[patch_to], idx);
    }

    // Security
    assert(BackwardMigration.size() == GP->PatchNumber);

    
    #ifdef DEBUG
        std::cout << "Exiting DispersalData::SetBackwardMigrationAndGetNextGenerationPatchSizes\n";
    #endif      

    if (SSP->fecundityForFitnessOfOne == -1)
    {
        return SSP->patchCapacity;
    } else
    {
        //std::cout << "nextGenerationPatchSizes[0] = " << nextGenerationPatchSizes[0] << "\n";
        return nextGenerationPatchSizes;
    }
        
}


std::vector<std::vector<double>> DispersalData::FromProbaLineToFullFormForwardMigration(
                                       std::vector<double>& probs,
                                       int center,
                                       int CurrentPatchNumber
                                       )
{
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

            double Migration = input.GetNextElementDouble();
            if (Migration < 0.0 || Migration > 1.0)
            {
                std::cout << "In '--DispMat', mode 'island' the received probability of dispersing from one patch to any other patch is " << Migration << ". Sorry, a probability must be bourned between 0 and 1.\n";
                abort();   
            }

            double noMigration = 1 - (Migration * (CurrentPatchNumber - 1));

            if (noMigration < 0.0 || noMigration > 1.0)
            {
                std::cout << "In '--DispMat', mode 'island' the received probability of dispersing from one patch to any other patch is " << Migration << ". Given that there are "<< CurrentPatchNumber <<" from generation " << generation << ", this means that the probability of not dispersing is " << noMigration << " ( 'noMigration = 1 - (Migration * (CurrentPatchNumber - 1))' ). Sorry, a probability must be bourned between 0 and 1.\n";
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
        __FullFormForwardMigration.push_back(FFFM);
    }
#ifdef DEBUG
    std::cout << "readDispMat is finished" << std::endl;
#endif
}


/*void DispersalData::FromFullFormForwardSetReducedFormForward(std::vector<std::vector<double>> FullFormForwardMigration)
{
    std::vector<std::vector<double>> ToSetForwardMigration;
    std::vector<std::vector<int>> ToSetForwardMigrationIndex;

    for (int patch_to = 0 ; patch_to < FullFormForwardMigration.size() ; patch_to++)
    {
        assert(FullFormForwardMigration.size() == FullFormForwardMigration[patch_to].size());
        std::vector<double> ForwardMigrationTmp;
        std::vector<int> ForwardMigrationIndexTmp;
        for (int patch_from = 0 ; patch_from < FullFormForwardMigration.size() ; patch_from++)
        {
            if (FullFormForwardMigration[patch_from][patch_to] != 0.0)
            {
                assert(FullFormForwardMigration[patch_from][patch_to] > 0.0 && FullFormForwardMigration[patch_from][patch_to] <= 1.0);
                ForwardMigrationTmp.push_back(FullFormForwardMigration[patch_from][patch_to]);
                ForwardMigrationIndexTmp.push_back(patch_from);
            }
        }
        ToSetForwardMigration.push_back(ForwardMigrationTmp);
        ToSetForwardMigrationIndex.push_back(ForwardMigrationIndexTmp);
    }

    // Security
    assert(ToSetForwardMigration.size() == FullFormForwardMigration.size());
    for (int patch_to = 0 ; patch_to < FullFormForwardMigration.size() ; patch_to++)
    {
        double sumOfProb = 0.0;
        for (int fake_patch_from = 0 ; fake_patch_from < ToSetForwardMigration[patch_to].size() ; fake_patch_from++)
        {
            sumOfProb += ToSetForwardMigration[patch_to][fake_patch_from];
        }
        if (std::abs(sumOfProb - 1.0) > 0.00001)
        {
            std::cout << "Internal error! In DispersalData::FromFullFormForwardSetReducedFormForward(), the sum of probabilities for patch_to = " << patch_to << " is " << sumOfProb << " while it should be one.\n";
            abort(); 
        }
    }

    // set value
    this->__ForwardMigration.push_back(ToSetForwardMigration);
    this->__ForwardMigrationIndex.push_back(ToSetForwardMigrationIndex);
}*/
