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

double DispersalData::ComputeEffectOfSpeciesEcologicalInteractions(int patch_index)
{
    if (GP->nbSpecies == 1)
    {
        return 0.0;
    } else
    {   
        int focal_speciesIndex = SSP->speciesIndex;

        double otherSpeciesEffect = 0.0;
        //std::cout << "\n";
        for (int speciesIndex = 0 ; speciesIndex < GP->nbSpecies; speciesIndex++)
        {
            /*
            std::cout << "patch size of focal_speciesIndex (" <<  SSP->speciesName<<") = " << GP->allSpeciesPatchSizes[patch_index][focal_speciesIndex] << "\n";
            std::cout << "patch size of speciesIndex = " << GP->allSpeciesPatchSizes[patch_index][speciesIndex] << "\n";
            std::cout << "GP->speciesEcoRel_type["<<focal_speciesIndex<<"]["<<speciesIndex<<"] = " << GP->speciesEcoRel_type[focal_speciesIndex][speciesIndex] << "\n";
            std::cout << "GP->speciesEcoRel_effect["<<focal_speciesIndex<<"]["<<speciesIndex<<"] = " << GP->speciesEcoRel_effect[focal_speciesIndex][speciesIndex] << "\n";
            */

            double effectOnEntirePatch = GP->allSpeciesPatchSizes[patch_index][focal_speciesIndex];

            if (speciesIndex == focal_speciesIndex)
            {
                assert(GP->speciesEcoRel_type[focal_speciesIndex][speciesIndex] == '0');
                assert(GP->speciesEcoRel_effect[focal_speciesIndex][speciesIndex] == 0.0);
                effectOnEntirePatch *= SSP->growthK[patch_index];
            } else
            {
                effectOnEntirePatch *= GP->speciesEcoRel_effect[focal_speciesIndex][speciesIndex];

                if (GP->speciesEcoRel_type[focal_speciesIndex][speciesIndex] == 'A')
                {
                    otherSpeciesEffect += 
                        effectOnEntirePatch *
                        GP->allSpeciesPatchSizes[patch_index][speciesIndex];   
                } else if (GP->speciesEcoRel_type[focal_speciesIndex][speciesIndex] == 'B')
                {
                    otherSpeciesEffect += 
                        effectOnEntirePatch *
                        GP->allSpeciesPatchSizes[patch_index][speciesIndex] /
                        GP->allSpeciesPatchSizes[patch_index][focal_speciesIndex];
                } else if (GP->speciesEcoRel_type[focal_speciesIndex][speciesIndex] == 'C')
                {
                    otherSpeciesEffect += 
                        effectOnEntirePatch /
                        GP->allSpeciesPatchSizes[patch_index][focal_speciesIndex];
                } else if (GP->speciesEcoRel_type[focal_speciesIndex][speciesIndex] == 'D')
                {
                    otherSpeciesEffect += 
                        effectOnEntirePatch *
                        GP->allSpeciesPatchSizes[patch_index][speciesIndex] *
                        GP->allSpeciesPatchSizes[patch_index][focal_speciesIndex];
                } else if (GP->speciesEcoRel_type[focal_speciesIndex][speciesIndex] != '0')
                {
                    std::cout << "Internal error (but you might want to check your input). Unkown Ecological relationship type " <<  GP->speciesEcoRel_type[focal_speciesIndex][speciesIndex] << " detected in DispersalData::ComputeEffectOfSpeciesEcologicalInteractions.\n";
                    abort();
                }
            }
        }
        assert(GP->allSpeciesPatchSizes[patch_index][focal_speciesIndex] == SSP->patchSize[patch_index]);

        //std::cout << "effect on species '" << SSP->speciesName << "' is "<< otherSpeciesEffect << "\n";

        return otherSpeciesEffect;
    }
}

double DispersalData::ComputeEffectiveNumberOfMigrants(
    std::vector<std::vector<std::vector<double>>>& CumSumFits,
    bool computingOriginal,
    int patch_to,
    int patch_from
)
{
#ifdef DEBUG
    std::cout << "Enters in DispersalData::ComputeEffectiveNumberOfMigrants\n";
#endif        

    double effectiveNumberOfMigrants;
    //std::cout << "computingOriginal = " << computingOriginal << "\n";
    //std::cout << "from = " << patch_from << " to " << patch_to <<"\n";
    

    if (computingOriginal)
    {
        effectiveNumberOfMigrants = (double) SSP->patchCapacity[patch_from] * this->FullFormForwardMigration[patch_from][patch_to];
    } else
    {
        //std::cout << "SSP->patchSize[patch_from] = " << SSP->patchSize[patch_from]  << "\n";
        assert(CumSumFits.size() == GP->PatchNumber);
        assert(patch_from < GP->PatchNumber);
        /*
        std::cout << "SSP->malesAndFemales = " << SSP->malesAndFemales << "\n";
        std::cout << "SSP->DispWeightByFitness = " << SSP->DispWeightByFitness << "\n";
        std::cout << "SSP->fecundityForFitnessOfOne = " << SSP->fecundityForFitnessOfOne << "\n";
        */
        if (SSP->malesAndFemales)
        {
            assert(CumSumFits[patch_from].size() == 2);
            //std::cout << "CumSumFits[patch_from][0].size() = " << CumSumFits[patch_from][0].size()  << "\n";
            //std::cout << "CumSumFits[patch_from][1].size() = " << CumSumFits[patch_from][1].size()  << "\n";
            assert(SSP->patchSize[patch_from] == CumSumFits[patch_from][0].size() + CumSumFits[patch_from][1].size());
        } else
        {
            assert(CumSumFits[patch_from].size() == 1);
            //std::cout << "CumSumFits[patch_from][0].size() = " << CumSumFits[patch_from][0].size()  << "\n";
            assert(SSP->patchSize[patch_from] == CumSumFits[patch_from][0].size());
        }

        if (CumSumFits[patch_from][0].size() == 0 || (SSP->malesAndFemales && CumSumFits[patch_from][1].size() == 0))
        {
            effectiveNumberOfMigrants = 0.0;
        } else
        {
            //std::cout << "line 150\n";
            if (SSP->DispWeightByFitness)
            {
                //std::cout << "line 153\n";
                if (SSP->fecundityForFitnessOfOne != -1.0)   
                {
                    //std::cout << "line 156\n";
                    effectiveNumberOfMigrants =
                        CumSumFits[patch_from][0].back() * // There is in here patchSize within CumSumFits[patch_from][0].back(). Note btw that CumSumFits[patch_from][0] is of length of the patchSize of patch_from
                        SSP->fecundityForFitnessOfOne *
                        this->FullFormForwardMigration[patch_from][patch_to]
                        ; 
                    //std::cout << "line 162 effectiveNumberOfMigrants = " << effectiveNumberOfMigrants << "\n";
                    /*
                    if (effectiveNumberOfMigrants < 0)
                    {
                        std::cout << "line 164 effectiveNumberOfMigrants = " << effectiveNumberOfMigrants << "\n";
                        std::cout << "CumSumFits[patch_from][0].back() = " << CumSumFits[patch_from][0].back() << "\n";
                        std::cout << "SSP->fecundityForFitnessOfOne = " << SSP->fecundityForFitnessOfOne << "\n";
                        std::cout << "this->FullFormForwardMigration["<<patch_from<<"]["<<patch_to<<"] = " << this->FullFormForwardMigration[patch_from][patch_to] << "\n";
                    }
                    */
                } else
                {
                    //std::cout << "line 171\n";
                    effectiveNumberOfMigrants =
                        CumSumFits[patch_from][0].back() * // There is in here patchSize within CumSumFits[patch_from].back(). Note btw that CumSumFits[patch_from][0] is of length of the patchSize of patch_from
                        this->FullFormForwardMigration[patch_from][patch_to]
                        ; 
                }
            } else
            {
                //std::cout << "line 179\n";
                if (SSP->fecundityForFitnessOfOne != -1.0)
                {
                    //std::cout << "line 182\n";
                    // actually if SSP->fecundityForFitnessOfOne dispersal is necessarily weigthed by fitness
                    effectiveNumberOfMigrants =
                        CumSumFits[patch_from][0].back() * // There is in here patchSize within CumSumFits[patch_from].back(). Note btw that CumSumFits[patch_from][0] is of length of the patchSize of patch_from
                        SSP->fecundityForFitnessOfOne *
                        this->FullFormForwardMigration[patch_from][patch_to]
                        ; 
                } else
                {
                    //std::cout << "line 190\n";
                    effectiveNumberOfMigrants =
                        (double) SSP->patchCapacity[patch_from] * 
                        this->FullFormForwardMigration[patch_from][patch_to]
                        ; 
                }
            }
            if (GP->nbSpecies > 1)
            {
                effectiveNumberOfMigrants += ComputeEffectOfSpeciesEcologicalInteractions(patch_to);
                if (effectiveNumberOfMigrants < 0) effectiveNumberOfMigrants = 0;
            }
        }
    }
    
   /* // Get value into integer by watching for for potential overflow!
    double limitToAvoidOverflow = DBL_MAX / GP->PatchNumber - 100;
    if (effectiveNumberOfMigrants > limitToAvoidOverflow)
    {
        effectiveNumberOfMigrants = limitToAvoidOverflow;
    }
    if (effectiveNumberOfMigrants < 0)
    {
        effectiveNumberOfMigrants = 0;
    }*/

    //std::cout << effectiveNumberOfMigrants << " eff. migants from " << patch_from << " to " << patch_to<< "\n";
    //std::cout << "SSP->patchSize[patch_from] = " << SSP->patchSize[patch_from]  << "\n";
    //std::cout << "effectiveNumberOfMigrants = " << effectiveNumberOfMigrants << "\n";
    //std::cout << "line 216 effectiveNumberOfMigrants = " << effectiveNumberOfMigrants << "\n";
    assert(effectiveNumberOfMigrants >= 0.0);
    
    if (!computingOriginal && effectiveNumberOfMigrants > 0.0)
    {
        assert(SSP->patchSize[patch_from] > 0);
    }

#ifdef DEBUG
    std::cout << "Exits in DispersalData::ComputeEffectiveNumberOfMigrants\n";
#endif   
    return effectiveNumberOfMigrants;
}

std::vector<int> DispersalData::SetBackwardMigrationAndGetNextGenerationPatchSizes(
    std::vector<std::vector<std::vector<double>>>& CumSumFits,
    bool computingOriginal
)
{
#ifdef DEBUG
    std::cout << "Enters in DispersalData::SetBackwardMigrationAndGetNextGenerationPatchSizes\n";
#endif
    //std::cout << "computingOriginal = " << computingOriginal << "\n";
    // Fastest solution
    if (!computingOriginal && !SSP->DispWeightByFitness && SSP->fecundityForFitnessOfOne == -1)
    {
        // Nothing to do as BackwardMigration and BackwardMigrationIndex must already have been assigned to the right constant values from when computingOriginal was set.
        return SSP->patchCapacity; // The next generation patch size is simply the patchCapacity
    }


    BackwardMigration.resize(GP->PatchNumber);
    BackwardMigrationIndex.resize(GP->PatchNumber);

    // CumSumFits assertions
    if (!computingOriginal && SSP->DispWeightByFitness)
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
        



    // return object (although the data type will be modified)
    std::vector<double> nextGenerationPatchSizes; // use doubles to avoid overflow trouble
    nextGenerationPatchSizes.resize(GP->PatchNumber,0.0);

    // Loop thorugh each patch_to
    assert(this->FullFormForwardMigration.size() == GP->PatchNumber);
    for (int patch_to = 0 ; patch_to < GP->PatchNumber ; patch_to++)
    {   
        assert(this->FullFormForwardMigration[patch_to].size() == GP->PatchNumber); // Here patch_to is used as patch_from but just for an assertion

        std::vector<double> effectiveNumberOfMigrants; // use doubles to avoid overflow trouble
        effectiveNumberOfMigrants.resize(GP->PatchNumber);

        
        BackwardMigration[patch_to].resize(0);
        BackwardMigrationIndex[patch_to].resize(0);

        // Sum of migrants
        for (int patch_from = 0 ; patch_from < GP->PatchNumber ; patch_from++)
        {   

            if (
                SSP->malesAndFemales &&
                !computingOriginal &&
                (  
                    CumSumFits[patch_from][0].size() == 0.0 || 
                    CumSumFits[patch_from][1].size() == 0.0
                )
            )
            {
                effectiveNumberOfMigrants[patch_from] = 0;
            } else
            {
                if (this->FullFormForwardMigration[patch_from][patch_to] > 0.0)
                {
                    effectiveNumberOfMigrants[patch_from] = this->ComputeEffectiveNumberOfMigrants(
                                        CumSumFits,
                                        computingOriginal,
                                        patch_to,
                                        patch_from
                                    );
                    assert(effectiveNumberOfMigrants[patch_from] >= 0);
                    //std::cout << "effectiveNumberOfMigrants from "<<patch_from<<" to " << patch_to <<": " << effectiveNumberOfMigrants[patch_from] << "\n";
                    nextGenerationPatchSizes[patch_to] += effectiveNumberOfMigrants[patch_from];
                } else
                {
                    assert(this->FullFormForwardMigration[patch_from][patch_to] == 0);
                    effectiveNumberOfMigrants[patch_from] = 0; // non-sense number
                }
            }
        }

        // Fraction of migrants from patch
        assert(nextGenerationPatchSizes[patch_to] >= 0);
        if (nextGenerationPatchSizes[patch_to] > 0)
        {
            double sumOfBMR = 0.0;
            for (int patch_from = 0 ; patch_from < GP->PatchNumber ; patch_from++)
            {
                assert( this->FullFormForwardMigration[patch_from].size() == GP->PatchNumber );
                if (this->FullFormForwardMigration[patch_from][patch_to] > 0.0)
                {
                    //std::cout << "nextGenerationPatchSizes["<<patch_to<<"] = " << nextGenerationPatchSizes[patch_to] << "\n";
                    double backwardMigrationRate = effectiveNumberOfMigrants[patch_from] / nextGenerationPatchSizes[patch_to];
                    //std::cout << "backwardMigrationRate = " << backwardMigrationRate << "\n";
                    assert( backwardMigrationRate >= 0.0 && backwardMigrationRate <= 1.0 );
                    BackwardMigration[patch_to].push_back(backwardMigrationRate);
                    BackwardMigrationIndex[patch_to].push_back(patch_from);
                    sumOfBMR += backwardMigrationRate;
                } else
                {
                    assert(effectiveNumberOfMigrants[patch_from] == 0);
                }
            }

            if (sumOfBMR != 1.0)
            {
                //std::cout << "sumOfBM = " << sumOfBM << "\n";
                for (int patch_from = 0 ; patch_from < BackwardMigration[patch_to].size() ; patch_from++)
                {
                    BackwardMigration[patch_to][patch_from] /= sumOfBMR;
                }
            } 

           /* for (int patch_from = 0 ; patch_from < BackwardMigration[patch_to].size() ; patch_from++)
            {
                std::cout << "to " << patch_to << " from " << BackwardMigrationIndex[patch_to][patch_from] << ": " << BackwardMigration[patch_to][patch_from] << "\n";
            }*/

            // print for debugging purpose
            /*for (int patch_from = 0 ; patch_from < SSP->dispersalData.BackwardMigration[patch_to].size() ; ++patch_from)
            {
                double prob    = SSP->dispersalData.BackwardMigration[patch_to][patch_from];
                assert(SSP->dispersalData.BackwardMigrationIndex[patch_to][patch_from] == patch_from);
                std::cout << "before reordering... from " << patch_from << " to " << patch_to << ": " << prob << "\n";
            }  */
                
            // Reorder
            auto idx = reverse_sort_indexes(BackwardMigration[patch_to]);
            reorder(BackwardMigration[patch_to], idx);
            reorder(BackwardMigrationIndex[patch_to], idx);
        }    
    }

    // print for debugging purpose
    /*for (int patch_to = 0 ; patch_to < SSP->dispersalData.BackwardMigration.size() ; ++patch_to)
    {
        for (int fake_patch_from = 0 ; fake_patch_from < SSP->dispersalData.BackwardMigration[patch_to].size() ; ++fake_patch_from)
        {
            double prob    = SSP->dispersalData.BackwardMigration[patch_to][fake_patch_from];
            int patch_from = SSP->dispersalData.BackwardMigrationIndex[patch_to][fake_patch_from];
            std::cout << "after reordering... from " << patch_from << " to " << patch_to << ": " << prob << "\n";
        }  
    }*/

    // Security
    assert(BackwardMigration.size() == GP->PatchNumber);

  /*  std::cout << "nextGenerationPatchSizes: ";
    for (int p = 0 ; p < GP->PatchNumber ; p++)
        std::cout << nextGenerationPatchSizes[p] << " ";
    std::cout << std::endl;*/



    // copy nextGenerationPatchSizes into r to go from double to int
    std::vector<int> r;
    if (SSP->fecundityForFitnessOfOne == -1.0)
    {
        r = SSP->patchCapacity;
    } else
    {
        r.resize(GP->PatchNumber);
        assert(nextGenerationPatchSizes.size() == GP->PatchNumber);
        for (int p = 0 ; p < GP->PatchNumber ; p++)
        {       
            // Stochasticity in patchSize. I would assume std::poisson_distribution<int> is good at dealing with potential overflow
            std::poisson_distribution<int> dist(nextGenerationPatchSizes[p]);
            r[p] = dist(GP->mt);
            

            // Bounds at Carrying Capacity
            assert(SSP->patchCapacity[p] >= 0);
            assert(r[p] >= 0);
            if (r[p] > SSP->patchCapacity[p])
            {
                r[p] = SSP->patchCapacity[p];
            }
        }
    } 
#ifdef DEBUG
    std::cout << "Exits in DispersalData::SetBackwardMigrationAndGetNextGenerationPatchSizes\n";
#endif      

    /*
    std::cout << "r: ";
    for (int p = 0 ; p < GP->PatchNumber ; p++)
        std::cout << r[p] << " ";
    std::cout << std::endl;
    */
    return r;
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
            
        } else if (Mode.compare("OnePatch") == 0 || Mode.compare("NA") == 0) // A single population, no dispersal
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
