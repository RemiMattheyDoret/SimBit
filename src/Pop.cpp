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


Pop Pop::operator=(Pop&& p)
{
    //std::cout << "Pop move assignment\n";
    patches = std::move(p.patches);
    T2LociToCorrect = p.T2LociToCorrect;
    indexFirstMale = p.indexFirstMale;
    maleCurrentSum = p.maleCurrentSum;
    femaleCurrentSum = p.femaleCurrentSum;
    CumSumFits = p.CumSumFits;
    return *this;
}

Pop::Pop(Pop&& p)
{
    //std::cout << "Pop move constructor\n";
    patches = std::move(p.patches);
    T2LociToCorrect = p.T2LociToCorrect;
    indexFirstMale = p.indexFirstMale;
    maleCurrentSum = p.maleCurrentSum;
    femaleCurrentSum = p.femaleCurrentSum;
    CumSumFits = p.CumSumFits;
}

Pop::Pop(const Pop& p)
{
    //std::cout << "Pop copy constructor\n";
    patches = std::move(p.patches);
    T2LociToCorrect = p.T2LociToCorrect;
    indexFirstMale = p.indexFirstMale;
    maleCurrentSum = p.maleCurrentSum;
    femaleCurrentSum = p.femaleCurrentSum;
    CumSumFits = p.CumSumFits;
}

Pop Pop::operator=(const Pop& p)
{
    //std::cout << "Pop copy assignment\n";
    patches = p.patches;
    T2LociToCorrect = p.T2LociToCorrect;
    indexFirstMale = p.indexFirstMale;
    maleCurrentSum = p.maleCurrentSum;
    femaleCurrentSum = p.femaleCurrentSum;
    CumSumFits = p.CumSumFits;
    return *this;
}

Pop::Pop(){}

Pop::Pop(bool ShouldReadPopFromBinary)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'Pop::Pop'\n";
#endif
    if (ShouldReadPopFromBinary)
    {
        // What is below sounds misleading but it is used to indicate to SimulationTracker::Initialization that there is (potentially) genetic variance.
        SSP->T1_Initial_AlleleFreqs_AllZeros = false;
        SSP->T1_Initial_AlleleFreqs_AllOnes = false;


        GP->BinaryFileToRead.open(SSP->readPopFromBinaryPath.c_str(), std::ios::in | std::ios::binary);
        if (!GP->BinaryFileToRead.is_open())
        {
            std::cout << "GP->BinaryFileToRead (" << SSP->readPopFromBinaryPath << ") failed to open! Please note that 'readPopFromBinary' should NOT be relative to 'GeneralPath' but should be a root path.\n";
            abort();
        }
            

        //.read each patch (each patch starts with its PatchSize)
        for (int patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
        {
            patches.push_back(Patch(patch_index, ShouldReadPopFromBinary));
        }

        // Test that we've reached the end
        {
            if (GP->BinaryFileToRead.eof())
            {
                std::cout << "There seems to have not enough data in the binary file '" << SSP->readPopFromBinaryPath << "', given the input parameters.\n";
                abort();   
            }

            char nothing;
            GP->BinaryFileToRead >> nothing; // sets EOF to true if there was nothing left to read

            if (!GP->BinaryFileToRead.eof())
            {
                std::cout << "There seems to have more data in the binary file '" << SSP->readPopFromBinaryPath << "', than expected from the input parameters.\n";
                int nbExtraBytes = 1;
                while (GP->BinaryFileToRead >> nothing)
                {
                    nbExtraBytes++;
                    if (nbExtraBytes > 99999)
                    {
                        std::cout << "There were more than " << nbExtraBytes << " extra bytes!\n";
                        abort();   
                    }
                }
                std::cout << "There were " << nbExtraBytes << " extra bytes!\n";
                abort();
            }
        }
            

        // close
        GP->BinaryFileToRead.close();
    } else
    {
        for (int patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
        {
            patches.push_back(Patch(patch_index, 'A'));
        }
    }
    GP->saveSSPPatchSize_toGP();
    assert(patches.size() == GP->PatchNumber);

    // I won't need these anymore
    std::vector<std::vector<int>>().swap(SSP->funkyMathForQuasiRandomT1AllFreqInitialization);
    std::vector<std::vector<int>>().swap(SSP->funkyMathForQuasiRandomT56AllFreqInitialization);
    // std::map<std::string, Individual>().swap(SSP->individualTypes); // Will need these for resetGenetics eventB
    std::vector<std::vector<std::string>>().swap(SSP->individualTypesForInitialization);
    std::vector<std::vector<double>>().swap(SSP->T1_Initial_AlleleFreqs);
}


void Pop::PrintBinaryFile()
{

    if (outputWriter.isFile(SaveBinaryFile))
    {
        auto& file = outputWriter.get_OutputFiles(SaveBinaryFile)[0];
        if (file.isTime())
        {
            ////// Seed
            if (SSP->speciesIndex == (GP->nbSpecies-1))
            {
                file.openForSeed(); // opens a different path only for the seed
                file.writeBinary(GP->mt);
                file.close();
            }

            ////// Population

            // open
            file.open(); // reopen a different path for saving the population

            // Write each patch (each patch will start with PatchSize info)
            for (int patch_index = 0 ; patch_index < GP->PatchNumber; patch_index++)
            {
                this->getPatch(patch_index).PrintBinaryFile(file, patch_index);
            }

            // close
            file.close();
        }
    }
}

void Pop::AddPatch(Patch& newPatch)
{
    patches.push_back(newPatch); // copy
}
void Pop::RemoveLastPatch()
{
    patches.pop_back(); // Will call destructor of Patch
    patches.shrink_to_fit();
}

int Pop::getNbPatches()
{
    return patches.size();
}

Patch& Pop::getPatch(const int& patch_index)
{
    //std::cout << "patch_index = " << patch_index << "   | patches.size() = " << patches.size() << std::endl;
    //assert(patch_index < patches.size());
    /*if (patch_index >= patches.size())
    {
        std::cout << "patch_index = " << patch_index << "  patches.size() = " << patches.size() << "\n";
    }*/
    return patches[patch_index];
}


void Pop::prepareNextGenerationAndIndexFirstMale()
{
    // Used in CalculateFitnessForNextGeneration
    maleCurrentSum   = 0.0; 
    femaleCurrentSum = 0.0; 

    // resize
    CumSumFits.resize(GP->PatchNumber); // This should be true as nothing has change recently
    if (SSP->malesAndFemales) indexFirstMale.resize(GP->PatchNumber);


    // More resizing + indexFirstMale
    for ( int patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index )
    {
        // CumSumFitsNextGeneration
        if (SSP->malesAndFemales)
        {
            CumSumFits[patch_index].resize(2);
            CumSumFits[patch_index][0].resize(0);
            CumSumFits[patch_index][1].resize(0);
        } else
        {
            CumSumFits[patch_index].resize(1);
            CumSumFits[patch_index][0].resize(0);
            CumSumFits[patch_index][0].reserve(SSP->TotalpatchSize);
        }


        // indexFirstMale
        if (SSP->malesAndFemales)
        {
            indexFirstMale[patch_index] = (int) (SSP->sexRatio * (double) SSP->patchSize[patch_index] + 0.5);
            assert(indexFirstMale[patch_index] >= 0 && indexFirstMale[patch_index] <= SSP->patchSize[patch_index]);
        }
    }
}




double Pop::CalculateFitnessForNextGeneration(Individual& Offspring, int patch_index, int ind_index)
{
    assert(CumSumFits.size() > patch_index);
    assert(CumSumFits[0].size() > 0);
    
    double w = Offspring.CalculateFitness(patch_index);
    
    if (SSP->malesAndFemales && ind_index >= indexFirstMale[patch_index])
    {
        maleCurrentSum += w;
        CumSumFits[patch_index][1].push_back(maleCurrentSum);
    } else
    {
        femaleCurrentSum += w;
        CumSumFits[patch_index][0].push_back(femaleCurrentSum);
    }

    return w;
}


void Pop::checkIfCumSumFitsIsNotTooSmall(int patch_index)
{
    if (SSP->fecundityForFitnessOfOne == -1.0)
    {
        if (SSP->malesAndFemales)
        {
            int nbFemales = indexFirstMale[patch_index];
            int nbMales   = SSP->patchSize[patch_index] - indexFirstMale[patch_index];

            if (nbFemales <= 0)
            {
                std::cout << "In patch " << patch_index << ", there are " << nbFemales << " females (or hermaphrodites). When fecundity (option --fec) is set to -1, SimBit enforces that there is at least one female and one male per patch (nbMales = "<<nbMales<<" nbFemales = "<<nbFemales<<"). (this error message was written at a time when sex ratio is deterministic. If at some point SimBit evovles to have stochasticity in the sex ratio, then I should remove this error message and deal with cases where the sex ratio is 0 or 1).\n";
                abort();
            }

            if (nbMales <= 0)
            {
                std::cout << "In patch " << patch_index << ", there are " << nbMales << " males. When fecundity (option --fec) is set to -1, SimBit enforces that there is at least one female and one male per patch (nbMales = "<<nbMales<<" nbFemales = "<<nbFemales<<"). (this error message was written at a time when sex ratio is deterministic. If at some point SimBit evovles to have stochasticity in the sex ratio, then I should remove this error message and deal with cases where the sex ratio is 0 or 1).\n";
                abort();
            }
        }
        
        assert(CumSumFits.size() > patch_index);
        assert(CumSumFits[patch_index].size() > 0);
        assert(CumSumFits[patch_index][0].size() > 0);
        if ((CumSumFits[patch_index][0].back() / CumSumFits[patch_index][0].size()) < 0.00000000000001)
        {
            std::cout << "The average fitness of the females (or of the hermaphrodites) in patch " << patch_index << " is " << CumSumFits[patch_index][0].back() / CumSumFits[patch_index][0].size() << ". Selection appears too strong and round off error could become non-negligible! You might want to check your fitness related parameters as well as the mutation rates. As a reminder the fitness effect among loci is multiplicative. If the fecundity was not set to -1.0, this message would not have appeared and the patch size would have simply dropped to zero." << std::endl;
            if (GP->CurrentGeneration == 1)
            {
                std::cout << "As the error message poped up at the first generation, you should probably check the initial conditions. Maybe you start with 'AllOnes' with selection against the the '1' variant (selection is against the '1' variant when making the 'Multiplicity' assumption). It is also possible that you start the simulation with a fixed lethal variant." << std::endl;
            }
            abort();
        }
    }
}  


void Pop::CalculateFitnesses()
{
#ifdef CALLENTRANCEFUNCTIONS
    std::cout << "Enters in Pop::CalculateFitnesses\n";
#endif      

    ///////////////////////////////////////////////////
    ///// Exit if nothing needs to be calculated ////// // Should be a common case
    ///////////////////////////////////////////////////

    // Start with common special case where we know for sure that all fitnesses are already correct and that there is no need 
    if (!SSP->malesAndFemales) // no males index to set
    {
        if (GP->startAtGeneration != GP->CurrentGeneration && (GP->__GenerationChange.size() == 1 || std::find(GP->__GenerationChange.begin(), GP->__GenerationChange.end(), GP->CurrentGeneration) != GP->__GenerationChange.end())) // no recent updates in population
        {
            if (!this->hasCrazyResettingHappened) // No recent weird genetic change
            {
                // No need to compute fitness as we should already know them all!
                // CumSumFitsNextGeneration must get ready to get news values
                
                assert(CumSumFits.size() == GP->PatchNumber); // security
                for ( int patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index ) checkIfCumSumFitsIsNotTooSmall(patch_index); // To avoid round off errors
                return;
            }
        }
    }
    


    //////////////////////////////////////////////////////////
    ///// Compute CumSumFits for the current generation ////// // Should be a rare case
    //////////////////////////////////////////////////////////
    CumSumFits.resize(GP->PatchNumber);
    
    for ( int patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index )
    {
        // Resizing (necessary for the first generation at least)        
        if (SSP->malesAndFemales)
        {
            CumSumFits[patch_index].resize(2);
            CumSumFits[patch_index][0].resize(0);
            CumSumFits[patch_index][1].resize(0);
        } else
        {
            CumSumFits[patch_index].resize(1);
            CumSumFits[patch_index][0].resize(0);
            CumSumFits[patch_index][0].reserve(SSP->TotalpatchSize);
        }


        // Compute fitnesses
        femaleCurrentSum = 0.0; // hermaphrodite are also called females here
        maleCurrentSum = 0.0;
        Patch& patch = this->getPatch(patch_index);
        for ( int ind_index=0 ; ind_index < SSP->patchSize[patch_index] ; ++ind_index )
        {
            double w = patch.getInd(ind_index).CalculateFitness(patch_index);
            //std::cout << "SSP->malesAndFemales " << SSP->malesAndFemales << "\n";
            //std::cout << "ind_index " << ind_index << " indexFirstMale[patch_index] " << indexFirstMale[patch_index] << "\n";
            if (SSP->malesAndFemales && ind_index >= indexFirstMale[patch_index])
            {
                maleCurrentSum += w;
                CumSumFits[patch_index][1].push_back(maleCurrentSum);
            } else
            {
                femaleCurrentSum += w;
                CumSumFits[patch_index][0].push_back(femaleCurrentSum);
            }
        }

        // Set them back to zero (I don't think it is needed though)
        femaleCurrentSum = 0.0;
        maleCurrentSum = 0.0;

        // Security to avoid round off error in extreme parameter sets
        checkIfCumSumFitsIsNotTooSmall(patch_index);

    }
    
    //std::cout << "sum of fitness = " << CumSumFits[0][0].back() << "\n";
#ifdef CALLENTRANCEFUNCTIONS
    std::cout << "Exits in Pop::CalculateFitnesses\n";
#endif      
}

int Pop::SelectionParent(int patch_from, int sex)
{
#ifdef CALLENTRANCEFUNCTIONS
    std::cout << "Enters in Pop::SelectionParent\n";
#endif        
    assert(this->CumSumFits[patch_from].size() > sex);
    
    int parent_index;
    if (SSP->isAnySelection && SSP->selectionOn != 1) // selection on fertility (maybe on both)
    {

        //std::cout << "this->CumSumFits["<<patch_from<<"]["<<sex<<"].size() = " << this->CumSumFits[patch_from][sex].size() << "\n";
        //assert(this->CumSumFits[patch_from][sex].back() > 0.0);
        std::uniform_real_distribution<double> runiform_double_0andSumOfFit(0.0, this->CumSumFits[patch_from][sex].back()); // CumSumFits must be of patchSize length, not of patchCapacity length
        double rnd = runiform_double_0andSumOfFit(GP->mt);
        // binary search
        std::vector<double>::iterator high = std::upper_bound(CumSumFits[patch_from][sex].begin(), CumSumFits[patch_from][sex].end(), rnd);
        
        // Get the index from the iterator
        if (SSP->malesAndFemales)
        {
            int parent_index_inSex = distance(CumSumFits[patch_from][sex].begin(), high);
            //std::cout << "parent_index_inSex = " << parent_index_inSex << "\n";
            if (sex == 0)
            {
                parent_index = parent_index_inSex;
            } else
            {
                parent_index = this->indexFirstMale[patch_from] + parent_index_inSex;
            }
        } else
        {
            parent_index = distance(CumSumFits[patch_from][sex].begin(), high);   
        }
    } else
    {
        // no selection on fertility
        
        if (SSP->malesAndFemales)
        {
            if (sex == 0)
            {
                std::uniform_int_distribution<int> runiform_int_0andNbInds(0,indexFirstMale[patch_from]-1);
                parent_index = runiform_int_0andNbInds(GP->mt);
            } else
            {
                std::uniform_int_distribution<int> runiform_int_0andNbInds(indexFirstMale[patch_from],SSP->patchSize[patch_from]-1);
                parent_index = runiform_int_0andNbInds(GP->mt);
            }
        } else
        {
            std::uniform_int_distribution<int> runiform_int_0andNbInds(0,SSP->patchSize[patch_from]-1);
            parent_index = runiform_int_0andNbInds(GP->mt);   
        }
    }

    
    //std::cout << "parent_index = " << parent_index << "\n";
    //std::cout << "SSP->patchSize[patch_from] = " << SSP->patchSize[patch_from] << "\n";
    assert(parent_index >= 0 && parent_index < SSP->patchSize[patch_from]);
#ifdef CALLENTRANCEFUNCTIONS
    std::cout << "Exits in Pop::SelectionParent\n";
#endif       
    return parent_index;
}


int Pop::SelectionOriginPatch(size_t patch_to)
{
#ifdef CALLENTRANCEFUNCTIONS
    std::cout << "Enters in Pop::SelectionOriginPatch\n";
#endif    
    if (GP->PatchNumber == 1)
    {
        return 0;
    }
    int patch_from = -1;
    
    double rnd = GP->random_0and1(GP->mt); // random between 0 and 1

    for (int fake_patch_from = 0 ; fake_patch_from < SSP->dispersalData.BackwardMigration[patch_to].size() ; ++fake_patch_from)
    {
        double probability = SSP->dispersalData.BackwardMigration[patch_to][fake_patch_from];
        
        if (rnd < probability)
        {
            patch_from = SSP->dispersalData.BackwardMigrationIndex[patch_to][fake_patch_from];
            break;
        }
        rnd -= probability;
    }   
 
    
    assert(patch_from >= 0 && patch_from < GP->PatchNumber);
    //std::cout << "from patch " << patch_from << " to patch " << patch_to << "\n";
#ifdef CALLENTRANCEFUNCTIONS
    std::cout << "Exits in Pop::SelectionOriginPatch\n";
#endif  
    return patch_from;
}


void Pop::updatePops(Pop& pop1, Pop& pop2, int speciesIndex, std::vector<int> previousPatchSizes)
{
    pop1.hasCrazyResettingHappened = true;
    pop2.hasCrazyResettingHappened = true;

    //std::cout << "Enters in Pop::updatePops\n";
    // Change number of patches
    // If there are more patches than before, then the first patch is simply copied several times.
    // If there are fewer patches than before, then last patches are just removed.
    //std::cout << "GP->CurrentGeneration = " << GP->CurrentGeneration << "\n";
    //std::cout << "GP->allSpeciesPatchSizes[0][0] = " << GP->allSpeciesPatchSizes[0][0] << "\n";
    assert(GP->PatchNumber > 0);
    //pop1.CumSumFits.resize(GP->PatchNumber); // I don't think this is needed though
    //pop2.CumSumFits.resize(GP->PatchNumber); // I don't think this is needed though
    if (pop1.getNbPatches() > GP->PatchNumber) // If it needs to remove patches
    {
        int NbPatchesToRemove = pop1.getNbPatches() - GP->PatchNumber;
        for (int i = 0 ; i < NbPatchesToRemove ; ++i)
        {
            pop1.RemoveLastPatch();    // Will call destructor of Patch
            pop2.RemoveLastPatch(); // Will call destructor of Patch
        }
    } else if (pop1.getNbPatches() < GP->PatchNumber) // if it needs to add patches
    {
        assert(GP->PatchNumber - pop1.getNbPatches() > 0);
        int nbPatchesToAdd = GP->PatchNumber - pop1.getNbPatches();
        Patch NewPatch;  // empty patch
        for (int newPatch_index = 0 ; newPatch_index < nbPatchesToAdd ; ++newPatch_index)
        {
            pop1.AddPatch(NewPatch);  // Will copy in 'AddPatch'
            pop2.AddPatch(NewPatch);  // Will copy in 'AddPatch'
        }
    }
    assert(pop1.getNbPatches() == GP->PatchNumber);
    assert(pop2.getNbPatches() == GP->PatchNumber);


    // Change patch carrying capacity and patchSize
    // If there are more individuals than before in a given patch, then the individuals are sampled (with replacement of course) from the same patch (which mifght well be just a copy of patch 0) and cloned until new carrying capacity is reached.
    // If there are fewer individuals than before in a given patch, random individuals are removed.
    // As a reminder, the actual number of individuals per patch is always the carrying capacity. The patch size may differ from it but not the number of individuals kept in memory. This is to avoid freeing and reallcoation of memory.
    assert(SSP->patchCapacity.size() == GP->PatchNumber);
    assert(SSP->patchSize.size() == GP->PatchNumber);
    assert(previousPatchSizes.size() == GP->PatchNumber);
    assert(GP->allSpeciesPatchSizes.size() == GP->PatchNumber);

    for (int patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
    {
        assert(pop1.getPatch(patch_index).getpatchCapacity() == pop2.getPatch(patch_index).getpatchCapacity());
        assert(GP->allSpeciesPatchSizes[patch_index].size() == GP->nbSpecies);
        //std::cout << "In Pop::updatePops, SSP->patchCapacity["<<patch_index<<"] = "<<SSP->patchCapacity[patch_index]<<"\n";
        //std::cout << "In Pop::updatePops, pop1.getPatch(patch_index).getpatchCapacity() = "<<pop1.getPatch(patch_index).getpatchCapacity()<<"\n";
            

        if (pop1.getPatch(patch_index).getpatchCapacity() < SSP->patchCapacity[patch_index]) // If the carrying capacity increased.
        {        
            // even if fec == -1, I still have to add individuals (even if I won't make sense of them)
            // Figure out what patch the individuals should be sampled from            
            int patch_index_ToSampleFrom = SSP->selectNonEmptyPatch(patch_index, previousPatchSizes, true);
            assert(patch_index_ToSampleFrom >= 0 && patch_index_ToSampleFrom < GP->PatchNumber);
            
            if (patch_index_ToSampleFrom != -1)
            {
                assert(previousPatchSizes[patch_index_ToSampleFrom] > 0);
                //std::cout << "patch_index_ToSampleFrom = " << patch_index_ToSampleFrom << "\n";

                int NbIndsToAdd = SSP->patchCapacity[patch_index] - pop1.getPatch(patch_index).getpatchCapacity();
                for (int fake_ind_index = 0 ; fake_ind_index < NbIndsToAdd ; ++fake_ind_index)
                {
                    int ind_index = fake_ind_index % previousPatchSizes[patch_index_ToSampleFrom];
                    
                    assert(ind_index >= 0 && ind_index < previousPatchSizes[patch_index_ToSampleFrom]);
                    //std::cout << "RandomIndividual = " << RandomIndividual << "\n";
                    //assert(RandomIndividual >= 0 && RandomIndividual < pop1.getPatch(patch_index_ToSampleFrom).getpatchCapacity());
                    Individual NewIndividual = pop1.getPatch(patch_index_ToSampleFrom).getInd(ind_index); // Get the individual from this patch. DO NOT TAKE A REFERENCE ONLY. APPARENTLY AddIndividual DOES NOT COPY AS I EXPECTED IT TO.
                    pop1.getPatch(patch_index).AddIndividual(NewIndividual); // Is copied in 'AddIndividual'
                    pop2.getPatch(patch_index).AddIndividual(NewIndividual); // Is copied in 'AddIndividual'
                }
            }

        } else if (pop1.getPatch(patch_index).getpatchCapacity() > SSP->patchCapacity[patch_index]) // If the carrying capacity decreased.
        {
            int NbIndsToRemove = pop1.getPatch(patch_index).getpatchCapacity() - SSP->patchCapacity[patch_index];
            if (SSP->patchSize[patch_index] > SSP->patchCapacity[patch_index])
            {
                SSP->patchSize[patch_index] = SSP->patchCapacity[patch_index];
            }
            for (int ind_index = 0 ; ind_index < NbIndsToRemove ; ++ind_index)
            {
                pop1.getPatch(patch_index).removeLastIndividual(); // call destructor of Individual
                pop2.getPatch(patch_index).removeLastIndividual(); // call destructor of Individual
            }
        }
    }

    // Set all fitnesses to -1 to force recalculation (Some other alternative would have better performance but I am assuming that a user won't ask for many parameters changes through time)$
    for (int patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
    {
        for (int ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ++ind_index)
        {
            for (int haplo_index = 0 ; haplo_index < SSP->ploidy ; ++haplo_index)
            {
                if (SSP->T1_nbBits)
                {
                    if (SSP->T1_isMultiplicitySelection)
                    {
                        pop1.getPatch(patch_index).getInd(ind_index, patch_index).getHaplo(haplo_index).setAllW_T1(-1);
                        pop2.getPatch(patch_index).getInd(ind_index, patch_index).getHaplo(haplo_index).setAllW_T1(-1);
                    } 
                }
                    
                if (SSP->T2_nbChars)
                {
                    if (SSP->T2_isSelection)
                    {
                        pop1.getPatch(patch_index).getInd(ind_index, patch_index).getHaplo(haplo_index).setAllW_T2(-1);
                        pop2.getPatch(patch_index).getInd(ind_index, patch_index).getHaplo(haplo_index).setAllW_T2(-1);
                    }
                }

                if (SSP->T56sel_nbBits)
                {
                    assert(SSP->T56_isSelection);
                    if (SSP->T56_isMultiplicitySelection)
                    {
                        pop1.getPatch(patch_index).getInd(ind_index, patch_index).getHaplo(haplo_index).setAllW_T56(-1);
                        pop2.getPatch(patch_index).getInd(ind_index, patch_index).getHaplo(haplo_index).setAllW_T56(-1);
                    }
                }
            }
        }
    }
    //std::cout << "Exits in Pop::updatePops\n";
}




void Pop::addT2LocusToCorrect(int T2Locus) // is a static
{
    assert(SSP != nullptr);
    assert(T2Locus >= 0 && T2Locus < SSP->T2_nbChars);
    if
    (
        std::find(
            T2LociToCorrect.begin(),
            T2LociToCorrect.end(),
            T2Locus
        ) == T2LociToCorrect.end()
    )
    {
        T2LociToCorrect.push_back(T2Locus);
        std::sort(
            T2LociToCorrect.begin(),
            T2LociToCorrect.end()
        );
    }
        

}

int Pop::correctT2Loci()
{
    if (T2LociToCorrect.size() > 0)
    {
        assert(SSP != nullptr);
        if (SSP->fecundityForFitnessOfOne != -1)
        {
            std::cout << "A T2 block reached more than 250 mutations. SimBit cannot reset the values of the entire population without affecting the fecundity of the individuals of this population. Hence, this correction can only be carried out if fecundity is set to -1 (infinite fecundity). Fecundity was set to " << SSP->fecundityForFitnessOfOne << " for species " << SSP->speciesName << ". SimBit will abort.\n";
            abort();
        }

        // Find out who has the lowest number of mutations for all loci that need correction
        std::vector<unsigned char> lowestNbMuts;
        lowestNbMuts.resize(T2LociToCorrect.size(), UCHAR_MAX);
        assert(lowestNbMuts.size() == T2LociToCorrect.size());

        for (int patch_index = 0 ; patch_index < GP->PatchNumber ; patch_index++)
        {
            assert(SSP->patchSize.size() == GP->PatchNumber);
            Patch& patch = this->getPatch(patch_index);
            for (int ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ind_index++)
            {
                Individual& individual = patch.getInd(ind_index);
                for (int haplo_index = 0 ; haplo_index < 2 ; haplo_index++)
                {
                    // Loop through loci that need to be corrected
                    for (int i = 0 ; i < T2LociToCorrect.size() ; i++)
                    {
                        if (T2LociToCorrect[i] < 0 || T2LociToCorrect[i] >= SSP->T2_nbChars)
                        {
                            printf("Internal error in Pop::correctT2Loci. Supposed to look for T2 block index %d but there are only %d T2 blocks in the whole genome.\n", T2LociToCorrect[i], SSP->T2_nbChars);
                            abort();
                        }
                        auto T2Allele = individual.getHaplo(haplo_index).getT2_Allele(T2LociToCorrect[i]);
                        if (T2Allele < 0 || T2Allele > UCHAR_MAX)
                        {
                            printf("Internal error in Pop::correctT2Loci. T2Allele of the T2 block index %d is %d.\n", T2LociToCorrect[i], T2Allele);
                            abort();
                        }
                        lowestNbMuts[i] = std::min(lowestNbMuts[i],T2Allele);
                    }
                }
            }
        }

        // Make sure the correction will be helpful
        for (int i = 0 ; i < T2LociToCorrect.size() ; i++)
        {
            if (lowestNbMuts[i] < 3)
            {
                assert(lowestNbMuts[i] >= 0);
                std::cout << "In Pop::correctT2Loci, the locus " << T2LociToCorrect[i] << " for species " << SSP->speciesName << "reached too a hihg number of mutations and needed to be corrected (already a slow process). But when correcting it, it appears that there is a haplotype on which the number of mutation at this locus is lower than 3. With such variance, the correction is not very helpful and the simulations become very slow. SimBit will therefore abort. Please consider making sure that the number of mutations per T2 block is lower and that the variance in the number of mutations in T2 block is lower too because a T2 block cannot keep track of that much information (a T2 block is only 8 bits).\n ";
                abort();
            }
        }

        // Correct
        for (int patch_index = 0 ; patch_index < GP->PatchNumber ; patch_index++)
        {
            assert(SSP->patchSize.size() == GP->PatchNumber);
            Patch& patch = this->getPatch(patch_index);
            for (int ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ind_index++)
            {
                Individual& individual = patch.getInd(ind_index);
                for (int haplo_index = 0 ; haplo_index < 2 ; haplo_index++)
                {
                    // Loop through loci that need to be corrected
                    for (int i = 0 ; i < T2LociToCorrect.size() ; i++)
                    {
                        auto oldValue = individual.getHaplo(haplo_index).getT2_Allele(T2LociToCorrect[i]);

                        auto newValue = oldValue - lowestNbMuts[i];
                        if (newValue < 0 || newValue > UCHAR_MAX)
                        {
                            printf("Internal error in Pop::correctT2Loci. newValue = %d  oldValue = %d lowestNbMuts[%d] = %d. Also note that UCHAR_MAX = %d\n", newValue, oldValue, i, lowestNbMuts[i], UCHAR_MAX);
                            abort();
                        }
        
                        individual.getHaplo(haplo_index).setT2_Allele(T2LociToCorrect[i], newValue);
                    }
                }
            }
        }

        // Clear things up
        T2LociToCorrect.clear();

        
        return 1;
    }
    return 0;
}

std::vector<int> Pop::findWhatMustBeToggledAndUpdateFlipped(const std::vector<unsigned>& popFreqs, std::vector<unsigned int>& flipped, const int& nbLoci,  const double freqThreshold)
{
    std::vector<int> lociToToggle;

    // look for what is to toggle
    assert(flipped.size() <= nbLoci);
    //std::cout << "toggle: ";

    std::vector<bool> isFlipped = inverseWhich(flipped, nbLoci);
    assert(isFlipped.size() == nbLoci);

    for (size_t locus = 0 ; locus < nbLoci ; ++locus)
    {
        double relFreq = (double) popFreqs[locus] / (2.0 * (double) SSP->TotalpatchSize);
        assert(relFreq >= 0.0 && relFreq <= 1.0);

        // Wow hard part because frequency is already reversed in computeT56ntrlFrequencies or computeT56selFrequencies
        if (isFlipped[locus] && relFreq < (1-freqThreshold))
        {
            isFlipped[locus] = !isFlipped[locus];
            lociToToggle.push_back(locus);
            //std::cout << locus << "(unflip) ";
        } else if (!isFlipped[locus] && relFreq > freqThreshold)
        {
            isFlipped[locus] = !isFlipped[locus];
            lociToToggle.push_back(locus);
            //std::cout << locus << "(flip) ";
        }
    }
    //std::cout << "\n";
    flipped = whichUnsignedInt(isFlipped); // Update the flipped information

    return lociToToggle;                    // return what loci must be toggled

    /*
    if (T56ntrlLociToToggle.size() > 0)
    {
        std::cout << "T56ntrlLociToToggle :";
        for (auto& elem : T56ntrlLociToToggle)
            std::cout << elem << " ";
        std::cout << "\n";
    }
    */
}

void Pop::toggleT56FixedMutations()
{
    //////////////
    /// T5ntrl ///
    //////////////
    std::vector<int> T56ntrlLociToToggle;
    if (SSP->T5ntrl_nbBits)
    {
        std::vector<unsigned> popFreqs = this->computeT56ntrlFrequencies();
        T56ntrlLociToToggle = findWhatMustBeToggledAndUpdateFlipped(popFreqs, SSP->T5ntrl_flipped, SSP->T5ntrl_nbBits, SSP->T56ntrl_frequencyThresholdForFlippingMeaning);
    } else if (SSP->T6ntrl_nbBits)
    {
        std::vector<unsigned> popFreqs = this->computeT56ntrlFrequencies();
        auto flipped = SSP->T6ntrl_flipped.toVector();
        T56ntrlLociToToggle = findWhatMustBeToggledAndUpdateFlipped(popFreqs, flipped, SSP->T6ntrl_nbBits, SSP->T56ntrl_frequencyThresholdForFlippingMeaning);
        SSP->T6ntrl_flipped = CompressedSortedDeque(flipped, SSP->T6ntrl_nbBits);
    }


    //////////////
    /// T5sel ////
    //////////////
    /*std::vector<int> T56selLociToToggle;

    if (SSP->T5sel_nbBits)
    {
        assert(!SSP->T56sel_compress);
        std::vector<unsigned> popFreqs = this->computeT56selFrequencies();
        T56selLociToToggle = findWhatMustBeToggledAndUpdateFlipped(popFreqs, SSP->T5sel_flipped, SSP->T5sel_nbBits, SSP->T56sel_frequencyThresholdForFlippingMeaning);
    } else if (SSP->T6sel_nbBits)
    {
        assert(SSP->T56sel_compress);
        std::vector<unsigned> popFreqs = this->computeT56selFrequencies();
        auto flipped = SSP->T6sel_flipped.toVector();
        T56selLociToToggle = findWhatMustBeToggledAndUpdateFlipped(popFreqs, flipped, SSP->T6sel_nbBits, SSP->T56sel_frequencyThresholdForFlippingMeaning);
        SSP->T6sel_flipped = CompressedSortedDeque(flipped, SSP->T6sel_nbBits);
    }*/


    ///////////////
    /// Toggle ////
    ///////////////
    if ((T56ntrlLociToToggle.size()) > 0)
    {
        std::vector<int> empty;
        this->toggleT56LociFromEveryone(T56ntrlLociToToggle, empty);
        //                                                    ^ This spot was reserve for sel loci but this version does not toggle sel loci anymore!
    }

    //std::cout << "Toggled " << T56ntrlLociToToggle.size() << "+" << T56selLociToToggle.size() << " loci\n";
}

void Pop::toggleT56LociFromEveryone(std::vector<int>& T56ntrlLociToToggle, std::vector<int>& T56selLociToToggle)
{
    assert(SSP->Habitats.size() == patches.size() );
    for (size_t patchIndex = 0 ; patchIndex < patches.size() ; ++patchIndex)
        patches[patchIndex].toggleT56LociFromEveryone(T56ntrlLociToToggle, T56selLociToToggle, SSP->Habitats[patchIndex]);
}


void Pop::toggleT56MutationsIfNeeded()
{
    //std::cout << "Enters in Pop::toggleT56FixedNtrlMutationsIfNeeded\n";
    if (SSP->T56_toggleMutsEveryNGeneration_nextGeneration == GP->CurrentGeneration && SSP->TotalpatchSize > 0)
    {
        SSP->T56_toggleMutsEveryNGeneration_nextGeneration += SSP->T56_toggleMutsEveryNGeneration;
        
        this->toggleT56FixedMutations();
    }
    //std::cout << "Exits Pop::toggleT56FixedNtrlMutationsIfNeeded\n";
}



std::vector<std::vector<unsigned>> Pop::computePatchSpecificT56ntrlFrequencies()
{
    // std::cout << "Enters in Pop::computePatchSpecificT5ntrlFrequencies()\n";
    // Build obsFreqs
    std::vector<std::vector<unsigned>> obsFreqs(GP->PatchNumber); // obsFreqs[patch_index][locus]

    if (SSP->T56ntrl_nbBits)
    {
        // count
        for (size_t patch_index = 0 ; patch_index < GP->PatchNumber ; patch_index++)
        {        
            obsFreqs[patch_index].resize(SSP->T56ntrl_nbBits, 0.0);
            for (size_t ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ind_index++)
            {
                for (size_t haplo_index = 0 ; haplo_index < 2 ; haplo_index++)
                {
                    Haplotype& haplo = this->getPatch(patch_index).getInd(ind_index).getHaplo(haplo_index);

                    if (SSP->T56ntrl_compress)
                    {
                        auto it = haplo.T6ntrl_AllelesBegin();
                        auto itEnd = haplo.T6ntrl_AllelesEnd();
                        for (; it != itEnd ; ++it)
                        {
                            auto value = *it;
                            assert(value >= 0 && value < SSP->T6ntrl_nbBits);
                            obsFreqs[patch_index][value]++;
                        }
                    } else
                    {
                        auto it = haplo.T5ntrl_AllelesBegin();
                        auto itEnd = haplo.T5ntrl_AllelesEnd();
                        for (; it != itEnd ; ++it)
                        {
                            auto value = *it;
                            assert(value >= 0 && value < SSP->T5ntrl_nbBits);
                            obsFreqs[patch_index][value]++;
                        }
                    }
                }
            }
        }
    }
        
    //std::cout << "About to exit in Pop::computePatchSpecificT5ntrlFrequencies()\n";
    return obsFreqs;
}


std::vector<std::vector<unsigned>> Pop::computePatchSpecificT56selFrequencies()
{    
    //std::cout << "Enters in Pop::computePatchSpecificT5selFrequencies()\n";
    // Build obsFreqs
    std::vector<std::vector<unsigned>> obsFreqs(GP->PatchNumber); // obsFreqs[patch_index][locus]

    if (SSP->T56sel_nbBits)
    {
        for (size_t patch_index = 0 ; patch_index < GP->PatchNumber ; patch_index++)
        {
            // count
            obsFreqs[patch_index].resize(SSP->T56sel_nbBits, 0.0);
            for (size_t ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ind_index++)
            {
                for (size_t haplo_index = 0 ; haplo_index < 2 ; haplo_index++)
                {
                    Haplotype& haplo = this->getPatch(patch_index).getInd(ind_index).getHaplo(haplo_index);


                    if (SSP->T56sel_compress)
                    {
                        auto it = haplo.T6sel_AllelesBegin();
                        auto itEnd = haplo.T6sel_AllelesEnd();
                        for (; it != itEnd ; ++it)
                        {
                            auto value = *it;
                            assert(value >= 0 && value < SSP->T6sel_nbBits);
                            obsFreqs[patch_index][value]++;
                        }
                    } else
                    {
                        auto it = haplo.T5sel_AllelesBegin();
                        auto itEnd = haplo.T5sel_AllelesEnd();
                        for (; it != itEnd ; ++it)
                        {
                            auto value = *it;
                            assert(value >= 0 && value < SSP->T5sel_nbBits);
                            obsFreqs[patch_index][value]++;
                        }
                    }
                }
            }
        }
    }

    //std::cout << "About to exit in Pop::computePatchSpecificT5selFrequencies()\n";
    return obsFreqs;
}

std::vector<std::vector<unsigned>> Pop::computePatchSpecificT56Frequencies()
{
    //std::cout << "enters Pop::computePatchSpecificT56Frequencies\n";

    assert(SSP->T56_nbBits);

    std::vector<std::vector<unsigned>> ntrlFreqs  = this->computePatchSpecificT56ntrlFrequencies();
    std::vector<std::vector<unsigned>> selFreqs = this->computePatchSpecificT56selFrequencies();
    
    assert(SSP->FromT56selLocusToLocus.size() == SSP->T56sel_nbBits);
    assert(SSP->FromT56ntrlLocusToLocus.size() == SSP->T56ntrl_nbBits);
    assert(SSP->FromT56LocusToT56genderLocus.size() == SSP->T56_nbBits);

    std::vector<std::vector<unsigned>> r(GP->PatchNumber);
    for (size_t patch_index = 0 ; patch_index < GP->PatchNumber ; patch_index++)
    {
        size_t selFreqsIndex = 0;
        size_t ntrlFreqsIndex = 0;
        r[patch_index].resize(SSP->T56_nbBits);
        for (size_t T56locus = 0 ; T56locus < SSP->T56_nbBits ; T56locus++)
        {
            //std::cout << "SSP->FromT56LocusToT56genderLocus["<<T5locus<<"].first = " << SSP->FromT56LocusToT56genderLocus[T5locus].first << "\n";
            if (SSP->FromT56LocusToT56genderLocus[T56locus].first)
            {
                r[patch_index][T56locus] = ntrlFreqs[0][ntrlFreqsIndex]; // Why '0' and not 'patch_index'? Because the previous patch index has been removed to keep RAM as low as possible
                ntrlFreqsIndex++;
            } else
            {
                r[patch_index][T56locus] = selFreqs[0][selFreqsIndex]; // Why '0' and not 'patch_index'? Because the previous patch index has been removed to keep RAM as low as possible
                selFreqsIndex++;
            }
        }

        ntrlFreqs.erase(ntrlFreqs.begin());  // Why remove the first element? Because we always refer to the zeroth element above. It helps keep teh RAM usage as low as possible
        selFreqs.erase(selFreqs.begin());    // Why remove the first element? Because we always refer to the zeroth element above. It helps keep teh RAM usage as low as possible
        
        assert(selFreqsIndex == SSP->T56sel_nbBits);
        assert(ntrlFreqsIndex == SSP->T56ntrl_nbBits);
    }
    assert(ntrlFreqs.size() == 0);
    assert(selFreqs.size() == 0);
    assert(r.size() == GP->PatchNumber);
    assert(r[0].size() == SSP->T56_nbBits);

    //std::cout << "About to exit Pop::computePatchSpecificT56Frequencies\n";
    return r;

}




std::vector<unsigned> Pop::computeT56ntrlFrequencies()
{
    assert(SSP->T56ntrl_nbBits == SSP->T5ntrl_nbBits + SSP->T6ntrl_nbBits);
    std::vector<unsigned> obsFreqs(SSP->T56ntrl_nbBits,0); // obsFreqs[patch_index][locus]

    if (SSP->T56ntrl_nbBits)
    {
        // count
        for (size_t patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
        {        
            for (size_t ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ++ind_index)
            {
                for (size_t haplo_index = 0 ; haplo_index < 2 ; ++haplo_index)
                {
                    Haplotype& haplo = this->getPatch(patch_index).getInd(ind_index).getHaplo(haplo_index);

                    if (SSP->T56ntrl_compress)
                    {
                        auto it    = haplo.T6ntrl_AllelesBegin();
                        auto itEnd = haplo.T6ntrl_AllelesEnd();
                        for (; it != itEnd ; ++it)
                        {
                            auto value = *it;
                            assert(value >= 0 && value < SSP->T6ntrl_nbBits);
                            ++obsFreqs[value];
                        }
                    } else
                    {
                        auto it = haplo.T5ntrl_AllelesBegin();
                        auto itEnd = haplo.T5ntrl_AllelesEnd();
                        for (; it != itEnd ; ++it)
                        {
                            auto value = *it;
                            assert(value >= 0 && value < SSP->T5ntrl_nbBits);
                            ++obsFreqs[value];
                        }
                    }
                }
            }
        }
    }

    /*
    assert(obsFreqs.size()== SSP->T56ntrl_nbBits);
    std::cout << "what is polymorphic: ";
    for (size_t locus = 0 ; locus < obsFreqs.size() ; ++locus)
    {
        if (obsFreqs[locus] != 0 && obsFreqs[locus] != 2 * SSP->TotalpatchSize)
            std::cout << locus << " ";
    }
    std::cout << "\n";
    */
        
    //std::cout << "About to exit in Pop::computePatchSpecificT5ntrlFrequencies()\n";
    return obsFreqs;
}


std::vector<unsigned> Pop::computeT56selFrequencies()
{
    // Build obsFreqs
    assert(SSP->T56sel_nbBits == SSP->T5sel_nbBits + SSP->T6sel_nbBits);
    std::vector<unsigned> obsFreqs(SSP->T56sel_nbBits,0); // obsFreqs[patch_index][locus]
    if (SSP->T56sel_nbBits)
    {
        // count
        size_t TotalpatchSize = 0;
        for (size_t patch_index = 0 ; patch_index < GP->PatchNumber ; patch_index++)
        {
            TotalpatchSize += SSP->patchSize[patch_index];
        
            for (size_t ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ind_index++)
            {
                for (size_t haplo_index = 0 ; haplo_index < 2 ; haplo_index++)
                {
                    Haplotype& haplo = this->getPatch(patch_index).getInd(ind_index).getHaplo(haplo_index);

                    if (SSP->T56sel_compress)
                    {
                        auto it = haplo.T6sel_AllelesBegin();
                        auto itEnd = haplo.T6sel_AllelesEnd();
                        for (; it != itEnd ; ++it)
                        {
                            auto value = *it;
                            assert(value >= 0 && value < SSP->T6sel_nbBits);
                            ++obsFreqs[value];
                        }
                    } else
                    {
                        auto it = haplo.T5sel_AllelesBegin();
                        auto itEnd = haplo.T5sel_AllelesEnd();
                        for (; it != itEnd ; ++it)
                        {
                            auto value = *it;
                            assert(value >= 0 && value < SSP->T5sel_nbBits);
                            ++obsFreqs[value];
                        }
                    }
                }
            }
        }  
    }


    return obsFreqs;
}


std::vector<unsigned> Pop::computeT56Frequencies()
{
    std::vector<unsigned> selFreqs = this->computeT56selFrequencies();
    std::vector<unsigned> ntrlFreqs = this->computeT56ntrlFrequencies();

    std::vector<unsigned> r(SSP->T56_nbBits);
    
    size_t selFreqsIndex = 0;
    size_t ntrlFreqsIndex = 0;

    for (size_t locus = 0 ; locus < SSP->T56_nbBits ; locus++)
    {
        if (SSP->FromT56LocusToT56genderLocus[locus].first)
        {
            r[locus] = ntrlFreqs[ntrlFreqsIndex];
            ntrlFreqsIndex++;
        } else
        {
            r[locus] = selFreqs[selFreqsIndex];
            selFreqsIndex++;
        }
        //std::cout << "freq of locus " << T5locus << " = " <<  r[T5locus] << "\n";
    }

    assert(selFreqsIndex == SSP->T56sel_nbBits);
    assert(ntrlFreqsIndex == SSP->T56ntrl_nbBits);


    assert(r.size() == SSP->T56_nbBits);

    return r;
}



std::vector<T1_locusDescription> Pop::listT1PolymorphicLoci()
{
    std::vector<bool> isPoly(SSP->T1_nbBits, false);

    for ( int patch_index = 0 ; patch_index < GP->maxEverPatchNumber ; ++patch_index )
    {
        for (int ind_index=0;ind_index<SSP->patchSize[patch_index];++ind_index)
        {
            for (int haplo_index=0;haplo_index<SSP->ploidy;haplo_index++)
            { 
                auto& haplo = getPatch(patch_index).getInd(ind_index).getHaplo(haplo_index);

                unsigned locus = 0;
                // First bytes
                for (unsigned byte = 0 ; byte < (SSP->T1_nbChars-1) ; ++byte )
                {
                    if (haplo.getT1_char(byte))
                    {
                        for (unsigned bit = 0 ; bit < 8 ; ++bit )
                        {
                            isPoly[locus] = haplo.getT1_Allele(byte, bit);
                            ++locus;
                        }
                    } else
                    {
                        locus+=8;
                    }
                }
                assert(locus == (SSP->T1_nbChars-1) * 8);

                // last byte
                auto lastByte = SSP->T1_nbChars-1;
                if (haplo.getT1_char(lastByte))
                {
                    for (unsigned bit = 0 ; bit < SSP->T1_nbBitsLastByte ; ++bit )
                    {
                        isPoly[locus] = haplo.getT1_Allele(lastByte, bit);
                        ++locus;
                    }
                }
            }
        }
    }
    
    return whichT1_locusDescription(isPoly);
}