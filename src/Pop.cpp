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

Pop::Pop(){
    std::cout << "Internal error: This Pop::Pop() constructor should not be used\n";
    abort();
}

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
            std::cout << "GP->BinaryFileToRead (" << SSP->readPopFromBinaryPath << ") failed to open! Please note that .readPopFromBinary' should NOT be relative to 'GeneralPath' but should be a root path.\n";
            abort();
        }
            

        //.read each patch (each patch starts with its PatchSize)
        for (int patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
        {
            Patch patch(patch_index, ShouldReadPopFromBinary);
            patches.push_back(std::move(patch));
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
                    if (nbExtraBytes > 999999)
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
            Patch patch(patch_index, 'A');
            patches.push_back(std::move(patch));
        }    
    }
    GP->saveSSPPatchSize_toGP();
    assert(patches.size() == GP->PatchNumber);
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
            file.open();

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


void Pop::CalculateFitnesses()
{
#ifdef CALLENTRANCEFUNCTIONS
    std::cout << "Enters in Pop::CalculateFitnesses\n";
#endif      
    // Calculate fitnesses. Should only be called when there is a temporal change
    CumSumFits.resize(GP->PatchNumber);

    indexFirstMale.resize(0);

    //#pragma omp parallel for
    for ( int patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index )
    {
        //std::cout << "SSP->sexRatio = " << SSP->sexRatio << " SSP->patchSize[patch_index] = " << SSP->patchSize[patch_index] << "\n";
        if (SSP->malesAndFemales)
        {
            indexFirstMale.push_back((int) (SSP->sexRatio * (double) SSP->patchSize[patch_index] + 0.5) );
            assert(indexFirstMale[patch_index] >= 0 && indexFirstMale[patch_index] <= SSP->patchSize[patch_index]);
        }
        

        double femaleCurrentSum = 0.0; // hermaphrodite are also called females here
        double maleCurrentSum = 0.0;
        //printf("Hello from thread %d, nthreads %d\n", omp_get_thread_num(), omp_get_num_threads());
        if (SSP->malesAndFemales)
        {
            CumSumFits[patch_index].resize(2);
            CumSumFits[patch_index][0].resize(0);
            CumSumFits[patch_index][1].resize(0);
        } else
        {
            CumSumFits[patch_index].resize(1);
            CumSumFits[patch_index][0].resize(0);
        }
            
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

            if ((femaleCurrentSum / CumSumFits[patch_index][0].size()) < 0.00000000000001)
            {
                std::cout << "The average fitness of the females (or of the hermaphrodites) in patch " << patch_index << " is " << femaleCurrentSum / CumSumFits[patch_index][0].size() << ". Selection appears too strong and round off error could become non-negligible! You might want to check your fitness related parameters as well as the mutation rates. As a reminder the fitness effect among loci is multiplicative. If the fecundity was not set to -1.0, this message would not have appeared and the patch size would have simply dropped to zero." << std::endl;
                if (GP->CurrentGeneration == 1)
                {
                    std::cout << "As the error message poped up at the first generation, you should probably check the initial conditions. Maybe you start with 'AllOnes' with selection against the the '1' variant (selection is against the '1' variant when making the 'Multiplicity' assumption). It is also possible that you start the simulation with a fixed lethal variant." << std::endl;
                }
                abort();
            }
        }
        
    }
#ifdef CALLENTRANCEFUNCTIONS
    std::cout << "Exits in Pop::CalculateFitnesses\n";
#endif      
}

int Pop::SelectionParent(int& patch_from, int sex)
{
#ifdef CALLENTRANCEFUNCTIONS
    std::cout << "Enters in Pop::SelectionParent\n";
#endif        
    assert(this->CumSumFits[patch_from].size() > sex);
    
    int parent_index;
    if (SSP->selectionOn != 1 && (SSP->T1_isSelection || SSP->T1_isEpistasis || SSP->T2_isSelection || SSP->T3_isSelection || SSP->T5_isSelection))
    {
        // selection on fertility

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


int Pop::SelectionOriginPatch(int& patch_to)
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
    //std::cout << "Enters in Pop::updatePops\n";
    // Change number of patches
    // If there are more patches than before, then the first patch is simply copied several times.
    // If there are fewer patches than before, then last patches are just removed.
    //std::cout << "GP->CurrentGeneration = " << GP->CurrentGeneration << "\n";
    //std::cout << "GP->allSpeciesPatchSizes[0][0] = " << GP->allSpeciesPatchSizes[0][0] << "\n";
    assert(GP->PatchNumber > 0);
    pop1.CumSumFits.resize(GP->PatchNumber);
    pop2.CumSumFits.resize(GP->PatchNumber);
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
                if (SSP->T1_isSelection)
                {
                    pop1.getPatch(patch_index).getInd(ind_index, patch_index).getHaplo(haplo_index).setAllW_T1(-1);
                    pop2.getPatch(patch_index).getInd(ind_index, patch_index).getHaplo(haplo_index).setAllW_T1(-1);
                } else
                {
                    for (int fitnessMapIndex = 0 ; fitnessMapIndex < SSP->NbElementsInFitnessMap ; fitnessMapIndex++)
                    {
                        assert(pop1.getPatch(patch_index).getInd(ind_index, patch_index).getHaplo(haplo_index).getW_T1(fitnessMapIndex) == 1.0);
                        assert(pop2.getPatch(patch_index).getInd(ind_index, patch_index).getHaplo(haplo_index).getW_T1(fitnessMapIndex) == 1.0);
                    }
                }

                if (SSP->T2_isSelection)
                {
                    pop1.getPatch(patch_index).getInd(ind_index, patch_index).getHaplo(haplo_index).setAllW_T2(-1);
                    pop2.getPatch(patch_index).getInd(ind_index, patch_index).getHaplo(haplo_index).setAllW_T2(-1);
                } else
                {
                    for (int fitnessMapIndex = 0 ; fitnessMapIndex < SSP->NbElementsInFitnessMap ; fitnessMapIndex++)
                    {
                        assert(pop1.getPatch(patch_index).getInd(ind_index, patch_index).getHaplo(haplo_index).getW_T2(fitnessMapIndex) == 1.0);
                        assert(pop2.getPatch(patch_index).getInd(ind_index, patch_index).getHaplo(haplo_index).getW_T2(fitnessMapIndex) == 1.0);
                    }
                }

                if (SSP->T5_isSelection)
                {
                    pop1.getPatch(patch_index).getInd(ind_index, patch_index).getHaplo(haplo_index).setAllW_T5(-1);
                    pop2.getPatch(patch_index).getInd(ind_index, patch_index).getHaplo(haplo_index).setAllW_T5(-1);
                } else
                {
                    for (int fitnessMapIndex = 0 ; fitnessMapIndex < SSP->NbElementsInFitnessMap ; fitnessMapIndex++)
                    {
                        assert(pop1.getPatch(patch_index).getInd(ind_index, patch_index).getHaplo(haplo_index).getW_T5(fitnessMapIndex) == 1.0);
                        assert(pop2.getPatch(patch_index).getInd(ind_index, patch_index).getHaplo(haplo_index).getW_T5(fitnessMapIndex) == 1.0);
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




