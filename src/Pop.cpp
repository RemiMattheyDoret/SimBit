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

/*
Pop Pop::operator=(Pop&& p)
{
    //std::cout << "Pop move assignment\n";
    patches = std::move(p.patches);
    T2LociToCorrect = p.T2LociToCorrect;
    indexFirstMale = p.indexFirstMale;
    //maleCurrentSum = p.maleCurrentSum;
    //femaleCurrentSum = p.femaleCurrentSum;
    CumSumFits = p.CumSumFits;
    return *this;
}

Pop::Pop(Pop&& p)
{
    //std::cout << "Pop move constructor\n";
    patches = std::move(p.patches);
    T2LociToCorrect = p.T2LociToCorrect;
    indexFirstMale = p.indexFirstMale;
    //maleCurrentSum = p.maleCurrentSum;
    //femaleCurrentSum = p.femaleCurrentSum;
    CumSumFits = p.CumSumFits;
}

Pop::Pop(const Pop& p)
{
    //std::cout << "Pop copy constructor\n";
    patches = p.patches;
    T2LociToCorrect = p.T2LociToCorrect;
    indexFirstMale = p.indexFirstMale;
    //maleCurrentSum = p.maleCurrentSum;
    //femaleCurrentSum = p.femaleCurrentSum;
    CumSumFits = p.CumSumFits;
}

Pop Pop::operator=(const Pop& p)
{
    //std::cout << "Pop copy assignment\n";
    patches = p.patches;
    T2LociToCorrect = p.T2LociToCorrect;
    indexFirstMale = p.indexFirstMale;
    //maleCurrentSum = p.maleCurrentSum;
    //femaleCurrentSum = p.femaleCurrentSum;
    CumSumFits = p.CumSumFits;
    return *this;
}

Pop::Pop(){}*/

Pop::Pop(bool ShouldReadPopFromBinary)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'Pop::Pop'\n";
#endif
    if (ShouldReadPopFromBinary)
    {
        std::cout << "Binary files have not been extensively tested (esp. reusing seeds from biniary file and T4 loci). Please pay attention and ensure your results are meaningful.\n";
        if (GP->nbSpecies > 1)
        {
            std::cout << "Trying to read several species population from a binary file. Sorry, in the current version, binary files only work for a single species.\n";
            abort();
        }

        // What is below sounds misleading but it is used to indicate to SimulationTracker::Initialization that there is (potentially) genetic variance.
        SSP->T1_Initial_AlleleFreqs_AllZeros = false;
        SSP->T1_Initial_AlleleFreqs_AllOnes = false;


        GP->binaryFileToRead.open(SSP->readPopFromBinaryPath);

        // Read each patch (each patch starts with its PatchSize)
        {
            int PN;
            GP->binaryFileToRead.read(PN);
            if (GP->PatchNumber != PN)
            {
                std::cout << "The binary file contained " << PN << " patches while this simulation is set to start with " << GP->PatchNumber << " patches.\n";
                abort();
            }
        }
        
       	patches.reserve(GP->PatchNumber);
        for (int patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
        {
            patches.push_back(Patch(patch_index, ShouldReadPopFromBinary));
        }

        // Read T4tree
        if (SSP->Gmap.T4_nbLoci)
            SSP->T4Tree.readFromBinaryFile(GP->binaryFileToRead);

        // Test everything is read
        GP->binaryFileToRead.testReachedEnd();

        // close
        GP->binaryFileToRead.close();
    } else
    {
    	patches.reserve(GP->PatchNumber);
        for (int patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
        {
            patches.push_back(Patch(patch_index, 'A'));
        }
    }
    
    patches.shrink_to_fit();
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
            std::cout << "Binary files have not been extensively tested (esp. reusing seeds from biniary file and T4 loci). Please pay attention and ensure your results are meaningful.\n";

            ////// Seed
            if (SSP->speciesIndex == (GP->nbSpecies-1))
            {
                file.openForSeed(); // opens a different path only for the seed
                file.writeBinary(GP->rngw.getRNG());
                file.close();
            }

            ////// Population

            // open
            file.open(); // reopen a different path for saving the population

            file.writeBinary(GP->PatchNumber);

            // Write each patch (each patch will start with PatchSize info)
            for (int patch_index = 0 ; patch_index < GP->PatchNumber; patch_index++)
            {
                this->getPatch(patch_index).PrintBinaryFile(file, patch_index);
            }

            // Write T4 tree
            if (SSP->Gmap.T4_nbLoci)
            {
                SSP->T4Tree.PrintBinaryFile(file);
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
void Pop::AddPatch()
{
    Patch p;
    patches.push_back(p);
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


/*void Pop::prepareNextGenerationAndIndexFirstMale(const int patch_index, const std::vector<int>& patchSizeNextGeneration)
{
    // Used in CalculateFitnessForNextGeneration
    maleCurrentSum   = 0.0; 
    femaleCurrentSum = 0.0; 

    CumSumFits.resize(GP->PatchNumber);

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
        indexFirstMale.resize(GP->PatchNumber);
        indexFirstMale[patch_index] = (int) (SSP->sexRatio * (double) patchSizeNextGeneration[patch_index] + 0.5);
        assert(indexFirstMale[patch_index] >= 0 && indexFirstMale[patch_index] <= patchSizeNextGeneration[patch_index]);
    }
}*/




/*double Pop::CalculateFitnessForNextGeneration(Individual& Offspring, int patch_index, int ind_index)
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
}*/


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
        
        if (!SSP->malesAndFemales) 
        {
            assert(CumSumFits[patch_index][0].size() == SSP->patchSize[patch_index]);
            if (SSP->patchSize[patch_index] > 0)
            {
                if ((CumSumFits[patch_index][0].back() / CumSumFits[patch_index][0].size()) < 0.00000000000001)
                {
                    std::cout << "The average fitness of the females (or of the hermaphrodites) in patch " << patch_index << " is " << CumSumFits[patch_index][0].back() / CumSumFits[patch_index][0].size() << ". Selection appears too strong and round off error could become non-negligible! You might want to check your fitness related parameters as well as the mutation rates. As a reminder the fitness effect among loci is multiplicative. If the fecundity was not set to -1.0, this message would not have appeared and the patch size would have simply dropped to zero. One nice way to figure out what is wrong is to print fitness either for each individual ('--fitness_file') or only the mean and variance in fitness ('--fitnessStats_file'). You might also to specify the seed in order to know exactly for what generation(s) you need fitness information." << std::endl;
                    if (GP->CurrentGeneration == 1)
                    {
                        std::cout << "As the error message poped up at the first generation, you should probably check the initial conditions. Maybe you start with 'AllOnes' with selection against the the '1' variant (selection is against the '1' variant when making the 'Multiplicity' assumption). It is also possible that you start the simulation with a fixed lethal variant." << std::endl;
                    }
                    abort();
                }
            }
        }
    } else // SSP->fecundityForFitnessOfOne != -1.0
    {
        if (!SSP->fecundityDependentOfFitness)
        {
            assert(CumSumFits[patch_index][0].size() == SSP->patchSize[patch_index]);
            if (SSP->patchSize[patch_index] > 0)
            {
                if ((CumSumFits[patch_index][0].back() / CumSumFits[patch_index][0].size()) < 0.00000000000001)
                {
                    std::cout << "The average fitness of the females (or of the hermaphrodites) in patch " << patch_index << " is " << CumSumFits[patch_index][0].back() / CumSumFits[patch_index][0].size() << ". Selection appears too strong and round off error could become non-negligible! You might want to check your fitness related parameters as well as the mutation rates." << std::endl;
                    if (GP->CurrentGeneration == 1)
                    {
                        std::cout << "As the error message poped up at the first generation, you should probably check the initial conditions. Maybe you start with 'AllOnes' with selection against the the '1' variant (selection is against the '1' variant when making the 'Multiplicity' assumption). It is also possible that you start the simulation with a fixed lethal variant." << std::endl;
                    }
                    abort();
                }
            }
        }
    }



    if (CumSumFits[patch_index][0].size() && SSP->isAnySelection)
    {
        if (CumSumFits[patch_index][0].back() < 0.0)
        {
            std::cout << "The average fitness of hermaphrodites in patch " << patch_index << " is negative (it is " << CumSumFits[patch_index][0].back() / CumSumFits[patch_index][0].size() << "). This makes no sense. It must be caused by an internal bug!" << std::endl;
            abort();
        }

        if (CumSumFits[patch_index][0].back() > 1e300)
        {
            std::cout << "The sum of fitnesses of hermaphrodites in patch " << patch_index << " is gigantic (it is " << CumSumFits[patch_index][0].back() / CumSumFits[patch_index][0].size() << "). This is likely to cause round-off errors! Please ensure to keep fitness values within reasonable values." << std::endl;
            abort();
        }
    }
}  


void Pop::CalculateFitnesses()
{
#ifdef CALLENTRANCEFUNCTIONS
    std::cout << "Enters in Pop::CalculateFitnesses\n";
#endif 

    /*
    std::cout << "SSP->malesAndFemales = " << SSP->malesAndFemales << "\n";
    std::cout << " GP->CurrentGeneration = " <<  GP->CurrentGeneration << "\n";
    std::cout << " GP->startAtGeneration = " <<  GP->startAtGeneration << "\n";
    std::cout << " GP->__GenerationChange.size() = " <<  GP->__GenerationChange.size() << "\n";
    bool bFind = std::find(GP->__GenerationChange.begin(), GP->__GenerationChange.end(), GP->CurrentGeneration) != GP->__GenerationChange.end();
    std::cout << " std::find(GP->__GenerationChange.begin(), GP->__GenerationChange.end(), GP->CurrentGeneration) != GP->__GenerationChange.end() = " <<  bFind << "\n";
    std::cout << " hasCrazyResettingHappened = " <<  hasCrazyResettingHappened << "\n";
    */


    // resize index first male
    if (SSP->malesAndFemales && indexFirstMale.size() != GP->PatchNumber ) indexFirstMale.resize(GP->PatchNumber);
    


    //////////////////////////////////////////////////////////
    ///// Compute CumSumFits for the current generation //////
    //////////////////////////////////////////////////////////
    if (CumSumFits.size() != GP->PatchNumber) CumSumFits.resize(GP->PatchNumber);
    if (SSP->individualSampling_withWalker && walkers.size() != GP->PatchNumber) walkers.resize(GP->PatchNumber);
    
    //std::cout << "CumSumFits has been resized to " << GP->PatchNumber << "\n";
    
    for ( int patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index )
    {
        // Resizing (necessary for the first generation at least)        
        if (SSP->malesAndFemales)
        {
            CumSumFits[patch_index].resize(2);
            CumSumFits[patch_index][0].resize(0);
            CumSumFits[patch_index][1].resize(0);
            if (SSP->individualSampling_withWalker) walkers[patch_index].resize(2);
        } else
        {
            CumSumFits[patch_index].resize(1);
            CumSumFits[patch_index][0].resize(0);
            CumSumFits[patch_index][0].reserve(SSP->patchSize[patch_index]);
            if (SSP->individualSampling_withWalker) walkers[patch_index].resize(1);
        }


        // index first male
        if (SSP->malesAndFemales) indexFirstMale[patch_index] = (int) (SSP->sexRatio * (double) SSP->patchSize[patch_index] + 0.5);


        // Compute fitnesses
        double femaleCurrentSum = 0.0; // hermaphrodite are also called females here
        double maleCurrentSum = 0.0;
        Patch& patch = this->getPatch(patch_index);
        for ( int ind_index=0 ; ind_index < SSP->patchSize[patch_index] ; ++ind_index )
        {
            double w = patch.getInd(ind_index).CalculateFitness(patch_index);
            if (!SSP->isAnySelection) assert(w == 1.0);
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

        // Security to avoid round off error in extreme parameter sets
        if (SSP->isAnySelection) {checkIfCumSumFitsIsNotTooSmall(patch_index);}

        assert(CumSumFits.size() == GP->PatchNumber);

        // Walker
        if (SSP->individualSampling_withWalker)
        {
            if (SSP->malesAndFemales)
            {
                //walkers[patch_index] = std::pair<Walker, Walker>(Walker(CumSumFits[patch_index][0]), Walker(CumSumFits[patch_index][1]));
                walkers[patch_index][0] = Walker(CumSumFits[patch_index][0]);
                walkers[patch_index][1] = Walker(CumSumFits[patch_index][1]);
            } else
            {
                //walkers[patch_index] = std::pair<Walker, Walker>(Walker(CumSumFits[patch_index][0]), {});
                walkers[patch_index][0] = Walker(CumSumFits[patch_index][0]);
            }
        }
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
    
    assert(SSP->patchSize[patch_from]);
    int parent_index;
    if (SSP->isAnySelection && SSP->selectionOn != 1) // selection on fertility (maybe on both)
    {

        //std::cout << "this->CumSumFits["<<patch_from<<"]["<<sex<<"].size() = " << this->CumSumFits[patch_from][sex].size() << "\n";
        //assert(this->CumSumFits[patch_from][sex].back() > 0.0);
        /*
        std::uniform_real_distribution<double> runiform_double_0andSumOfFit(0.0, this->CumSumFits[patch_from][sex].back()); // CumSumFits must be of patchSize length, not of patchCapacity length
        double rnd = runiform_double_0andSumOfFit(GP->rngw.getRNG());
        */


        //sample
        int index; // might be sex specific or not for the moment
        if (SSP->individualSampling_withWalker)
        {
            index = walkers[patch_from][sex](GP->rngw);
        } else
        {
            assert(this->CumSumFits[patch_from].size() > sex);
            double rnd = GP->rngw.uniform_real_distribution(this->CumSumFits[patch_from][sex].back());
            
            //std::cout << "this->CumSumFits["<<patch_from<<"]["<<sex<<"].back() = " << this->CumSumFits[patch_from][sex].back() << "\n";
            //std::cout << "rnd = " << rnd << "\n";
            std::vector<double>::iterator high = std::upper_bound(CumSumFits[patch_from][sex].begin(), CumSumFits[patch_from][sex].end(), rnd);
            index = distance(CumSumFits[patch_from][sex].begin(), high);
        }

        
        
        // Get the right index depending on the sex
        if (SSP->malesAndFemales)
        {
            if (sex == 0)
            {
                parent_index = index;
            } else
            {
                parent_index = this->indexFirstMale[patch_from] + index;
            }
        } else
        {
            parent_index = index;
        }
    } else
    {
        // no selection on fertility
        
        if (SSP->malesAndFemales)
        {
            if (sex == 0)
            {
                parent_index = GP->rngw.uniform_int_distribution(indexFirstMale[patch_from]);
            } else
            {
                parent_index = GP->rngw.uniform_int_distribution(indexFirstMale[patch_from],SSP->patchSize[patch_from]);
            }
        } else
        {
            parent_index = GP->rngw.uniform_int_distribution(SSP->patchSize[patch_from]);
        }
    }

    /*
    if (!(parent_index >= 0 && parent_index < SSP->patchSize[patch_from]))
    {
        std::cout << "selected parent " << parent_index << " from patch " << patch_from <<  "\n";
        std::cout << "SSP->patchSize["<<patch_from<<"] = " << SSP->patchSize[patch_from] << "\n";
        std::cout << "this->CumSumFits["<<patch_from<<"]["<<sex<<"].back() = " << this->CumSumFits[patch_from][sex].back() << "\n";
        std::cout << "rndSave = " << rndSave << "\n";
    }
    */

    assert(parent_index >= 0 && parent_index < SSP->patchSize[patch_from]);
#ifdef CALLENTRANCEFUNCTIONS
    std::cout << "Exits in Pop::SelectionParent\n";
#endif       
    return parent_index;
}

int Pop::SelectionOriginPatch(uint32_t patch_to, double rndIfNotStochastic)
{
#ifdef CALLENTRANCEFUNCTIONS
    std::cout << "Enters in Pop::SelectionOriginPatch\n";
#endif    
    if (GP->PatchNumber == 1)
    {
        //std::cout << "from patch " << 0 << " to patch " << patch_to << " because only one patch\n";
        return 0;
    }

    if (SSP->dispersalData.BackwardMigration[patch_to][0] == 1.0) // if deterministic
    {
        //std::cout << "from patch " << SSP->dispersalData.BackwardMigrationIndex[patch_to][0] << " to patch " << patch_to << " (because deterministic)\n";
        return SSP->dispersalData.BackwardMigrationIndex[patch_to][0];
    }


    int patch_from = -1;

    double rnd;
    if (SSP->isStochasticMigration)
    {
        rnd = GP->rngw.uniform_real_distribution(1.0); // random between 0 and 1
    } else
    {
        rnd = rndIfNotStochastic;

        if (rnd == 1.0)
        {
            return SSP->dispersalData.BackwardMigrationIndex[patch_to].back();
        } else if (rnd == 0.0)
        {
            return SSP->dispersalData.BackwardMigrationIndex[patch_to].front();
        }
        //std::cout << "bef rnd = " << rnd << "\n";
        assert(rnd >= 0.0 && rnd <= 1.0);
    }
    
        
    /*std::cout << "SSP->dispersalData.BackwardMigration[patch_to].size() = " << SSP->dispersalData.BackwardMigration[patch_to].size() << "\n";
    for (int fake_patch_from = 0 ; fake_patch_from < SSP->dispersalData.BackwardMigration[patch_to].size() ; ++fake_patch_from)
     std::cout << SSP->dispersalData.BackwardMigration[patch_to][fake_patch_from] << " ";
     std::cout << "\n";   
    */
     
    //std::cout << "-----\n";
    for (int fake_patch_from = 0 ; fake_patch_from < SSP->dispersalData.BackwardMigration[patch_to].size() ; ++fake_patch_from)
    {
        double probability = SSP->dispersalData.BackwardMigration[patch_to][fake_patch_from];
        
        if (rnd <= probability)
        {
            patch_from = SSP->dispersalData.BackwardMigrationIndex[patch_to][fake_patch_from];
            break;
        }
        rnd -= probability;
    }
 
    
    //if (!(patch_from >= 0 && patch_from < GP->PatchNumber))
    //{
    //std::cout << "from patch " << patch_from << " to patch " << patch_to << "\n";
    //std::cout << "rnd = " << rnd << "\n";
    //}
    
   
    assert(patch_from >= 0 && patch_from < GP->PatchNumber);
    assert(SSP->patchSize[patch_from] > 0);
        
    
#ifdef CALLENTRANCEFUNCTIONS
    std::cout << "Exits in Pop::SelectionOriginPatch\n";
#endif  
    //std::cout << "----\n";
    return patch_from;
}


void Pop::updatePops(Pop& pop1, Pop& pop2, int speciesIndex, int oldNbPatches, std::vector<int>& previousPatchSizes)
{
    //pop1.hasCrazyResettingHappened = true;
    //pop2.hasCrazyResettingHappened = true;

    
    //////////////////////////////
    // Change number of patches //
    //////////////////////////////

    assert(GP->PatchNumber > 0);

    if (pop1.getNbPatches() > GP->PatchNumber) // If it needs to remove patches
    {
        assert(pop1.getNbPatches() == oldNbPatches);
        int NbPatchesToRemove = pop1.getNbPatches() - GP->PatchNumber;
        for (int i = 0 ; i < NbPatchesToRemove ; ++i)
        {
            pop1.RemoveLastPatch();    // Will call destructor of Patch
            pop2.RemoveLastPatch();    // Will call destructor of Patch
        }
    } else if (pop1.getNbPatches() < GP->PatchNumber) // if it needs to add patches
    {
        assert(pop1.getNbPatches() == oldNbPatches);
        assert(GP->PatchNumber - pop1.getNbPatches() > 0);

        // Loop through new patches
        
        for ( auto patch_index = pop1.getNbPatches(); patch_index < GP->PatchNumber ; ++patch_index)
        {
            pop1.AddPatch();
            pop2.AddPatch();
        }
    }
    assert(pop1.getNbPatches() == GP->PatchNumber);
    assert(pop2.getNbPatches() == GP->PatchNumber);

    ////////////////////////////////////
    // Change patch carrying capacity //
    ////////////////////////////////////
    
    assert(SSP->patchCapacity.size() == GP->PatchNumber);
    assert(SSP->patchSize.size() == GP->PatchNumber);
    assert(previousPatchSizes.size() == GP->PatchNumber); // Yes it is the size of the current number of patches
    

    for (int patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
    {
        assert(pop1.getPatch(patch_index).getpatchCapacity() == pop2.getPatch(patch_index).getpatchCapacity());
            
        if (pop1.getPatch(patch_index).getpatchCapacity() < SSP->patchCapacity[patch_index]) // If the carrying capacity increased.
        {   
       
            // Figure which patch to copy individuals from (if needed to copy)
            int patchToSample = -1;
            if (SSP->fecundityForFitnessOfOne == -1)
            {
                if (pop1.getPatch(patch_index).getpatchCapacity())
                    patchToSample = patch_index;
                else
                {
                
                    double probMigr_patchToSample = -1.0;
                    assert(SSP->dispersalData.forwardMigration.size() == GP->PatchNumber);
                    assert(SSP->dispersalData.forwardMigrationIndex.size() == GP->PatchNumber);
                
                    for (size_t patch_from = 0 ; patch_from < GP->PatchNumber ; ++patch_from )
                    {
                        assert(SSP->dispersalData.forwardMigrationIndex[patch_from].size() == SSP->dispersalData.forwardMigration[patch_from].size());
                        for (size_t fake_patch_to = 0 ; fake_patch_to < SSP->dispersalData.forwardMigrationIndex[patch_from].size() ; ++fake_patch_to)
                        {
                            if (patch_index == SSP->dispersalData.forwardMigrationIndex[patch_from][fake_patch_to])
                            {
                                if (SSP->dispersalData.forwardMigration[patch_from][fake_patch_to] > probMigr_patchToSample && previousPatchSizes[patch_from])
                                {
                                
                                    probMigr_patchToSample = SSP->dispersalData.forwardMigration[patch_from][fake_patch_to];
                                    patchToSample = patch_from;
                                }
                            }
                        }
                    }
                
                    if (!(patchToSample >= 0 && patchToSample < oldNbPatches))
                    {
                        std::cout << "At generation " << GP->CurrentGeneration << " SimBit tried to add individual in the patch index " << patch_index << ". When the fecundityForFitnessOfOne is set to -1 (as it is the case; this is default setting), SimBit will simply copy individuals already present in the patch to reach carrying capacity. If the patch was empty, then it looks at the patch that has the highest forward migration rate toward the focal patch. If no patch has any migration rate to this focal patch, then SimBit does not know how to add individuals into this new patch. If you wanted the patch to remain empty despite the carrying capacity being high, then just set the fecundity to something else than -1 (which is the default). Hopefully, one day I will make this security check before the beginning of the simulation!\n";
                        abort();
                    }
                    assert(probMigr_patchToSample > 0.0 && probMigr_patchToSample <= 1.0);
                }
            }

            // add individuals
            int NbIndsToAdd = SSP->patchCapacity[patch_index] - pop1.getPatch(patch_index).getpatchCapacity();
            for (int ind_index = 0 ; ind_index < NbIndsToAdd ; ++ind_index)
            {
                Individual newInd;
                if (SSP->fecundityForFitnessOfOne == -1)
                {
                    assert(patchToSample != -1);
                    assert(pop1.getPatch(patchToSample).getpatchCapacity());
                    assert(previousPatchSizes[patchToSample]);
                    int indToSample = ind_index % previousPatchSizes[patchToSample];
                
                    newInd = pop1.getPatch(patchToSample).getInd(indToSample);
                }
                
                pop1.getPatch(patch_index).AddIndividual(newInd);
                pop2.getPatch(patch_index).AddIndividual(newInd);
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
        
    /////////////////////////////
    // Set all fitnesses to -1 //
    /////////////////////////////
    // Some other alternative would have better performance but I am assuming that a user won't ask for many parameters changes through time
    for (int patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
    {
        for (int ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ++ind_index)
        {
            for (int haplo_index = 0 ; haplo_index < SSP->ploidy ; ++haplo_index)
            {
                if (SSP->Gmap.T1_nbLoci)
                {
                    if (SSP->T1_isMultiplicitySelection)
                    {
                        pop1.getPatch(patch_index).getInd(ind_index, patch_index).getHaplo(haplo_index).setAllW_T1(-1);
                        pop2.getPatch(patch_index).getInd(ind_index, patch_index).getHaplo(haplo_index).setAllW_T1(-1);
                    } 
                }
                    
                if (SSP->Gmap.T2_nbLoci)
                {
                    if (SSP->T2_isSelection)
                    {
                        pop1.getPatch(patch_index).getInd(ind_index, patch_index).getHaplo(haplo_index).setAllW_T2(-1);
                        pop2.getPatch(patch_index).getInd(ind_index, patch_index).getHaplo(haplo_index).setAllW_T2(-1);
                    }
                }

                if (SSP->Gmap.T56sel_nbLoci)
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
    assert(T2Locus >= 0 && T2Locus < SSP->Gmap.T2_nbLoci);
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
                        if (T2LociToCorrect[i] < 0 || T2LociToCorrect[i] >= SSP->Gmap.T2_nbLoci)
                        {
                            printf("Internal error in Pop::correctT2Loci. Supposed to look for T2 block index %d but there are only %u T2 blocks in the whole genome.\n", T2LociToCorrect[i], SSP->Gmap.T2_nbLoci);
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

std::vector<int> Pop::findWhatMustBeToggledAndUpdateFlipped(const std::vector<unsigned>& popFreqs, std::vector<uint32_t>& flipped, const int& nbLoci,  const double freqThreshold)
{
    std::vector<int> lociToToggle;

    // look for what is to toggle
    assert(flipped.size() <= nbLoci);
    //std::cout << "toggle: ";

    std::vector<bool> isFlipped = inverseWhich(flipped, nbLoci);
    assert(isFlipped.size() == nbLoci);

    for (uint32_t locus = 0 ; locus < nbLoci ; ++locus)
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
    if (SSP->Gmap.T5ntrl_nbLoci)
    {
        std::vector<unsigned> popFreqs = this->computeT56ntrlFrequencies();
        T56ntrlLociToToggle = findWhatMustBeToggledAndUpdateFlipped(popFreqs, SSP->T5ntrl_flipped, SSP->Gmap.T5ntrl_nbLoci, SSP->T56ntrl_frequencyThresholdForFlippingMeaning);
    } else if (SSP->Gmap.T6ntrl_nbLoci)
    {
        std::vector<unsigned> popFreqs = this->computeT56ntrlFrequencies();
        auto flipped = SSP->T6ntrl_flipped.toVector();
        T56ntrlLociToToggle = findWhatMustBeToggledAndUpdateFlipped(popFreqs, flipped, SSP->Gmap.T6ntrl_nbLoci, SSP->T56ntrl_frequencyThresholdForFlippingMeaning);
        SSP->T6ntrl_flipped = CompressedSortedDeque(flipped, SSP->Gmap.T6ntrl_nbLoci);
    }


    //////////////
    /// T5sel ////
    //////////////
    /*std::vector<int> T56selLociToToggle;

    if (SSP->Gmap.T5sel_nbLoci)
    {
        assert(!SSP->Gmap.isT56selCompress);
        std::vector<unsigned> popFreqs = this->computeT56selFrequencies();
        T56selLociToToggle = findWhatMustBeToggledAndUpdateFlipped(popFreqs, SSP->T5sel_flipped, SSP->Gmap.T5sel_nbLoci, SSP->T56sel_frequencyThresholdForFlippingMeaning);
    } else if (SSP->Gmap.T6sel_nbLoci)
    {
        assert(SSP->Gmap.isT56selCompress);
        std::vector<unsigned> popFreqs = this->computeT56selFrequencies();
        auto flipped = SSP->T6sel_flipped.toVector();
        T56selLociToToggle = findWhatMustBeToggledAndUpdateFlipped(popFreqs, flipped, SSP->Gmap.T6sel_nbLoci, SSP->T56sel_frequencyThresholdForFlippingMeaning);
        SSP->T6sel_flipped = CompressedSortedDeque(flipped, SSP->Gmap.T6sel_nbLoci);
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
    for (uint32_t patchIndex = 0 ; patchIndex < patches.size() ; ++patchIndex)
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

    if (SSP->Gmap.T56ntrl_nbLoci)
    {
        // count
        for (uint32_t patch_index = 0 ; patch_index < GP->PatchNumber ; patch_index++)
        {        
            obsFreqs[patch_index].resize(SSP->Gmap.T56ntrl_nbLoci, 0.0);
            for (uint32_t ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ind_index++)
            {
                for (uint32_t haplo_index = 0 ; haplo_index < 2 ; haplo_index++)
                {
                    Haplotype& haplo = this->getPatch(patch_index).getInd(ind_index).getHaplo(haplo_index);

                    if (SSP->Gmap.isT56ntrlCompress)
                    {
                        auto it = haplo.T6ntrl_AllelesBegin();
                        auto itEnd = haplo.T6ntrl_AllelesEnd();
                        for (; it != itEnd ; ++it)
                        {
                            auto value = *it;
                            assert(value >= 0 && value < SSP->Gmap.T6ntrl_nbLoci);
                            obsFreqs[patch_index][value]++;
                        }
                    } else
                    {
                        auto it = haplo.T5ntrl_AllelesBegin();
                        auto itEnd = haplo.T5ntrl_AllelesEnd();
                        for (; it != itEnd ; ++it)
                        {
                            auto value = *it;
                            assert(value >= 0 && value < SSP->Gmap.T5ntrl_nbLoci);
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

    if (SSP->Gmap.T56sel_nbLoci)
    {
        for (uint32_t patch_index = 0 ; patch_index < GP->PatchNumber ; patch_index++)
        {
            // count
            obsFreqs[patch_index].resize(SSP->Gmap.T56sel_nbLoci, 0.0);
            for (uint32_t ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ind_index++)
            {
                for (uint32_t haplo_index = 0 ; haplo_index < 2 ; haplo_index++)
                {
                    Haplotype& haplo = this->getPatch(patch_index).getInd(ind_index).getHaplo(haplo_index);


                    if (SSP->Gmap.isT56selCompress)
                    {
                        auto it = haplo.T6sel_AllelesBegin();
                        auto itEnd = haplo.T6sel_AllelesEnd();
                        for (; it != itEnd ; ++it)
                        {
                            auto value = *it;
                            assert(value >= 0 && value < SSP->Gmap.T6sel_nbLoci);
                            obsFreqs[patch_index][value]++;
                        }
                    } else
                    {
                        auto it = haplo.T5sel_AllelesBegin();
                        auto itEnd = haplo.T5sel_AllelesEnd();
                        for (; it != itEnd ; ++it)
                        {
                            auto value = *it;
                            assert(value >= 0 && value < SSP->Gmap.T5sel_nbLoci);
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

    assert(SSP->Gmap.T56_nbLoci);

    std::vector<std::vector<unsigned>> ntrlFreqs  = this->computePatchSpecificT56ntrlFrequencies();
    std::vector<std::vector<unsigned>> selFreqs = this->computePatchSpecificT56selFrequencies();
    

    std::vector<std::vector<unsigned>> r(GP->PatchNumber);
    for (uint32_t patch_index = 0 ; patch_index < GP->PatchNumber ; patch_index++)
    {
        uint32_t selFreqsIndex = 0;
        uint32_t ntrlFreqsIndex = 0;
        r[patch_index].resize(SSP->Gmap.T56_nbLoci);
        for (uint32_t T56locus = 0 ; T56locus < SSP->Gmap.T56_nbLoci ; T56locus++)
        {
            //std::cout << "SSP->FromT56LocusToT56genderLocus["<<T5locus<<"].first = " << SSP->FromT56LocusToT56genderLocus[T5locus].first << "\n";
            if (SSP->Gmap.FromT56LocusToT56genderLocus(T56locus).isNtrl)
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
        
        assert(selFreqsIndex == SSP->Gmap.T56sel_nbLoci);
        assert(ntrlFreqsIndex == SSP->Gmap.T56ntrl_nbLoci);
    }
    assert(ntrlFreqs.size() == 0);
    assert(selFreqs.size() == 0);
    assert(r.size() == GP->PatchNumber);
    assert(r[0].size() == SSP->Gmap.T56_nbLoci);

    //std::cout << "About to exit Pop::computePatchSpecificT56Frequencies\n";
    return r;

}




std::vector<unsigned> Pop::computeT56ntrlFrequencies()
{
    assert(SSP->Gmap.T56ntrl_nbLoci == SSP->Gmap.T5ntrl_nbLoci + SSP->Gmap.T6ntrl_nbLoci);
    std::vector<unsigned> obsFreqs(SSP->Gmap.T56ntrl_nbLoci,0); // obsFreqs[patch_index][locus]

    if (SSP->Gmap.T56ntrl_nbLoci)
    {
        // count
        for (uint32_t patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
        {        
            for (uint32_t ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ++ind_index)
            {
                for (uint32_t haplo_index = 0 ; haplo_index < 2 ; ++haplo_index)
                {
                    Haplotype& haplo = this->getPatch(patch_index).getInd(ind_index).getHaplo(haplo_index);

                    if (SSP->Gmap.isT56ntrlCompress)
                    {
                        auto it    = haplo.T6ntrl_AllelesBegin();
                        auto itEnd = haplo.T6ntrl_AllelesEnd();
                        for (; it != itEnd ; ++it)
                        {
                            auto value = *it;
                            assert(value >= 0 && value < SSP->Gmap.T6ntrl_nbLoci);
                            ++obsFreqs[value];
                        }
                    } else
                    {
                        auto it = haplo.T5ntrl_AllelesBegin();
                        auto itEnd = haplo.T5ntrl_AllelesEnd();
                        for (; it != itEnd ; ++it)
                        {
                            auto value = *it;
                            assert(value >= 0 && value < SSP->Gmap.T5ntrl_nbLoci);
                            ++obsFreqs[value];
                        }
                    }
                }
            }
        }
    }

    /*
    assert(obsFreqs.size()== SSP->Gmap.T56ntrl_nbLoci);
    std::cout << "what is polymorphic: ";
    for (uint32_t locus = 0 ; locus < obsFreqs.size() ; ++locus)
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
    assert(SSP->Gmap.T56sel_nbLoci == SSP->Gmap.T5sel_nbLoci + SSP->Gmap.T6sel_nbLoci);
    std::vector<unsigned> obsFreqs(SSP->Gmap.T56sel_nbLoci,0); // obsFreqs[patch_index][locus]
    if (SSP->Gmap.T56sel_nbLoci)
    {
        // count
        uint32_t TotalpatchSize = 0;
        for (uint32_t patch_index = 0 ; patch_index < GP->PatchNumber ; patch_index++)
        {
            TotalpatchSize += SSP->patchSize[patch_index];
        
            for (uint32_t ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ind_index++)
            {
                for (uint32_t haplo_index = 0 ; haplo_index < 2 ; haplo_index++)
                {
                    Haplotype& haplo = this->getPatch(patch_index).getInd(ind_index).getHaplo(haplo_index);

                    if (SSP->Gmap.isT56selCompress)
                    {
                        auto it = haplo.T6sel_AllelesBegin();
                        auto itEnd = haplo.T6sel_AllelesEnd();
                        for (; it != itEnd ; ++it)
                        {
                            auto value = *it;
                            assert(value >= 0 && value < SSP->Gmap.T6sel_nbLoci);
                            ++obsFreqs[value];
                        }
                    } else
                    {
                        auto it = haplo.T5sel_AllelesBegin();
                        auto itEnd = haplo.T5sel_AllelesEnd();
                        for (; it != itEnd ; ++it)
                        {
                            auto value = *it;
                            assert(value >= 0 && value < SSP->Gmap.T5sel_nbLoci);
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

    std::vector<unsigned> r(SSP->Gmap.T56_nbLoci);
    
    uint32_t selFreqsIndex = 0;
    uint32_t ntrlFreqsIndex = 0;

    for (uint32_t locus = 0 ; locus < SSP->Gmap.T56_nbLoci ; locus++)
    {
        if (SSP->Gmap.FromT56LocusToT56genderLocus(locus).isNtrl)
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

    assert(selFreqsIndex == SSP->Gmap.T56sel_nbLoci);
    assert(ntrlFreqsIndex == SSP->Gmap.T56ntrl_nbLoci);


    assert(r.size() == SSP->Gmap.T56_nbLoci);

    return r;
}

std::vector<double> Pop::computeT56RelativeFrequencies()
{
    auto freqs = computeT56Frequencies();
    std::vector<double> relFreqs(freqs.size());
    for (uint32_t locus = 0 ; locus < SSP->Gmap.T56_nbLoci ; locus++)
    {
        relFreqs[locus] = (double)freqs[locus] / (2.0*SSP->TotalpatchSize);
        assert(relFreqs[locus]>=0.0 && relFreqs[locus] <= 1.0);
    }
    return relFreqs;
}

std::vector<std::vector<double>> Pop::computeT1RelativeFrequenciesPerPatch()
{
    std::vector<std::vector<double>> r(GP->PatchNumber);
    for (uint32_t patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
    {
        r[patch_index] = this->getPatch(patch_index).computeT1RelativeFrequencies(patch_index);
    }
    return r;
}




std::vector<T1_locusDescription> Pop::listT1PolymorphicLoci()
{
    // Gather the genetics of one "random" individual
    std::vector<bool> oneAllele(SSP->Gmap.T1_nbLoci);
    {
        assert(SSP->patchSize.size() == GP->PatchNumber);
        auto& oneHaplo = getPatch(GP->PatchNumber-1).getInd(SSP->patchSize[GP->PatchNumber-1]-1).getHaplo(0);


        uint32_t locus = 0;

        // First bytes
        for (uint32_t byte = 0 ; byte < (SSP->Gmap.T1_nbChars-1) ; ++byte )
        {
            if (oneHaplo.getT1_char(byte))
            {
                for (uint32_t bit = 0 ; bit < 8 ; ++bit )
                {
                    oneAllele[locus] = oneHaplo.getT1_Allele(byte, bit);
                    ++locus;
                }
            } else
            {
                locus+=8;
            }
        }
        assert(locus == (SSP->Gmap.T1_nbChars-1) * 8);

        // last byte
        auto lastByte = SSP->Gmap.T1_nbChars-1;
        if (oneHaplo.getT1_char(lastByte))
        {
            for (uint32_t bit = 0 ; bit < SSP->Gmap.T1_nbLociLastByte ; ++bit )
            {
                oneAllele[locus] = oneHaplo.getT1_Allele(lastByte, bit);
                ++locus;
            }
        }
    }
        



    // Figure out what is polymorphic
    std::vector<bool> isPoly(SSP->Gmap.T1_nbLoci, false);

    for ( int patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index )
    {
        for (int ind_index=0;ind_index<SSP->patchSize[patch_index];++ind_index)
        {
            assert(SSP->ploidy == 2);
            for ( int haplo_index=0 ; haplo_index < SSP->ploidy ; haplo_index++ )
            { 
                auto& haplo = getPatch(patch_index).getInd(ind_index).getHaplo(haplo_index);

                unsigned locus = 0;

                // First bytes
                for (unsigned byte = 0 ; byte < (SSP->Gmap.T1_nbChars-1) ; ++byte )
                {
                    for (unsigned bit = 0 ; bit < 8 ; ++bit )
                    {
                        if (haplo.getT1_Allele(byte, bit) != oneAllele[locus])
                            isPoly[locus] = true;
                        ++locus;
                    }
                }
                assert(locus == (SSP->Gmap.T1_nbChars-1) * 8);

                // last byte
                auto lastByte = SSP->Gmap.T1_nbChars-1;
                
                for (uint32_t bit = 0 ; bit < SSP->Gmap.T1_nbLociLastByte ; ++bit )
                {
                    if (haplo.getT1_Allele(lastByte, bit) != oneAllele[locus])
                            isPoly[locus] = true;
                    ++locus;
                }
            }
        }
    }
    
    return whichT1_locusDescription(isPoly);
}


void Pop::freeMemory()
{
    freeT56Memory(); // should not be needed
    patches.clear();
    patches.shrink_to_fit();

    T2LociToCorrect.clear();
    T2LociToCorrect.shrink_to_fit();

    indexFirstMale.clear();
    indexFirstMale.shrink_to_fit();
    
    walkers.clear();
    walkers.shrink_to_fit();
}


void Pop::freeT56Memory() // reset all T56 to zero
{
    for (auto& patch : patches)
        patch.freeT56Memory();
}

void Pop::shrink_to_fitT56() // reset all T56 to zero
{
    for (auto& patch : patches)
        patch.shrink_to_fitT56();
}


