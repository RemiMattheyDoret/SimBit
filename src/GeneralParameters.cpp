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


void GeneralParameters::saveSSPPatchSize_toGP()
{
    //std::cout << "Enters in GeneralParameters::saveSSPPatchSize_toGP. this->PatchNumber = " << this->PatchNumber << "\n";
    assert(this->allSpeciesPatchSizes.size() == this->PatchNumber);
    for (int patch_index = 0 ; patch_index < this->PatchNumber ; patch_index++)
    {
        assert(this->allSpeciesPatchSizes[patch_index].size() == this->nbSpecies);
        assert(SSP->speciesIndex >= 0 && SSP->speciesIndex < this->nbSpecies);
        assert(SSP->patchSize[patch_index] >= 0);
        assert(SSP->patchSize[patch_index] <= SSP->patchCapacity[patch_index]);

        this->allSpeciesPatchSizes[patch_index][SSP->speciesIndex] = SSP->patchSize[patch_index];
    }
}
/*
void GeneralParameters::saveSSPPatchSize_toGP_lowerSecurity()
{
    //std::cout << "Enters in GeneralParameters::saveSSPPatchSize_toGP. this->PatchNumber = " << this->PatchNumber << "\n";
    assert(this->allSpeciesPatchSizes.size() == this->PatchNumber);
    for (int patch_index = 0 ; patch_index < this->PatchNumber ; patch_index++)
    {
        assert(this->allSpeciesPatchSizes[patch_index].size() == this->nbSpecies);
        assert(SSP->speciesIndex >= 0 && SSP->speciesIndex < this->nbSpecies);
        
        this->allSpeciesPatchSizes[patch_index][SSP->speciesIndex] = SSP->patchSize[patch_index];
    }
}
*/


void GeneralParameters::readPatchNumber(InputReader input)
{
#ifdef DEBUG
    std::cout << "For option '--PN (--PatchNumber)', the std::string that is read is: " << input.print() << std::endl;
#endif
    maxEverPatchNumber = 0;
    for (int generation_index = 0 ; generation_index < this->__GenerationChange.size() ; generation_index++)
    {
        (void) input.GetNextGenerationMarker(generation_index);
        
        int PN = input.GetNextElementInt();
        if (PN < 1)
        {
            std::cout << "patchNumber must be greater than 0 (receveived " << PN << ")\n";
            abort();
        }

        this->__PatchNumber.push_back(PN);

        maxEverPatchNumber = std::max(maxEverPatchNumber, PN);
        
    }
    assert(this->__PatchNumber.size() == this->__GenerationChange.size());
    input.workDone();

    this->PatchNumber = this->__PatchNumber[0];
    assert(maxEverPatchNumber >= this->PatchNumber);
    assert(maxEverPatchNumber > 0);
}

void GeneralParameters::initializeAllSpeciesPatchSizes()
{
    assert(this->PatchNumber == this->__PatchNumber[0]);
    assert(this->PatchNumber > 0);
    GP->allSpeciesPatchSizes.resize(this->PatchNumber);
    assert(this->nbSpecies > 0);
    for (int patch_index = 0 ; patch_index < this->PatchNumber ; patch_index++)
    {
        GP->allSpeciesPatchSizes[patch_index].resize(this->nbSpecies);
    }
    assert(GP->allSpeciesPatchSizes.size() == GP->PatchNumber);
}



void GeneralParameters::readTemporalChanges(InputReader input)
{
#ifdef DEBUG
    std::cout << "For option '--T (--TemporalChanges)', the std::string that is read is: " << input.print() << std::endl;
#endif
    
    //assert(this->__GenerationChange.size() == 1); // as it should already contain 0
    
    this->__GenerationChange.push_back(0); // must always contain at least zero
    while (input.IsThereMoreToRead())
    {
        int generation = input.GetNextElementInt();
        if (generation != 0)
        {
            assert(generation > 0);
            assert(generation <= this->nbGenerations);
            this->__GenerationChange.push_back(generation);
        }
    }
    // Sort and remove duplicates
    std::sort(this->__GenerationChange.begin(), this->__GenerationChange.end());
    this->__GenerationChange.erase(std::unique(this->__GenerationChange.begin(), this->__GenerationChange.end()), this->__GenerationChange.end());
/*
    std::cout << "__GenerationChange:\n";
    for (auto& x : this->__GenerationChange) std::cout << x << " ";
    std::cout << "\n";
*/
    if (this->__GenerationChange[0] != 0)
    {
        this->__GenerationChange.insert(this->__GenerationChange.begin(), 0);
    }

    input.workDone();
}


void GeneralParameters::UpdateParametersPatchNumber(int generation_index)
{
    assert(generation_index < GP->__PatchNumber.size());
    GP->PatchNumber = GP->__PatchNumber[generation_index];

    // Resize allSpeciesPatchSizes
    allSpeciesPatchSizes.resize(this->PatchNumber);
    for (int patch_index = 0 ; patch_index < this->PatchNumber ; patch_index++)
    {
        allSpeciesPatchSizes[patch_index].resize(this->nbSpecies);
    }
}

void GeneralParameters::readSeed(InputReader input)
{
    if (input.PeakNextElementString() == "binfile" || input.PeakNextElementString() == "f")
    {
        input.skipElement();
        std::string seedBinaryFilePath = input.GetNextElementString();
        std::ifstream file;
        file.open(seedBinaryFilePath, std::ios::out | std::ios::binary);
        file >> mt;
    } else 
    {
        if (input.PeakNextElementString() == "default")
        {
            GP->random_seed = time(NULL);
            
        } else
        {
            random_seed = input.GetNextElementInt();
        }
        std::mt19937 tmp(GP->random_seed);
        GP->mt = tmp;
    }

}


/*void GeneralParameters::UpdateParametersallSpeciesPatchSizes()
{
    // This function should only be used when updating parameters when there is a temporal change. allSpeciesPatchSizes should be updated each generation otherwise directly from the 'main' function and using the higher security version of saveSSPPatchSize_toGP_lowerSecurity

    // Set the values
    assert(SSP == nullptr);
    for (int speciesIndex = 0 ; speciesIndex < GP->nbSpecies ; speciesIndex++)
    {
        SSP = allParameters.getSSPaddress(speciesIndex);
        this->saveSSPPatchSize_toGP_lowerSecurity();
    }
    SSP = nullptr;
}

*/



