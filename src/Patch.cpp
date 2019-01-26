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


void Patch::PrintBinaryFile(OutputFile& file, int patch_index)
{
    assert(SSP->patchSize[patch_index] <= SSP->patchCapacity[patch_index]);
    file.writeBinary(SSP->patchSize[patch_index]);
    for (int ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ind_index++)
    {
        this->getInd(ind_index).PrintBinaryFile(file);
    }
}

int Patch::getpatchCapacity()
{
    return inds.size();
}

void Patch::AddIndividual(Individual& newInd)
{
    inds.push_back(newInd); // Copy individual received as a reference
}

void Patch::removeLastIndividual()
{
    inds.pop_back(); // call destructor of Individual
    inds.shrink_to_fit();
}

Individual& Patch::getInd(const int& ind_index)
{
    assert(ind_index < inds.size());
    return inds[ind_index];
}

Individual& Patch::getInd(const int& ind_index, const int& patch_index)
{
    assert(ind_index < SSP->patchSize[patch_index]);
    return this->getInd(ind_index);
}

Patch::Patch(){} // nothing to do

Patch::Patch(int patch_index, bool ShouldReadPopFromBinary)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'Patch::Patch(int patch_index, bool ShouldReadPopFromBinary)'\n";
#endif
    GP->BinaryFileToRead >> SSP->patchSize[patch_index];
    assert(SSP->patchSize[patch_index] >= 0);
    assert(SSP->patchSize[patch_index] <= SSP->patchCapacity[patch_index]);

    if (SSP->patchSize[patch_index] != SSP->patchCapacity[patch_index])
    {
        if (SSP->fecundityForFitnessOfOne == -1)
        {
            std::cout << "When reading the binary file for species "<< SSP->speciesName <<", it indicated a patch size of " << SSP->patchSize[patch_index] << " for patch index (zero based counting) " << patch_index << ". This value differs from the carrying capacity (" <<  SSP->patchCapacity[patch_index] << "), which should not be possible as you have set fecundityForFitnessOfOne to -1\n";
            abort();
        }
    }

    std::cout << "SSP->patchSize[patch_index] = " << SSP->patchSize[patch_index] << "\n";

    for (int ind_index=0;ind_index<SSP->patchSize[patch_index];++ind_index)
    {
        Individual ind(ShouldReadPopFromBinary);
        inds.push_back(std::move(ind));
    }

    std::cout << "inds.size() = " << inds.size() << "\n";
}

Patch::Patch(const int patch_index, char Abiogenesis)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'Patch::Patch(const int patch_index, char Abiogenesis)'\n";
#endif


    if (
        (SSP->T1_Initial_AlleleFreqs_AllZeros || SSP->T1_Initial_AlleleFreqs_AllOnes)
        &&
        (SSP->T5_Initial_AlleleFreqs_AllZeros || SSP->T5_Initial_AlleleFreqs_AllOnes)
        )
    {   
        Haplotype ConstantHaplotype(patch_index, Abiogenesis, -1.0);
        Individual ConstantInd(ConstantHaplotype, Abiogenesis);

        assert(patch_index < SSP->patchCapacity.size());
        for (int ind_index = 0 ; ind_index < SSP->patchCapacity[patch_index] ; ++ind_index)
        {
            inds.push_back(ConstantInd);
        }
    } else
    {
        for (int ind_index = 0 ; ind_index < SSP->patchCapacity[patch_index] ; ++ind_index)
        {
            Individual ind(patch_index, Abiogenesis, ind_index);
            inds.push_back(std::move(ind));
        }
    }
}

