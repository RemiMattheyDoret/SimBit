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

void Patch::swapIndividuals(size_t i, size_t j)
{
    std::swap(inds[j],inds[i]);
}

int Patch::getpatchCapacity()
{
    return inds.size();
}

void Patch::AddIndividual(Individual& newInd)
{
    inds.push_back(newInd); // Copy individual received as a reference
   /* if (SSP->Gmap.T4_nbLoci > 0)
    {
        std::cout << "newInd.getHaplo(0).T4ID = " << newInd.getHaplo(0).T4ID << "\n";
        std::cout << "newInd.getHaplo(1).T4ID = " << newInd.getHaplo(1).T4ID << "\n";
        assert(newInd.getHaplo(0).T4ID == inds.back().getHaplo(0).T4ID);
        assert(newInd.getHaplo(1).T4ID == inds.back().getHaplo(1).T4ID);
    }*/
}

void Patch::removeLastIndividual()
{
    inds.pop_back(); // call destructor of Individual
}

void Patch::shrink_to_fit_inds()
{
    inds.shrink_to_fit();
}

void Patch::setInd(Individual& ind, const int& ind_index)
{
    assert(inds.size() > ind_index);
    inds[ind_index] = ind;
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

Patch::Patch(){}

/*Patch::Patch(const Patch& p)
{
    inds = p.inds;
}

Patch Patch::operator=(const Patch& p)
{
    inds = p.inds;
    return *this;
}

Patch::Patch(Patch&& p)
{
    inds.swap(p.inds);
}

Patch Patch::operator=(Patch&& p)
{
    inds.swap(p.inds);
    return *this;
}*/

Patch::Patch(int patch_index, bool ShouldReadPopFromBinary)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'Patch::Patch(int patch_index, bool ShouldReadPopFromBinary)'\n";
#endif
    GP->binaryFileToRead.read(SSP->patchSize[patch_index]);

    if (SSP->patchSize[patch_index] < 0)
    {
        std::cout << "For patch index " << patch_index << " (zero-based counting), the binary file indicate a negative patch size. Patch size read is " << SSP->patchSize[patch_index] << "\n";
        abort();   
    }

    if (SSP->patchSize[patch_index] > SSP->patchCapacity[patch_index])
    {
        std::cout << "For patch index " << patch_index << " (zero-based counting), you set a patch capacity of " << SSP->patchCapacity[patch_index] << " but it received a patch size that is greater than that. Patch size from the binary file is " << SSP->patchSize[patch_index] << "\n";
        abort();
    }

    if (SSP->patchSize[patch_index] != SSP->patchCapacity[patch_index])
    {
        if (SSP->fecundityForFitnessOfOne == -1)
        {
            std::cout << "When reading the binary file for species "<< SSP->speciesName <<", it indicated a patch size of " << SSP->patchSize[patch_index] << " for patch index (zero based counting) " << patch_index << ". This value differs from the carrying capacity (" <<  SSP->patchCapacity[patch_index] << "), which should not be possible as you have set fecundityForFitnessOfOne to -1\n";
            abort();
        }
    }

    //std::cout << "SSP->patchSize["<<patch_index<<"] = " << SSP->patchSize[patch_index] << "\n";

    inds.reserve(SSP->patchSize[patch_index]);
    for (int ind_index=0 ; ind_index < SSP->patchSize[patch_index] ; ++ind_index)
    {
        //std::cout << "About to read individual index " << ind_index << "\n";
        inds.push_back(Individual(ShouldReadPopFromBinary));
        //std::cout << "Individual read\n";
    }

    //std::cout << "inds.size() = " << inds.size() << "\n";
}

Patch::Patch(const int patch_index, char Abiogenesis)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'Patch::Patch(const int patch_index, char Abiogenesis)'\n";
#endif

    inds.reserve(SSP->patchCapacity[patch_index]);

    //////////////////////////////////////
    /// Initialize with indiviualTypes ///
    //////////////////////////////////////
    if (SSP->isIndividualInitialization)
    {
        assert(SSP->individualTypesForInitialization[patch_index].size() == SSP->patchSize[patch_index]);
        for (auto& individualTypeName : SSP->individualTypesForInitialization[patch_index])
        {
            inds.push_back(SSP->individualTypes[individualTypeName]); // copy
        }        
    }

    //////////////////////////////////////////
    /// Initialize with default parameters ///
    //////////////////////////////////////////
    else if (
        (SSP->T1_Initial_AlleleFreqs_AllZeros || SSP->T1_Initial_AlleleFreqs_AllOnes)
        &&
        (SSP->T56_Initial_AlleleFreqs_AllZeros || SSP->T56_Initial_AlleleFreqs_AllOnes)
        )
    {   
        Haplotype ConstantHaplotype(patch_index, Abiogenesis, -1.0);
        Individual ConstantInd(ConstantHaplotype, Abiogenesis);

        assert(patch_index < SSP->patchCapacity.size());
        for (int ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ++ind_index)
        {
            inds.push_back(ConstantInd);
        }
    }


    ////////////////////////////////////////////////////////
    /// Initialize with T1_ini and other stuff like that /// // That's an old method. Might need to be deprecated
    ///////////////////////////////////////////////////////
    else
    {
        for (int ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ++ind_index)
        {
            Individual ind(patch_index, Abiogenesis, ind_index);
            inds.push_back(ind);
        }
    }


    //// Ensures that they are at carrying capacity in case the initial patch size differs from the carrying capacity
    if (inds.size() != SSP->patchCapacity[patch_index])
    {
        Haplotype ConstantHaplotype(patch_index, Abiogenesis, -1.0);
        Individual ConstantInd(ConstantHaplotype, Abiogenesis);
        for (int ind_index = inds.size() ; ind_index < SSP->patchCapacity[patch_index] ; ++ind_index)
        {
            inds.push_back(ConstantInd); // Those individuals won't be used be we still need to allocate their memory
        }

        assert(inds.size() == SSP->patchCapacity[patch_index]);
    }

    inds.shrink_to_fit();
}


void Patch::toggleT56LociFromEveryone(std::vector<int>& T5ntrlLociToToggle, std::vector<int>& T5selLociToToggle, int Habitat)
{
    for (auto& ind : inds)
        ind.toggleT56LociFromHaplotypes(T5ntrlLociToToggle, T5selLociToToggle, Habitat);
}


std::vector<double> Patch::computeT1RelativeFrequencies(uint32_t self_patch_index)
{
    std::vector<double> r(SSP->Gmap.T1_nbLoci, 0.0);   
    for (uint32_t ind_index = 0 ; ind_index < SSP->patchSize[self_patch_index] ; ++ind_index)
    {
        for (uint32_t haplo_index = 0 ; haplo_index < 2 ; ++haplo_index)
        {
            for (uint32_t locus = 0 ; locus < SSP->Gmap.T1_nbLoci ; ++locus)
            {
                r[locus] += getInd(ind_index).getHaplo(haplo_index).getT1_Allele(locus);
            }
        }
    }

    for (uint32_t locus = 0 ; locus < SSP->Gmap.T1_nbLoci ; ++locus)
    {
        r[locus] /= 2*SSP->patchSize[self_patch_index];
    }

    return r;
}

void Patch::shrink_to_fitT56()
{
    for (auto& ind : inds)
        ind.shrink_to_fitT56();
}


void Patch::freeT56Memory()
{
    for (auto& ind : inds)
        ind.freeT56Memory();
}
