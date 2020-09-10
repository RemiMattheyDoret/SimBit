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


bool KillOnDemand::mustKill(Pop& pop)
{
    if (functionToCall.size() > 0) // no kill on demand asked by the user
    {
        if (functionToCall == "isT1LocusFixedAfterGeneration")
        {
            justAboutToKill = isT1LocusFixedAfterGeneration(pop);
        } else if (functionToCall == "isT1HaplotypeConfigurationFixed")
        {
            justAboutToKill = isT1HaplotypeConfigurationFixed(pop);
        } else if (functionToCall == "isTotalPatchSizeGreaterThan")
        {
            justAboutToKill = isTotalPatchSizeGreaterThan(pop);
        } else if (functionToCall == "isT4LocusFixed")
        {
            justAboutToKill = isT4LocusFixed(pop);
        } else
        {
            std::cout << "Internal error in KillOnDemand::mustKill(Pop& pop). Unknown functionToCall '" << functionToCall << "'.\n";
            abort();
        }        
    }

    return justAboutToKill;
}

bool KillOnDemand::isT4LocusFixed(Pop& pop)
{
    SSP->T4Tree.setLocusForWhichFixationMustBeComputedAtTheNextSimplify(isT4LocusFixed_T4Locus);
    SSP->T4Tree.simplify(pop);
    return SSP->T4Tree.isLocusForWhichFixationHadToBeComputedFixed();
}

bool KillOnDemand::isT1LocusFixedAfterGeneration(Pop& pop)
{
    if (GP->CurrentGeneration < isT1LocusFixedAfterGeneration_generation)
    {
        return false;
    }

    bool isAnyOne  = false;
    bool isAnyZero = false;
    for (unsigned patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
    {
        Patch& patch = pop.getPatch(patch_index);
        for (unsigned ind_index = 0; ind_index < SSP->patchSize[patch_index] ; ++ind_index)
        {
            Individual& ind = patch.getInd(ind_index);
            //std::cout << ind.getHaplo(0).getT1_Allele(isT1LocusFixedAfterGeneration_T1locus) << " ";
            if (ind.getHaplo(0).getT1_Allele(isT1LocusFixedAfterGeneration_T1locus))
            {
                isAnyOne = true;
            } else
            {
                isAnyZero = true;
            }

            if (isAnyZero && isAnyOne)
            {
                return false;
            }
        }
    }

    return true;
}

bool KillOnDemand::isT1HaplotypeConfigurationFixed(Pop& pop)
{
    // Set a configuration
    std::vector<bool> configuration;
    configuration.reserve(killWhenT1HaplotypeConfigurationIsFixed_T1loci.size());
    {
        for (uint32_t patch_index = 0 ; patch_index < GP->PatchNumber; ++patch_index )
        {
            if (SSP->patchSize[patch_index])
            {
                auto& haplo = pop.getPatch(GP->PatchNumber-1).getInd(SSP->patchSize[patch_index]-1).getHaplo(0);
                for (uint32_t locus_index = 0 ; locus_index < killWhenT1HaplotypeConfigurationIsFixed_T1loci.size() ; ++locus_index)
                {
                    auto& locus = killWhenT1HaplotypeConfigurationIsFixed_T1loci[locus_index];
                    configuration.push_back(haplo.getT1_Allele(locus));
                }
                break;
            }
        }        
    }

    if (configuration.size() == 0)
    {
        return true;
    }


        
    // Test they all have the same configuration
    for (uint32_t patch_index = 0 ; patch_index < GP->PatchNumber; ++patch_index )
    {
        auto& patch = pop.getPatch(patch_index);
        for (uint32_t ind_index = 0 ; ind_index < SSP->patchSize[patch_index]; ++ind_index )
        {
            auto& ind = patch.getInd(ind_index);
            for (uint32_t haplo_index = 0 ; haplo_index < 2; ++haplo_index )
            {
                auto& haplo = ind.getHaplo(haplo_index);
                for (uint32_t locus_index = 0 ; locus_index < killWhenT1HaplotypeConfigurationIsFixed_T1loci.size() ; ++locus_index)
                {
                    auto& locus = killWhenT1HaplotypeConfigurationIsFixed_T1loci[locus_index];
                    auto state = configuration[locus_index];

                    if (haplo.getT1_Allele(locus) != state)
                    {
                        return false;
                    }
                }
            }   
        }
    }
    return true;
}


bool KillOnDemand::isTotalPatchSizeGreaterThan(Pop& pop)
{
    // pop is actually not even used
    int sum = 0; // int because patchSize is an int (which is stupid I think, it should be unsigned)
    for (auto& patch : isTotalPatchSizeGreaterThan_patches)
    {
        if (SSP->patchSize.size() > patch)
        {
            sum += SSP->patchSize[patch];
        }
    }
    return sum > isTotalPatchSizeGreaterThan_N;
}



void KillOnDemand::readUserInput(InputReader& input)
{
    functionToCall = input.GetNextElementString();
    if (functionToCall == "isT1LocusFixedAfterGeneration")
    {
        isT1LocusFixedAfterGeneration_T1locus = input.GetNextElementInt();
        if (isT1LocusFixedAfterGeneration_T1locus >= SSP->Gmap.T1_nbLoci)
        {
            std::cout << "In option --killOnDemand, with function " << functionToCall << ", T1 locus received ("<< isT1LocusFixedAfterGeneration_T1locus <<" is larger or equal to the total number of T1 loci ("<<SSP->Gmap.T1_nbLoci<<")\n";
            abort();       
        }
        isT1LocusFixedAfterGeneration_generation = input.GetNextElementInt();

        if (isT1LocusFixedAfterGeneration_generation >= GP->nbGenerations)
        {
            std::cout << "In option --killOnDemand, with function " << functionToCall << ", generation received ("<< isT1LocusFixedAfterGeneration_generation <<" is larger or equal to the number of generations to simulate ("<<GP->nbGenerations<<")\n";
            abort();       
        }
    } else if (functionToCall == "killWhenT1HaplotypeConfigurationIsFixed")
    {
        while (input.IsThereMoreToRead())
        {
            killWhenT1HaplotypeConfigurationIsFixed_T1loci.push_back(input.GetNextElementInt());
        }
    } else if (functionToCall == "isTotalPatchSizeGreaterThan")
    {
        isTotalPatchSizeGreaterThan_N = input.GetNextElementInt();
        while (input.IsThereMoreToRead())
        {
            auto patch = input.GetNextElementInt();

            if (patch < 0 || patch >= GP->maxEverPatchNumber)
            {
                std::cout << "In --KillOnDemand with functionToCall" << functionToCall << ", received a patch that is either a negative value or a value that is greater or equal to the maximum number of patches ever in the simulation\n";
                abort();
            }
            isTotalPatchSizeGreaterThan_patches.push_back(patch);
        }
        isTotalPatchSizeGreaterThan_patches.shrink_to_fit();

        // assertion condition is possible
        int sum = 0; // int because patchSize is an int (which is stupid I think, it should be unsigned)
        //std::string s;
        for (auto& patch : isTotalPatchSizeGreaterThan_patches)
        {
            assert(SSP->maxEverpatchCapacity.size() > patch);
            sum += SSP->maxEverpatchCapacity[patch];
            //s += std::to_string(SSP->maxEverpatchCapacity[patch]) + "+";
        }
        //s.pop_back();
        //s += " = " + std::to_string(sum) + "\n";
        //std::cout << s;
        //std::cout << "isTotalPatchSizeGreaterThan_N = " << isTotalPatchSizeGreaterThan_N << "\n";
        if (sum <= isTotalPatchSizeGreaterThan_N)
        {
            std::cout << "\tWarning: given the carrying capacities at the patches indicated, the function isTotalPatchSizeGreaterThan will never kill the simulation\n";
            functionToCall = "";
        }
    } else if (functionToCall == "isT4LocusFixed")
    {
        isT4LocusFixed_T4Locus = input.GetNextElementInt();
        if (isT4LocusFixed_T4Locus < 0 || isT4LocusFixed_T4Locus >= SSP->Gmap.T4_nbLoci)
        {
            std::cout << "In --KillOnDemand with functionToCall" << functionToCall << ", received a T4 locus patch that is either a negative value or a value that is greater or equal to the total number of loci (received '"<<isT4LocusFixed_T4Locus<<"')\n";
                abort();
        }
    } else if (functionToCall != "default")
    {
        std::cout << "In option --killOnDemand, unknown 'functionToCall' received (received " << functionToCall << "\n";
        abort();
    }
}
