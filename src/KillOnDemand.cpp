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
    if (functionToCall.size() == 0) // no kill on demand asked by the user
    {
        return false; 
    }


    if (functionToCall == "isT1LocusFixedAfterGeneration")
    {
        return isT1LocusFixedAfterGeneration(pop);
    } else
    {
        std::cout << "Internal error in KillOnDemand::mustKill(Pop& pop). Unknown functionToCall '" << functionToCall << "'.\n";
        abort();
    }
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



void KillOnDemand::readUserInput(InputReader& input)
{
    functionToCall = input.GetNextElementString();
    if (functionToCall == "isT1LocusFixedAfterGeneration")
    {
        isT1LocusFixedAfterGeneration_T1locus = input.GetNextElementInt();
        if (isT1LocusFixedAfterGeneration_T1locus >= SSP->T1_nbLoci)
        {
            std::cout << "In option --killOnDemand, with function " << functionToCall << ", T1 locus received ("<< isT1LocusFixedAfterGeneration_T1locus <<" is larger or equal to the total number of T1 loci ("<<SSP->T1_nbLoci<<")\n";
            abort();       
        }
        isT1LocusFixedAfterGeneration_generation = input.GetNextElementInt();

        if (isT1LocusFixedAfterGeneration_generation >= GP->nbGenerations)
        {
            std::cout << "In option --killOnDemand, with function " << functionToCall << ", generation received ("<< isT1LocusFixedAfterGeneration_generation <<" is larger or equal to the number of generations to simulate ("<<GP->nbGenerations<<")\n";
            abort();       
        }
    } else if (functionToCall != "default")
    {
        std::cout << "In option --killOnDemand, unknown 'functionToCall' received (received " << functionToCall << "\n";
        abort();
    }
}
