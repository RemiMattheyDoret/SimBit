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


GeneralParameters::GeneralParameters(){}

/*void GeneralParameters::saveSSPPatchSize_toGP()
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
}*/
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

void GeneralParameters::readT1_FST_info(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option '--T1_FST_info', the std::string that is read is: " << input.print() << std::endl;
#endif
    if (!input.IsThereMoreToRead())
    {
        std::cout << "For option '--T1_FST_info', an empty string has been received!\n";
        abort();
    }
    if (input.PeakNextElementString() == "default")
    {
        input.skipElement();
        if (input.IsThereMoreToRead())
        {
            std::cout << "For option '--T1_FST_info', received 'default' followed by some other stuff. The entire input received for this option is: "<<input.print()<<" \n";
            abort();       
        }
        if (outputWriter.isFile(T1_FST))
        {
            if (GP->maxEverPatchNumber < 2)
            {
                std::cout << "You seem to be asking for FST outputs but you never have more than 1 patch in the population as indicated in --PatchNumber (--PN).\n";
                abort();           
            }
            output_FST_nbPatchesToConsider.push_back( 2 );
        }
    } else
    {
        if (!outputWriter.isFile(T1_FST))
        {
            std::cout << "For option '--T1_FST_info', received an input different than 'default' but you have not indicated a file with the option '--T1_FST_file'.\n";
            abort();
        }
        while (input.IsThereMoreToRead())
        {
            int nbPatchesToConsider = input.GetNextElementInt();
            if (nbPatchesToConsider < 2)
            {
                std::cout << "For option '--T1_FST_info', received the value '" << nbPatchesToConsider << "'. This number is lower than 2 and no FST can hence be computed.\n";
                abort();
            }
            if (nbPatchesToConsider > GP->maxEverPatchNumber)
            {
                std::cout << "For option '--T1_FST_info', received the value '" << nbPatchesToConsider << "'. This number is greater than the 'GP->maxEverPatchNumber' (the greatest number of patches at any point during the simulation as indiciated in '--PatchNumber (--PN)'; GP->maxEverPatchNumber = "<<GP->maxEverPatchNumber<<").\n";
                abort();
            }

            output_FST_nbPatchesToConsider.push_back( nbPatchesToConsider );
        }
        sortAndRemoveDuplicates(output_FST_nbPatchesToConsider);
    }
        

    input.workDone();
}

void GeneralParameters::readSpeciesEcologicalRelationships(InputReader& input)
{
    // Allocate memory first
    assert(GP->nbSpecies > 0);
    speciesInteraction.resize(GP->nbSpecies);
    speciesCompetition.resize(GP->nbSpecies);
    for (int speciesIndex = 0 ; speciesIndex < GP->nbSpecies ; ++speciesIndex)
    {
        speciesInteraction[speciesIndex].resize(GP->nbSpecies);
        speciesCompetition[speciesIndex].resize(GP->nbSpecies);
    }
    __speciesInteraction.reserve(this->__GenerationChange.size());
    __speciesCompetition.reserve(this->__GenerationChange.size());
    //__typeOfSpeciesInteraction.reserve(this->__GenerationChange.size());
        



    // read input
    for (int generation_index = 0 ; generation_index < this->__GenerationChange.size() ; generation_index++)
    {
        (void) input.GetNextGenerationMarker(generation_index);
        


        /////////////////
        // interaction //
        /////////////////
        std::string s_interaction_exp = "predation";
        std::string s_interaction = input.GetNextElementString();
        if (s_interaction != s_interaction_exp)
        {
            std::cout << "In --eco, expected the keyword '" << s_interaction_exp << "' but got '" << s_interaction << "' instead\n";
            abort();
        }


        /*
        Actually I'll only do additive
        std::string s_mulOrAdd = input.GetNextElementString();
        if (s_mulOrAdd == "multiply")
        {
            typeOfSpeciesInteraction = 'M';
        } else if (s_mulOrAdd == "add")
        {
            typeOfSpeciesInteraction = 'A';
        } else if (s_mulOrAdd == "no" || s_mulOrAdd == "default")
        {
            typeOfSpeciesInteraction = '0';
        } else
        {
            std::cout << "In --eco, expected either keyword 'multiply', 'add' or 'no' but instead received " << s_mulOrAdd << "\n";
            abort();
        }*/

        if (input.PeakNextElementString() == "default")
        {
            input.skipElement();
            
            for (int speciesIndexFrom = 0 ; speciesIndexFrom < GP->nbSpecies ; ++speciesIndexFrom)
            {
                for (int speciesIndexTo = 0 ; speciesIndexTo < GP->nbSpecies ; ++speciesIndexTo)
                {
                    speciesInteraction[speciesIndexTo][speciesIndexFrom] = {'A',0.0};
                }
            }
        } else
        {
            for (int speciesIndexFrom = 0 ; speciesIndexFrom < GP->nbSpecies ; ++speciesIndexFrom)
            {
                for (int speciesIndexTo = 0 ; speciesIndexTo < GP->nbSpecies ; ++speciesIndexTo)
                {
                    if (speciesIndexFrom == speciesIndexTo)
                    {
                        if (input.GetNextElementString() != "self")
                        {
                            std::cout << "For option --eco, received something else than 'self' for the predation matrix for the effect of species index " << speciesIndexFrom << " onto itself.\n";
                            abort();
                        }
                        speciesInteraction[speciesIndexTo][speciesIndexFrom] = {'A',0.0}; 
                        
                    } else
                    {
                        std::string stype = input.GetNextElementString(); // I don't direclty cast to avoid problems with string that has more than one character
                        char type;
                        double magnitude = input.GetNextElementDouble();
                        if (stype == "A")
                        {
                            type = 'A';
                        } else if (stype == "B")
                        {
                            type = 'B';
                        } else if (stype == "C")
                        {
                            type = 'C';
                        } else if (stype == "D")
                        {
                            type = 'D';
                        } else if (stype == "0")
                        {
                            type = '0';
                        } else
                        {
                            std::cout << "For option --eco, for the predation matrix, received the type " << stype <<" for the effect of species index "<<speciesIndexFrom << " on species index "<<speciesIndexTo << ". Sorry only types 'A', 'B', 'C' and 'D' are recognized.\n";
                            abort();
                        }
                        speciesInteraction[speciesIndexTo][speciesIndexFrom] = {type,magnitude};
                    }
                }
            }
        }




        /////////////////
        // competition //
        /////////////////
        std::string s_competition_exp = "competition";
        std::string s_competition = input.GetNextElementString();
        if (s_competition != s_competition_exp)
        {
            std::cout << "In --eco, expected the keyword '" << s_competition_exp << "' but got '" << s_competition << "' instead\n";
            abort();
        }


        if (input.PeakNextElementString() == "default")
        {
            input.skipElement();

            for (int speciesIndexFrom = 0 ; speciesIndexFrom < GP->nbSpecies ; ++speciesIndexFrom)
            {
                for (int speciesIndexTo = 0 ; speciesIndexTo < GP->nbSpecies ; ++speciesIndexTo)
                {
                    if (speciesIndexFrom == speciesIndexTo)
                    {
                        speciesCompetition[speciesIndexTo][speciesIndexFrom] = 1.0;
                    } else
                    {
                        speciesCompetition[speciesIndexTo][speciesIndexFrom] = 0.0;
                    }
                }
            }
        } else
        {
            for (int speciesIndexFrom = 0 ; speciesIndexFrom < GP->nbSpecies ; ++speciesIndexFrom)
            {
                for (int speciesIndexTo = 0 ; speciesIndexTo < GP->nbSpecies ; ++speciesIndexTo)
                {
                    if (speciesIndexFrom == speciesIndexTo)
                    {
                        if (input.GetNextElementString() != "self")
                        {
                            std::cout << "For option --eco, received something else than 'self' for the competition matrix for the effect of species index " << speciesIndexFrom << " onto itself.\n";
                            abort();
                        }
                        speciesCompetition[speciesIndexTo][speciesIndexFrom] = 1.0;
                    } else
                    {
                        double magnitude = input.GetNextElementDouble();
                        if (magnitude < 0.0)
                        {
                            std::cout << "For option --eco, received a negative value for the magnitude (received "<<magnitude<< ") for the competition matrix for the effect of species index " << speciesIndexFrom << " onto species index "<<speciesIndexTo << ".\n";
                            abort();
                        }
                        speciesCompetition[speciesIndexTo][speciesIndexFrom] = magnitude;
                    }
                }
            }
        }

        ///////////////////////////
        // Assign to __variables //
        ///////////////////////////
        __speciesInteraction.push_back(speciesInteraction);
        __speciesCompetition.push_back(speciesCompetition);
        //__typeOfSpeciesInteraction.push_back(typeOfSpeciesInteraction);
        

    }
}


void GeneralParameters::readPatchNumber(InputReader& input)
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

void GeneralParameters::setAllPatchSizePreviousGenerationIfNeeded()
{
    if (nbSpecies > 1) // I could add if at least one fec is different from -1 but I am lazy and it does not matter
    {
        assert(allSpeciesPatchSizePreviousGeneration.size() == PatchNumber);
        for (int patch_index = 0 ; patch_index < PatchNumber ; patch_index++)
        {
            assert(allSpeciesPatchSizePreviousGeneration[patch_index].size() == nbSpecies);
            for (int speciesIndex = 0 ; speciesIndex < nbSpecies ; ++speciesIndex)
            {
                allSpeciesPatchSizePreviousGeneration[patch_index][speciesIndex] = allParameters.SSPs[speciesIndex].patchSize[patch_index];
            }
        }
    }
}

/*void GeneralParameters::initializeAllSpeciesPatchSizes()
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
}*/



/*void GeneralParameters::readTemporalChanges(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option '--T (--TemporalChanges)', the std::string that is read is: " << input.print() << std::endl;
#endif

    std::cout << "Internal error. This function should not be used anymore!\n";
    abort();
    
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

    //std::cout << "__GenerationChange:\n";
    //for (auto& x : this->__GenerationChange) std::cout << x << " ";
    //std::cout << "\n";

    if (this->__GenerationChange[0] != 0)
    {
        this->__GenerationChange.insert(this->__GenerationChange.begin(), 0);
    }

    input.workDone();
}
*/



void GeneralParameters::readTemporalChanges(std::vector<int>& T)
{
#ifdef DEBUG
    std::cout << "Entering in 'readTemporalChanges'" << std::endl;
#endif
    
    //assert(this->__GenerationChange.size() == 1); // as it should already contain 0
    
    this->__GenerationChange.push_back(0); // must always contain at least zero
    for (auto& t : T)
    {
        if (t < 0)
        {
            std::cout << "Received a generation specific marker (@G) followed by a negative number (followed by the number "<< t <<"). Generations must be zero or positive integers.\n";
            abort();
        }
        this->__GenerationChange.push_back(t);
    }

    // Sort and remove duplicates
    std::sort(this->__GenerationChange.begin(), this->__GenerationChange.end());
    this->__GenerationChange.erase(std::unique(this->__GenerationChange.begin(), this->__GenerationChange.end()), this->__GenerationChange.end());
/*
    std::cout << "__GenerationChange:\n";
    for (auto& x : this->__GenerationChange) std::cout << x << " ";
    std::cout << "\n";
*/
/*
    if (this->__GenerationChange[0] != 0)
    {
        this->__GenerationChange.insert(this->__GenerationChange.begin(), 0);
    }
*/
}


void GeneralParameters::update(int generation_index)
{
    assert(generation_index < __PatchNumber.size());
    //assert(generation_index < __typeOfSpeciesInteraction.size());
    
    GP->PatchNumber = GP->__PatchNumber[generation_index];
    //typeOfSpeciesInteraction = __typeOfSpeciesInteraction[generation_index];
    speciesInteraction = __speciesInteraction[generation_index];
    speciesCompetition = __speciesCompetition[generation_index];

    // Resize allSpeciesPatchSizePreviousGeneration
    allSpeciesPatchSizePreviousGeneration.resize(this->PatchNumber);
    for (int patch_index = 0 ; patch_index < this->PatchNumber ; patch_index++)
    {
        allSpeciesPatchSizePreviousGeneration[patch_index].resize(this->nbSpecies);
    }
}

void GeneralParameters::readSeed(InputReader& input)
{
    //std::cout << "in readSeed: VIndex = " << input.getVIndex() << "\n";
    //std::cout << "in readSeed: V.size() = " << input.getSizeOfV() << "\n";
    int random_seed;
    if (input.PeakNextElementString() == "binfile" || input.PeakNextElementString() == "f")
    {
        input.skipElement();
        std::string seedBinaryFilePath = input.GetNextElementString();
        std::ifstream file;
        file.open(seedBinaryFilePath, std::ios::out | std::ios::binary);
        file >> rngw.getRNG();
        random_seed = std::numeric_limits<int>::quiet_NaN();
    } else 
    {
        if (input.PeakNextElementString() == "default")
        {
            input.skipElement();
            unsigned int time_ui = time(NULL) ;
            srand(time_ui);
            random_seed = rand();
        } else
        {
            random_seed = input.GetNextElementInt();
        }
        RNG_wrapper tmp(random_seed);
        rngw = tmp;
    }
    //std::cout << "GP->random_seed = " << GP->random_seed << "\n";

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



