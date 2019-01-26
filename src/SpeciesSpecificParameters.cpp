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


int SpeciesSpecificParameters::selectNonEmptyPatch(int firstPatchToLookAt, std::vector<int>& PSs, bool increasingOrder)
{
    int r = -1;
    if (firstPatchToLookAt != -1 && PSs[firstPatchToLookAt] > 0)
    {
        r = firstPatchToLookAt;
    } else
    {
        if (increasingOrder)
        {
            for (int patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
            {
                if (PSs[patch_index] > 0)
                {
                    r = patch_index;
                    break;
                }
            }
        } else
        {
            for (int patch_index = GP->PatchNumber - 1 ; patch_index >= 0 ; --patch_index)
            {
                if (PSs[patch_index] > 0)
                {
                    r = patch_index;
                    break;
                }
            }
        }
            
    }
    return r;
}

void SpeciesSpecificParameters::setRandomDistributions()
{
    // The distributions always include the bounds
    std::poisson_distribution<int>             tmp1(TotalRecombinationRate);
    this->rpois_nbRecombination                   = tmp1;
    std::poisson_distribution<int>             tmp2(T1_Total_Mutation_rate);
    this->T1_rpois_nbMut                          = tmp2;
    std::poisson_distribution<int>             tmp3(T2_Total_Mutation_rate);
    this->T2_rpois_nbMut                          = tmp3;
    std::poisson_distribution<int>             tmp3h(T3_Total_Mutation_rate);
    this->T3_rpois_nbMut                          = tmp3h;
    std::poisson_distribution<int>             tmp3i(T5_Total_Mutation_rate);
    this->T5_rpois_nbMut                          = tmp3i;
    //std::cout << "T1_Total_Mutation_rate = " << T1_Total_Mutation_rate << "\n";
    //std::cout << "T5_Total_Mutation_rate = " << T5_Total_Mutation_rate << "\n";
    
    std::uniform_int_distribution<int>         tmp4(0, TotalNbLoci-2);
    this->runiform_int_ForRecPos                  = tmp4;
    std::uniform_real_distribution<double>     tmp5(0, TotalRecombinationRate);
    this->runiform_double_ForRecPos               = tmp5;
    
    std::uniform_int_distribution<int>         tmp6(0, T1_nbBits-1);
    this->T1_runiform_int_ForMutPos               = tmp6;
    std::uniform_int_distribution<int>         tmp7(0, T2_nbChars-1);
    this->T2_runiform_int_ForMutPos               = tmp7;
    std::uniform_int_distribution<int>         tmp7h(0, T3_nbChars-1);
    this->T3_runiform_int_ForMutPos               = tmp7h;
    std::uniform_int_distribution<int>         tmp7i(0, T5_nbBits-1);
    this->T5_runiform_int_ForMutPos               = tmp7i;
    //std::cout << "T5_runiform_int_ForMutPos initialized with boundaries ["<<0<<";"<<T5_nbBits-1<<std::endl;
    
    std::uniform_real_distribution<double>     tmp8(0, T1_Total_Mutation_rate);
    this->T1_runiform_double_ForMutPos            = tmp8;
    std::uniform_real_distribution<double>     tmp9(0, T2_Total_Mutation_rate);
    this->T2_runiform_double_ForMutPos            = tmp9;
    std::uniform_real_distribution<double>     tmp10(0, T3_Total_Mutation_rate);
    this->T3_runiform_double_ForMutPos            = tmp10;
    std::uniform_real_distribution<double>     tmp11(0, T5_Total_Mutation_rate);
    this->T5_runiform_double_ForMutPos            = tmp11;
    
    
    std::uniform_real_distribution<double> tmp12(0,1);
    GP->random_0and1                       = tmp12;
    std::uniform_int_distribution<int>     tmp13(0,1);
    GP->random_0or1                        = tmp13;

#ifdef DEBUG
    std::cout << "Exiting 'setRandomDistributions' " << std::endl;
#endif
}


void SpeciesSpecificParameters::readDispMat(InputReader& input)
{
    dispersalData.readDispMat(input);
}

void SpeciesSpecificParameters::readGameteDispersal(InputReader& input)
{
    std::string s = input.GetNextElementString();
    if (s == "no" || s == "n" || s == "N" || s == "No" || s == "NO" || s == "0")
    {
        gameteDispersal = false;
    } else if (s == "yes" || s == "y" || s == "Y" || s == "Yes" || s == "YES" || s == "1")
    {
        gameteDispersal = true;
    } else
    {
        std::cout << "In option --gameteDispersal, expected either 'yes' (or 'y','Y', '1' and a few others) or 'no' (or 'n', 'N', '0' and a few others) but instead received " << s << "\n";
        abort();
    }
    
}

void SpeciesSpecificParameters::readOutputSFSbinSize(InputReader& input)
{
    while (input.IsThereMoreToRead())
    {
        if (input.PeakNextElementString() == "default")
        {
            input.skipElement();
            outputSFSbinSizes.push_back(0.01);
        } else
        {
            double x = input.GetNextElementDouble();            
            if (x <= 0 || x > 1)
            {
                std::cout << "While reading option --outputSFSbinSize, a bin size of " << x << " has been received. A bin size must be great than zero and must be lower (or equal) to one.\n";
                abort();
            }
            outputSFSbinSizes.push_back(x);
        }
    }
}


SpeciesSpecificParameters::SpeciesSpecificParameters(std::string sN, int sI)
: speciesName(sN), speciesIndex(sI)
{
    //std::cout << "Constructor of SpeciesSpecificParameters\n";
}

SpeciesSpecificParameters::~SpeciesSpecificParameters()
{
    //std::cout << "Destructor of SpeciesSpecificParameters\n";
}

/*
template <typename T>
T log(T x, double base)
{
    return std::log(x) / std::log(base);
}
*/

bool SpeciesSpecificParameters::setFromLocusToFitnessMapIndex_DecidingFunction(double sumOfProb, int nbLoci)
{
    assert(this->FitnessMapCoefficient == -9.0 || this->FitnessMapProbOfEvent == -9.0);
    assert(nbLoci > 0);
    
    if (this->FitnessMapCoefficient > 0.0)
    {
        return (nbLoci > 80 && (sumOfProb * nbLoci / 500) > this->FitnessMapCoefficient);
    } else if (this->FitnessMapProbOfEvent > 0.0 && nbLoci > this->FitnessMapMinimNbLoci)
    {
        return sumOfProb > this->FitnessMapProbOfEvent;
    }
    return false;
}

void SpeciesSpecificParameters::setFromLocusToFitnessMapIndex()
{
    assert(this->FromLocusToFitnessMapIndex.size() == 0);
    if (this->FitnessMapProbOfEvent < 0.0)
    {
        std::cout << "In 'setFromLocusToFitnessMapIndex', it appears that the inputted parameter for 'FitnessMapProbOfEvent' is negative. Value received is " << this->FitnessMapProbOfEvent << "\n";
        abort();
    }

    bool ShouldThereBeSeveralFitnessBlocks = (this->T1_isSelection && this->T1_isMultiplicitySelection) || (this->T5_isSelection && this->T5_isMultiplicitySelection) || this->T2_isSelection;

    double sumOfProb = 0.0;
    int FitnessMapIndex = 0;
    int nbLociInPiece = 0;
    for (int interlocus = 0 ; interlocus < (this->TotalNbLoci - 1) ; interlocus++)
    {
        nbLociInPiece++;

        // Get Locus info
        int T1_locus = this->FromLocusToTXLocus[interlocus].T1;
        int T2_locus = this->FromLocusToTXLocus[interlocus].T2;
        int T3_locus = this->FromLocusToTXLocus[interlocus].T3;
        int T4_locus = this->FromLocusToTXLocus[interlocus].T4;
        int T5_locus = this->FromLocusToTXLocus[interlocus].T5;
        int locusType;


        if (ShouldThereBeSeveralFitnessBlocks)
        {
            // Long security
            if (interlocus == 0)
            {
                assert(T1_locus + T2_locus + T3_locus + T4_locus + T5_locus == 1);
                if (T1_locus == 1)
                {
                    locusType = 1;
                } else if (T2_locus == 1)
                {
                    locusType = 2;
                } else if (T3_locus == 1)
                {
                    locusType = 3;
                } else if (T4_locus == 1)
                {
                    locusType = 4;
                } else if (T5_locus == 1)
                {
                    locusType = 5;
                } else
                {
                    std::cout << "Internal error in void SpeciesSpecificParameters::setFromLocusToFitnessMapIndex(): unknown locusType received.\n";
                    abort();
                }
            } else
            {
                assert(T1_locus + T2_locus + T3_locus + T4_locus + T5_locus - this->FromLocusToTXLocus[interlocus - 1].T1 - this->FromLocusToTXLocus[interlocus - 1].T2 - this->FromLocusToTXLocus[interlocus - 1].T3 - this->FromLocusToTXLocus[interlocus - 1].T4 - this->FromLocusToTXLocus[interlocus - 1].T5 == 1);
                if (T1_locus - this->FromLocusToTXLocus[interlocus - 1].T1 == 1)
                {
                    locusType = 1;
                } else if (T2_locus - this->FromLocusToTXLocus[interlocus - 1].T2 == 1)
                {
                    locusType = 2;
                } else if (T3_locus - this->FromLocusToTXLocus[interlocus - 1].T3 == 1)
                {
                    locusType = 3;
                } else if (T4_locus - this->FromLocusToTXLocus[interlocus - 1].T4 == 1)
                {
                    locusType = 4;
                } else if (T5_locus - this->FromLocusToTXLocus[interlocus - 1].T5 == 1)
                {
                    locusType = 5;
                }else
                {
                    std::cout << "Internal error in 'void SpeciesSpecificParameters::setFromLocusToFitnessMapIndex()'.\n";
                    abort();
                }
            }
            assert(locusType == this->FromLocusToTXLocus[interlocus].TType);

            // Get Recombination Rate info
            double r;
            if (locusType == 1 || locusType == 2 || locusType == 5)
            {
                if (this->RecombinationRate.size() == 1) // which is the case if the rate is constant between any two loci.
                {
                    r = this->RecombinationRate[0];
                } else
                {
                    if (interlocus == 0)
                    {
                        r = this->RecombinationRate[0];
                    } else
                    {
                        assert(this->RecombinationRate.size() > interlocus);
                        r = this->RecombinationRate[interlocus] - this->RecombinationRate[interlocus - 1];
                    }
                }
            } else
            {
                r = 0.0;
            }
            if (locusType == 5)
            {
                r *= this->FitnessMapT5WeightProbOfEvent;
            }
            

            // Get Mutation rate info
            /*double mu;
            bool ShouldIConsiderMutations = false;
            if (ShouldIConsiderMutations)
            {
                if (locusType == 1) // T1
                {
                    if (T1_locus == 0 || this->T1_MutationRate.size()==1)
                    {
                        mu = this->T1_MutationRate[0];
                    } else
                    {
                        assert(this->T1_MutationRate.size() > T1_locus);
                        mu = this->T1_MutationRate[T1_locus] - this->T1_MutationRate[T1_locus - 1];
                    }
                } else if (locusType == 2) // T2
                {
                    if (T2_locus == 0 || this->T2_MutationRate.size()==1)
                    {
                        mu = this->T2_MutationRate[0];
                    } else
                    {
                        assert(this->T2_MutationRate.size() > T2_locus);
                        mu = this->T2_MutationRate[T2_locus] - this->T2_MutationRate[T2_locus - 1];
                    }
                } else if (locusType == 3) // T3
                {
                    if (T3_locus == 0 || this->T3_MutationRate.size()==1)
                    {
                        mu = this->T3_MutationRate[0];
                    } else
                    {
                        assert(this->T3_MutationRate.size() > T3_locus);
                        mu = this->T3_MutationRate[T3_locus] - this->T3_MutationRate[T3_locus - 1];
                    }
                } else
                {
                    std::cout << "Internal error in 'void SpeciesSpecificParameters::setFromLocusToFitnessMapIndex()'.\n";
                    abort();
                }
            } else
            {
                mu = 0.0;
            }*/

            // Sum probabilities of an event
            sumOfProb += r;

            // test if needs to make a new boundary
            if (setFromLocusToFitnessMapIndex_DecidingFunction(sumOfProb, nbLociInPiece))
            {
                FitnessMapIndex++;
                sumOfProb = 0.0;
                FromLocusToTXLocusElement E(T1_locus, T2_locus, T3_locus, T4_locus, T5_locus, locusType);
                FromLocusToFitnessMapBoundaries.push_back(E);
                nbLociInPiece = 0;
            }
        } 
    }


    assert(this->FromLocusToTXLocus.size() == this->TotalNbLoci);
    NbElementsInFitnessMap = FitnessMapIndex + 1;

    // Add last element to boundaries
    FromLocusToTXLocusElement E(this->T1_nbBits, this->T2_nbChars, this->T3_nbChars, this->T4_nbBits, this->T5_nbBits, this->FromLocusToTXLocus[this->TotalNbLoci-1].TType);
    FromLocusToFitnessMapBoundaries.push_back(E);
    assert(FromLocusToFitnessMapBoundaries.size() == NbElementsInFitnessMap);

    // Create FromLocusToFitnessMapIndex from FromLocusToFitnessMapBoundaries
    int previous_b_interlocus = 0;
    for (int i = 0 ; i < FromLocusToFitnessMapBoundaries.size() ; i++)
    {
        auto& b = FromLocusToFitnessMapBoundaries[i];
        assert(b.T1 + b.T2 + b.T3 + b.T4 + b.T5 <= this->TotalNbLoci);
        int b_interlocus = b.T1 + b.T2 + b.T3 + b.T4 + b.T5;

        for (int j = previous_b_interlocus  ; j < b_interlocus ; j++)
        {
            FromLocusToFitnessMapIndex.push_back(i);        
        }
        previous_b_interlocus = b_interlocus;
    }
    //std::cout << "FromLocusToFitnessMapIndex.size() = " << FromLocusToFitnessMapIndex.size() << " this->TotalNbLoci = " << this->TotalNbLoci << "\n";
    assert(FromLocusToFitnessMapIndex.size() == this->TotalNbLoci);
    //std::cout << "NbElementsInFitnessMap = " << NbElementsInFitnessMap << "\n";
    //std::cout << "FromLocusToFitnessMapIndex.back() = " << FromLocusToFitnessMapIndex.back() << "\n";
    assert(NbElementsInFitnessMap == FromLocusToFitnessMapIndex.back()+1);
    
    // print to console
    if ( ShouldThereBeSeveralFitnessBlocks )
    {
        if (NbElementsInFitnessMap == 1)
        {
            std::cout << "\t\tThe fitnessMap of species '"<<this->speciesName<<"' contains " << NbElementsInFitnessMap << " element\n\n";
        } else
        {
            assert(NbElementsInFitnessMap > 1);
            std::cout << "\t\tThe fitnessMap of species '"<<this->speciesName<<"' contains " << NbElementsInFitnessMap << " elements\n\n";
        }
    }

    /*
    std::cout << "FromLocusToFitnessMapBoundaries";
    std::cout << "("<<FromLocusToFitnessMapBoundaries.size()<<"):\n";
    for (auto& e : FromLocusToFitnessMapBoundaries)
        std::cout << e.T1 + e.T2 + e.T3 << " ";
    std::cout << "\n";

    std::cout << "FromLocusToFitnessMapIndex" << std::flush;
    std::cout << "("<<FromLocusToFitnessMapIndex.size()<<"):\n" << std::flush;
    for (int i = 0 ; i < FromLocusToFitnessMapIndex.size() ; i++)
        std::cout << "FromLocusToFitnessMapIndex["<<i<<"] = " << FromLocusToFitnessMapIndex[i] << "\n" << std::flush;
    */
    

}

void SpeciesSpecificParameters::readT4_printTree(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option '--T4_printTree', the std::string that is read is: " << input.print() << std::endl;
#endif
    OutputFile file(input.GetNextElementString(), T4_printTree);
    T4Tree.indicateOutputFile(&file);
}

void SpeciesSpecificParameters::readLoci(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option '--Loci', the std::string that is read is: " << input.print() << std::endl;
#endif

    if (input.PeakNextElementString() == "NA" || input.PeakNextElementString() == "0")
    {
        std::cout << "In option --L (--Loci), received 'NA' or '0' at the beginning of input. I suppose you meant that you don't want any genetics. Do you really want that? SimBit is a population genetics software but I can understand that you might be interested in demographic and ecological processes only. It is probably easy to allow SimBit to run simulations without genetic architecture but I did not made sure that would not cause bugs. So either, you really need it and you ask Remi to implement the possiblity to run simulations without genetics or you just add a single locus without selection on it (do '--L 2 1 --T2_mu unif 0 --T2_fit 1' for example) and you just don't ask for outputs concerning the genetics.\n";
        abort();
    }
    
    this->T1_nbChars  = 0;
    this->T2_nbChars  = 0;
    this->T3_nbChars  = 0;
    this->T4_nbBits   = 0;
    this->T5_nbBits   = 0;
    this->T1_nbBits   = 0;
    this->TotalNbLoci = 0;

    while( input.IsThereMoreToRead() )
    {
        std::string Type = input.GetNextElementString();
        if (Type == "1" || Type == "T1" || Type == "t1")
        {
            int nbBits = input.GetNextElementInt();
            
            if (nbBits >= 1)
            {
                for (int i = 0 ; i < nbBits ; i++)
                {

                    this->FromT1LocusToLocus.push_back(this->TotalNbLoci);
                    this->T1_nbBits++;
                    this->TotalNbLoci++;
                    FromLocusToTXLocusElement CME(this->T1_nbBits, this->T2_nbChars, this->T3_nbChars, this->T4_nbBits, this->T5_nbBits, 1);
                    this->FromLocusToTXLocus.push_back(CME);
                }
            }
        } else if (Type == "2" || Type == "T2" || Type == "t2")
        {
            int nbBytes = input.GetNextElementInt();
            if (nbBytes >= 1)
            {
                for (int byte = 0 ; byte < nbBytes ; byte++)
                {
                    this->FromT2LocusToLocus.push_back(this->TotalNbLoci);
                    this->T2_nbChars++;
                    this->TotalNbLoci++;
                    FromLocusToTXLocusElement CME(this->T1_nbBits, this->T2_nbChars, this->T3_nbChars, this->T4_nbBits, this->T5_nbBits, 2);
                    this->FromLocusToTXLocus.push_back(CME);                
                }
            }
        } else if (Type == "3" || Type == "T3" || Type == "t3")
        {
            int nbBytes = input.GetNextElementInt();
            
            if (nbBytes >= 1)
            {
                for (int byte = 0 ; byte < nbBytes ; byte++)
                {
                    this->FromT3LocusToLocus.push_back(this->TotalNbLoci);
                    this->T3_nbChars++;
                    this->TotalNbLoci++;
                    FromLocusToTXLocusElement CME(this->T1_nbBits, this->T2_nbChars, this->T3_nbChars, this->T4_nbBits, this->T5_nbBits, 3);
                    this->FromLocusToTXLocus.push_back(CME);            
                }
            }
        } else if (Type == "4" || Type == "T4" || Type == "t4")
        {
            int nbBits = input.GetNextElementInt();
            
            if (nbBits >= 1)
            {
                for (int i = 0 ; i < nbBits ; i++)
                {

                    this->FromT4LocusToLocus.push_back(this->TotalNbLoci);
                    this->T4_nbBits++;
                    this->TotalNbLoci++;
                    FromLocusToTXLocusElement CME(this->T1_nbBits, this->T2_nbChars, this->T3_nbChars, this->T4_nbBits, this->T5_nbBits, 4);
                    this->FromLocusToTXLocus.push_back(CME);
                }
            }
        } else if (Type == "5" || Type == "T5" || Type == "t5")
        {
            int nbBits = input.GetNextElementInt();
            
            if (nbBits >= 1)
            {
                for (int i = 0 ; i < nbBits ; i++)
                {

                    this->FromT5LocusToLocus.push_back(this->TotalNbLoci);
                    this->T5_nbBits++;
                    this->TotalNbLoci++;
                    FromLocusToTXLocusElement CME(this->T1_nbBits, this->T2_nbChars, this->T3_nbChars, this->T4_nbBits, this->T5_nbBits, 5);
                    this->FromLocusToTXLocus.push_back(CME);
                }
            }
        } else
        {
            std::cout << "Sorry! Unknown Type if '--Loci' (Received Mode " << Type << ". Only three types of traits are recognized for the moment, type '1' (bitwise), type '2' (track nb mutations) and type '3' multidimensional fitness representation";
            abort();
        }
    }

    assert(this->T1_nbBits == this->FromT1LocusToLocus.size());
    assert(this->T2_nbChars == this->FromT2LocusToLocus.size());
    assert(this->T3_nbChars == this->FromT3LocusToLocus.size());
    assert(this->T4_nbBits == this->FromT4LocusToLocus.size());
    assert(this->T5_nbBits == this->FromT5LocusToLocus.size());
    assert(this->TotalNbLoci == this->FromLocusToTXLocus.size());
    assert(this->TotalNbLoci == T1_nbBits + T2_nbChars + T3_nbChars + T4_nbBits + T5_nbBits );

    this->T1_nbChars  = ceil(this->T1_nbBits/ (double) EIGHT);

    if (this->T1_nbBits > this->T1_nbChars * EIGHT || this->T1_nbBits < (this->T1_nbChars-1) * EIGHT)
    {
        std::cout << "Internal error in option '--Loci'. this->T1_nbBits > this->T1_nbChars * EIGHT. \nthis->T1_nbBits = "<< this->T1_nbBits << "\nthis->T1_nbChars = " << this->T1_nbChars << "\n";
        abort();
    }
    this->T1_nbBitsLastByte = this->T1_nbBits % EIGHT;
    if (this->T1_nbBitsLastByte == 0) {this->T1_nbBitsLastByte = EIGHT;}
    assert(this->T1_nbBitsLastByte >= 0 && this->T1_nbBitsLastByte <= EIGHT);
    assert(this->TotalNbLoci == this->T1_nbBits + this->T2_nbChars + this->T3_nbChars + this->T4_nbBits + this->T5_nbBits);
    if (this->TotalNbLoci <= 0)
    {
        std::cout << "In option '--Loci' you seem to have indicated zero loci. The platform requires to have at least one locus (of either Mode) to run.\n";
        abort();
    }

    if (this->T4_nbBits > 0)
    {
        if (this->nbSubGenerationsPerGeneration != 1)
        {
            std::cout << "You have asked for loci of type T4 via option --Loci (--L). You have also asked for having a number of subGeneration per generation different than 1 via option --nbSubGenerations (--nbSubGens). Sorry, the coalescent tree used to track T4 loci assumes on subgeneration per generation. It would be easy for Remi to get rid of this assumption. Please let me know!\n";
            abort();
        }
    }
    
    // Assert FromLocusToTXLocus is good
    assert(this->FromLocusToTXLocus.size() == this->TotalNbLoci);

    int previous = 0;
    int current;
    for (int i = 0 ; i < this->FromLocusToTXLocus.size() ; i++)
    {
        current = this->FromLocusToTXLocus[i].T1 + this->FromLocusToTXLocus[i].T2 + this->FromLocusToTXLocus[i].T3 + this->FromLocusToTXLocus[i].T4 + this->FromLocusToTXLocus[i].T5;
        //std::cout << "i = " << i << " current = " << current << " previous = " << previous << std::endl;
        if (previous + 1 != current)
        {
            std::cout << "Internal Error 1 when building FromLocusToTXLocus found for element " << i << ". current = " << current << " previous = " << previous << std::endl;
            abort();
        }
        if (current != i + 1 )
        {
            std::cout << "Internal Error 2 when building FromLocusToTXLocus found for element " << i << ". current = " << current << std::endl;
            abort();
        }     
        previous = current;
    }

    //std::cout << "this->T1_nbBits = " << this->T1_nbBits << "\nthis->T1_nbChars = " << this->T1_nbChars << "\nthis->T2_nbChars = " << this->T2_nbChars << "\nthis->T1_nbBitsLastByte = " << this->T1_nbBitsLastByte << "\n";
    //std::cout << "TotalNbLoci = " << TotalNbLoci << " FromLocusToTXLocus.size() = " << FromLocusToTXLocus.size() << "\n";

/*    std::cout << "FromT1LocusToLocus:\n";
    for (auto& e : FromT1LocusToLocus)
    {
        std::cout << e << " ";
    }
    std::cout << std::endl;*/
    
}

void SpeciesSpecificParameters::readCloningRate(InputReader& input)
{
    this->cloningRate = input.GetNextElementDouble();
    if (this->cloningRate > 1.0 || this->cloningRate < 0.0)
    {
        std::cout << "In --cloningRate, received a cloning rate for species '"<<this->speciesName<<"' that is either lower than zero or greater than one. cloningRate receveid is " << this->cloningRate << "\n";
        abort();
    }
}

void SpeciesSpecificParameters::readSelfingRate(InputReader& input)
{
    this->selfingRate = input.GetNextElementDouble();
    if ((this->selfingRate > 1.0 || this->selfingRate < 0.0) || this->selfingRate != -1.0)
    {
        std::cout << "In --selfingRate, received a selfing rate that is either lower than zero or greater than one. selfing rate receveid is " << this->selfingRate << ". Note that a selfingRate of -1 (default) refers to a simple Wright-Fisher model where the selfing rate is at 1/2N (a little bit lower with migration)\n";
        abort();
    }
}

void SpeciesSpecificParameters::readMatingSystem(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option '--matingSystem', the std::string that is read is: " << input.print() << std::endl;
#endif
    std::string s = input.GetNextElementString();
    if (s == "H" || s == "h")
    {
        this->malesAndFemales = false;
        this->sexRatio = 0.0;
    } else if (s == "MF" || s == "FM" || s == "fm"  || s == "mf")
    {
        this->malesAndFemales = true;
        this->sexRatio = input.GetNextElementDouble();

        if (this->sexRatio >= 1.0 || this->sexRatio <= 0.0)
        {
            std::cout << "For option '--matingSystem', the mating system chosen is 'fm' (males and females). The sex ratio received is " << this->sexRatio << ". A sex ratio must be between 0 and 1 (non-bounded). If set at zero, all individuals are female and no reproduction is possible. If set at one, all individuals are males and no reproduction is possible. If you want hermaphrodites, please indicate 'h' (or 'H') and no sex-ratio\n";
            abort();
        }


        assert(SSP->__patchCapacity.size() == GP->__GenerationChange.size());
        assert(GP->__PatchNumber.size() == GP->__GenerationChange.size());

        if (SSP->fecundityForFitnessOfOne == -1.0)
        {
            for (int generation_index = 0 ; generation_index < GP->__GenerationChange.size() ; generation_index++)
            {
                assert(SSP->__patchCapacity[generation_index].size() == GP->__PatchNumber[generation_index]);
                for (double N : SSP->__patchCapacity[generation_index])
                {
                    int nbMales = (int) (sexRatio * N + 0.5);
                    int nbFemales = N - nbMales;
                    //std::cout << "nbMales = "<<nbMales<<" nbFemales = "<<nbFemales<<"\n";
                    if (nbMales <= 0 || nbFemales <= 0)
                    {
                        std::cout << "For option '--matingSystem', the mating system chosen is 'fm' (males and females). The sex ratio received for species "<< this->speciesName<< " is " << this->sexRatio << ". At the " << generation_index << "th generation index (after generation " << GP->__GenerationChange[generation_index] << "), the patchCapacity multiplied by the sex ratio lead to either the number of males or the number of females to be zero (nbMales = "<<nbMales<<" nbFemales = "<<nbFemales<<"). This error message would not appear if the fecundity (from option --fec) was set to -1.0 (fecundityForFitnessOfOne = "<<SSP->fecundityForFitnessOfOne<<"). If fecundity was not set to -1.0, then the patch size would just reach 0. SimBit enforces that there is at least one female and one male per patch when fecundity is different from -1.0. (this error message was written at a time when sex ratio is deterministic. If at some point SimBit evovles to have stochasticity in the sex ratio, then I should remove this error message and deal with cases where the sex ratio is 0 or 1).\n";
                        abort(); 
                    }
                }
            }
        }
    } else
    {
        std::cout << "For option '--matingSystem', the mating system received is " << s << ". Only types 'fm' (or 'mf', 'FM','MF') and 'h' (or 'H') are recognized. Sorry.\n";
        abort();
    }
}

void SpeciesSpecificParameters::readDispWeightByFitness(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option '--DispWeightByFitness', the std::string that is read is: " << input.print() << std::endl;
#endif
    if (input.PeakNextElementString() == "NA")
    {
        if (GP->PatchNumber != 1)
        {
            std::cout << "For option '--DispWeightByFitness', received 'NA' but the number of patches is different from 1 (GP->PatchNumber = "<<GP->PatchNumber<<"). Please just write '0' if you don't want disperal to be weighted by patch mean fitness.\n";
            abort();
        }
        assert(GP->PatchNumber == 1);
        this->DispWeightByFitness = false;
        input.skipElement();
    } else
    {
        int value = input.GetNextElementInt();
        if (value == 0)
        {
            this->DispWeightByFitness = false;
        } else if (value == 1)
        {
            this->DispWeightByFitness = true;
        } else
        {
            std::cout << "'--DispWeightByFitness' must be either 0 (not weighted by fitness) or 1 (weighted by fitness). Received " <<  value << std::endl;
            abort();
        }
    }
}

void SpeciesSpecificParameters::readPloidy(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option '--ploidy', the std::string that is read is: " << input.print() << std::endl;
#endif    
    this->ploidy = input.GetNextElementInt();
    if (this->ploidy!=2)
    {
        std::cout << "Sorry, for the moment ploidy can only be equal to 2" << std::endl;
        abort();
    }
}

void SpeciesSpecificParameters::readResetGenetics(InputReader& input)
{
    #ifdef DEBUG
    std::cout << "For option '--resetGenetics', the std::string that is read is: " << input.print() << std::endl;
#endif 

    // input example --resetGenetics @S0 event 100 T1 setTo0 lociList 2 4 6 8 haplo both patch 0 allInds patch 2 indsList 1 2 3 4 5 
    //                                          ^generation
    
    int eventIndex = 0;
    while (input.IsThereMoreToRead())
    {
        std::string eventString = input.GetNextElementString();
        if (eventString != "event")
        {
            std::cout << "For option '--resetGenetics', for the "<<eventIndex<<"th event (species "<<SSP->speciesName<<"), expected the keyword 'event' but got "<<eventString<<"\n";
            abort();
        }

        // generation
        int generation = input.GetNextElementInt();
        if (generation < 0)
        {
            std::cout << "For option '--resetGenetics', the generation must be 0 or larger. Received generation = " << generation << ".\n";
            abort();
        }
        int generation_index = std::upper_bound(GP->__GenerationChange.begin(), GP->__GenerationChange.end(), generation) - GP->__GenerationChange.begin() - 1;

        // loci type
        std::string locusType = input.GetNextElementString();
        if (locusType != "T1" && locusType != "T2" && locusType != "T3"  && locusType != "T5")
        {
            std::cout << "For option '--resetGenetics', " << eventIndex << "th event, for the "<<eventIndex<<"th event (generation "<<generation<<"; species "<<SSP->speciesName<<"), the locus type is  '" << locusType <<"'. Sorry only mutation types 'T1', 'T2', 'T3' and 'T5' (but not T4) are accepted.\n";
            abort();
        }

        // mutation type
        char mutationType;
        if (locusType == "T1" || locusType == "T5")
        {
            std::string string_mutationType = input.GetNextElementString();
            if (string_mutationType == "setTo0")
            {
                mutationType = 0;
            } else if (string_mutationType == "setTo1")
            {
                mutationType = 1;
            } else if (string_mutationType == "toggle")
            {
                mutationType = 2;
            } else
            {
                std::cout << "For option '--resetGenetics', " << eventIndex << "th event, for the "<<eventIndex<<"th event (generation "<<generation<<"; species "<<SSP->speciesName<<") for loci of type "<<locusType<<", the mutation type received is '" << string_mutationType <<"'. Sorry only mutation types 'setTo0', 'setTo1' and 'toggle' are accepted.\n";
                abort();
            }
        } else
        {
            mutationType = -1;
        }

        // Loci indices
        std::string lociKeyword = input.GetNextElementString();
        if (lociKeyword != "lociList" && lociKeyword != "allLoci")
        {
            std::cout << "For option '--resetGenetics', " << eventIndex << "th event, for the "<<eventIndex<<"th event (generation "<<generation<<"; species "<<SSP->speciesName<<") expected the keyword 'lociList' or 'allLoci' but instead received "<<lociKeyword<<".\n";
                abort();
        }

        std::vector<int> T1loci;
        std::vector<int> T2loci;
        std::vector<int> T3loci;
        std::vector<int> T5loci;
        if (lociKeyword == "lociList")
        {
            while (input.PeakNextElementString() != "haplo")
            {
                if (locusType == "T1")
                {
                    int T1locus = input.GetNextElementInt();
                    if (T1locus < 0)
                    {
                        std::cout << "For option --resetGenetics received T1locus index of "<<T1locus<<" for an event happening at generation "<<generation<<".\n";
                        abort();
                    }
                    if (T1locus >= SSP->T1_nbBits)
                    {
                        std::cout << "For option --resetGenetics for an event happening at generation "<<generation<<", received T1locus index of "<<T1locus<<" while there are only "<<SSP->T1_nbBits<<" T1 loci (in bits) in total. As a reminder, the first locus has index 0.\n";
                        abort();
                    }
                    T1loci.push_back(T1locus);
                } else if (locusType == "T2")
                {
                    int T2locus = input.GetNextElementInt();
                    if (T2locus < 0)
                    {
                        std::cout << "For option --resetGenetics received T2 locus index of "<<T2locus<<" for an event happening at generation "<<generation<<".\n";
                        abort();
                    }
                    if (T2locus >= SSP->T2_nbChars)
                    {
                        std::cout << "For option --resetGenetics for an event happening at generation "<<generation<<", received T2 locus index of "<<T2locus<<" while there are only "<<SSP->T2_nbChars<<" T2 loci in total. As a reminder, the first locus has index 0.\n";
                        abort();
                    }
                    T2loci.push_back(T2locus);
                } else if (locusType == "T3")
                {
                    int T3locus = input.GetNextElementInt();
                    if (T3locus < 0)
                    {
                        std::cout << "For option --resetGenetics received T3locus index of "<<T3locus<<" for an event happening at generation "<<generation<<".\n";
                        abort();
                    }
                    if (T3locus >= SSP->T3_nbChars)
                    {
                        std::cout << "For option --resetGenetics for an event happening at generation "<<generation<<", received T3 locus index of "<<T3locus<<" while there are only "<<SSP->T3_nbChars<<" T3 loci in total. As a reminder, the first locus has index 0.\n";
                        abort();
                    }
                    T3loci.push_back(T3locus);
                } else if (locusType == "T5")
                {
                    int T5locus = input.GetNextElementInt();
                    if (T5locus < 0)
                    {
                        std::cout << "For option --resetGenetics received T5locus index of "<<T5locus<<" for an event happening at generation "<<generation<<".\n";
                        abort();
                    }
                    if (T5locus >= SSP->T5_nbBits)
                    {
                        std::cout << "For option --resetGenetics for an event happening at generation "<<generation<<", received T5locus index of "<<T5locus<<" while there are only "<<SSP->T5_nbBits<<" T5 loci in total. As a reminder, the first locus has index 0.\n";
                        abort();
                    }
                    T5loci.push_back(T5locus);
                }
            }
        } else // lociKeyword == "allLoci"
        {
            if (locusType == "T1")
            {
                for (int T1locus = 0 ; T1locus < SSP->T1_nbBits; T1locus++)
                {
                    T1loci.push_back(T1locus);
                }
            } else if (locusType == "T2")
            {
                for (int T2locus = 0 ; T2locus < SSP->T2_nbChars; T2locus++)
                {
                    T2loci.push_back(T2locus);
                }
            } else if (locusType == "T3")
            {
                for (int T3locus = 0 ; T3locus < SSP->T3_nbChars; T3locus++)
                {
                    T3loci.push_back(T3locus);
                }
            } else if (locusType == "T5")
            {
                for (int T5locus = 0 ; T5locus < SSP->T5_nbBits; T5locus++)
                {
                    T5loci.push_back(T5locus);
                }
            }
        }

        // haplotypes concerned
        std::string haploString = input.GetNextElementString();
        if (haploString != "haplo")
        {
            std::cout << "For option '--resetGenetics', " << eventIndex << "th event, for the "<<eventIndex<<"th event (generation "<<generation<<"; species "<<SSP->speciesName<<"), expected the keyword 'haplo' but instead received '" << haploString <<"'.\n";
            abort();
        }
        std::string haploDescription = input.GetNextElementString();
        std::vector<int> haplotypes;
        if (haploDescription == "0" || haploDescription == "f" || haploDescription == "F")
        {
            haplotypes.push_back(0);
        } else if (haploDescription == "1" || haploDescription == "m" || haploDescription == "M")
        {
            haplotypes.push_back(1);
        } else if (haploDescription == "both")
        {
            haplotypes.push_back(0);
            haplotypes.push_back(1);
        }



        // patches and individuals
        std::vector<std::vector<int>> individuals;
        std::vector<int> patches;
        while (input.IsThereMoreToRead() && input.PeakNextElementString() != "event")
        {
            // patch
            std::string patchString = input.GetNextElementString();
            if (patchString != "patch")
            {
                std::cout << "For option '--resetGenetics', " << eventIndex << "th event, for the "<<eventIndex<<"th event (generation "<<generation<<"; species "<<SSP->speciesName<<"), expected the keyword 'patch' but instead received '" << patchString <<"'.\n";
                abort();
            }

            

            int patch_index = input.GetNextElementInt();
            if (patch_index < 0 || patch_index > GP->__PatchNumber[generation_index])
            {
                std::cout << "For option '--resetGenetics', " << eventIndex << "th event, for the "<<eventIndex<<"th event (generation "<<generation<<"; species "<<SSP->speciesName<<"), after the keyword 'patch' received the patch index '" << patch_index <<"'. At the generation "<<generation<<" there are "<<GP->__PatchNumber[generation_index]<< "patches. The patch index must be lower than the number of patches (because it is a zero based counting) and positive.\n";
                abort();
            }
            patches.push_back(patch_index);


            // individuals
            std::string modeOfEntryOfIndividuals = input.GetNextElementString();

            if (modeOfEntryOfIndividuals == "allInds")
            {

                assert(SSP->__patchCapacity[generation_index].size() == GP->__PatchNumber[generation_index]);

                std::vector<int> onePatchIndviduals;
                onePatchIndviduals.reserve(SSP->__patchCapacity[generation_index][patch_index]);
                for (int ind_index = 0 ; ind_index < SSP->__patchCapacity[generation_index][patch_index]; ind_index++)
                {
                    onePatchIndviduals.push_back(ind_index);
                }
                individuals.push_back(onePatchIndviduals);

            } else if (modeOfEntryOfIndividuals == "indsList")
            {
                std::vector<int> onePatchIndviduals;
                while (input.IsThereMoreToRead() && input.PeakNextElementString() != "event" && input.PeakNextElementString() != "patch")
                {
                    onePatchIndviduals.push_back(input.GetNextElementInt());
                }
                if (onePatchIndviduals.size() == 0)
                {
                    std::cout << "For option '--resetGenetics', " << eventIndex << "th event, for the "<<eventIndex<<"th event (generation "<<generation<<"; species "<<SSP->speciesName<<"), for patch index "<< patch_index <<" there seems to have no individuals indicated.\n";
                    abort();
                }
                individuals.push_back(onePatchIndviduals);


            } else
            {
                std::cout << "For option '--resetGenetics', " << eventIndex << "th event, for the "<<eventIndex<<"th event (generation "<<generation<<"; species "<<SSP->speciesName<<"), for patch index "<< patch_index <<" expected mode of entry of individuals (either 'indsList' or 'allInds') but instead received '" << modeOfEntryOfIndividuals <<"'.\n";
                abort();
            }
        }
        if (patches.size() == 0)
        {
            std::cout << "For option '--resetGenetics', " << eventIndex << "th event, for the "<<eventIndex<<"th event (generation "<<generation<<"; species "<<SSP->speciesName<<"), there seems to have no patch indicated.\n";
            abort();
        }


        resetGenetics.addEvent(
            ResetGeneticsEvent(
                generation,
                mutationType,
                T1loci,
                T2loci,
                T3loci,
                T5loci,
                patches,
                haplotypes,
                individuals
            )
        );
        eventIndex++;
    }
}

void SpeciesSpecificParameters::readfecundityForFitnessOfOne(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option '--fec', the std::string that is read is: " << input.print() << std::endl;
#endif
    this->fecundityForFitnessOfOne = input.GetNextElementDouble();
    //std::cout << "fecundity set to " << this->fecundityForFitnessOfOne << "\n";

    if (this->fecundityForFitnessOfOne > 10000.0)
    {
        std::cout << "this->fecundityForFitnessOfOne received is greater than " << 10000.0 << " (is " << this->fecundityForFitnessOfOne << "). If you want the patch to always be at carrying capacity, just use fecundityForFitnessOfOne = -1. It will make everything faster. Change the code if you really want to use such high fecundity\n";
        abort();
    } else if (this->fecundityForFitnessOfOne < 0.0)
    {
        if (this->fecundityForFitnessOfOne != -1.0)
        {
            std::cout << "this->fecundityForFitnessOfOne received is negative but different from -1.0 (which would mean patch size is always at carrying capacity).\n";
            abort();
        }
    }
}

void SpeciesSpecificParameters::readResetTrackedT1Muts(InputReader& input)
{
    IsThereSelection();
    if (input.PeakNextElementString() == "default")
    {
        input.skipElement();
        if (T1_isSelection)
        {
            if (T1_isMultiplicitySelection)
            {
                recomputeLociOverWhichFitnessMustBeComputedEveryHowManyGenerations = -1;
            } else
            {
                
                double averageMu = T1_MutationRate.back() / T1_nbBits;
                
                double a = 0.0;
                for (int i = 1 ; i < TotalpatchCapacity ; i++)
                {
                    a += 1/i;
                }
                //std::cout << "TotalpatchCapacity = " << TotalpatchCapacity << " averageMu = " << averageMu << " a = " <<  a << "\n";
                //std::cout << "TotalpatchCapacity * averageMu * 4 * a = " << TotalpatchCapacity * averageMu * 4 * a << "\n";
                if (TotalpatchCapacity * averageMu * 4 * a < 1.0)
                {
                    recomputeLociOverWhichFitnessMustBeComputedEveryHowManyGenerations = 20;
                } else
                {
                    recomputeLociOverWhichFitnessMustBeComputedEveryHowManyGenerations = -1;
                }
            }
        } else
        {
            recomputeLociOverWhichFitnessMustBeComputedEveryHowManyGenerations = -1;
        }  
    } else
    {
        recomputeLociOverWhichFitnessMustBeComputedEveryHowManyGenerations = input.GetNextElementInt();
        if (recomputeLociOverWhichFitnessMustBeComputedEveryHowManyGenerations == 0 || recomputeLociOverWhichFitnessMustBeComputedEveryHowManyGenerations < -1)
        {
            std::cout << "In --T1_resetTrackedT1Muts, the number of generations received is " << recomputeLociOverWhichFitnessMustBeComputedEveryHowManyGenerations << ". Please use a number that is strictly positive or -1 (-1 means, to loop through all the sites systematically).\n";
            abort();
        }
    }
        

   /* std::string s = input.GetNextElementString();
    if (s == "n" || s == "N" || s == "no" || s == "0")
    {
        allowToCorrectRecomputeLociOverWhichFitnessMustBeComputedEveryHowManyGenerations = false;
    } else if (s == "y" || s == "Y" || s == "yes" || s == "1")
    {
        std::cout << "In --T1_resetTrackedT1Muts, the string received for the second argument is " << s << ". This string is valid. However, the current version of SimBit is not able to modify the generations at which the listing of tracked mutations is done. Sorry! The only possible entry for the second argument is therefore 'no' (or 'n' or '0' or a few other equivalents).\n";
        abort();
        allowToCorrectRecomputeLociOverWhichFitnessMustBeComputedEveryHowManyGenerations = true;
    } else
    {
        std::cout << "In --T1_resetTrackedT1Muts, expected yes or no (or a few equivalents) as second argument but instead received " << s << ".\n";
        abort();
    }
*/
    
    //std::cout << "recomputeLociOverWhichFitnessMustBeComputedEveryHowManyGenerations = "<< recomputeLociOverWhichFitnessMustBeComputedEveryHowManyGenerations << "\n";
}

void SpeciesSpecificParameters::readFitnessMapInfo(InputReader& input)
{
    this->FitnessMapMinimNbLoci = 0; // This is set to a different value only if used default
    std::string mode;
    if (!this->T1_isMultiplicitySelection && !this->T5_isMultiplicitySelection &&  this->T2_nbChars == 0)
    {
        mode = "prob";
        this->FitnessMapProbOfEvent = DBL_MAX;
        while (input.IsThereMoreToRead())
        {
            input.skipElement();
        }
        this->FitnessMapT5WeightProbOfEvent = 0.1;
    } else
    {
        mode = input.GetNextElementString();
        if (mode == "prob")
        {
            this->FitnessMapProbOfEvent = input.GetNextElementDouble();
            if (this->FitnessMapProbOfEvent <= 0.0)
            {
                std::cout << "For species " << this->speciesName << " when reading --FitnessMapInfo ( mode = " << mode << ") received not strictly positive value for 'FitnessMapProbOfEvent' (received "<<this->FitnessMapProbOfEvent<<").\n";
                abort();        
            }

            this->FitnessMapT5WeightProbOfEvent = input.GetNextElementDouble();
            if (this->FitnessMapT5WeightProbOfEvent <= 0.0)
            {
                std::cout << "For species " << this->speciesName << " when reading --FitnessMapInfo ( mode = " << mode << ") received not strictly positive value for 'FitnessMapT5WeightProbOfEvent'(received "<<this->FitnessMapT5WeightProbOfEvent<<").\n";
                abort();        
            }

            this->FitnessMapCoefficient = -9.0;
        } else if (mode == "coef")
        {
            std::cout << "Sorry for --FitnessMapInfo, mode coef is not available on this version.\n";
            abort();
            this->FitnessMapProbOfEvent = -9.0;
            this->FitnessMapCoefficient = input.GetNextElementDouble();
            if (this->FitnessMapCoefficient <= 0.0)
            {
                std::cout << "For species " << this->speciesName << " when reading --FitnessMapInfo ( mode = " << mode << ") received a not strictly positive value (received "<<this->FitnessMapCoefficient<<").\n";
                abort();        
            }

            this->FitnessMapT5WeightProbOfEvent = input.GetNextElementDouble();
            if (this->FitnessMapT5WeightProbOfEvent <= 0.0)
            {
                std::cout << "For species " << this->speciesName << " when reading --FitnessMapInfo ( mode = " << mode << ") received not strictly positive value for 'FitnessMapT5WeightProbOfEvent'(received "<<this->FitnessMapT5WeightProbOfEvent<<").\n";
                abort();        
            }
        } else
        {
            std::cout << "For species " << this->speciesName << " when reading --FitnessMapInfo, received mode = " << mode << ". Sorry only modes 'prob' and 'coef' are allowed\n";
            std::cout << "Hum.... actually even mode 'coef' is not allowed on this version due to potential misfunctioning and slow down. Sorry! Please use 'prob'.\n";
            abort(); 
        }
    }
    //std::cout << "\n\n\n\n\n\nin SSP : FitnessMapT5WeightProbOfEvent = " << FitnessMapT5WeightProbOfEvent << "\n\n\n\n\n\n\n";
}

void SpeciesSpecificParameters::readInitialpatchSize(InputReader& input)
{
    int TotalNbIndsAtStart = 0;

    std::string Mode = input.GetNextElementString();

    if (Mode == "A")
    {
        for (int patch_index = 0 ; patch_index < GP->__PatchNumber[0] ; patch_index++)
        {
            int x = input.GetNextElementInt();
            int ps;

            if (x < 0)
            {
                ps = this->__patchCapacity[0][patch_index];
            } else
            {
                ps = x;
            }

            if (ps > this->__patchCapacity[0][patch_index])
            {
                std::cout << "For patch " << patch_index << " (zero-based counting), the capacity is " << this->__patchCapacity[0][patch_index] << " but the initial patch size is " << ps << "\n";
                abort();
            }
            if (this->fecundityForFitnessOfOne == -1 && ps != this->__patchCapacity[0][patch_index])
            {
                std::cout << "For patch " << patch_index << " (zero-based counting), the capacity is " << this->__patchCapacity[0][patch_index] << " and the initial patch size is " << ps << ". When the fecundityForFitnessOfOne differs from -1 (which is the default and means infinite), the initial patch size must be at carrying capacity. This error message should only appear if you set both --InitialpatchSize and --fec yourself and made a mismatch between the initial patch size and the carrying capacity at the zeroth generation (indicated by option --N)\n";
                abort();
            }
            this->patchSize.push_back(ps);
            TotalNbIndsAtStart+=ps;
        }
    } else if (Mode == "unif")
    {
        int x = input.GetNextElementInt();
        
        for (int patch_index = 0 ; patch_index < GP->__PatchNumber[0] ; patch_index++)
        {
            int ps; // stands for patchSize
            if (x < 0)
            {
                ps = this->__patchCapacity[0][patch_index];
            } else
            {
                ps = x;
            }

            if (ps > this->__patchCapacity[0][patch_index])
            {
                std::cout << "For patch " << patch_index << " (zero-based counting), the capacity is " << this->__patchCapacity[0][patch_index] << " but the initial patch size is " << ps << "\n";
            }
            if (this->fecundityForFitnessOfOne == -1 && ps != this->__patchCapacity[0][patch_index])
            {
                std::cout << "For patch " << patch_index << " (zero-based counting), the capacity is " << this->__patchCapacity[0][patch_index] << " and the initial patch size is " << ps << ". When the fecundityForFitnessOfOne differs from -1 (which is the default and means infinite), the initial patch size must be at carrying capacity. This error message should only appear if you set both --InitialpatchSize and --fec yourself and made a mismatch between the initial patch size and the carrying capacity at the zeroth generation (indicated by option --N)\n";
                abort();
            }
            this->patchSize.push_back(ps);
            TotalNbIndsAtStart+=ps;
        }

    } else
    {
        std::cout << "For option --InitialpatchSize, mode received is " << Mode << ". Sorry only modes 'A' and 'unif' are recognized.\n";
        abort();
    }

    assert(this->patchSize.size() == GP->__PatchNumber[0]);

    if (TotalNbIndsAtStart == 0)
    {
        std::cout << "You asked for a total of exactly 0 individual (over all patches) at initialization (via option '--InitialpatchSize')! Nothing can be simulated for this species!\n" << std::endl;
        abort();
    }
    assert(TotalNbIndsAtStart > 0);
}


void SpeciesSpecificParameters::readT1_Initial_AlleleFreqs(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option --T1_Initial_AlleleFreqs, the std::string that is read is: " << input.print() << std::endl;
#endif
    if (input.PeakNextElementString() == "NA")
    {
        if (this->T1_nbBits != 0)
        {
            std::cout << "For option --T1_Initial_AlleleFreqs, received 'NA' however there are T1 loci as indiciated by --L (--Loci) (this->T1_nbBits = "<<this->T1_nbBits<<")" << "\n";
            abort();
        }
        input.skipElement();
    } else
    {
        std::string Mode;

        Mode = input.GetNextElementString();

        this->T1_Initial_AlleleFreqs_AllZeros = false; // Set to true laster if received "AllZeros"
        this->T1_Initial_AlleleFreqs_AllOnes = false;  // Set to true laster if received "AllOnes"
        
        if (Mode.compare("A") == 0)
        {
            for (int patch_index = 0 ; patch_index < GP->__PatchNumber[0] ; ++patch_index)
            {
                int nbRepeatsLeft = 0;
                double freq;
                std::vector<double> OneLine;
                for (int T1Locus = 0 ; T1Locus < this->T1_nbBits ; T1Locus++)
                {
                    if (nbRepeatsLeft == 0 && input.PeakNextElementString() == "R")
                    {
                        input.skipElement();
                        freq = input.GetNextElementDouble();
                        nbRepeatsLeft = input.GetNextElementInt();
                        if (nbRepeatsLeft < 0)
                        {
                            std::cout << "In '--T1_Initial_AlleleFreqs', received a number of repeats lower than zero. Number of repeats received is " << nbRepeatsLeft << "\n";
                            abort();
                        }
                    }
                    if (nbRepeatsLeft == 0)
                    {
                        freq = input.GetNextElementDouble();
                    } else
                    {
                        nbRepeatsLeft--;
                    }
                    
                    if (freq < 0.0 || freq > 1.0)
                    {
                        std::cout << "In '--T1_Initial_AlleleFreqs', received an impossible allele frequency. Received '" << freq << "'\n";
                        abort();
                    }
                    //std::cout << "T1Locus = " << T1Locus << " x = " << x << " nbRepeatsLeft = "<<nbRepeatsLeft<<"\n";
                    OneLine.push_back(freq);
                }
                assert(nbRepeatsLeft >= 0);
                if (nbRepeatsLeft != 0)
                {
                    std::cout << "In '--T1_Initial_AlleleFreqs', too many values received!'\n";
                    abort();
                }
                assert(OneLine.size() == this->T1_nbBits);
                this->T1_Initial_AlleleFreqs.push_back(OneLine);
            }
        } else if (Mode.compare("Shift")==0)
        {
            int shiftPosition = input.GetNextElementInt();
            if (shiftPosition < 0 || shiftPosition > GP->__PatchNumber[0])
            {
                std::cout << "In '--T1_Initial_AlleleFreqs', received an impossible shiftPosition. Received '" << shiftPosition << ". For information, the simulation is starting with " << GP->__PatchNumber[0] << " as indicated in option '--PN (--PatchNumber)'\n";
                abort();
            }
            
            for (int patch_index = 0 ; patch_index < GP->__PatchNumber[0] ; ++patch_index)
            {
                std::vector<double> line;
                for (int T1Locus = 0 ; T1Locus < this->T1_nbBits ; T1Locus++)
                {
                    if (patch_index < shiftPosition)
                    {
                        line.push_back(0.0);
                    } else
                    {
                        line.push_back(1.0);
                    }
                }
                assert(line.size() == this->T1_nbBits);
                this->T1_Initial_AlleleFreqs.push_back(line);
            }
        } else if (Mode.compare("AllZeros") == 0)
        {
            this->T1_Initial_AlleleFreqs_AllZeros = true;
        } else if (Mode.compare("AllOnes") == 0)
        {
            this->T1_Initial_AlleleFreqs_AllOnes = true;
        } else
        {
            std::cout << "Sorry, for '--T1_Initial_AlleleFreqs', the Mode " << Mode << " has not been implemented yet. Only Modes 'A', 'AllZeros', 'AllOnes' and 'Shift' are accepted for the moment.";
            abort();
        }
    }

    // Security check
    if (!this->T1_Initial_AlleleFreqs_AllZeros && !this->T1_Initial_AlleleFreqs_AllOnes)
    {
        assert(this->T1_Initial_AlleleFreqs.size() == GP->__PatchNumber[0]);
        for (int patch_index = 0 ; patch_index < GP->__PatchNumber[0] ; ++patch_index)
        {
            assert(this->T1_Initial_AlleleFreqs[patch_index].size() == this->T1_nbBits);
        }
    }
}


void SpeciesSpecificParameters::readT5_Initial_AlleleFreqs(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option --T5_Initial_AlleleFreqs, the std::string that is read is: " << input.print() << std::endl;
#endif
    if (input.PeakNextElementString() == "NA")
    {
        if (this->T5_nbBits != 0)
        {
            std::cout << "For option --T5_Initial_AlleleFreqs, received 'NA' however there are T5 loci as indiciated by --L (--Loci) (this->T5_nbBits = "<<this->T5_nbBits<<")" << "\n";
            abort();
        }
        input.skipElement();
    } else
    {
        std::string Mode;

        Mode = input.GetNextElementString();

        this->T5_Initial_AlleleFreqs_AllZeros = false; // Set to true laster if received "AllZeros"
        this->T5_Initial_AlleleFreqs_AllOnes = false;  // Set to true laster if received "AllOnes"
        
        if (Mode.compare("A") == 0)
        {
            for (int patch_index = 0 ; patch_index < GP->__PatchNumber[0] ; ++patch_index)
            {
                int nbRepeatsLeft = 0;
                double freq;
                std::vector<double> OneLine;
                for (int T5Locus = 0 ; T5Locus < this->T5_nbBits ; T5Locus++)
                {
                    if (nbRepeatsLeft == 0 && input.PeakNextElementString() == "R")
                    {
                        input.skipElement();
                        freq = input.GetNextElementDouble();
                        nbRepeatsLeft = input.GetNextElementInt();
                        if (nbRepeatsLeft < 0)
                        {
                            std::cout << "In '--T5_Initial_AlleleFreqs', received a number of repeats lower than zero. Number of repeats received is " << nbRepeatsLeft << "\n";
                            abort();
                        }
                    }
                    if (nbRepeatsLeft == 0)
                    {
                        freq = input.GetNextElementDouble();
                    } else
                    {
                        nbRepeatsLeft--;
                    }
                    
                    if (freq < 0.0 || freq > 1.0)
                    {
                        std::cout << "In '--T5_Initial_AlleleFreqs', received an impossible allele frequency. Received '" << freq << "'\n";
                        abort();
                    }
                    //std::cout << "T1Locus = " << T1Locus << " x = " << x << " nbRepeatsLeft = "<<nbRepeatsLeft<<"\n";
                    OneLine.push_back(freq);
                }
                assert(nbRepeatsLeft >= 0);
                if (nbRepeatsLeft != 0)
                {
                    std::cout << "In '--T1_Initial_AlleleFreqs', too many values received!'\n";
                    abort();
                }
                assert(OneLine.size() == this->T5_nbBits);
                this->T5_Initial_AlleleFreqs.push_back(OneLine);
            }
        } else if (Mode.compare("Shift")==0)
        {
            int shiftPosition = input.GetNextElementInt();
            if (shiftPosition < 0 || shiftPosition > GP->__PatchNumber[0])
            {
                std::cout << "In '--T5_Initial_AlleleFreqs', received an impossible shiftPosition. Received '" << shiftPosition << ". For information, the simulation is starting with " << GP->__PatchNumber[0] << " as indicated in option '--PN (--PatchNumber)'\n";
                abort();
            }
            
            for (int patch_index = 0 ; patch_index < GP->__PatchNumber[0] ; ++patch_index)
            {
                std::vector<double> line;
                for (int T5Locus = 0 ; T5Locus < this->T5_nbBits ; T5Locus++)
                {
                    if (patch_index < shiftPosition)
                    {
                        line.push_back(0.0);
                    } else
                    {
                        line.push_back(1.0);
                    }
                }
                assert(line.size() == this->T5_nbBits);
                this->T5_Initial_AlleleFreqs.push_back(line);
            }
        } else if (Mode.compare("AllZeros") == 0)
        {
            this->T5_Initial_AlleleFreqs_AllZeros = true;
        } else if (Mode.compare("AllOnes") == 0)
        {
            this->T5_Initial_AlleleFreqs_AllOnes = true;
        } else
        {
            std::cout << "Sorry, for '--T5_Initial_AlleleFreqs', the Mode " << Mode << " has not been implemented yet. Only Modes 'A', 'AllZeros', 'AllOnes' and 'Shift' are accepted for the moment.";
            abort();
        }
    }

    // Security check
    if (!this->T5_Initial_AlleleFreqs_AllZeros && !this->T5_Initial_AlleleFreqs_AllOnes)
    {
        assert(this->T5_Initial_AlleleFreqs.size() == GP->__PatchNumber[0]);
        for (int patch_index = 0 ; patch_index < GP->__PatchNumber[0] ; ++patch_index)
        {
            assert(this->T5_Initial_AlleleFreqs[patch_index].size() == this->T5_nbBits);
        }
    }
}

void SpeciesSpecificParameters::readSelectionOn(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option '--selectionOn', the std::string that is read is: " << input.print() << std::endl;
#endif    

    std::string s = input.GetNextElementString();

    if (s == "fertility")
    {
        selectionOn = 0;
    } else if (s == "viability")
    {
        selectionOn = 1;
    } else if (s == "both" || s == "fertilityAndViability")
    {
        selectionOn = 2;
    } else
    {
        std::cout << "For option '--selectionOn', expected either 'fertility' (which should be the default), 'viability' or 'fertilityAndViability' (aka. 'both'). Instead it received '" << s << "'.\n";
        abort();
    }
}

void SpeciesSpecificParameters::readGrowthK(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option '--GrowthK', the std::string that is read is: " << input.print() << std::endl;
#endif

    this->__growthK.resize(GP->__GenerationChange.size());
    for ( int generation_index = 0; generation_index < GP->__GenerationChange.size() ; generation_index++)
    {
        (void) input.GetNextGenerationMarker(generation_index);
        int currentPatchNumber = GP->__PatchNumber[generation_index];
        
        this->__growthK[generation_index].resize(currentPatchNumber);

        std::string Mode = input.GetNextElementString();

        if (Mode == "A")
        {
            for (int patch_index = 0 ; patch_index < currentPatchNumber ; patch_index++)
            {
                double x;
                std::string s = input.PeakNextElementString();
                if (s == "def" || s == "Def")
                {
                    input.skipElement();
                    assert(this->__patchCapacity[generation_index][patch_index] > 0);
                    x = -1.0;
                } else
                {
                    x = input.GetNextElementDouble();
                }
                
                this->__growthK[generation_index][patch_index] = x;
            }
        } else if (Mode == "unif")
        {
            std::string s = input.PeakNextElementString();
            double x;
            if (s == "def" || s == "Def")
            {
                x = -1;
                input.skipElement();
            } else
            {
                x = input.GetNextElementDouble();
            }

            for (int patch_index = 0 ; patch_index < currentPatchNumber ; patch_index++)
            {
                this->__growthK[generation_index][patch_index] = x;
            }
        } else
        {
            std::cout << "In --growthK, received " << Mode << " for mode. Only modes 'A' and 'unif' are recognized. Sorry\n";
            abort(); 
        }
        
            
    }
    this->growthK = this->__growthK[0];
}

void SpeciesSpecificParameters::readHabitats(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option '--H (--Habitats)', the std::string that is read is: " << input.print() << std::endl;
#endif

    for ( int generation_index = 0; generation_index < GP->__GenerationChange.size() ; generation_index++)
    {
        int Generation = input.GetNextGenerationMarker(generation_index);
        //std::cout << "Generation = " << Generation << "\n";
        std::string Mode = input.GetNextElementString();
        //std::cout << "Mode = " << Mode << "\n";

        std::vector<int> ForASingleGeneration;

        if (Mode.compare("unif")==0)
        {
            int UniqueHabitat = input.GetNextElementInt();
            assert(GP->__PatchNumber.size() > generation_index);
            for (int patch_index = 0 ; patch_index < GP->__PatchNumber[generation_index] ; ++patch_index)
            {
                ForASingleGeneration.push_back(UniqueHabitat);
            }
            this->__MaxHabitat.push_back(UniqueHabitat);
        } else if (Mode.compare("A")==0)
        {
            assert(GP->__PatchNumber[generation_index] > 0);
            int max = 0;
            assert(GP->__PatchNumber.size() > generation_index);
            for (int patch_index = 0 ; patch_index < GP->__PatchNumber[generation_index] ; ++patch_index)
            {
                int h = input.GetNextElementInt();
                if (max < h) {max = h;}                
                ForASingleGeneration.push_back(h);
            }
            this->__MaxHabitat.push_back(max);
        } else
        {
            std::cout << "In '--H (--Habitats)', at generation " << Generation << ", Mode received is " << Mode << ". Sorry, only Modes 'unif' and 'A' are implemented for the moment.\n";
            abort();
        }
        assert(ForASingleGeneration.size() == GP->__PatchNumber[generation_index]);
        
        // Assign to __Habitats
        this->__Habitats.push_back(ForASingleGeneration);
    }
    assert(this->__Habitats.size() == GP->__GenerationChange.size());
    assert(this->__Habitats.size() == this->__MaxHabitat.size());
    

    // flatten, sort and remove duplicates from the matrix to 1) get the maximum number and 2) make sure no habitat index is missing.
    std::vector<int> flat;
    for (int generation_index = 0 ; generation_index < GP->__GenerationChange.size() ; generation_index++)
    {
        auto& ForASingleGeneration = this->__Habitats[generation_index];
        assert(ForASingleGeneration.size() == GP->__PatchNumber[generation_index]);
        for (int patch_index = 0 ; patch_index < GP->__PatchNumber[generation_index] ; patch_index++)
        {
            flat.push_back(ForASingleGeneration[patch_index]);
        }
    }
    std::sort(flat.begin(), flat.end());
    flat.erase(std::unique(flat.begin(), flat.end()), flat.end());
    assert(this->__MaxHabitat.size() > 0);
    assert(this->__Habitats.size() > 0);
    this->MaxHabitat = this->__MaxHabitat[0];
    this->Habitats = this->__Habitats[0];

    if (flat[0] != 0)
    {
        std::cout << "In option '--H (--Habitats)', the lowest element received is " << flat[0] << ". The lowest element has to be 0 (zero based counting). Please also note that one cannot input an habitat index of n if all [0,n-1] have have not been included.\n";
        abort();
    }

    // 2) make sure no habitat index is missing.
    int previous = flat[0] - 1;
    for (int current : flat)
    {
        if (previous + 1 != current)
        {
            std::cout << "In option '--H (--Habitats)', you included the index " << current << " while the previous index has not been included. The last index included was " << previous << ".  Please also note that you used. Please note --Habitats use zero based counting.\n";
            abort();
        }
        previous = current;
    }

    
}

void SpeciesSpecificParameters::readpatchCapacity(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option '--N (--patchCapacity)', the std::string that is read is: " << input.print() << std::endl;
#endif

    this->TotalpatchCapacity = 0;
    maxEverpatchCapacity.resize(GP->maxEverPatchNumber,-1);

    for (int generation_index = 0 ; generation_index < GP->__GenerationChange.size() ; generation_index++)
    {
        (void) input.GetNextGenerationMarker(generation_index);
        int CurrentPatchNumber = GP->__PatchNumber[generation_index];
        assert(CurrentPatchNumber <= GP->maxEverPatchNumber);

        std::string Mode = input.GetNextElementString();
        std::vector<int> patchCapacity_line; // Will contain a single entry (for one interval of time) that will be pushed to '__patchCapacity'
        if (Mode.compare("A")==0)
        {
            for (int patch_index = 0 ; patch_index < CurrentPatchNumber; patch_index++)
            {
                int PC = input.GetNextElementInt();
                if (PC < 0)
                {
                    std::cout << "For option 'patchCapacity', received a patch capacity that is negative for the patch index "<<patch_index<<". Value received is "<<PC<<"."<< std::endl;
                    abort();
                }
                if (generation_index == 0) this->TotalpatchCapacity += PC;
                patchCapacity_line.push_back(PC);
            }
        } else if (Mode.compare("unif")==0)
        {
            int uniquepatchCapacity = input.GetNextElementInt();
            if (uniquepatchCapacity < 0)
            {
                std::cout << "For option 'patchCapacity', received a unique patch capacity that is negative. Value received is "<<uniquepatchCapacity<<"."<< std::endl;
                abort();
            }
            for (int patch_index=0 ; patch_index < CurrentPatchNumber ; ++patch_index)
            {
                if (generation_index == 0) this->TotalpatchCapacity += uniquepatchCapacity; 
                patchCapacity_line.push_back(uniquepatchCapacity);
            }
            
        } else
        {
            std::cout << "Sorry, for option 'patchCapacity', only Modes 'A' and 'unif' are implemented for the moment. Mode received is " << Mode << std::endl;
            abort();
        }

        for (int patch_index=0 ; patch_index < CurrentPatchNumber ; ++patch_index)
        {
            maxEverpatchCapacity[patch_index] = std::max(maxEverpatchCapacity[patch_index], patchCapacity_line[patch_index]);
        }

        this->__patchCapacity.push_back(patchCapacity_line);
            

        if (generation_index == 0) 
        {
            if (this->TotalpatchCapacity <= 0)
            {
                std::cout << "The total patch capacity over all patches for species " << this->speciesName << " is " << this->TotalpatchCapacity << ". Nothing can be simulated.\n";
                abort();
            }
        }
    }
    
    assert(this->__patchCapacity.size() == GP->__GenerationChange.size());
    this->patchCapacity = this->__patchCapacity[0];
    

    for (int patch_index=0 ; patch_index < GP->maxEverPatchNumber ; ++patch_index)
    {
        if (maxEverpatchCapacity[patch_index] == 0)
        {
            std::cout << "patch_index " << patch_index << " never gets a carrying capacity greater than " << 0 << ". This is a weird unseless patch!\n";
            abort();
        }
    }
}

void SpeciesSpecificParameters::readnbSubGenerations(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option '--nbSubGens (--nbSubGenerations)', the std::string that is read is: " << input.print() << std::endl;
#endif

    this->nbSubGenerationsPerGeneration = input.GetNextElementInt();
    if (this->nbSubGenerationsPerGeneration < 1)
    {
        std::cout << "In option '--nbSubGens (--nbSubGenerations)', for species " << this->speciesName << " received " << this->nbSubGenerationsPerGeneration << ". The number of subGeneration must be equal or greater than 1.\n";
        abort();
    }
    
    
}


void SpeciesSpecificParameters::readT1_EpistaticFitnessEffects(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option '--T1_epistasis (--T1_EpistaticFitnessEffects)', the std::string that is read is: " << input.print() << std::endl;
#endif 

    if (input.PeakNextElementString() == "NA")
    {
        if (this->T1_nbBits != 0)
        {
            std::cout << "In option '--T1_epistasis (--T1_EpistaticFitnessEffects)', for species "<<this->speciesName<<", received 'NA' but there are T1 loci as indiciated in --L (--Loci) (this->T1_nbBits = "<<this->T1_nbBits<<").\n";
            abort();   
        }
        input.skipElement();
    } else
    {
        for (int habitat = 0; habitat <= this->MaxHabitat ; habitat++)
        {
            input.GetNextHabitatMarker(habitat);

            // Get number of groups of loci for this habitat
            int nbGroupsOfLoci = input.GetNextElementInt();
            if (nbGroupsOfLoci < 0)
            {
                std::cout << "In option '--T1_epistasis (--T1_EpistaticFitnessEffects)', for species "<<this->speciesName<<", receive a value lower than 0 (received "<<nbGroupsOfLoci<<")\n";
                abort();
            }

            // Get info for this habitat
            std::vector<std::vector<T1_locusDescription>> LociIndices_ForASingleHabitat;
            std::vector<std::vector<double>>              FitnessEffects_ForASingleHabitat;

            for (int groupOfLociIndex = 0 ; groupOfLociIndex < nbGroupsOfLoci ; groupOfLociIndex++)
            {
                //std::cout << "groupOfLociIndex = " << groupOfLociIndex << "\n";
                std::vector<T1_locusDescription> LociIndices_ForASingleGroupOfLoci;
                std::vector<double>              FitnessEffects_ForASingleGroupOfLoci;

                // Get 'nbLociUnderEpistasis'
                int nbLociUnderEpistasis = input.GetNextElementInt();
                if (nbLociUnderEpistasis >= 20)
                {
                    std::cout << "In option '--T1_epistasis (--T1_EpistaticFitnessEffects)', the nbLociUnderEpistasis indicated for habitat " << habitat << " and groupOfLociIndex "<< groupOfLociIndex <<" is " << nbLociUnderEpistasis << ". Such a 'nbLociUnderEpistasis' will require you to indicate" << pow(3,nbLociUnderEpistasis) << " fitness values! Is this really what you want to do? If yes comment out this security check in file SpeciesSpecificParameters.cpp! Please be well aware that all the fitness values that you will input will be stored in the RAM!" << std::endl;
                    abort();
                }
                if (nbLociUnderEpistasis < 2)
                {
                    std::cout << "In option '--T1_epistasis (--T1_EpistaticFitnessEffects)', the nbLociUnderEpistasis indicated for habitat " << habitat << " and groupOfLociIndex "<< groupOfLociIndex <<" is " << nbLociUnderEpistasis << ". To have an epistatic interaction, you need at least two loci in a group of loci." << std::endl;
                    abort();
                }

                for (int i = 0 ; i < nbLociUnderEpistasis ; i++)
                {
                    int LocusPosition = input.GetNextElementInt();
                    if ( LocusPosition < 0 || LocusPosition >= this->T1_nbBits )
                    {
                      std::cout << "In option '--T1_epistasis (--T1_EpistaticFitnessEffects)', the " << i << "th locus position indicated for habitat " << habitat << " and groupOfLociIndex "<< groupOfLociIndex <<" is either lower than zero or greater or equal to the number of T1 sites (" << this->T1_nbBits <<"). As a reminder the first locus is the zeroth locus (zero based counting). The locus position received is " << LocusPosition << ".\n";
                        abort();  
                    }
                    //std::cout << "LocusPosition = " << LocusPosition << "\n";
                    T1_locusDescription T1_locus(LocusPosition / EIGHT, LocusPosition % EIGHT, LocusPosition);
                    LociIndices_ForASingleGroupOfLoci.push_back(T1_locus);
                }

                /*
                SSP->T1_Epistasis_LociIndices[Habitat][0].
                locus 1:    00 00 00 00 00 00 00 00 00 01 01 01 01 01 01 01 01 01 11 11 11 11 11 11 11 11 11
                locus 2:    00 00 00 01 01 01 11 11 11 00 00 00 01 01 01 11 11 11 00 00 00 01 01 01 11 11 11
                locus 3:    00 01 11 00 01 11 00 01 11 00 01 11 00 01 11 00 01 11 00 01 11 00 01 11 00 01 11
                Index:      0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 
                */

                int TotalNbFitnessesToIndicate = pow(3,nbLociUnderEpistasis);
                for (int i = 0 ; i < TotalNbFitnessesToIndicate; i++)
                {
                    double fit = input.GetNextElementDouble();
                    if ( fit < 0.0 || fit > 1.0 )
                    {
                      std::cout << "In option '--T1_epistasis (--T1_EpistaticFitnessEffects)', the " << i << "th fitness value indicated after habitat " << habitat << "is either lower than zero or greater than one. The fitness value indicated is  " << fit << ".\n";
                        abort();  
                    }
                    FitnessEffects_ForASingleGroupOfLoci.push_back(fit);               
                }

                LociIndices_ForASingleHabitat.push_back(LociIndices_ForASingleGroupOfLoci);
                FitnessEffects_ForASingleHabitat.push_back(FitnessEffects_ForASingleGroupOfLoci);   

            } // end of for groupLociIndex
            this->T1_Epistasis_LociIndices.push_back(LociIndices_ForASingleHabitat);
            this->T1_Epistasis_FitnessEffects.push_back(FitnessEffects_ForASingleHabitat); 
        } // end of for habitat
        
        assert(this->T1_Epistasis_LociIndices.size() == this->MaxHabitat + 1);
        assert(this->T1_Epistasis_FitnessEffects.size() == this->MaxHabitat + 1);

        //get fitness value with this->FitnessEffects_ForASingleHabitat[habitat][groupOfLoci][fitnessValueIndex]
    } // end of if "NA"

    
}

void SpeciesSpecificParameters::readT5_FitnessEffects(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option '--T5_FitnessEffects', the std::string that is read is: " << input.print() << std::endl;
#endif
    assert(SSP != nullptr);
    if (input.PeakNextElementString() == "NA")
    {
        if (this->T5_nbBits != 0)
        {
            std::cout << "In option '--T5_FitnessEffects', for species "<<this->speciesName<<", received 'NA' but there are T5 loci as indiciated in --L (--Loci) (this->T5_nbBits = "<<this->T5_nbBits<<").\n";
            abort();   
        }
        input.skipElement();
    } else
    {
        for (int habitat = 0; habitat <= this->MaxHabitat ; habitat++)
        {
            (void) input.GetNextHabitatMarker(habitat);
            std::string Mode = input.GetNextElementString();
            std::string firstMode;
            if (habitat == 0)
            {
                firstMode = Mode;
            } else
            {
                if  (
                        (
                            (firstMode.compare("A")==0 || firstMode.compare("unif")==0)
                            &&
                            (Mode.compare("MultiplicityA")==0 || Mode.compare("MultiplicityAGamma")==0 || Mode.compare("MultiplicityUnif")==0)
                        )
                        ||
                        (
                            (Mode.compare("A")==0 || Mode.compare("unif")==0)
                            &&
                            (firstMode.compare("MultiplicityA")==0 || Mode.compare("MultiplicityAGamma")==0 || firstMode.compare("MultiplicityUnif")==0  )
                        )
                    )
                {
                    std::cout << "The Mode received for the zeroth habitat is '" << firstMode << "'. The Mode received for the " << habitat << "th habitat is '" << Mode << "'. In the current version, if one patch make the assumption of 'multiplicity', then all patches must make this assumption.\n";
                    abort();
                }
            }
       
            std::vector<double> ForASingleHabitat;
                   
            if (Mode.compare("A")==0)
            {
                this->T5_isMultiplicitySelection = false;
                for (int entry_index = 0 ; entry_index < (2 * this->T5_nbBits) ; entry_index++)    
                {
                    double fit = input.GetNextElementDouble();
                    if (fit < 0 || fit > 1 )
                    {
                        std::cout << "In option '--T5_FitnessEffects', habitat " << habitat << ", Mode 'A', one 'fit' value is either greater than 1 or lower than 0 ( is " << fit << " ).\n";
                        abort();
                    }
                    ForASingleHabitat.push_back(fit);
                }
                if (ForASingleHabitat.size() != 2 * this->T5_nbBits)
                {
                    std::cout << "In option '--T5_FitnessEffects', habitat "  << habitat << ", Mode 'A' " << ForASingleHabitat.size() << " have been received while " << 2 * (this->T5_nbBits) << " were expected\n";
                    abort();
                }
            } else if (Mode.compare("MultiplicityA")==0)
            {
                this->T5_isMultiplicitySelection = true;
                for (int entry_index = 0 ; entry_index < this->T5_nbBits ; )    
                {
                    double fit01 = input.GetNextElementDouble();
                    if (fit01 < 0 || fit01 > 1 )
                    {
                        std::cout << "In option '--T5_FitnessEffects', habitat " << habitat << ", Mode 'MultiplicityA', one 'fit01' value is either greater than 1 or lower than 0 ( is " << fit01 << " ).\n";
                        abort();
                    }
                    ForASingleHabitat.push_back(fit01);
                    entry_index++;
                }
                if (ForASingleHabitat.size() != this->T5_nbBits)
                {
                    std::cout << "In option '--T5_FitnessEffects', habitat " << habitat << ", Mode 'MultiplicityA' " << ForASingleHabitat.size() << " have been received while " << this->T5_nbBits << "  were expected\n";
                    abort();
                }
            } else if (Mode.compare("unif")==0)
            {
                this->T5_isMultiplicitySelection = false;
                
                double fitHet = input.GetNextElementDouble();
                double fitHom = input.GetNextElementDouble();
                
                if (fitHet<0 || fitHom<0 || fitHet>1 || fitHom>1)
                {
                    std::cout << "In option '--T5_FitnessEffects', habitat " << habitat << ", Mode 'unif', value for fitness are negative or great than 1 (they are " << fitHet << " " << fitHom << " )\n";
                    abort();
                }
                
                for (int locus = 0 ; locus < this->T5_nbBits ; ++locus)
                {
                    ForASingleHabitat.push_back(fitHet);
                    ForASingleHabitat.push_back(fitHom);
                }
                if (ForASingleHabitat.size() != 2 * this->T5_nbBits)
                {
                    std::cout << "In option '--T5_FitnessEffects', habitat " << habitat << ", Mode 'unif' " << ForASingleHabitat.size() << " have been received while " << 2 * (this->T5_nbBits) << " were expected\n";
                    abort();
                }
            } else if (Mode.compare("MultiplicityUnif")==0)
            {
                this->T5_isMultiplicitySelection = true;
                
                double fit01 = input.GetNextElementDouble();
                if (fit01<0 || fit01>1)
                {
                    std::cout << "In option '--T5_FitnessEffects' Mode 'MultiplicityUnif', habitat " << habitat << ", value for 'fit01' is negative or greater than zero (is " << fit01 << ")\n";
                    abort();
                }
                for (int locus = 0 ; locus < this->T5_nbBits ; locus++)
                {
                    ForASingleHabitat.push_back(fit01);
                }
                if (ForASingleHabitat.size() != this->T5_nbBits)
                {
                    std::cout << "In option '--T5_FitnessEffects', habitat " << habitat << ", Mode 'MultiplicityUnif' " << ForASingleHabitat.size() << " have been received while " << this->T5_nbBits << "  were expected (likely an internal error)\n";
                    abort();
                }
            } else if (Mode.compare("MultiplicityAGamma")==0)
            {
                double alpha = input.GetNextElementDouble();
                double beta = input.GetNextElementDouble();
                if (alpha <= 0 || beta <= 0)
                {
                    std::cout << "Either alpha ("<<alpha<<") or beta ("<<beta<<") have non-positive values.\n";
                    abort();
                }
                std::gamma_distribution<double> dist(alpha, beta);
                for (int locus = 0 ; locus < this->T5_nbBits ; locus++)
                {
                    double fit = 1 - dist(GP->mt);
                    if (fit < 0.0) fit = 0.0;
                    if (fit > 1.0) fit = 1.0;
                    ForASingleHabitat.push_back(fit);
                }
            } else
            {
                std::cout << "Sorry, for option '--T5_fit (--T5_fitnessEffects)', only Modes 'A', 'unif', MultiplicityA', 'MultiplicityUnif' and 'MultiplicityAGamma' are implemented for the moment (received Mode " << Mode << ")" << std::endl;
                abort();
            } // end of ifelse Mode
            this->T5_FitnessEffects.push_back(ForASingleHabitat);
        } // end of for habitat
    
        assert(this->T5_FitnessEffects.size() == this->MaxHabitat + 1);
    }
        
}

void SpeciesSpecificParameters::readT1_FitnessEffects(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option '--T1_FitnessEffects', the std::string that is read is: " << input.print() << std::endl;
#endif
    assert(SSP != nullptr);
    if (input.PeakNextElementString() == "NA")
    {
        if (this->T1_nbBits != 0)
        {
            std::cout << "In option '--T1_FitnessEffects', for species "<<this->speciesName<<", received 'NA' but there are T1 loci as indiciated in --L (--Loci) (this->T1_nbBits = "<<this->T1_nbBits<<").\n";
            abort();   
        }
        input.skipElement();
    } else
    {
        for (int habitat = 0; habitat <= this->MaxHabitat ; habitat++)
        {
            (void) input.GetNextHabitatMarker(habitat);
            std::string Mode = input.GetNextElementString();
            std::string firstMode;
            if (habitat == 0)
            {
                firstMode = Mode;
            } else
            {
                if  (
                        (
                            (firstMode.compare("A")==0 || firstMode.compare("unif")==0)
                            &&
                            (Mode.compare("MultiplicityA")==0 || Mode.compare("MultiplicityAGamma")==0 || Mode.compare("MultiplicityUnif")==0)
                        )
                        ||
                        (
                            (Mode.compare("A")==0 || Mode.compare("unif")==0)
                            &&
                            (firstMode.compare("MultiplicityA")==0 || Mode.compare("MultiplicityAGamma")==0 || firstMode.compare("MultiplicityUnif")==0  )
                        )
                    )
                {
                    std::cout << "The Mode received for the zeroth habitat is '" << firstMode << "'. The Mode received for the " << habitat << "th habitat is '" << Mode << "'. In the current version, if one patch make the assumption of 'multiplicity', then all patches must make this assumption.\n";
                    abort();
                }
            }
       
            std::vector<double> ForASingleHabitat;
            if (Mode.compare("DomA")==0 || Mode.compare("domA")==0)
            {
            
                this->T1_isMultiplicitySelection = false;
                std::string cst_or_fun = input.GetNextElementString();
                if (cst_or_fun.compare("cst") != 0 && cst_or_fun.compare("fun") != 0)
                {
                    std::cout << "In option '--T1_FitnessEffects', habitat " << habitat << ", Mode 'domA', you had to indicate 'cst' or 'fun' for how dominance coefficient is defined but instead, SimBit received " << cst_or_fun << ".\n";
                    abort();
                }
            

                double dominance = input.GetNextElementDouble(); // either cst or mean dominance
                if (dominance < 0.0 || dominance > 0.5)
                {
                    std::cout << "In option '--T1_FitnessEffects', habitat " << habitat << ", Mode 'domA', the dominance value is either smaller than 0.0 or larger than 1.0. dominance recevied = " << dominance << std::endl;
                    abort();
                }
                double mean_s = 0; // used only if cst_or_fun is 'fun'
            
                for (int entry_index = 0 ; entry_index < this->T1_nbBits ; )    
                {
                
                    if (input.isNextRKeyword())
                    {
                        double fit11 = input.GetNextElementDouble();
                        if (fit11 < 0 || fit11 > 1 )
                        {
                            std::cout << "In option '--T1_FitnessEffects', habitat " << habitat << ", Mode 'domA', one 'fit' value (directly after keyword 'R') is either greater than 1 or lower than 0 ( is " << fit11 << " ).\n";
                            abort();
                        }
                        double fit01;
                        if (cst_or_fun.compare("cst") == 0)
                        {
                            fit01 = 1 - ((1-fit11) * dominance);
                        } else if (cst_or_fun.compare("fun") == 0)
                        {
                            fit01 = -1.0; // this is temporary because we need to read all s, before being able to compute the locus specific dominance
                            mean_s += 1 - fit11;
                        } else
                        {
                            std::cout << "Internal error in Parameters::readT1_FitnessEffects.\n";
                            abort();
                        }
                            
                        int nbRepeats = input.GetNextElementInt();
                        for (int i = 0 ; i < nbRepeats ; i++)
                        {
                            ForASingleHabitat.push_back(1.0);
                            ForASingleHabitat.push_back(fit01);
                            ForASingleHabitat.push_back(fit11);
                            entry_index++;
                            if (entry_index > this->T1_nbBits)
                            {
                                std::cout << "In 'T1_FitnessEffects', mode 'DomA', when reading after keyword 'R', you asked to repeat the value " << fit11 << " for " << nbRepeats << " times. This (on top with previous entries) results of a number of entry greater than " << this->T1_nbBits << " as requested given the mode here ('CstDomA') and the arguments given to the option '--Loci'\n";
                                abort();
                            }
                        }
                    } else
                    {
                    
                        double fit11 = input.GetNextElementDouble();
                        if (fit11 < 0 || fit11 > 1 )
                        {
                            std::cout << "In option '--T1_FitnessEffects', habitat " << habitat << ", Mode 'A', one 'fit' value is either greater than 1 or lower than 0 ( is " << fit11 << " ).\n";
                            abort();
                        }
                    
                        double fit01;
                        if (cst_or_fun.compare("cst") == 0)
                        {
                            fit01 = 1 - ((1-fit11) * dominance);
                        } else if (cst_or_fun.compare("fun") == 0)
                        {
                            fit01 = -1.0; // this is temporary because we need to read all s, before being able to compute the locus specific dominance
                            mean_s += 1 - fit11;
                        } else
                        {
                            std::cout << "Internal error in Parameters::readT1_FitnessEffects.\n";
                            abort();
                        }

                        ForASingleHabitat.push_back(1.0);
                        ForASingleHabitat.push_back(fit01);
                        ForASingleHabitat.push_back(fit11);
                        entry_index++;
                    
                    }
                
                }  
                if (ForASingleHabitat.size() != THREE * this->T1_nbBits)
                {
                    std::cout << "In option '--T1_FitnessEffects', habitat "  << habitat << ", Mode 'DomA' " << ForASingleHabitat.size() / 3 << " have been received while " << this->T1_nbBits << " (= 8 * ( " << this->T1_nbChars << "))  were expected\n";
                    abort();
                }

                // Now if used the function, compute locus-specific dominance and insert them in place of the current -1.0
                if (cst_or_fun.compare("fun") == 0)
                {
                    mean_s /= this->T1_nbBits;
                    assert(mean_s >= 0.0 && mean_s <= 1.0);
                    if (mean_s == 0.0) {std::cout << "All T1 are neutral!\n";}
                    double k = -1.0;
                    if (mean_s > 0.0)
                    {
                        k = log(2 * dominance) / mean_s; // See Nemo manual for definition
                    }
                    int securityCount = 0;
                    for (int i = 2 ; i < ForASingleHabitat.size() ; i += 3)
                    {
                        double fit11 = ForASingleHabitat[i];
                        assert(fit11 >= 0.0 && fit11 <= 1.0);
                        double s11 = 1 - fit11;
                        assert(s11 >= 0.0 && s11 <= 1.0);
                        double dom = -1.0;
                        double fit01 = 1.0;
                        if (mean_s > 0.0)
                        {
                            dom = exp(k * s11) / 2;   
                            assert(dom >= 0.0 && dom <= 1.0);
                            fit01 = 1 - (s11 * dom);
                        }
                        assert(fit01 >= 0.0 && fit01 <= 1.0);
                        assert(ForASingleHabitat[i-2] == 1.0);  // homozygote wildtype
                        assert(ForASingleHabitat[i-1] == -1.0); // heterozygote should currently be set to -1.0
                        ForASingleHabitat[i-1] = fit01;
                        securityCount++;
                    }
                    assert(securityCount == this->T1_nbBits);
                }
                    
            } else if (Mode.compare("A")==0)
            {
                this->T1_isMultiplicitySelection = false;
                for (int entry_index = 0 ; entry_index < THREE * this->T1_nbBits ; )    
                {
                    if (input.isNextRKeyword())
                    {
                        double fit00 = input.GetNextElementDouble();
                        double fit01 = input.GetNextElementDouble();
                        double fit11 = input.GetNextElementDouble();
                        int nbRepeats = input.GetNextElementInt();
                        for (int i = 0 ; i < nbRepeats ; i++)
                        {
                            ForASingleHabitat.push_back(fit00);
                            ForASingleHabitat.push_back(fit01);
                            ForASingleHabitat.push_back(fit11);
                            entry_index += 3;
                            if (entry_index > (THREE * this->T1_nbBits))
                            {
                                std::cout << "In 'T1_FitnessEffects', mode 'A', when reading after keyword 'R', you asked to repeat the values " << fit00 << " " << fit01 << " " << fit11 << " for " << nbRepeats << " times. This (on top with previous entries) results of a number of entry greater than " << THREE * this->T1_nbBits << " as requested given the mode here ('A') and the arguments given to the optinn '--Loci'\n";
                                abort();
                            }
                        }
                    } else
                    {
                        double fit = input.GetNextElementDouble();
                        if (fit < 0 || fit > 1 )
                        {
                            std::cout << "In option '--T1_FitnessEffects', habitat " << habitat << ", Mode 'A', one 'fit' value is either greater than 1 or lower than 0 ( is " << fit << " ).\n";
                            abort();
                        }
                        ForASingleHabitat.push_back(fit);
                        entry_index++;
                    }
                }
                if (ForASingleHabitat.size() != THREE * this->T1_nbBits)
                {
                    std::cout << "In option '--T1_FitnessEffects', habitat "  << habitat << ", Mode 'A' " << ForASingleHabitat.size() << " have been received while " << THREE * (this->T1_nbBits) << " (= 3 * 8 * ( " << this->T1_nbChars << "))  were expected\n";
                    abort();
                }
            } else if (Mode.compare("MultiplicityA")==0)
            {
                this->T1_isMultiplicitySelection = true;
                for (int entry_index = 0 ; entry_index < this->T1_nbBits ; )    
                {
                    if (input.isNextRKeyword())
                    {
                        double fit01 = input.GetNextElementDouble();
                        int nbRepeats = input.GetNextElementInt();
                        if (fit01 < 0 || fit01 > 1 )
                        {
                            std::cout << "In option '--T1_FitnessEffects', habitat " << habitat << ", Mode 'MultiplicityA', one 'fit01' value (coming after keyword 'R') is either greater than 1 or lower than 0 ( is " << fit01 << " ).\n";
                            abort();
                        }
                        for (int i = 0 ; i < nbRepeats; i++)
                        {
                            ForASingleHabitat.push_back(fit01);
                            entry_index++;
                            if (entry_index > this->T1_nbBits)
                            {
                                std::cout << "In 'T1_FitnessEffects', mode 'MultiplicityA', when reading after keyword 'R', you asked to repeat the value " << fit01 << " for " << nbRepeats << " times. This (on top with previous entries) results of a number of entry greater than " << this->T1_nbBits << " as requested given the mode here (MultiplicityA) and the arguments given to the option '--Loci'\n";
                                abort();
                            }
                        }
                    } else
                    {
                        double fit01 = input.GetNextElementDouble();
                        if (fit01 < 0 || fit01 > 1 )
                        {
                            std::cout << "In option '--T1_FitnessEffects', habitat " << habitat << ", Mode 'MultiplicityA', one 'fit01' value is either greater than 1 or lower than 0 ( is " << fit01 << " ).\n";
                            abort();
                        }
                        ForASingleHabitat.push_back(fit01);
                        entry_index++;
                    }
                }
                if (ForASingleHabitat.size() != this->T1_nbBits)
                {
                    std::cout << "In option '--T1_FitnessEffects', habitat " << habitat << ", Mode 'MultiplicityA' " << ForASingleHabitat.size() << " have been received while " << this->T1_nbBits << " (= 8 * ( " << this->T1_nbChars << "))  were expected\n";
                    abort();
                }
            } else if (Mode.compare("unif")==0)
            {
                this->T1_isMultiplicitySelection = false;
                
                double fit00 = input.GetNextElementDouble();
                double fit01 = input.GetNextElementDouble();
                double fit11 = input.GetNextElementDouble();
                if (fit00<0 || fit01<0 || fit11<0 || fit00>1 || fit01>1 || fit11>1)
                {
                    std::cout << "In option '--T1_FitnessEffects', habitat " << habitat << ", Mode 'unif', value for fitness are negative or great than 1 (they are " << fit00 << " " << fit01 << " " << fit11 << " )\n";
                    abort();
                }
                
                for (int locus = 0 ; locus < this->T1_nbBits ; ++locus)
                {
                    ForASingleHabitat.push_back(fit00);
                    ForASingleHabitat.push_back(fit01);
                    ForASingleHabitat.push_back(fit11);
                }
                if (ForASingleHabitat.size() != THREE * this->T1_nbBits)
                {
                    std::cout << "In option '--T1_FitnessEffects', habitat " << habitat << ", Mode 'unif' " << ForASingleHabitat.size() << " have been received while " << THREE * (this->T1_nbBits) << " (= 3 * 8 * ( " << this->T1_nbChars << "))  were expected\n";
                    abort();
                }
            } else if (Mode.compare("MultiplicityUnif")==0)
            {
                this->T1_isMultiplicitySelection = true;
                
                double fit01 = input.GetNextElementDouble();
                if (fit01<0 || fit01>1)
                {
                    std::cout << "In option '--T1_FitnessEffects' Mode 'MultiplicityUnif', habitat " << habitat << ", value for 'fit01' is negative or greater than zero (is " << fit01 << ")\n";
                    abort();
                }
                for (int locus = 0 ; locus < this->T1_nbBits ; locus++)
                {
                    ForASingleHabitat.push_back(fit01);
                }
                if (ForASingleHabitat.size() != this->T1_nbBits)
                {
                    std::cout << "In option '--T1_FitnessEffects', habitat " << habitat << ", Mode 'MultiplicityUnif' " << ForASingleHabitat.size() << " have been received while " << this->T1_nbBits << " (= 8 * ( " << this->T1_nbChars << "))  were expected\n";
                    abort();
                }
            } else if (Mode.compare("MultiplicityAGamma")==0)
            {
                double alpha = input.GetNextElementDouble();
                double beta = input.GetNextElementDouble();
                if (alpha <= 0 || beta <= 0)
                {
                    std::cout << "Either alpha ("<<alpha<<") or beta ("<<beta<<") have non-positive values.\n";
                    abort();
                }
                std::gamma_distribution<double> dist(alpha, beta);
                for (int locus = 0 ; locus < this->T1_nbBits ; locus++)
                {
                    double fit = 1 - dist(GP->mt);
                    if (fit < 0.0) fit = 0.0;
                    if (fit > 1.0) fit = 1.0;
                    ForASingleHabitat.push_back(fit);
                }
            } else
            {
                std::cout << "Sorry, for option '--T1_fit (--T1_fitnessEffects)', only Modes 'A', 'unif', 'domA' (or 'DomA'), MultiplicityA', 'MultiplicityUnif' and 'MultiplicityAGamma' are implemented for the moment (received Mode " << Mode << ")" << std::endl;
                abort();
            } // end of ifelse Mode
            this->T1_FitnessEffects.push_back(ForASingleHabitat);
        } // end of for habitat
    
        assert(this->T1_FitnessEffects.size() == this->MaxHabitat + 1);
    }
        
    
}


void SpeciesSpecificParameters::readT2_FitnessEffects(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option '--T2_FitnessEffects', the std::string that is read is: " << input.print() << std::endl;
#endif

    if (input.PeakNextElementString() == "NA")
    {
        if (this->T2_nbChars != 0)
        {
            std::cout << "In option '--T2_FitnessEffects', for species "<<this->speciesName<<", received 'NA' but there are T2 loci as indiciated in --L (--Loci) (this->T2_nbChars = "<<this->T2_nbChars<<").\n";
            abort();   
        }
        input.skipElement();
    } else
    {
        for (int habitat = 0 ; habitat <= this->MaxHabitat ; habitat++)
        {
            (void) input.GetNextHabitatMarker(habitat);

            std::vector<double> ForASingleHabitat;
            
            std::string Mode = input.GetNextElementString();
            if (Mode.compare("unif") == 0)
            {
                double fit = input.GetNextElementDouble();
                if (fit < 0.0 || fit > 1.0)
                {
                    std::cout <<  "In option '--T2_FitnessEffects', habitat " << habitat << ", Mode 'unif', 'fit' value received is " << fit << " which is either lower than zero or greater than 1.\n";
                    abort();
                }
                for (int char_index = 0 ; char_index < this->T2_nbChars ; char_index++)
                {
                    ForASingleHabitat.push_back(fit);
                }
            } else if (Mode.compare("A") == 0)
            {
                for (int entry_index = 0 ; entry_index < this->T2_nbChars ; ++entry_index)
                {
                    double fit = input.GetNextElementDouble();
                    if (fit < 0.0 || fit > 1.0)
                    {
                        std::cout <<  "In option '--T2_FitnessEffects' Mode 'A', 'fit' value received is " << fit << " which is either lower than zero or greater than 1.\n";
                        abort();
                    }
                    ForASingleHabitat.push_back(fit);
                }
                if (ForASingleHabitat.size() != this->T2_nbChars)
                {
                    std::cout << "In option '--T2_FitnessEffects', habitat " << habitat << ", Mode 'A', " << ForASingleHabitat.size() << " values received while " << this->T2_nbChars << " were expected." << std::endl;
                    abort();
                }
            } else if (Mode.compare("Gamma") == 0)
            {
                double alpha = input.GetNextElementDouble();
                double beta = input.GetNextElementDouble();
                if (alpha <= 0 || beta <= 0)
                {
                    std::cout << "Either alpha ("<<alpha<<") or beta ("<<beta<<") have non-positive values.\n";
                    abort();
                }
                std::gamma_distribution<double> dist(alpha, beta);
                for (int locus = 0 ; locus < this->T2_nbChars ; locus++)
                {
                    double fit = 1 - dist(GP->mt);
                    if (fit < 0.0) fit = 0.0;
                    if (fit > 1.0) fit = 1.0;
                    ForASingleHabitat.push_back(fit);
                }
            } else
            {
                std::cout << "In option '--T2_FitnessEffects', habitat " << habitat << ", Mode received is " << Mode << ". Sorry only Modes 'A' and 'unif' are recognized for the moment!\n";
                abort();
            } // end of ifelse Mode
            assert(this->T2_FitnessEffects.size() == habitat);
            this->T2_FitnessEffects.push_back(ForASingleHabitat);
        } // end of Habitat loop

        assert(this->T2_FitnessEffects.size() == this->MaxHabitat + 1);
    }

    
}

void SpeciesSpecificParameters::readT5_MutationRate(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option '--T5_MutationRate', the std::string that is read is: " << input.print() << std::endl;
#endif

    if (input.PeakNextElementString() == "NA")
    {
        if (this->T5_nbBits != 0)
        {
            std::cout << "In option '--T5_MutationRate', for species "<<this->speciesName<<", received 'NA' but there are T5 loci as indiciated in --L (--Loci) (this->T5_nbBits = "<<this->T5_nbBits<<").\n";
            abort();   
        }
        input.skipElement();
    } else
    {
        std::string Mode = input.GetNextElementString();
        if (Mode.compare("A")==0)
        {
            double currentSum = 0;
            while (input.IsThereMoreToRead())
            {
                double inputValue = input.GetNextElementDouble();
                currentSum += inputValue;
                this->T5_MutationRate.push_back(currentSum);

                // security
                if (inputValue < 0)
                {
                    std::cout << "In 'T5_MutationRate' of Mode 'A' the value received (" << inputValue << ") is negative." << std::endl;
                    abort();
                }
                if (inputValue >= 0.1)
                {
                    std::cout << "In 'T5_MutationRate' of Mode 'A' the value received (" << inputValue << ") is greater than 0.1. Is it really what you want to do?" << std::endl;
                    abort();
                }          
            }

            // Did I get the info I expected?
            if (this->T5_MutationRate.size() != T5_nbBits)
            {
                std::cout << "In '--T5_MutationRate', expected " << this->T5_nbBits << " elements. But " << this->T5_MutationRate.size() << " elements were received.\n";
                abort();
            }

        } else if (Mode.compare("unif")==0)
        {
            this->T5_MutationRate.push_back(input.GetNextElementDouble() * this->T5_nbBits);
            if (this->T5_MutationRate.back() < 0)
            {
                std::cout << "T5_MutationRate of Mode 'unif' is lower than zero" << std::endl;
                abort();
            }
        } else
        {
            std::cout << "Sorry, for 'T5_MutationRate' only Mode 'unif' and 'A' are recognized so far. Mode received (" << Mode << ") is unknown!" << std::endl;
            abort();
        }
        
        assert(this->T5_MutationRate.size() == this->T5_nbBits || this->T5_MutationRate.size() == 1);
        this->T5_Total_Mutation_rate = this->T5_MutationRate.back();
    }    
}

void SpeciesSpecificParameters::readT1_MutationRate(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option '--T1_MutationRate', the std::string that is read is: " << input.print() << std::endl;
#endif

    if (input.PeakNextElementString() == "NA")
    {
        if (this->T1_nbBits != 0)
        {
            std::cout << "In option '--T1_MutationRate', for species "<<this->speciesName<<", received 'NA' but there are T1 loci as indiciated in --L (--Loci) (this->T1_nbBits = "<<this->T1_nbBits<<").\n";
            abort();   
        }
        input.skipElement();
    } else
    {
        std::string Mode = input.GetNextElementString();
        if (Mode.compare("A")==0)
        {
            double currentSum = 0;
            int T1_char_index = 0;
            int bit_index = 0;
            while (input.IsThereMoreToRead())
            {
                if (input.isNextRKeyword())
                {
                    double inputValue = input.GetNextElementDouble();
                    int NbRepeats  = input.GetNextElementInt();

                    for (int repeat = 0 ; repeat < NbRepeats ; ++repeat)
                    {
                        currentSum += inputValue;
                        this->T1_MutationRate.push_back(currentSum); 

                        // Adjust the bit_index and T1_char_index
                        if (bit_index == 7)
                        {
                            bit_index = 0;
                            ++T1_char_index;
                        } else
                        {
                            ++bit_index;
                        }   
                    }

                    // security
                    if (inputValue < 0)
                    {
                        std::cout << "In 'T1_MutationRate' of Mode 'A' the value for the " << T1_char_index << "th T1_char_index is negative (is " << inputValue << ")" << std::endl;
                        abort();
                    }
                    if (inputValue >= 0.1)
                    {
                        std::cout << "In 'T1_MutationRate' of Mode 'A' the value for the " << T1_char_index << "th T1_char_index is greater or equal to 0.1 (is " << inputValue << "). Is it really what you want to do?" << std::endl;
                        abort();
                    }

                } else
                {
                    double inputValue = input.GetNextElementDouble();
                    currentSum += inputValue;
                    this->T1_MutationRate.push_back(currentSum);

                    // security
                    if (inputValue < 0)
                    {
                        std::cout << "In 'T1_MutationRate' of Mode 'A' the value for the " << T1_char_index << "th T1_char_index is negative (is " << inputValue << ")" << std::endl;
                        abort();
                    }
                    if (inputValue >= 0.1)
                    {
                        std::cout << "In 'T1_MutationRate' of Mode 'A' the value for the " << T1_char_index << "th T1_char_index is greater than 0.1 (is " << inputValue << "). Is it really what you want to do?" << std::endl;
                        abort();
                    }

                    // Adjust the bit_index and T1_char_index
                    if (bit_index == 7)
                    {
                        bit_index = 0;
                        ++T1_char_index;
                    } else
                    {
                        ++bit_index;
                    }             
                }
            }
            // Did I get the info I expected?
            if ((T1_char_index * EIGHT + bit_index) != this->T1_nbBits)
            {
                std::cout << "In '--T1_MutationRate', expected " << this->T1_nbBits << " elements. But " << T1_char_index * EIGHT + bit_index << " elements were received.\n";
                abort();
            }

        } else if (Mode.compare("unif")==0)
        {
            this->T1_MutationRate.push_back(input.GetNextElementDouble() * this->T1_nbBits);
            if (this->T1_MutationRate.back() < 0)
            {
                std::cout << "T1_MutationRate of Mode 'unif' is lower than zero" << std::endl;
                abort();
            }
        } else
        {
            std::cout << "Sorry, for 'T1_MutationRate' only Mode 'unif' and 'A' are recognized so far. Mode received (" << Mode << ") is unknown!" << std::endl;
            abort();
        }
        
        assert(this->T1_MutationRate.size() == this->T1_nbBits || this->T1_MutationRate.size() == 1);
        this->T1_Total_Mutation_rate = this->T1_MutationRate.back();
    }

    
}

void SpeciesSpecificParameters::readT2_MutationRate(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option '--T2_MutationRate', the std::string that is read is: " << input.print() << std::endl;
#endif

    if (input.PeakNextElementString() == "NA")
    {
        if (this->T2_nbChars != 0)
        {
            std::cout << "In option '--T2_MutationRate', for species "<<this->speciesName<<", received 'NA' but there are T2 loci as indiciated in --L (--Loci) (this->T2_nbChars = "<<this->T2_nbChars<<").\n";
            abort();   
        }
        input.skipElement();
    } else
    {
        std::string Mode = input.GetNextElementString();
        if (Mode.compare("A")==0)
        {
            double currentSum = 0;
            int T2_char_index = 0;
            while (input.IsThereMoreToRead())
            {
                if (input.isNextRKeyword())
                {
                    double inputValue = input.GetNextElementDouble();
                    int nbRepeats     = input.GetNextElementInt();
                    T2_char_index += nbRepeats;
                    for (int repeat = 0 ; repeat < nbRepeats ; ++repeat)
                    {
                        currentSum += inputValue;
                        this->T2_MutationRate.push_back(currentSum);
                    }

                    // security
                    if (inputValue < 0)
                    {
                        std::cout << "In 'T2_MutationRate' of Mode 'A' the value for the " << T2_char_index << "th T2_char_index is negative (is " << inputValue << ")" << std::endl;
                        abort();
                    }
                    if (inputValue > 0.5)
                    {
                        std::cout << "In 'T2_MutationRate' of Mode 'A' the value for the " << T2_char_index << "th T2_char_index is greater than 0.5 (is " << inputValue << "). Is it really what you want to do?" << std::endl;
                        abort();
                    }
                } else
                {
                    double inputValue = input.GetNextElementDouble();
                    currentSum += inputValue;
                    this->T2_MutationRate.push_back(currentSum);
                    ++T2_char_index;

                    // security
                    if (inputValue < 0)
                    {
                        std::cout << "In 'T2_MutationRate' of Mode 'A' the value for the " << T2_char_index << "th T2_char_index is negative (is " << this->T2_MutationRate.back() << ")" << std::endl;
                        abort();
                    }
                    if (inputValue > 0.5)
                    {
                        std::cout << "In 'T2_MutationRate' of Mode 'A' the value for the " << T2_char_index << "th T2_char_index is greater than 0.5 (is " << this->T2_MutationRate.back() << "). Is it really what you want to do?" << std::endl;
                        abort();
                    }
                }
            }
        } else if (Mode.compare("unif")==0)
        {
            this->T2_MutationRate.push_back(input.GetNextElementDouble() * this->T2_nbChars);
            if (this->T2_MutationRate.back() < 0)
            {
                std::cout << "T2_MutationRate of Mode 'unif' is lower than zero" << std::endl;
                abort();
            }
        } else
        {
            std::cout << "Sorry, for 'T2_MutationRate' only Mode 'unif' and 'A' are recognized so far. Mode received (" << Mode << ") is unknown!" << std::endl;
            abort();
        }
        this->T2_Total_Mutation_rate = this->T2_MutationRate.back();
        assert(this->T2_MutationRate.size() == this->T2_nbChars || this->T2_MutationRate.size() == 1);
        
    }
        
    
}

void SpeciesSpecificParameters::readT3_MutationRate(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option 'T3_MutationRate', the std::string that is read is: " << input.print() << std::endl;
#endif

    if (input.PeakNextElementString() == "NA")
    {
        if (this->T3_nbChars != 0)
        {
            std::cout << "In option '--T3_MutationRate', for species "<<this->speciesName<<", received 'NA' but there are T3 loci as indiciated in --L (--Loci) (this->T3_nbChars = "<<this->T3_nbChars<<").\n";
            abort();   
        }
        input.skipElement();
    } else
    {
        std::string Mode = input.GetNextElementString();
        
        if (Mode.compare("A")==0)
        {
            double currentSum = 0;
            int T3_char_index = 0;
            while (input.IsThereMoreToRead())
            {
                if (input.isNextRKeyword())
                {
                    double inputValue = input.GetNextElementDouble();
                    int nbRepeats     = input.GetNextElementInt();
                    T3_char_index += nbRepeats;
                    for (int repeat = 0 ; repeat < nbRepeats ; ++repeat)
                    {
                        currentSum += inputValue;
                        this->T3_MutationRate.push_back(currentSum);
                    }

                    // security
                    if (inputValue < 0)
                    {
                        std::cout << "In 'T3_MutationRate' of Mode 'A' the value for the " << T3_char_index << "th T3_char_index is negative (is " << inputValue << ")" << std::endl;
                        abort();
                    }
                    if (inputValue > 0.5)
                    {
                        std::cout << "In 'T3_MutationRate' of Mode 'A' the value for the " << T3_char_index << "th T3_char_index is greater than 0.5 (is " << inputValue << "). Is it really what you want to do?" << std::endl;
                        abort();
                    }
                } else
                {
                    double inputValue = input.GetNextElementDouble();
                    currentSum += inputValue;
                    this->T3_MutationRate.push_back(currentSum);
                    ++T3_char_index;

                    // security
                    if (inputValue < 0)
                    {
                        std::cout << "In 'T3_MutationRate' of Mode 'A' the value for the " << T3_char_index << "th T3_char_index is negative (is " << this->T3_MutationRate.back() << ")" << std::endl;
                        abort();
                    }
                    if (inputValue > 0.5)
                    {
                        std::cout << "In 'T3_MutationRate' of Mode 'A' the value for the " << T3_char_index << "th T3_char_index is greater than 0.5 (is " << this->T3_MutationRate.back() << "). Is it really what you want to do?" << std::endl;
                        abort();
                    }
                }
            }
        } else if (Mode.compare("unif")==0)
        {
            this->T3_MutationRate.push_back(input.GetNextElementDouble() * this->T3_nbChars);
            if (this->T3_MutationRate.back() < 0)
            {
                std::cout << "T3_MutationRate of Mode 'unif' is lower than zero" << std::endl;
                abort();
            }
        } else
        {
            std::cout << "Sorry, for 'T3_MutationRate' only Mode 'unif' and 'A' are recognized so far. Mode received (" << Mode << ") is unknown!" << std::endl;
            abort();
        }
        this->T3_Total_Mutation_rate = this->T3_MutationRate.back();
        assert(this->T3_MutationRate.size() == this->T3_nbChars || this->T3_MutationRate.size() == 1);
        
    }
    
    
}


void SpeciesSpecificParameters::readT3_PhenotypicEffects(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option 'T3_PhenotypicEffects', the std::string that is read is: " << input.print() << std::endl;
#endif

    if (input.PeakNextElementString() == "NA")
    {
        if (this->T3_nbChars != 0)
        {
            std::cout << "In option 'T3_PhenotypicEffects', for species "<<this->speciesName<<", received 'NA' but there are T3 loci as indiciated in --L (--Loci) (this->T3_nbChars = "<<this->T3_nbChars<<").\n";
            abort();   
        }
        input.skipElement();
    } else
    {
        this->T3_PhenoNbDimensions = input.GetNextElementInt();
        if (this->T3_PhenoNbDimensions < 1)
        {
            std::cout << "In option 'T3_PhenotypicEffects', the first entry should be the number of dimensions of the fitness space. The minimal possible value is 1. Received " << this->T3_PhenoNbDimensions  << "\n";
            abort();
        }
        if (this->T3_PhenoNbDimensions > this->T3_nbChars)
        {
            std::cout << "In option 'T3_PhenotypicEffects', the first entry should be the number of dimensions of the fitness space. The number of dimensions received ("<<this->T3_PhenoNbDimensions <<") is higher than the number of T3 loci ("<<this->T3_nbChars<<"). Is it really what you want? If yes, turn this security check off.\n";
            abort();
        }
        
        for (int i = 0 ; i < this->T3_PhenoNbDimensions; i++)
        {
            Individual::T3_IndPhenotype.push_back(0.0); // This is a static variable that is used by all inds to calculate fitness on T3 trait
        }

        for (int habitat = 0; habitat <= this->MaxHabitat ; habitat++)
        {
            (void) input.GetNextHabitatMarker(habitat);

            std::vector<double> T3_PhenotypicEffects_OneHabitat;
            
            std::string Mode = input.GetNextElementString();

            if (Mode.compare("A") == 0)
            {
                for (int T3_locus = 0 ; T3_locus < this->T3_nbChars ; ++T3_locus)
                {
                    for (int dim = 0 ; dim < this->T3_PhenoNbDimensions ; dim++)
                    {
                        T3_PhenotypicEffects_OneHabitat.push_back(input.GetNextElementDouble());
                    }
                }  
            } else if (Mode.compare("unif") == 0)
            {
                std::vector<double> v;
                for (int dim = 0 ; dim < this->T3_PhenoNbDimensions ; dim++)
                {
                    v.push_back(input.GetNextElementDouble());
                }
                assert(v.size() == this->T3_PhenoNbDimensions );
                for (int T3_locus = 0 ; T3_locus < this->T3_nbChars ; ++T3_locus)
                {
                    for (int dim = 0 ; dim < this->T3_PhenoNbDimensions ; dim++)
                    {
                        T3_PhenotypicEffects_OneHabitat.push_back(v[dim]);
                    }
                }  
            } else
            {
                std::cout << "In option 'T3_PhenotypicEffects', received Mode "<< Mode << ". Sorry only modes 'A' and 'unif' are recognized for the moment. \n";
                abort();
            }
            this->T3_PhenotypicEffects.push_back(T3_PhenotypicEffects_OneHabitat);
        } // end of habitat for loop
    }
 
   
    
}

void SpeciesSpecificParameters::readCentralT1LocusForExtraGeneticInfo(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option '--centralT1LocusForExtraGeneticInfo', the std::string that is read is: " << input.print() << std::endl;
#endif
    this->centralT1LocusForExtraGeneticInfo = input.GetNextElementInt();
      
}


void SpeciesSpecificParameters::readT3_FitnessLandscape(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option '--T3_fit (--T3_fitnessLandscape)', the std::string that is read is: " << input.print() << std::endl;
#endif

    if (input.PeakNextElementString() == "NA")
    {
        if (this->T3_nbChars != 0)
        {
            std::cout << "In option 'T3_FitnessLandscape', for species "<<this->speciesName<<", received 'NA' but there are T3 loci as indiciated in --L (--Loci) (this->T3_nbChars = "<<this->T3_nbChars<<").\n";
            abort();   
        }
        input.skipElement();
    } else
    {
        std::string SelectionMode = input.GetNextElementString();
        for (int habitat = 0; habitat <= this->MaxHabitat ; habitat++)
        {
            (void) input.GetNextHabitatMarker(habitat);
        
            std::string EntryMode = input.GetNextElementString();

            std::vector<double> T3_fitnessLandscapeOptimum_OneHabitat;
            std::vector<double> T3_fitnessLandscapeLinearGradient_OneHabitat;
            std::vector<double> T3_fitnessLandscapeGaussStrength_OneHabitat;

            if (SelectionMode.compare("linear")==0)
            {
                this->T3_fitnessLandscapeType='L';
                if (EntryMode.compare("A")==0)
                {
                    for (int dim = 0 ; dim < this->T3_PhenoNbDimensions  ; dim++)
                    {
                        double mean = input.GetNextElementDouble();
                        double gradient = input.GetNextElementDouble();

                        if (mean < CHAR_MIN || mean > CHAR_MAX)
                        {
                            std::cout << "In '--T3_FitnessLandscape', SelectionMode 'linear', EntryMode 'A', for dimension "<<dim<<" (counting zero-based) received a mean of "<<mean<<" which is either lower than the minimal possible mean (" << CHAR_MIN << ") or higher that the maximal possible mean ("<<CHAR_MAX<<")\n";
                            abort();
                        }
                        if (gradient < 0.0)
                        {
                            std::cout << "In '--T3_FitnessLandscape', SelectionMode 'linear', EntryMode 'A', for dimension "<<dim<<" (counting zero-based) received a gradient of "<<gradient<<". Gradient must be non-negative.\n";
                            abort();
                        }
                        if (gradient > 1.0)
                        {
                            std::cout << "In '--T3_FitnessLandscape', SelectionMode 'linear', EntryMode 'A', for dimension "<<dim<<" (counting zero-based) received a gradient of "<<gradient<<". Gradient cannot be greater than 1. As 1 is the smallest possible deviation from the optimal phenotype, a gradient of 1 would mean that any deviation would be lethal. Can be worst than that!\n";
                            abort();
                        }

                        T3_fitnessLandscapeOptimum_OneHabitat.push_back(mean);
                        T3_fitnessLandscapeLinearGradient_OneHabitat.push_back(gradient);

                    }
                } else if (EntryMode.compare("unif")==0) 
                {
                    double mean = input.GetNextElementDouble();
                    double gradient = input.GetNextElementDouble();
                    if (mean < CHAR_MIN || mean > CHAR_MAX)
                    {
                        std::cout << "In '--T3_FitnessLandscape', SelectionMode 'linear', EntryMode 'unif', received a mean of "<<mean<<" which is either lower than the minimal possible mean (" << CHAR_MIN << ") or higher that the maximal possible mean ("<<CHAR_MAX<<")\n";
                        abort();
                    }
                    if (gradient < 0.0)
                    {
                        std::cout << "In '--T3_FitnessLandscape', SelectionMode 'linear', EntryMode 'unif', received a gradient of "<<gradient<<". Gradient must be non-negative.\n";
                        abort();
                    }
                    if (gradient > 1.0)
                        {
                            std::cout << "In '--T3_FitnessLandscape', SelectionMode 'linear', EntryMode 'unif', (counting zero-based) received a gradient of "<<gradient<<". Gradient cannot be greater than 1. As 1 is the smallest possible deviation from the optimal phenotype, a gradient of 1 would mean that any deviation would be lethal. Can be worst than that!\n";
                            abort();
                        }

                    for (int dim = 0 ; dim < this->T3_PhenoNbDimensions  ; dim++)
                    {
                        T3_fitnessLandscapeOptimum_OneHabitat.push_back(mean);
                        T3_fitnessLandscapeLinearGradient_OneHabitat.push_back(gradient);
                    }
                } else
                {
                    std::cout << "In option '--T3_fit (--T3_fitnessLandscape)', received EntryMode " << EntryMode << ". The only EntryMode accepted are 'A' and 'unif'. Sorry.\n";
                    abort();       
                }
            } else if (SelectionMode.compare("gaussian")==0 || SelectionMode.compare("gauss")==0)
            {
                this->T3_fitnessLandscapeType='G';
                if (EntryMode.compare("A")==0)
                {
                    for (int dim = 0 ; dim < this->T3_PhenoNbDimensions  ; dim++)
                    {
                        double mean = input.GetNextElementDouble();
                        double strength = input.GetNextElementDouble();

                        if (mean < CHAR_MIN || mean > CHAR_MAX)
                        {
                            std::cout << "In '--T3_FitnessLandscape', SelectionMode 'gauss', EntryMode 'A', for dimension "<<dim<<" (counting zero-based) received a mean of "<<mean<<" which is either lower than the minimal possible mean (" << CHAR_MIN << ") or higher that the maximal possible mean ("<<CHAR_MAX<<")\n";
                            abort();
                        }
                        if (strength <= 0.0)
                        {
                            std::cout << "In '--T3_FitnessLandscape', SelectionMode 'gauss', EntryMode 'A', for dimension "<<dim<<" (counting zero-based) received a strength of "<<strength<<". strength must be non-negative.\n";
                            abort();
                        }

                        T3_fitnessLandscapeOptimum_OneHabitat.push_back(mean);
                        T3_fitnessLandscapeGaussStrength_OneHabitat.push_back(strength);

                    }
                } else if (EntryMode.compare("unif")==0) 
                {
                    double mean = input.GetNextElementDouble();
                    double strength = input.GetNextElementDouble();
                    
                    if (mean < CHAR_MIN || mean > CHAR_MAX)
                    {
                        std::cout << "In '--T3_FitnessLandscape', SelectionMode 'gauss', EntryMode 'unif', received a mean of "<<mean<<" which is either lower than the minimal possible mean (" << CHAR_MIN << ") or higher that the maximal possible mean ("<<CHAR_MAX<<")\n";
                        abort();
                    }
                    if (strength <= 0.0)
                    {
                        std::cout << "In '--T3_FitnessLandscape', SelectionMode 'gauss', EntryMode 'unif', received a strength of "<<strength<<". strength must be non-negative.\n";
                        abort();
                    }

                    assert(this->T3_PhenoNbDimensions > 0);
                    for (int dim = 0 ; dim < this->T3_PhenoNbDimensions ; dim++)
                    {
                        T3_fitnessLandscapeOptimum_OneHabitat.push_back(mean);
                        T3_fitnessLandscapeGaussStrength_OneHabitat.push_back(strength);
                    }
                } else
                {
                    std::cout << "In option '--T3_fit (--T3_fitnessLandscape)', received EntryMode " << EntryMode << ". The only EntryMode accepted are 'A' and 'unif'.\n";
                    abort();       
                }
            } else
            {
                std::cout << "In option '--T3_fit (--T3_fitnessLandscape)', the only SelectionMode accepted are 'linear' and 'gaussian' (or 'gauss'). Mode received = '"<<SelectionMode<< "'.\n";
                abort();
            }
            this->T3_fitnessLandscapeOptimum.push_back(T3_fitnessLandscapeOptimum_OneHabitat);

            if (SelectionMode.compare("gaussian")==0 || SelectionMode.compare("gauss")==0)
            {
                this->T3_fitnessLandscapeGaussStrength.push_back(T3_fitnessLandscapeGaussStrength_OneHabitat);
            } else if (SelectionMode.compare("linear")==0)
            {
                this->T3_fitnessLandscapeLinearGradient.push_back(T3_fitnessLandscapeLinearGradient_OneHabitat);    
            } else if (SelectionMode != "NA")
            {
                std::cout << "Unknown SelectionMode in T3_FitnessLandscape. The error should have been detected earlier in the process. There is therefore an internal error but you might want to check your input parameters anyway.\n";
                abort();
            }
            
        } // end of for habitat
        assert(this->T3_fitnessLandscapeOptimum.size() == this->MaxHabitat+1);
        if (SelectionMode.compare("gaussian")==0 || SelectionMode.compare("gauss")==0)
        {
            assert(this->T3_fitnessLandscapeGaussStrength.size() == this->MaxHabitat+1);
        } else if (SelectionMode.compare("linear")==0)
        {
            assert(this->T3_fitnessLandscapeLinearGradient.size() == this->MaxHabitat+1);
        }
            
    }

    
}



void SpeciesSpecificParameters::readT3_DevelopmentalNoise(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option '--T3_DN (--T3_DevelopmentalNoise)', the std::string that is read is: " << input.print() << std::endl;
#endif

    if (this->T3_nbChars)
    {
        for (int habitat = 0; habitat <= this->MaxHabitat ; habitat++)
        {
            (void) input.GetNextHabitatMarker(habitat);
        
            std::string EntryMode = input.GetNextElementString();

            std::vector<double> DN_OneHabitat; // each value is for one dimension of the phenotypic space
            

            if (EntryMode.compare("A")==0)
            {
                for (int dim = 0 ; dim < this->T3_PhenoNbDimensions  ; dim++)
                {
                    double sd = input.GetNextElementDouble();

                    if (sd < 0)
                    {
                        std::cout << "In '--T3_DevelopmentalNoise', mode = A. Received a negative standard deviation. sd = "<< sd << "\n";
                        abort();
                    }

                    DN_OneHabitat.push_back(sd);
                }
            } else if (EntryMode.compare("unif")==0) 
            {
                double sd = input.GetNextElementDouble();
                if (sd < 0)
                {
                    std::cout << "In '--T3_DevelopmentalNoise', mode = unif. Received a negative standard deviation. sd = "<< sd << "\n";
                    abort();
                }
                for (int dim = 0 ; dim < this->T3_PhenoNbDimensions  ; dim++)
                {
                    DN_OneHabitat.push_back(sd);
                }
            } else
            {
                std::cout << "In option '--T3_fit (--T3_fitnessLandscape)', received EntryMode " << EntryMode << ". The only EntryMode accepted are 'A' and 'unif'. Sorry.\n";
                abort();       
            }

            this->T3_DevelopmentalNoiseStandardDeviation.push_back(DN_OneHabitat);
            
        } // end of for habitat

        assert(this->T3_DevelopmentalNoiseStandardDeviation.size() == this->MaxHabitat+1);
        
    } else
    {
        this->T3_PhenoNbDimensions = 0;
        input.consideredFullyRead();
    }
}


void SpeciesSpecificParameters::readT4_MutationRate(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option 'T4_MutationRate', the std::string that is read is: " << input.print() << std::endl;
#endif
    std::string mode = input.GetNextElementString();
    if (mode == "unif")
    {
        double constantRate = input.GetNextElementDouble();
        T4_MutationRate.push_back(constantRate); // yes, just one value. It will be taken care of. Not very clean though as a design
    } else if (mode == "A")
    {
        if (T4_nbBits > 0)
            T4_MutationRate.push_back(input.GetNextElementDouble());
        for (size_t locus = 1 ; locus < T4_nbBits ; locus++)
        {
            T4_MutationRate.push_back(T4_MutationRate.back() + input.GetNextElementDouble());
        }
    } else
    {
        std::cout << "For option 'T4_MutationRate', received the mode of entry '"<<mode<<"'. Sorry, only modes 'unif' and 'A' are accepted.\n";
        abort();
    }
}

void SpeciesSpecificParameters::readT4_maxAverageNbNodesPerHaplotype(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option 'T4_maxAverageNbNodesPerHaplotype', the std::string that is read is: " << input.print() << std::endl;
#endif 
    if (input.PeakNextElementString() == "default")
    {
        input.skipElement();
        T4_maxAverageNbNodesPerHaplotypeBeforeRecalculation = 100;
    } else
    {
        T4_maxAverageNbNodesPerHaplotypeBeforeRecalculation = input.GetNextElementDouble();
    }
}


void SpeciesSpecificParameters::readRecRateOnMismatch(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option 'recRateOnMismatch', the std::string that is read is: " << input.print() << std::endl;
#endif

    if (input.PeakNextElementString() == "NA")
    {
        recRateOnMismatch_bool = false;
        input.skipElement();
    } else
    {
        if (this->T1_nbBits != this->TotalNbLoci)
        {
            std::cout << "Sorry, recRateOnMismatch can only be used if you only use T1 loci.\n";
            abort();
        }

        recRateOnMismatch_bool = true;
        recRateOnMismatch_halfWindow = input.GetNextElementInt() / 2;
        recRateOnMismatch_factor = input.GetNextElementDouble();
        
        if (recRateOnMismatch_halfWindow*2 > this->T1_nbBits )
        {
            std::cout << "recRateOnMismatch_Window must be smaller or equal to the total number ot T1 loci. Received recRateOnMismatch_Window = " << 2 * recRateOnMismatch_halfWindow << " while there are " << this->T1_nbBits << " T1 loci." << std::endl;
            abort();   
        }
        if (recRateOnMismatch_factor >= 0.0)
        {
            std::cout << "recRateOnMismatch_factor must be negative. Received recRateOnMismatch_factor = " << recRateOnMismatch_factor << std::endl;
            abort();
        }
    }


    
}

void SpeciesSpecificParameters::readRecombinationRate(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option 'RecombinationRate', the std::string that is read is: " << input.print() << std::endl;
#endif

    if (input.PeakNextElementString() == "NA")
    {
        if (this->T1_nbBits != 0)
        {
            std::cout << "In option '--T1_FitnessEffects', for species "<<this->speciesName<<", received 'NA' but there are T1 loci as indiciated in --L (--Loci) (this->T1_nbBits = "<<this->T1_nbBits<<").\n";
            abort();   
        }
        input.skipElement();
    } else
    {

    }

    std::string unit = input.GetNextElementString();
    std::string Mode = input.GetNextElementString();
    
    if (Mode.compare("unif")==0)
    {
        double InputValue = input.GetNextElementDouble();
        double AddToParam;
        double cumsum = 0.0;
        
        if (InputValue == -1.0)
        {
            AddToParam = -1.0;
        } else
        {
            if (unit.compare("M")==0)
            {
                AddToParam = InputValue;
            } else if (unit.compare("cM")==0)
            {
                AddToParam = InputValue / 100;
            } else if (unit.compare("rate")==0)
            {
                AddToParam = - log(1 - 2 * InputValue)/2;
            } else
            {
                std::cout << "First element of option 'RecombinationRate' should be the unit (either 'M' (Morgan), 'cM' (centiMorgan) or 'c'(crossOver rate)) should be inputed. Value received is " << unit << ".\n";
                abort();
            }
        }

        if (AddToParam > 10)
        {
            std::cout << "In option '--RecombinationRate' received a value which equivalent 'M' is higher than 10 (received " << AddToParam << ". Is it really what you wanted? It is not going to be very efficient. Note you can have independent chromosomes with '-1' if this is what you want.\n";
            abort();
        }

        if (AddToParam == -1.0)
        {
            for (int pos=0 ; pos < (this->TotalNbLoci - 1); pos++)
            {
                this->ChromosomeBoundaries.push_back(pos);
            }
            cumsum += 0.5;
            this->RecombinationRate.push_back(cumsum);
        } else
        {
            assert(AddToParam >= 0.0);
            cumsum += AddToParam;
            this->RecombinationRate.push_back(cumsum);
        }
    } else if(Mode.compare("A")==0)
    {
        double cumsum = 0.0;
        int pos=0;
        
        while (input.IsThereMoreToRead())
        {
            int NbRepeats;
            double InputValue;
            double AddToParam;
            if (input.isNextRKeyword())
            {
                InputValue = input.GetNextElementDouble();
                NbRepeats = input.GetNextElementDouble();
                if (NbRepeats < 0)
                {
                    std::cout << "In option '--RecombinationRate', NbRepeats (not directly after indication 'R') is negative\n";
                    abort();
                }
            } else
            {
                InputValue = input.GetNextElementDouble();
                NbRepeats = 1;
            }
            if (unit.compare("M")==0)
            {
                if (InputValue==-1)
                {
                    AddToParam = -1;
                } else if (InputValue >= 0)
                {
                    AddToParam = InputValue;
                } else
                {
                    std::cout << "In option '--RecombinationRate', Unit 'M' received a value which is negative but different from -1. Value received is " << InputValue << ".\n";
                    abort();
                }
            } else if (unit.compare("cM")==0)
            {
                if (InputValue==-1)
                {
                    AddToParam = -1;
                } else if (InputValue >= 0)
                {
                    AddToParam = InputValue / 100;
                } else
                {
                    std::cout << "In option '--RecombinationRate', Unit 'cM' received a value which is negative but different from -1. Value received is " << InputValue << ".\n";
                    abort();
                }
            } else if (unit.compare("rate")==0)
            {
                if (InputValue==-1 || InputValue==0.5)
                {
                    AddToParam = -1;
                } else if (InputValue >= 0 && InputValue < 0.5) 
                {
                    AddToParam = - log(1 - 2 * InputValue)/2;
                } else
                {
                    std::cout << "In option '--RecombinationRate', Unit 'rate' received a value which is either negative but different from -1 or greater than 0.5. Value received is " << InputValue << ".\n";
                    abort();
                }
            } else
            {
                std::cout << "First element of option 'RecombinationRate' should be the unit (either 'M' (Morgan), 'cM' (centiMorgan) or 'c'(crossOver rate)) should be inputed. Value received is " << unit << ".\n";
                abort();
            }
            if (AddToParam > 10)
            {
                std::cout << "In option '--RecombinationRate' received a value which equivalent 'M' is higher than 10 (received " << AddToParam << ". Is it really what you wanted? It is not going to be very efficient!  Note you can have independent chromosomes with '-1' if this is what you want.\n";
                abort();
            }
            for (int rep = 0 ; rep < NbRepeats ; rep++)
            {
                if (AddToParam == -1)
                {
                    cumsum += 0.5;
                    this->ChromosomeBoundaries.push_back(pos);
                    this->RecombinationRate.push_back(cumsum);
                } else if (AddToParam >= 0)
                {
                    cumsum += AddToParam;
                    this->RecombinationRate.push_back(cumsum);
                } else
                {
                    std::cout << "Internal Error: Should not get here in function 'readRecombinationRate'\n";
                    abort();
                }
                pos++;
            }
        }
    } else
    {
        std::cout << "Sorry, for 'RecombinationRate' only Mode 'unif' and 'A' are recognized so far. Mode received (" << Mode << " and unit " << unit << ") is unknown!" << std::endl;
        abort();
    }
    if (this->RecombinationRate.size() == 1)
    {
        this->TotalRecombinationRate = this->RecombinationRate[0] * this->TotalNbLoci;
    } else if (this->RecombinationRate.size() > 1)
    {
        this->TotalRecombinationRate = this->RecombinationRate.back();
    } else
    {
        abort();
    }
        
    assert(this->TotalRecombinationRate >= 0);
    
    if (this->RecombinationRate.size() != this->TotalNbLoci - 1 && this->RecombinationRate.size() != 1)
    {
        std::cout << "\nthis->RecombinationRate.size() = " << this->RecombinationRate.size() << "    this->T1_nbBits = " << this->T1_nbBits << "    this->T2_nbChars = " << this->T2_nbChars  << "    this->T3_nbChars = " << this->T3_nbChars << "\n\n";
        std::cout << "In '--RecombinationRate' did not receive the expected number of elements!" << std::endl;
        abort();
    }
    assert(this->ChromosomeBoundaries.size() < this->TotalNbLoci);

    if (this->TotalNbLoci < 2)
    {
        this->TotalRecombinationRate = 0;
    }
    
    
}


/*void SpeciesSpecificParameters::readT1_vcfOutput_sequence(InputReader& input)
{
#ifdef DEBUG
    std::cout << "For option 'T1_vcfOutput_sequence', the std::string that is read is: " << input.print() << std::endl;
#endif

    std::string Mode = input.GetNextElementString();
    if (Mode.compare("range") == 0)
    {
        int from = input.GetNextElementInt();
        int to = input.GetNextElementInt();
        if (from > to)
        {
            std::cout << "In option 'T1_vcfOutput_sequence', 'from' is greater than 'to' from = " << from << ", to = " << to << ".\n";
            abort();
        }
        if (from < 0)
        {
            std::cout << "In option 'T1_vcfOutput_sequence', 'from' is lower than zero from = " << from << ", to = " << to << ".\n";
            abort();   
        }
        if (to >= this->T1_nbBits)
        {
            std::cout << "In option 'T1_vcfOutput_sequence', 'to' is greater or equal to the number of T1 loci from = " << from << ", to = " << to << ", number of T1 loci = " <<  this->T1_nbBits << ".\n";
            abort();   
        }
        for (int locus = from ; locus <= to ; locus++)
        {

            T1_locusDescription T1_locus(locus / EIGHT, locus % EIGHT, locus);
            T1_vcfOutput_sequenceList.push_back(T1_locus);
        }

        T1_vcfOutput_sequenceIsRange = true;
    } else if (Mode.compare("A") == 0)
    {
        std::string sub;
        T1_vcfOutput_sequenceIsRange = true; // will be set to false if loci dont follow in a range
        while (input.IsThereMoreToRead())
        {
            int locus = input.GetNextElementInt();
            if (locus < 0)
            {
                std::cout << "In option 'T1_vcfOutput_sequence', 'locus' is lower than zero. locus = " << locus << ".\n";
                abort();   
            }
            if (locus >= this->T1_nbBits)
            {
                std::cout << "In option 'T1_vcfOutput_sequence', 'locus' is equal or greater than the number the T1 loci. locus = " << locus << ". Number of T1 loci = " << this->T1_nbBits << ". Please note that the loci indices here must be zero-based counting. The first locus has index '0'. The last locus has index 'Number of T1 loci - 1'\n";
                abort();   
            }
            T1_locusDescription T1_locus(locus / EIGHT, locus % EIGHT, locus);
            // Test if it is in a range
            if (T1_vcfOutput_sequenceList.size() > 0)
            {
                if (T1_vcfOutput_sequenceList.back().locus !=  T1_locus.locus - 1)
                {
                    T1_vcfOutput_sequenceIsRange = false;
                }
            }
            // Add to T1_vcfOutput_sequenceList
            T1_vcfOutput_sequenceList.push_back(T1_locus);     
        }
    } else
    {
        std::cout << "In option '--T1_vcfOutput_sequence', the only two Modes accepted are 'range' and 'A'. Mode received is '" << Mode << "\n";
        abort();
    }
    std::sort(
        T1_vcfOutput_sequenceList.begin(),
        T1_vcfOutput_sequenceList.end(),
        [](const T1_locusDescription& left, const T1_locusDescription& right)
        {
          return left.locus < right.locus;
        }
    );

    if (T1_vcfOutput_sequenceList.size() > 1)
    {
        assert(T1_vcfOutput_sequenceList[0].locus < T1_vcfOutput_sequenceList[1].locus);
    }

    

}*/


void SpeciesSpecificParameters::IsThereSelection()
{
#ifdef CALLENTRANCEFUNCTIONS
    std::cout << "In 'IsThereSelection' start of T1_isSelection\n";
#endif        
    bool ShouldIbreak = false;
    T1_isSelection=0;
    
    if (T1_nbBits)
    {
        assert(this->T1_FitnessEffects.size() == this->MaxHabitat + 1);
        for (int Habitat = 0 ; Habitat <= this->MaxHabitat ; Habitat++)
        {
            for (int locus = 0 ; locus < this->T1_nbBits ; ++locus)
            {
                if (this->T1_isMultiplicitySelection)
                {
                    
                    assert(this->T1_FitnessEffects[Habitat].size() > locus);
                    assert(this->T1_FitnessEffects[Habitat][locus] >= 0.0);
                    if (this->T1_FitnessEffects[Habitat][locus] != 1.0)
                    {
                        this->T1_isSelection = 1;
                        ShouldIbreak = true;
                        break;
                    }   
                } else
                {
                    for (int geno = 0;geno < this->ploidy+1;geno++)
                    {
                        assert(this->T1_FitnessEffects[Habitat].size() > (THREE * locus + geno));
                        assert(this->T1_FitnessEffects[Habitat][THREE * locus + geno] >= 0.0);
                        if (this->T1_FitnessEffects[Habitat][THREE * locus + geno] != this->T1_FitnessEffects[Habitat][0])
                        {
                            this->T1_isSelection = 1;
                            ShouldIbreak = true;
                            break;
                        }
                    }
                } 
                if (ShouldIbreak) break;
            }
            if (ShouldIbreak) break;
        }
    }
    
  

#ifdef CALLENTRANCEFUNCTIONS
    std::cout << "In 'IsThereSelection' start of T1_isEpistasis\n";
#endif
    ShouldIbreak = false;
    this->T1_isEpistasis = false;
    if (this->T1_nbBits)
    {
        //std::cout << "this->T1_nbChars = " << this->T1_nbChars << "\n";
        int howManyWarningGiven = 0;
        //std::cout << "this->T1_Epistasis_LociIndices.size() = " << this->T1_Epistasis_LociIndices.size() << "\n";
        if (this->T1_Epistasis_LociIndices.size() != 0)
        {    
            //std::cout << "this->MaxHabitat = " << this->MaxHabitat << "\n";
            //std::cout << "this->T1_Epistasis_FitnessEffects.size() = " << this->T1_Epistasis_FitnessEffects.size() << "\n";
            assert(this->T1_Epistasis_FitnessEffects.size() == this->MaxHabitat + 1);
            for (int Habitat = 0 ; Habitat <= this->MaxHabitat ; Habitat++)
            {
                //std::cout << "this->T1_Epistasis_FitnessEffects[Habitat].size() = " << this->T1_Epistasis_FitnessEffects[Habitat].size() << "\n";
                for (int groupOfLociIndex = 0 ; groupOfLociIndex < this->T1_Epistasis_FitnessEffects[Habitat].size() ; groupOfLociIndex++)
                {
                    //std::cout << "groupOfLociIndex = " << groupOfLociIndex << "\n";
                    assert(this->T1_Epistasis_LociIndices[Habitat].size() > groupOfLociIndex);

                    // Send warnings if a locus is under both types of selection
                    for (auto& T1_locus : this->T1_Epistasis_LociIndices[Habitat][groupOfLociIndex])
                    {
                        //std::cout << "T1_locus.locus = " << T1_locus.locus << "\n";

                        
                        if (this->T1_FitnessEffects.size() != 0)
                        {
                            if (this->T1_isMultiplicitySelection)
                            {
                                assert(this->T1_FitnessEffects.size() > Habitat);
                                assert(this->T1_FitnessEffects[Habitat].size() > T1_locus.locus);
                                if (this->T1_FitnessEffects[Habitat][T1_locus.locus] != 1.0)
                                {
                                    if (howManyWarningGiven < 50)
                                    {
                                        std::cout << "\tWARNING: The " << T1_locus.locus << "th locus is under both epistatic selection (--T1_epistasis (--T1_EpistaticFitnessEffects)) and regular selection (--T1_fit (--T1_FitnessEffects)) under habitat " << Habitat << ".\n";
                                        howManyWarningGiven++;
                                    } else
                                    {
                                        if (howManyWarningGiven == 50)
                                        {
                                            std::cout << "You should have already received 50 warnings for having a locus under both types of selection. Further warnings will be silenced.\n";
                                            howManyWarningGiven++;
                                        }
                                    }
                                }
                            } else
                            {
                                bool ShouldIGiveWarning = false;
                                for (int geno=0 ; geno <= this->ploidy ; geno++)
                                {
                                    //std::cout << "geno = " << geno << "\n";
                                    assert(this->T1_FitnessEffects.size() > Habitat);
                                    assert(this->T1_FitnessEffects[Habitat].size() > THREE * T1_locus.locus + geno);
                                    //std::cout << "this->T1_FitnessEffects[Habitat][THREE * T1_locus.locus + geno] = " << this->T1_FitnessEffects[Habitat][THREE * T1_locus.locus + geno] << "\n";
                                    //std::cout << "this->T1_FitnessEffects[Habitat][THREE * T1_locus.locus] = " << this->T1_FitnessEffects[Habitat][THREE * T1_locus.locus] << "\n";
                                    if (this->T1_FitnessEffects[Habitat][THREE * T1_locus.locus + geno] != this->T1_FitnessEffects[Habitat][THREE * T1_locus.locus])
                                    {
                                        ShouldIGiveWarning = true;
                                        break;
                                    }
                                }
                                if (ShouldIGiveWarning)
                                {
                                    if (howManyWarningGiven < 50)
                                    {
                                        std::cout << "\tWARNING: The " << T1_locus.locus << "th locus is under both epistatic selection (--T1_epistasis (--T1_EpistaticFitnessEffects)) and regular selection (--T1_fit (--T1_FitnessEffects)) under habitat " << Habitat << ".\n";
                                        howManyWarningGiven++;
                                    } else
                                    {
                                        if (howManyWarningGiven == 50)
                                        {
                                            std::cout << "You should have already received 50 warnings for having a locus under both types of selection. Further warnings will be silenced.\n";
                                            howManyWarningGiven++;
                                        }
                                    }    
                                }
                            }
                        }
                    }

                    // Check if there is any selection
                    assert(this->T1_Epistasis_FitnessEffects.size() > Habitat);
                    assert(this->T1_Epistasis_FitnessEffects[Habitat].size() > groupOfLociIndex);
                    for (auto& fit : this->T1_Epistasis_FitnessEffects[Habitat][groupOfLociIndex])
                    {
                        if (fit != this->T1_Epistasis_FitnessEffects[Habitat][groupOfLociIndex][0])
                        {
                            this->T1_isEpistasis = true;
                        }
                    }
                } // end of groupOfLoci
            }
        } else
        {
            assert(this->T1_Epistasis_FitnessEffects.size() == 0);
        }
    }


#ifdef CALLENTRANCEFUNCTIONS
    std::cout << "In 'IsThereSelection' start of T2_isSelection\n";
#endif     

    ShouldIbreak = false;
    T2_isSelection=false;
    if (this->T2_nbChars)
    {
        assert(T2_FitnessEffects.size() == this->MaxHabitat + 1);
        for (int Habitat = 0 ; Habitat <= this->MaxHabitat ; Habitat++)
        {
            for (int T2_char_index=0;T2_char_index < this->T2_nbChars;T2_char_index++)
            {
                assert(this->T2_FitnessEffects[Habitat].size() > T2_char_index);
                if (this->T2_FitnessEffects[Habitat][T2_char_index] != 1.0)
                {
                    this->T2_isSelection = true;
                    ShouldIbreak = true;
                    break;
                }
            }
            if (ShouldIbreak) break;
        }
    }


#ifdef CALLENTRANCEFUNCTIONS
    std::cout << "In 'IsThereSelection' start of T3_isSelection\n";
#endif     

    ShouldIbreak = false;
    this->T3_isSelection=false;
    
    if (this->T3_nbChars)
    {
        if (this->T3_fitnessLandscapeType!='L' && this->T3_fitnessLandscapeType!='G')
        {
            std::cout << "Internal error: function 'IsThereSelection' only knows how to deal with T3_fitnessLandscapeType=='L' and T3_fitnessLandscapeType=='G' but it appears that T3_fitnessLandscapeType=="<<this->T3_fitnessLandscapeType<<"\n";
            abort();
        }
        if (this->T3_PhenoNbDimensions < 1)
        {
            std::cout << "Internal error: In 'IsThereSelection' the number of dimensions for T3 is lower than 1. It is "<<this->T3_PhenoNbDimensions<<".\n";
            abort();   
        }

        if (this->T3_fitnessLandscapeType=='L')
        {
            assert(this->T3_fitnessLandscapeLinearGradient.size() == this->T3_fitnessLandscapeOptimum.size());
            assert(this->T3_fitnessLandscapeLinearGradient.size() == this->MaxHabitat + 1);

            for (int Habitat = 0 ; Habitat <= this->MaxHabitat ; Habitat++)
            {
                assert(this->T3_fitnessLandscapeLinearGradient[Habitat].size() == this->T3_fitnessLandscapeOptimum[Habitat].size());
                assert(this->T3_fitnessLandscapeLinearGradient[Habitat].size() == this->T3_PhenoNbDimensions);
                for (int dim = 0 ; dim < this->T3_PhenoNbDimensions ; dim++)
                {
                    if (this->T3_fitnessLandscapeLinearGradient[Habitat][dim] != 0.0)
                    {
                        this->T3_isSelection = true;
                        ShouldIbreak = true;
                        break;
                    }
                }
                if (ShouldIbreak) break;
            }
        } else if (this->T3_fitnessLandscapeType=='G')
        {
            for (auto& elem1 : this->T3_fitnessLandscapeGaussStrength)
            {
                for (auto& elem2 : elem1)
                {
                    assert(elem2 > 0.0);
                }
            }      

            // There is necessarily some selection unless this->T3_fitnessLandscapeGaussian is infinity which is impossible.
            this->T3_isSelection = true;
        } else
        {
            std::cout << "Internal error in funciton IsThereSelection\n";
            abort();
        }
    }
    
#ifdef CALLENTRANCEFUNCTIONS
    std::cout << "In 'IsThereSelection' end of T3_isSelection\n";
#endif
    


#ifdef CALLENTRANCEFUNCTIONS
    std::cout << "In 'IsThereSelection' start of T5_isSelection\n";
#endif        
    ShouldIbreak = false;
    T5_isSelection=0;
    
    if (T5_nbBits)
    {
        assert(this->T5_FitnessEffects.size() == this->MaxHabitat + 1);
        for (int Habitat = 0 ; Habitat <= this->MaxHabitat ; Habitat++)
        {
            for (int locus = 0 ; locus < this->T5_nbBits ; ++locus)
            {
                if (this->T5_isMultiplicitySelection)
                {
                    
                    assert(this->T5_FitnessEffects[Habitat].size() > locus);
                    assert(this->T5_FitnessEffects[Habitat][locus] >= 0.0);
                    if (this->T5_FitnessEffects[Habitat][locus] != 1.0)
                    {
                        this->T5_isSelection = 1;
                        ShouldIbreak = true;
                        break;
                    }   
                } else
                {
                    for (int geno = 0; geno < 2; geno++) // only values for het and double mutant
                    {
                        assert(this->T5_FitnessEffects[Habitat].size() > (2 * locus + geno));
                        assert(this->T5_FitnessEffects[Habitat][2 * locus + geno] >= 0.0);
                        if (this->T5_FitnessEffects[Habitat][2 * locus + geno] != this->T5_FitnessEffects[Habitat][0])
                        {
                            this->T5_isSelection = 1;
                            ShouldIbreak = true;
                            break;
                        }
                    }
                } 
                if (ShouldIbreak) break;
            }
            if (ShouldIbreak) break;
        }
    }
#ifdef CALLENTRANCEFUNCTIONS
    std::cout << "In 'IsThereSelection' end of T5_isSelection\n";
#endif



    // Dispersal probability cannot be a function of fitness if there is no selection
    if (!this->T1_isSelection && !this->T1_isEpistasis && !this->T2_isSelection && !this->T3_isSelection  && !this->T5_isSelection )
    {
        this->DispWeightByFitness = false;
    }

    // if fecundity is different from -1, then dispersal rate necessarily depends on fitness
    if (this->fecundityForFitnessOfOne != -1.0)
    {
        this->DispWeightByFitness = true;   
    }
}

/*void SpeciesSpecificParameters::ClearT1_Initial_AlleleFreqs()
{
    for ( int patch_index = 0 ; patch_index < this->T1_Initial_AlleleFreqs.size() ; ++patch_index )
    {
        this->T1_Initial_AlleleFreqs[patch_index].clear();
    }
    this->T1_Initial_AlleleFreqs.clear();
}*/


/*void SpeciesSpecificParameters::setInitialT1_AlleleFreqTo(const int uniqueFreq)
{
    std::vector<double> line(this->T1_nbChars);
    for (int T1_char_index=0;T1_char_index<this->T1_nbChars;T1_char_index++)
    {
        line[T1_char_index] = uniqueFreq;
    }
    for (int patch_index=0;patch_index<GP->PatchNumber;++patch_index)
    {
        this->T1_Initial_AlleleFreqs.push_back(line);
    }
}*/

std::vector<int> SpeciesSpecificParameters::UpdateParameters(int generation_index)
{
#ifdef DEBUG
    std::cout << "Enters in SpeciesSpecificParameters::UpdateParameters\n";
#endif
    // Change patchCapacity Do not change patchSize yet as it will be used later
    assert(generation_index < this->__patchCapacity.size());
    this->patchCapacity = this->__patchCapacity[generation_index];
    this->patchCapacity.shrink_to_fit();
    assert(this->patchCapacity.size() == GP->PatchNumber);

    // Change Patch size
    assert(this->patchSize.size()>0);
    this->patchSize.resize(this->patchCapacity.size(), 0); // complete with zeros
    std::vector<int> previousPatchSizes = this->patchSize;
    this->patchSize.shrink_to_fit();
    assert(this->patchSize.size() == GP->PatchNumber);

    // Set patch size to carrying capacity if it is not allowed to differ from it
    if (this->fecundityForFitnessOfOne == -1.0)
    {
        for (int patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
        {
            this->patchSize[patch_index] = this->patchCapacity[patch_index];
        }
    }

    // Change GP->AllSpeciesPatchSizes
    GP->saveSSPPatchSize_toGP();

    // Change TotalpatchCapacity
    this->TotalpatchCapacity = 0;
    for ( auto& OnepatchCapacity : this->patchCapacity )
    {
        this->TotalpatchCapacity += OnepatchCapacity;
    }

    // change Habitat
    this->Habitats = this->__Habitats[generation_index];
    this->MaxHabitat = this->__MaxHabitat[generation_index];
    assert(this->Habitats.size() == GP->PatchNumber);

     // Change DispMat
    this->dispersalData.FullFormForwardMigration = this->dispersalData.__FullFormForwardMigration[generation_index];
    std::vector<std::vector<std::vector<double>>> emptyCumSumFits;
    (void) this->dispersalData.SetBackwardMigrationAndGetNextGenerationPatchSizes(
        emptyCumSumFits,
        true
    );

    assert(this->dispersalData.FullFormForwardMigration.size() == GP->PatchNumber);

    return previousPatchSizes; // This return value will be used to know what individuals to duplicate
}


void SpeciesSpecificParameters::readadditiveEffectAmongLoci(InputReader& input)
{
#ifdef DEBUG
    std::cout << "Enters in SpeciesSpecificParameters::readadditiveEffectAmongLoci\n";
#endif    
    auto s = input.GetNextElementString();
    if (s == "n"  || s == "no" || s == "No" || s == "NO" || s == "0")
    {
        //this->multiplicativeEffectAmongLoci = true;
    } else if (s == "y" || s == "yes" || s == "Yes" || s == "YES" || s == "1")
    {
        //this->multiplicativeEffectAmongLoci = false;
    } else
    {
        std::cout << "For option --additiveEffectAmongLoci, received '" << s << "'. Please indicate 'yes' or 'no' ('y' and 'n' and '1' and '0' and a few others work too)\n";
        abort();
    }
}

void SpeciesSpecificParameters::readReadPopFromBinary(InputReader& input)
{
#ifdef DEBUG
    std::cout << "Enters in SpeciesSpecificParameters::readReadPopFromBinary\n";
#endif    

    std::string yesNo = input.GetNextElementString();
    if (yesNo == "yes" || yesNo == "y" || yesNo == "Y" || yesNo == "YES" || yesNo == "true" || yesNo == "1")
    {
        readPopFromBinary = true;
        
        // read file path
        readPopFromBinaryPath = input.GetNextElementString();

        // test that the file exists
        std::ifstream f(readPopFromBinaryPath);
        if (!f.good())
        {
            std::cout << "For species " << speciesName << ", in function 'SpeciesSpecificParameters::readReadPopFromBinary', the file " << readPopFromBinaryPath << " was not found. Note that paths to read binary populations from are NOT relative to the general path (the general path only applies to outputs).\n";
            abort();
        }        
    } else if (yesNo == "no" || yesNo == "n" || yesNo == "N" || yesNo == "NO" || yesNo == "false" || yesNo == "0")
    {
        readPopFromBinary = false;
        readPopFromBinaryPath = "";
    } else
    {
        std::cout << "For option '--readReadPopFromBinary', the first element received is " << yesNo << ". Expected 'yes' (or 'y', 'Y', 'true', '1', etc...), 'no' (or 'n', 'NO', '0', 'false', etc...) or 'default' as to whether the random seed must be initialized based on the seed present in a binary file.\n";
        abort();
    }


    
}



void SpeciesSpecificParameters::readSubsetLociForfitnessSubsetLoci_file(InputReader& input)
{
    // Make sure it starts by LociSet!
    if (input.PeakNextElementString() != "LociSet")
    {
        std::cout << "For option --fitnessStats_file, in function SpeciesSpecificParameters::readSubsetLociForfitnessSubsetLoci_file(InputReader& input), the very first word expected (after the time information and species specific marker) was the keyword 'LociSet'. Instead, SimBit received received "<<input.PeakNextElementString()<<".\n";
        abort();
    }

    int nbLociSets = 0;
    std::string whichT;
    while (input.IsThereMoreToRead())
    {
        if (input.PeakNextElementString() == "LociSet")
        {
            input.skipElement();
            std::string whichTstring = input.PeakNextElementString();
            if (whichTstring != "T1" && whichTstring != "T2" && whichTstring != "T3" && whichTstring != "T1epistasis")
            {
                std::cout << "For option --fitnessStats_file, in function SpeciesSpecificParameters::readSubsetLociForfitnessSubsetLoci_file(InputReader& input), expected a description of the type of locus (T1, T2, T3 or T1epistasis) after the keyword 'LociSet' but instead it received "<<whichTstring<<".\n";
                abort();
            }
            nbLociSets++;

             // prepare vectors
            subsetT1LociForfitnessSubsetLoci_file.push_back({});
            subsetT2LociForfitnessSubsetLoci_file.push_back({});
            subsetT3LociForfitnessSubsetLoci_file.push_back({});
            subsetT5LociForfitnessSubsetLoci_file.push_back({});
            subsetT1epistasisLociForfitnessSubsetLoci_file.push_back({});
            assert(subsetT1LociForfitnessSubsetLoci_file.size() == nbLociSets);
            assert(subsetT2LociForfitnessSubsetLoci_file.size() == nbLociSets);
            assert(subsetT3LociForfitnessSubsetLoci_file.size() == nbLociSets);
            assert(subsetT5LociForfitnessSubsetLoci_file.size() == nbLociSets);
            assert(subsetT1epistasisLociForfitnessSubsetLoci_file.size() == nbLociSets);

        } else
        {
            std::string eventualTstring = input.PeakNextElementString();
            if (eventualTstring.at(0) == 'T')
            {
                input.skipElement();
                whichT = eventualTstring;
                
                if (whichT != "T1" && whichT != "T2" && whichT != "T3" && whichT != "T1epistasis")
                {
                    std::cout << "For option --fitnessStats_file, in function SpeciesSpecificParameters::readSubsetLociForfitnessSubsetLoci_file(InputReader& input), received a type of locus descriptor starting with the letter 'T', however it was not followed by '1', '2', '3' or '1epistasis'. The type of locus descriptor received is "<<whichT<<".\n";
                    abort();
                }
            }
            

            int locus = input.GetNextElementInt();
            if (whichT == "T1")
            {
                //std::cout << "locus = " << locus << "\n";
                //std::cout << "T1_nbBits = " << T1_nbBits << "\n";
                if (locus < 0 || locus >= T1_nbBits)
                {
                    std::cout << "For option --fitnessStats_file, in function SpeciesSpecificParameters::readSubsetLociForfitnessSubsetLoci_file(InputReader& input), received the locus index " << locus << " for trait of type "<<whichT<<". Only positive number lower than the total number of loci is accepted. As a reminder, the first locus has index 0 and the last locus has index nbLociForThisTrait-1 (T1_nbBits = "<<T1_nbBits<<", T2_nbChars = "<<T2_nbChars<<", T3_nbChars = "<<T3_nbChars<<", T5_nbBits = "<<T5_nbBits<<")\n";
                    abort();
                }
                if (
                    std::find(
                        subsetT1LociForfitnessSubsetLoci_file.back().begin(),
                        subsetT1LociForfitnessSubsetLoci_file.back().end(),
                        locus
                    ) != subsetT1LociForfitnessSubsetLoci_file.back().end()
                )
                {
                    std::cout << "For option --fitnessStats_file, in function SpeciesSpecificParameters::readSubsetLociForfitnessSubsetLoci_file(InputReader& input), received the locus index " << locus << " for trait type "<<whichT<<" several times within a single LociSet.\n";
                    abort();
                }
                subsetT1LociForfitnessSubsetLoci_file.back().push_back(locus);
            } else if (whichT == "T2")
            {
                if (locus < 0 || locus >= T2_nbChars)
                {
                    std::cout << "For option --fitnessStats_file, in function SpeciesSpecificParameters::readSubsetLociForfitnessSubsetLoci_file(InputReader& input), received the locus index " << locus << " for trait of type "<<whichT<<". Only positive number lower than the total number of loci is accepted. As a reminder, the first locus has index 0 and the last locus has index nbLociForThisTrait-1 (T1_nbBits = "<<T1_nbBits<<", T2_nbChars = "<<T2_nbChars<<", T3_nbChars = "<<T3_nbChars<<", T5_nbBits = "<<T5_nbBits<<")\n";
                    abort();
                }
                if (
                    std::find(
                        subsetT2LociForfitnessSubsetLoci_file.back().begin(),
                        subsetT2LociForfitnessSubsetLoci_file.back().end(),
                        locus
                    ) != subsetT2LociForfitnessSubsetLoci_file.back().end()
                )
                {
                    std::cout << "For option --fitnessStats_file, in function SpeciesSpecificParameters::readSubsetLociForfitnessSubsetLoci_file(InputReader& input), received the locus index " << locus << " for trait type "<<whichT<<" several times within a single LociSet.\n";
                    abort();
                }
                subsetT2LociForfitnessSubsetLoci_file.back().push_back(locus);
            } else if (whichT == "T3")
            {
                if (locus < 0 || locus >= T3_nbChars)
                {
                    std::cout << "For option --fitnessStats_file, in function SpeciesSpecificParameters::readSubsetLociForfitnessSubsetLoci_file(InputReader& input), received the locus index " << locus << " for trait of type "<<whichT<<". Only positive number lower than the total number of loci is accepted. As a reminder, the first locus has index 0 and the last locus has index nbLociForThisTrait-1 (T1_nbBits = "<<T1_nbBits<<", T2_nbChars = "<<T2_nbChars<<", T3_nbChars = "<<T3_nbChars<<", T5_nbBits = "<<T5_nbBits<<")\n";
                    abort();
                }  
                if (
                    std::find(
                        subsetT3LociForfitnessSubsetLoci_file.back().begin(),
                        subsetT3LociForfitnessSubsetLoci_file.back().end(),
                        locus
                    ) != subsetT3LociForfitnessSubsetLoci_file.back().end()
                )
                {
                    std::cout << "For option --fitnessStats_file, in function SpeciesSpecificParameters::readSubsetLociForfitnessSubsetLoci_file(InputReader& input), received the locus index " << locus << " for trait type "<<whichT<<" several times within a single LociSet.\n";
                    abort();
                }
                subsetT3LociForfitnessSubsetLoci_file.back().push_back(locus); 
            } else if (whichT == "T1epistasis")
            {
                if (locus < 0 || locus >= T1_nbBits)
                {
                    std::cout << "For option --fitnessStats_file, in function SpeciesSpecificParameters::readSubsetLociForfitnessSubsetLoci_file(InputReader& input), received the locus index " << locus << " for trait of type "<<whichT<<". Only positive number lower than the total number of loci is accepted. As a reminder, the first locus has index 0 and the last locus has index nbLociForThisTrait-1 (T1_nbBits = "<<T1_nbBits<<", T2_nbChars = "<<T2_nbChars<<", T3_nbChars = "<<T3_nbChars<<", T5_nbBits = "<<T5_nbBits<<")\n";
                    abort();
                }
                if (
                    std::find(
                        subsetT1epistasisLociForfitnessSubsetLoci_file.back().begin(),
                        subsetT1epistasisLociForfitnessSubsetLoci_file.back().end(),
                        locus
                    ) != subsetT1epistasisLociForfitnessSubsetLoci_file.back().end()
                )
                {
                    std::cout << "For option --fitnessStats_file, in function SpeciesSpecificParameters::readSubsetLociForfitnessSubsetLoci_file(InputReader& input), received the locus index " << locus << " for trait type "<<whichT<<" several times within a single LociSet.\n";
                    abort();
                }
                subsetT1epistasisLociForfitnessSubsetLoci_file.back().push_back(locus);
            } else if (whichT == "T5")
            {
                //std::cout << "locus = " << locus << "\n";
                //std::cout << "T1_nbBits = " << T1_nbBits << "\n";
                if (locus < 0 || locus >= T5_nbBits)
                {
                    std::cout << "For option --fitnessStats_file, in function SpeciesSpecificParameters::readSubsetLociForfitnessSubsetLoci_file(InputReader& input), received the locus index " << locus << " for trait of type "<<whichT<<". Only positive number lower than the total number of loci is accepted. As a reminder, the first locus has index 0 and the last locus has index nbLociForThisTrait-1 (T1_nbBits = "<<T1_nbBits<<", T2_nbChars = "<<T2_nbChars<<", T3_nbChars = "<<T3_nbChars<<", T5_nbBits = "<<T5_nbBits<<")\n";
                    abort();
                }
                if (
                    std::find(
                        subsetT5LociForfitnessSubsetLoci_file.back().begin(),
                        subsetT5LociForfitnessSubsetLoci_file.back().end(),
                        locus
                    ) != subsetT5LociForfitnessSubsetLoci_file.back().end()
                )
                {
                    std::cout << "For option --fitnessStats_file, in function SpeciesSpecificParameters::readSubsetLociForfitnessSubsetLoci_file(InputReader& input), received the locus index " << locus << " for trait type "<<whichT<<" several times within a single LociSet.\n";
                    abort();
                }
                subsetT5LociForfitnessSubsetLoci_file.back().push_back(locus);
            } else
            {
                std::cout << "For option --fitnessStats_file, in function SpeciesSpecificParameters::readSubsetLociForfitnessSubsetLoci_file(InputReader& input), received a type of locus descriptor starting with the letter 'T', however it was not followed by '1', '2', '3', '5' or '1epistasis'. The type of locus descriptor received is "<<whichT<<". Note that this error message comes at a bit of an unexpected time and might therefore be due to an internal error. You might want to check your input nevertheless.\n";
                abort();
            }
        }
    }

    // Make sure at least one set has been received
    assert(subsetT1LociForfitnessSubsetLoci_file.size() == nbLociSets);
    assert(subsetT2LociForfitnessSubsetLoci_file.size() == nbLociSets);
    assert(subsetT3LociForfitnessSubsetLoci_file.size() == nbLociSets);
    assert(subsetT5LociForfitnessSubsetLoci_file.size() == nbLociSets);
    assert(subsetT1epistasisLociForfitnessSubsetLoci_file.size() == nbLociSets);
    if (nbLociSets == 0)
    {
        std::cout << "For option --fitnessStats_file, in function SpeciesSpecificParameters::readSubsetLociForfitnessSubsetLoci_file(InputReader& input), received 0 LociSet.\n";
        abort();
    }

    // Make sure there is not empty set and also just sort them for fun!
    for (int lociSetIndex = 0 ; lociSetIndex < nbLociSets; lociSetIndex++)
    {
        if (
            subsetT1LociForfitnessSubsetLoci_file[lociSetIndex].size() == 0 &&
            subsetT2LociForfitnessSubsetLoci_file[lociSetIndex].size() == 0 &&
            subsetT3LociForfitnessSubsetLoci_file[lociSetIndex].size() == 0 &&
            subsetT5LociForfitnessSubsetLoci_file[lociSetIndex].size() == 0 &&
            subsetT1epistasisLociForfitnessSubsetLoci_file[lociSetIndex].size() == 0
        )
        {
            std::cout << "For option --fitnessStats_file, in function SpeciesSpecificParameters::readSubsetLociForfitnessSubsetLoci_file(InputReader& input), the "<<lociSetIndex<<"th LociSet seem to contain no locus of any type! Maybe the keyword 'LociSet' has been repeated twice (or more times) without indicating loci indices in between. Maybe the option argument ends with the keyword 'LociSet'.\n";
            abort();
        }
        std::sort(
            subsetT1LociForfitnessSubsetLoci_file[lociSetIndex].begin(),
            subsetT1LociForfitnessSubsetLoci_file[lociSetIndex].end()
        );
        std::sort(
            subsetT2LociForfitnessSubsetLoci_file[lociSetIndex].begin(),
            subsetT2LociForfitnessSubsetLoci_file[lociSetIndex].end()
        );
        std::sort(
            subsetT3LociForfitnessSubsetLoci_file[lociSetIndex].begin(),
            subsetT3LociForfitnessSubsetLoci_file[lociSetIndex].end()
        );
        std::sort(
            subsetT5LociForfitnessSubsetLoci_file[lociSetIndex].begin(),
            subsetT5LociForfitnessSubsetLoci_file[lociSetIndex].end()
        );
        std::sort(
            subsetT1epistasisLociForfitnessSubsetLoci_file[lociSetIndex].begin(),
            subsetT1epistasisLociForfitnessSubsetLoci_file[lociSetIndex].end()
        );
    }
}



