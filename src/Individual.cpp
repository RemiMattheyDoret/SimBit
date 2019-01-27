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


Haplotype& Individual::getHaplo(const int& haplo_index)
{
    if (haplo_index)
    {
        return haplo1;
    } else
    {
        return haplo0;
    }  
}

bool Individual::isFreeFromMutations()
{
    return haplo0.isFreeFromMutations() && haplo1.isFreeFromMutations();
}

bool Individual::isFreeFromMutations(int T1_locusFrom, int T1_locusTo)
{
    return haplo0.isFreeFromMutations(T1_locusFrom, T1_locusTo) && haplo1.isFreeFromMutations(T1_locusFrom, T1_locusTo);   
}


double Individual::CalculateT1FitnessMultiplicity(const int& Habitat, int fitnessMapIndex, int T1_locusFrom, int T1_locusTo)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'CalculateT1FitnessMultiplicity'\n";
#endif   

    assert(SSP->T1_isMultiplicitySelection);

    // Get the right fitness array
    auto& fits = SSP->T1_FitnessEffects[Habitat];
    assert(fits.size() == SSP->T1_nbBits);

    /*if (SSP->simTracker.shouldLoopOverAllSites)
    {
        // Find what bits to look at
        int T1_locusToIncluded = T1_locusTo - 1;  
        int T1_ByteIndexFrom = T1_locusFrom / 8;
        int T1_BitIndexFrom  = T1_locusFrom % 8;
        int T1_ByteIndexTo = T1_locusToIncluded / 8;
        int T1_BitIndexTo  = T1_locusToIncluded % 8;

        // calculate fitness
        if (haplo0.getW_T1(fitnessMapIndex) == -1.0 && haplo1.getW_T1(fitnessMapIndex) == -1.0)
        {
            double W_T1haplo0 = 1.0;
            double W_T1haplo1 = 1.0;
            for (int byte_index = T1_ByteIndexFrom ; byte_index <= T1_ByteIndexTo ; byte_index++ )
            {
                int bit_index_to;
                int bit_index_from;
                if (byte_index == T1_ByteIndexTo)
                {
                    bit_index_to = T1_BitIndexTo;
                } else
                {
                    bit_index_to = 7;
                }
                if (byte_index == T1_ByteIndexFrom)
                {
                    bit_index_from = T1_BitIndexFrom;
                } else
                {
                    bit_index_from = 0;
                }
                for (int bit_index = bit_index_from ; bit_index <= bit_index_to ; bit_index++ )
                {
                    //std::cout << byte_index*8+bit_index << " ";
                    if (haplo0.getT1_Allele(byte_index,bit_index))
                    {
                        //std::cout << "Mut at haplo0 " << byte_index << " " << bit_index << " (" << byte_index * 8 + bit_index << ").\n";
                        //std::cout << "SSP->T1_FitnessEffects[Habitat][locus] = " << SSP->T1_FitnessEffects[Habitat][locus] << "\n";
                        int locus  = byte_index * EIGHT + bit_index;
                        //std::cout << "l2085 fits["<<locus<<"] = " << fits[locus] << " ";
                        W_T1haplo0 *= fits[locus];
                    }
                    if (haplo1.getT1_Allele(byte_index,bit_index))
                    {
                        //std::cout << "Mut at haplo1 " << byte_index << " " << bit_index << " (" << byte_index * 8 + bit_index << ").\n";
                        int locus  = byte_index * EIGHT + bit_index;
                        //std::cout << "l2092 fits["<<locus<<"] = " << fits[locus] << " ";
                        W_T1haplo1 *= fits[locus];
                    }
                }
            }
            haplo0.setW_T1(W_T1haplo0, fitnessMapIndex);
            haplo1.setW_T1(W_T1haplo1, fitnessMapIndex);
            //std::cout << "W_T1haplo0 = " << W_T1haplo0 << std::endl;
            //std::cout << "W_T1haplo1 = " << W_T1haplo1 << std::endl;
        } else if (haplo0.getW_T1(fitnessMapIndex) == -1.0)
        {
            double W_T1haplo0 = 1.0;
            for (int byte_index = T1_ByteIndexFrom ; byte_index <= T1_ByteIndexTo ; byte_index++ )
            {
                int bit_index_to;
                int bit_index_from;
                if (byte_index == T1_ByteIndexTo)
                {
                    bit_index_to = T1_BitIndexTo;
                } else
                {
                    bit_index_to = 7;
                }
                if (byte_index == T1_ByteIndexFrom)
                {
                    bit_index_from = T1_BitIndexFrom;
                } else
                {
                    bit_index_from = 0;
                }
                for (int bit_index = bit_index_from ; bit_index <= bit_index_to ; bit_index++ )
                {
                    //std::cout << byte_index*8+bit_index << " ";
                    if (haplo0.getT1_Allele(byte_index,bit_index))
                    {
                        //std::cout << "Mut at haplo0 " << byte_index << " " << bit_index << " (" << byte_index * 8 + bit_index << ").\n";
                        int locus = byte_index * EIGHT + bit_index;
                        //std::cout << "l2118 fits["<<locus<<"] = " << fits[locus] << " ";
                        W_T1haplo0 *= fits[locus];
                    }
                }
            }
            haplo0.setW_T1(W_T1haplo0, fitnessMapIndex);
            //std::cout << "W_T1haplo0 = " << W_T1haplo0 << std::endl;
        } else if (haplo1.getW_T1(fitnessMapIndex) == -1.0)
        {
            double W_T1haplo1 = 1.0;
            for (int byte_index = T1_ByteIndexFrom ; byte_index <= T1_ByteIndexTo ; byte_index++ )
            {
                int bit_index_to;
                int bit_index_from;
                if (byte_index == T1_ByteIndexTo)
                {
                    bit_index_to = T1_BitIndexTo;
                } else
                {
                    bit_index_to = 7;
                }
                if (byte_index == T1_ByteIndexFrom)
                {
                    bit_index_from = T1_BitIndexFrom;
                } else
                {
                    bit_index_from = 0;
                }
                for (int bit_index = bit_index_from ; bit_index <= bit_index_to ; bit_index++ )
                {
                    //std::cout << byte_index*8+bit_index << " ";
                    if (haplo1.getT1_Allele(byte_index,bit_index))
                    {
                        //std::cout << "Mut at haplo1 " << byte_index << " " << bit_index << " (" << byte_index * 8 + bit_index << ").\n";
                        int locus = byte_index * EIGHT + bit_index;
                        //std::cout << "l2144 fits["<<locus<<"] = " << fits[locus] << " ";
                        W_T1haplo1 *= fits[locus];
                    }
                }
            }
            haplo1.setW_T1(W_T1haplo1, fitnessMapIndex);
            //std::cout << "W_T1haplo1 = " << W_T1haplo1 << std::endl;
        }
    } else
    {
        if (haplo0.getW_T1(fitnessMapIndex) == -1.0 && haplo1.getW_T1(fitnessMapIndex) == -1.0)
        {
            double W_T1haplo0 = 1.0;
            double W_T1haplo1 = 1.0;
            for (auto& polymorphicLocus : SSP->simTracker.T1SitesForFitnessMultiplicityCalculation[fitnessMapIndex])
            {
                //assert(polymorphicLocus.locus >= T1_locusFrom && polymorphicLocus.locus < T1_locusTo);
                if (haplo1.getT1_Allele(polymorphicLocus.byte_index,polymorphicLocus.bit_index))
                {
                    W_T1haplo0 *= fits[polymorphicLocus.locus];
                }
                if (haplo1.getT1_Allele(polymorphicLocus.byte_index,polymorphicLocus.bit_index))
                {
                    W_T1haplo0 *= fits[polymorphicLocus.locus];
                }
            }
            haplo0.setW_T1(W_T1haplo0, fitnessMapIndex);
            haplo1.setW_T1(W_T1haplo1, fitnessMapIndex);

        } else if (haplo0.getW_T1(fitnessMapIndex) == -1.0)
        {
            double W_T1haplo0 = 1.0;
            for (auto& polymorphicLocus : SSP->simTracker.T1SitesForFitnessMultiplicityCalculation[fitnessMapIndex])
            {
                //assert(polymorphicLocus.locus >= T1_locusFrom && polymorphicLocus.locus < T1_locusTo);
                if (haplo1.getT1_Allele(polymorphicLocus.byte_index,polymorphicLocus.bit_index))
                {
                    W_T1haplo0 *= fits[polymorphicLocus.locus];
                }
            }
            haplo0.setW_T1(W_T1haplo0, fitnessMapIndex);

        } else if (haplo1.getW_T1(fitnessMapIndex) == -1.0)
        {
            double W_T1haplo1 = 1.0;
            for (auto& polymorphicLocus : SSP->simTracker.T1SitesForFitnessMultiplicityCalculation[fitnessMapIndex])
            {
                //assert(polymorphicLocus.locus >= T1_locusFrom && polymorphicLocus.locus < T1_locusTo);
                if (haplo1.getT1_Allele(polymorphicLocus.byte_index,polymorphicLocus.bit_index))
                {
                    W_T1haplo1 *= fits[polymorphicLocus.locus];
                }
            }
            haplo1.setW_T1(W_T1haplo1, fitnessMapIndex);
        }
    }*/

    if (haplo0.getW_T1(fitnessMapIndex) == -1.0 && haplo1.getW_T1(fitnessMapIndex) == -1.0)
    {
        double W_T1haplo0 = 1.0;
        double W_T1haplo1 = 1.0;
        for (auto& polymorphicLocus : SSP->simTracker.T1SitesForFitnessMultiplicityCalculation[fitnessMapIndex])
        {
            //assert(polymorphicLocus.locus >= T1_locusFrom && polymorphicLocus.locus < T1_locusTo);
            if (haplo0.getT1_Allele(polymorphicLocus.byte_index,polymorphicLocus.bit_index))
            {
                W_T1haplo0 *= fits[polymorphicLocus.locus];
            }
            if (haplo1.getT1_Allele(polymorphicLocus.byte_index,polymorphicLocus.bit_index))
            {
                W_T1haplo1 *= fits[polymorphicLocus.locus];
            }
        }
        haplo0.setW_T1(W_T1haplo0, fitnessMapIndex);
        haplo1.setW_T1(W_T1haplo1, fitnessMapIndex);

    } else if (haplo0.getW_T1(fitnessMapIndex) == -1.0)
    {
        double W_T1haplo0 = 1.0;
        for (auto& polymorphicLocus : SSP->simTracker.T1SitesForFitnessMultiplicityCalculation[fitnessMapIndex])
        {
            //assert(polymorphicLocus.locus >= T1_locusFrom && polymorphicLocus.locus < T1_locusTo);
            if (haplo0.getT1_Allele(polymorphicLocus.byte_index,polymorphicLocus.bit_index))
            {
                W_T1haplo0 *= fits[polymorphicLocus.locus];
            }
        }
        haplo0.setW_T1(W_T1haplo0, fitnessMapIndex);

    } else if (haplo1.getW_T1(fitnessMapIndex) == -1.0)
    {
        double W_T1haplo1 = 1.0;
        for (auto& polymorphicLocus : SSP->simTracker.T1SitesForFitnessMultiplicityCalculation[fitnessMapIndex])
        {
            //assert(polymorphicLocus.locus >= T1_locusFrom && polymorphicLocus.locus < T1_locusTo);
            if (haplo1.getT1_Allele(polymorphicLocus.byte_index,polymorphicLocus.bit_index))
            {
                W_T1haplo1 *= fits[polymorphicLocus.locus];
            }
        }
        haplo1.setW_T1(W_T1haplo1, fitnessMapIndex);
    }
    
        
    double r = haplo0.getW_T1(fitnessMapIndex) * haplo1.getW_T1(fitnessMapIndex);  

    
    /*
    if (r != 1.0)
    {
        std::cout << "line 2187 r = "<< r << " ";
        assert (SSP->T1_isMultiplicitySelection);
        for (int locus = 0 ; locus < SSP->T1_nbBits ; locus++)
        {
            int byte_index = locus / EIGHT;
            int bit_index = locus % EIGHT;
            if (haplo0.getT1_Allele(byte_index,bit_index))
            {
                std::cout << "haplo0 : ";
                std::cout << SSP->T1_FitnessEffects[Habitat][locus];
                std::cout << " ";
            }
            
            if (haplo1.getT1_Allele(byte_index,bit_index))
            {
                std::cout << "haplo1 : ";
                std::cout << SSP->T1_FitnessEffects[Habitat][locus];
                std::cout << " ";
            }
        }
    }
    */
    //if (r!=1.0)  std::cout << "T1 fit: = " << r << " ";

    /*bool FFM = this->isFreeFromMutations(T1_locusFrom, T1_locusTo); // be careful, you might carry neutral mutations so this is a bad test
    if (r == 1.0 && !FFM)
    {
        std::cout<<"Prob 1\n";
        //abort();
    } else if (r != 1.0 && FFM)
    {
        std::cout<<"Prob 2\n";
        //abort();
    }*/

    assert(r<=1.0 && r>=0.0);

    return r;
}

double Individual::CalculateT1FitnessNoMultiplicity(const int& Habitat)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'CalculateT1FitnessNoMultiplicity'\n";
#endif   
    auto& fits = SSP->T1_FitnessEffects[Habitat];
    assert(fits.size() == SSP->T1_nbBits * THREE);


    double W_T1_WholeIndividual = 1.0;

    for (auto& polymorphicLocus : SSP->simTracker.T1SitesForFitnessNoMultiplicityCalculation)
    {
        W_T1_WholeIndividual *= 
            fits[
                THREE * polymorphicLocus.locus + 
                haplo0.getT1_Allele(polymorphicLocus.byte_index,polymorphicLocus.bit_index) +
                haplo1.getT1_Allele(polymorphicLocus.byte_index,polymorphicLocus.bit_index)
            ];
    }

    /*if (SSP->simTracker.shouldLoopOverAllSites)
    {
        for (int byte_index = 0 ; byte_index < SSP->T1_nbChars ; byte_index++)
        {
            int bit_index_to;
            if (byte_index == SSP->T1_nbChars - 1)
            {
                bit_index_to = SSP->T1_nbBitsLastByte;
            } else
            {
                bit_index_to = EIGHT;
            }
            for (int bit_index=0 ; bit_index < bit_index_to ; bit_index++ )
            {
                int locus = byte_index * EIGHT + bit_index;
                //std::cout << "SSP->T1_FitnessEffects.size() = "<< SSP->T1_FitnessEffects.size() << std::endl;
                //std::cout << "SSP->T1_FitnessEffects[Habitat].size() = "<< SSP->T1_FitnessEffects[Habitat].size() << std::endl;
                //assert(SSP->T1_FitnessEffects[Habitat].size() > THREE * locus + haplo0.getT1_Allele(byte_index,bit_index) + haplo1.getT1_Allele(byte_index,bit_index));
                //std::cout << SSP->T1_FitnessEffects[Habitat][THREE * locus + haplo0.getT1_Allele(byte_index,bit_index) + haplo1.getT1_Allele(byte_index,bit_index)] << std::endl;
                W_T1_WholeIndividual *= fits[THREE * locus + haplo0.getT1_Allele(byte_index,bit_index) + haplo1.getT1_Allele(byte_index,bit_index)];
                //std::cout << "bit_index = " << bit_index << " | byte_index = " << byte_index << " | fit Value = " << SSP->T1_FitnessEffects[Habitat][THREE * locus + haplo0.getT1_Allele(byte_index,bit_index) + haplo1.getT1_Allele(byte_index,bit_index)] << "\n";
            }
        }
    } else
    {
        for (auto& polymorphicLocus : SSP->simTracker.T1SitesForFitnessNoMultiplicityCalculation)
        {
            W_T1_WholeIndividual *= 
                fits[
                    THREE * polymorphicLocus.locus + 
                    haplo0.getT1_Allele(polymorphicLocus.byte_index,polymorphicLocus.bit_index) +
                    haplo1.getT1_Allele(polymorphicLocus.byte_index,polymorphicLocus.bit_index)
                ];
        }
    }*/

    
    return W_T1_WholeIndividual;
    //std::cout << "W_T1_WholeIndividual = " << W_T1_WholeIndividual << std::endl;
}

double Individual::CalculateT1EpistaticFitness(const int& Habitat)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'CalculateT1EpistaticFitness'\n";
#endif
    double r = 1.0;

    assert(SSP->T1_Epistasis_FitnessEffects.size() > Habitat);
    for (int groupOfLociIndex = 0 ; groupOfLociIndex < SSP->T1_Epistasis_FitnessEffects[Habitat].size() ; groupOfLociIndex++)
    {
        auto& fits = SSP->T1_Epistasis_FitnessEffects[Habitat][groupOfLociIndex];
        auto& loci = SSP->T1_Epistasis_LociIndices[Habitat][groupOfLociIndex];
        int FitnessValueIndex = 0;
        int nbLociUnderEpistasis = loci.size();
        for (int i = 0 ; i < nbLociUnderEpistasis ; i++)
        {
            /*
            Organization of SSP->T1_Epistasis_FitnessEffects[Habitat] example with 3 loci where the first line is the locus SSP->T1_Epistasis_LociIndices[Habitat][0].
            locus 1:    00 00 00 00 00 00 00 00 00 01 01 01 01 01 01 01 01 01 11 11 11 11 11 11 11 11 11
            locus 2:    00 00 00 01 01 01 11 11 11 00 00 00 01 01 01 11 11 11 00 00 00 01 01 01 11 11 11
            locus 3:    00 01 11 00 01 11 00 01 11 00 01 11 00 01 11 00 01 11 00 01 11 00 01 11 00 01 11
            Index:      0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 

            Tiny R code to create this kind of sequence
                cat(paste0(
                    rep(paste0("A",rep(c("0|0", "0|1","1|1"), each=3)), each=3),
                    "_",
                    rep(paste0("B",rep(c("0|0", "0|1","1|1"), 3)), each=3),
                    "_",
                    paste0("B",rep(c("0|0", "0|1","1|1"), 9))
                ))
            */
            
            T1_locusDescription& T1_locus = loci[i];
            FitnessValueIndex += (haplo0.getT1_Allele(T1_locus.byte_index, T1_locus.bit_index) + haplo1.getT1_Allele(T1_locus.byte_index, T1_locus.bit_index)) * pow(THREE,nbLociUnderEpistasis - i - 1);
        }
        assert(FitnessValueIndex >= 0 && FitnessValueIndex < pow(THREE,nbLociUnderEpistasis));
        assert(fits[FitnessValueIndex] >= 0.0 && fits[FitnessValueIndex] <= 1.0);
        //std::cout << "SSP->T1_Epistasis_FitnessEffects[Habitat][FitnessValueIndex] = " << SSP->T1_Epistasis_FitnessEffects[Habitat][FitnessValueIndex] << std::endl;
        r *= fits[FitnessValueIndex];
    }
        
    return r;
}

double Individual::CalculateT2Fitness(const int& Habitat, int fitnessMapIndex, int T2_locusFrom, int T2_locusTo)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'CalculateT2Fitness'\n";
#endif          
    assert(SSP->T2_FitnessEffects.size() > Habitat);
    auto& fits = SSP->T2_FitnessEffects[Habitat];
    assert(fits.size() == SSP->T2_nbChars);

    if (haplo0.getW_T2(fitnessMapIndex) == -1.0 && haplo1.getW_T2(fitnessMapIndex) == -1.0)
    {
        double W_T2haplo0 = 1.0;
        double W_T2haplo1 = 1.0;
        for (int byte_index = T2_locusFrom ; byte_index < T2_locusTo ; byte_index++ )
        {
            W_T2haplo0 *= pow(fits[byte_index], haplo0.getT2_Allele(byte_index));
            W_T2haplo1 *= pow(fits[byte_index], haplo1.getT2_Allele(byte_index));
        }
        haplo0.setW_T2(W_T2haplo0, fitnessMapIndex);
        haplo1.setW_T2(W_T2haplo1, fitnessMapIndex);
        //if (W_T2haplo0!=1.0) std::cout << "W_T2haplo0 = " << W_T2haplo0 << std::endl;
        //if (W_T2haplo1!=1.0) std::cout << "W_T2haplo1 = " << W_T2haplo1 << std::endl;
    } else
    {
        if (haplo0.getW_T2(fitnessMapIndex) == -1.0)
        {
            double W_T2haplo0 = 1.0;
            for (int byte_index = T2_locusFrom ; byte_index < T2_locusTo ; byte_index++ )
            {
                //assert(SSP->T2_FitnessEffects.size() > Habitat);
                //assert(SSP->T2_FitnessEffects[Habitat].size() > char_index);
                W_T2haplo0 *= pow(fits[byte_index], haplo0.getT2_Allele(byte_index));
            }
            haplo0.setW_T2(W_T2haplo0, fitnessMapIndex);
            //if (W_T2haplo0!=1.0) std::cout << "W_T2haplo0 = " << W_T2haplo0 << std::endl;
        }
        if (haplo1.getW_T2(fitnessMapIndex) == -1.0)
        {
            double W_T2haplo1 = 1.0;
            for (int byte_index = T2_locusFrom ; byte_index < T2_locusTo ; byte_index++ )
            {
                //assert(SSP->T2_FitnessEffects.size() > Habitat);
                //assert(SSP->T2_FitnessEffects[Habitat].size() > char_index);
                W_T2haplo1 *= pow(fits[byte_index], haplo1.getT2_Allele(byte_index)); 
            }
            haplo1.setW_T2(W_T2haplo1, fitnessMapIndex);
            //if (W_T2haplo1!=1.0) std::cout << "W_T2haplo1 = " << W_T2haplo1 << std::endl;
        }
    }
    double r = haplo0.getW_T2(fitnessMapIndex) * haplo1.getW_T2(fitnessMapIndex);


    //if (r!=1.0)  std::cout << "T2 fit: = " << r << " ";
    assert(r<=1.0 && r>=0.0);
    return r;
}

void Individual::CalculateT3Phenotype(const int& Habitat)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'CalculateT3Phenotype'\n";
#endif      
    // Directly write the phenotype in the stat 'T3_IndPhenotype'. 
    // Habitat is used here to allow plasticity. One must not confound plasticity on phenotypes and local selection!

    // set T3_IndPhenotype zero for all dimensions
    T3_IndPhenotype.resize(SSP->T3_PhenoNbDimensions);
    for (auto& DimPheno : T3_IndPhenotype)
    {
        DimPheno = 0.0;
    }

    // Calculate phenotype
    for (int byte_index = 0 ; byte_index < SSP->T3_nbChars; byte_index++)
    {
        double allele = (double) haplo0.getT3_Allele(byte_index);
        for (int dim = 0 ; dim < SSP->T3_PhenoNbDimensions; dim++)
        {
            std::normal_distribution<double> dist(0.0,SSP->T3_DevelopmentalNoiseStandardDeviation[Habitat][dim]);
            T3_IndPhenotype[dim] += SSP->T3_PhenotypicEffects[Habitat][byte_index * SSP->T3_PhenoNbDimensions + dim] * allele + dist(GP->mt);
        }
    }
}

double Individual::CalculateT3Fitness(const int& Habitat) // is a static function. The static specification must appear only in the class declaration
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'CalculateT3Fitness'\n";
#endif  
    // Calculate fitness of phenotype describe in the static 'T3_IndPhenotype'
    assert(Individual::T3_IndPhenotype.size() == SSP->T3_PhenoNbDimensions);

    // Calculate Fitness
    double W = 1.0;
    for (int dim = 0 ; dim < SSP->T3_PhenoNbDimensions; dim++)
    {
        double diffToOptimal = Individual::T3_IndPhenotype[dim] - SSP->T3_fitnessLandscapeOptimum[Habitat][dim];
        if (SSP->T3_fitnessLandscapeType == 'L')
        {
            W *= 1 - (std::abs(diffToOptimal) * SSP->T3_fitnessLandscapeLinearGradient[Habitat][dim]);
            if (W < 0.0) 
            {
                W = 0.0;
                break;
            }
        } else if (SSP->T3_fitnessLandscapeType == 'G')
        {
            W *= exp( - pow(diffToOptimal,2) / SSP->T3_fitnessLandscapeGaussStrength[Habitat][dim]);
        } else
        {
            std::cout << "Internal error in Individual::CalculateT3Fitness. Unkown SSP->T3_fitnessLandscapeType.  SSP->T3_fitnessLandscapeType = " << SSP->T3_fitnessLandscapeType << "\n";
        }
    }
    assert(W <= 1.0 && W >= 0.0);
    return W;
}



std::vector<double> Individual::CalculateFitnessComponents(const int& Habitat)
{
    double rT1 = 1.0;
    double rT2 = 1.0;
    double rT1epistasis = 1.0;
    double rT3 = 1.0;
    double rT5 = 1.0;

    // Trait 1 if multiplicity (excluding epistasis) and Trait 2
    if ((SSP->T1_isMultiplicitySelection && SSP->T1_isSelection) || (SSP->T5_isMultiplicitySelection && SSP->T5_isSelection) || SSP->T2_isSelection)
    {
        int T1_locusFrom = 0; // from included
        int T2_locusFrom = 0; // from included
        int T5_locusFrom = 0; // from included
        int T5_haplo0From = 0;
        int T5_haplo1From = 0;
        
        for (int fitnessMapIndex = 0 ; fitnessMapIndex < SSP->NbElementsInFitnessMap; fitnessMapIndex++)
        {
            //std::cout << "fitnessMapIndex = " << fitnessMapIndex << "\n";
            
            // Trait 1
            if (SSP->T1_isMultiplicitySelection && SSP->T1_isSelection)
            {
                int T1_locusTo = SSP->FromLocusToFitnessMapBoundaries[fitnessMapIndex].T1; // to excluded
                assert(T1_locusTo >= T1_locusFrom);

                //std::cout << "calculate T1 Fitness from " << T1_locusFrom << " to " << T1_locusTo << "\n";

                rT1 *= this->CalculateT1FitnessMultiplicity(Habitat, fitnessMapIndex, T1_locusFrom, T1_locusTo);

                // set from Locus
                T1_locusFrom = T1_locusTo;
            }

            // Trait 2
            if (SSP->T2_isSelection)
            {
                int T2_locusTo = SSP->FromLocusToFitnessMapBoundaries[fitnessMapIndex].T2; // to excluded
                assert(T2_locusTo >= T2_locusFrom);

                rT2 *= this->CalculateT2Fitness(Habitat, fitnessMapIndex, T2_locusFrom, T2_locusTo);

                // set from Locus
                T2_locusFrom = T2_locusTo;
            }

            // Trait 5
            if (SSP->T5_isMultiplicitySelection && SSP->T5_isSelection)
            {
                int T5_locusTo = SSP->FromLocusToFitnessMapBoundaries[fitnessMapIndex].T5; // to excluded
                assert(T5_locusTo >= T5_locusFrom);

                auto T5_haplo0To = haplo0.T5_AllelesPosition(T5_locusTo, T5_haplo0From);
                auto T5_haplo1To = haplo0.T5_AllelesPosition(T5_locusTo, T5_haplo1From);


                //std::cout << "calculate T5 Fitness from " << T5_locusFrom << " to " << T5_locusTo << "\n";

                rT5 *= this->CalculateT5FitnessMultiplicity(Habitat, fitnessMapIndex, T5_haplo0From, T5_haplo0To, T5_haplo1From, T5_haplo1To); // Will directly set the two "from" iterators to the "to" iterators for the next fitness block

                // set from Locus
                T5_locusFrom = T5_locusTo;
            }
        }

    }
      
    // Trait 1 (excluding epistasis) 
    if (!SSP->T1_isMultiplicitySelection && SSP->T1_isSelection)
    {
        assert(rT1 == 1.0);
        rT1 = this->CalculateT1FitnessNoMultiplicity(Habitat);
    }


    // Trait 1 epistasis
    if (SSP->T1_isEpistasis)
    {
        rT1epistasis = this->CalculateT1EpistaticFitness(Habitat);
    }

    // Trait 3
    if (SSP->T3_isSelection)
    {
        this->CalculateT3Phenotype(Habitat); // save phenotype in the static attribute 'T3_IndPhenotype'
        rT3 = Individual::CalculateT3Fitness(Habitat); // calculate fitness associated to teh phenotype saved in the static attribute 'T3_IndPhenotype'
    }

    // Trait 5 
    if (!SSP->T5_isMultiplicitySelection && SSP->T5_isSelection)
    {
        assert(rT5 == 1.0);
        rT5 = this->CalculateT5FitnessNoMultiplicity(Habitat);
    }

    return {rT1, rT2, rT1epistasis, rT3, rT5};
}

double Individual::CalculateFitness(const int& patch_index)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'CalculateFitness'\n";
#endif  

    // Figure out what Habitat it corresponds to
    const int& Habitat = SSP->Habitats[patch_index];

    // Get fitness components
    auto fitnessComponents = CalculateFitnessComponents(Habitat);

    // Get total fitness
    double r = fitnessComponents[0] * fitnessComponents[1] * fitnessComponents[2] * fitnessComponents[3] * fitnessComponents[4];
    
    assert(r >= 0 && r <= 1);  

    return  r;
}

void Individual::PrintBinaryFile(OutputFile& file)
{
    this->haplo0.PrintBinaryFile(file);
    this->haplo1.PrintBinaryFile(file);
}

void Individual::SetHaplo(int haplo_index, Haplotype& chrom)
{
    if (haplo_index)
    {
        haplo1 = chrom;
    } else
    {
        haplo0 = chrom;
    }
}

Individual::Individual(bool ShouldReadPopFromBinary)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'Individual::Individual(bool ShouldReadPopFromBinary)'\n";
#endif  
    Haplotype c0(ShouldReadPopFromBinary);
    Haplotype c1(ShouldReadPopFromBinary);
    haplo0 = c0;
    haplo1 = c1;
}

Individual::Individual(const int patch_index, char Abiogenesis, int ind_index)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'Individual::Individual(const int patch_index, char Abiogenesis)'\n";
#endif
    Haplotype c0(patch_index,Abiogenesis, (double) (ind_index*2));
    Haplotype c1(patch_index,Abiogenesis, (double) (ind_index*2+1));
    haplo0 = c0;
    haplo1 = c1;
}

Individual::Individual(Haplotype& knownHaplotype, char Abiogenesis)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'Individual::Individual(Haplotype knownHaplotype, char Abiogenesis)'\n";
#endif  
    (void) Abiogenesis;
    haplo0 = knownHaplotype;
    haplo1 = knownHaplotype;
}

Individual::Individual(const Haplotype& knownHaplotype)
{
    haplo0 = knownHaplotype;
    haplo1 = knownHaplotype;   
}














///////////// Fitness On Subset Of Loci

bool Individual::isLocusIsInSet(const int locus, const std::vector<int>& LociSet)
{
    return std::find(LociSet.begin(), LociSet.end(), locus) != LociSet.end();
}

std::vector<double> Individual::CalculateFitnessComponentsOnSubsetOfLoci(const int& Habitat, const int lociSetIndex)
{
    double rT1 = 1.0;
    double rT2 = 1.0;
    double rT1epistasis = 1.0;
    double rT3 = 1.0;
    double rT5 = 1.0;

    // Trait 1 if multiplicity (excluding epistasis) and Trait 2
    if ((SSP->T1_isMultiplicitySelection && SSP->T1_isSelection) || (SSP->T5_isMultiplicitySelection && SSP->T5_isSelection) || SSP->T2_isSelection)
    {    
        // Trait 1
        if (SSP->T1_isMultiplicitySelection && SSP->T1_isSelection)
        {
            rT1 *= this->CalculateT1FitnessMultiplicityOnSubsetOfLoci(
                Habitat,
                SSP->subsetT1LociForfitnessSubsetLoci_file[lociSetIndex]
            );
        }

        // Trait 2
        if (SSP->T2_isSelection)
        {
            rT2 *= this->CalculateT2FitnessOnSubsetOfLoci(
                Habitat,
                SSP->subsetT2LociForfitnessSubsetLoci_file[lociSetIndex]
            );
        }

        // Trait 5
        if (SSP->T5_isSelection)
        {
            rT5 *= this->CalculateT5FitnessMultiplicityOnSubsetOfLoci(
                Habitat,
                SSP->subsetT5LociForfitnessSubsetLoci_file[lociSetIndex]
            );
        }
    }
      
    // Trait 1 no multiplicity (excluding epistasis) 
    if (!SSP->T1_isMultiplicitySelection && SSP->T1_isSelection)
    {
        assert(rT1 == 1.0);
        rT1 = this->CalculateT1FitnessNoMultiplicityOnSubsetOfLoci(
            Habitat,
            SSP->subsetT1LociForfitnessSubsetLoci_file[lociSetIndex]
        );
    }


    // Trait 1 epistasis
    if (SSP->T1_isEpistasis)
    {
        rT1epistasis = this->CalculateT1EpistaticFitnessOnSubsetOfLoci(
            Habitat,
            SSP->subsetT1epistasisLociForfitnessSubsetLoci_file[lociSetIndex]
        );
    }

    // Trait 3
    if (SSP->T3_isSelection)
    {
        this->CalculateT3PhenotypeOnSubsetOfLoci(
            Habitat,
            SSP->subsetT3LociForfitnessSubsetLoci_file[lociSetIndex]
        ); // save phenotype in the static attribute 'T3_IndPhenotype'
        rT3 = Individual::CalculateT3Fitness(Habitat); // calculate fitness associated to teh phenotype saved in the static attribute 'T3_IndPhenotype'
    }
    //std::cout << "rT1 = "<<rT1<< " rT2 = " <<rT2<<" rT1epistasis = "<<rT1epistasis<< " rT3 = " <<rT3<<"\n";

    if (SSP->T5_isMultiplicitySelection && SSP->T5_isSelection)
    {
        assert(rT5 == 1.0);
        rT5 = this->CalculateT1FitnessNoMultiplicityOnSubsetOfLoci(
            Habitat,
            SSP->subsetT5LociForfitnessSubsetLoci_file[lociSetIndex]
        );
    }

    return {rT1, rT2, rT1epistasis, rT3, rT5};
}




double Individual::CalculateT1FitnessMultiplicityOnSubsetOfLoci(const int& Habitat, const std::vector<int>& LociSet)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'CalculateT1FitnessMultiplicityOnSubsetOfLoci'\n";
#endif   

    assert(SSP->T1_isMultiplicitySelection);

    // Get the right fitness array
    auto& fits = SSP->T1_FitnessEffects[Habitat];
    assert(fits.size() == SSP->T1_nbBits);


    double r = 1.0;
    for (const int& locus : LociSet)
    {
        int byte_index = locus / 8;
        int bit_index = locus % 8;

        //assert(polymorphicLocus.locus >= T1_locusFrom && polymorphicLocus.locus < T1_locusTo);
        if (haplo0.getT1_Allele(byte_index,bit_index))
        {
            r *= fits[locus];
        }
        if (haplo1.getT1_Allele(byte_index,bit_index))
        {
            r *= fits[locus];
        }
    }

    assert(r<=1.0 && r>=0.0);

    return r;
}

double Individual::CalculateT1FitnessNoMultiplicityOnSubsetOfLoci(const int& Habitat, const std::vector<int>& LociSet)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'CalculateT1FitnessNoMultiplicityOnSubsetOfLoci'\n";
#endif   
    auto& fits = SSP->T1_FitnessEffects[Habitat];
    assert(fits.size() == SSP->T1_nbBits * THREE);


    double W_T1_WholeIndividual = 1.0;

    for (const int locus : LociSet)
    {
        int byte_index = locus / 8;
        int bit_index = locus % 8;
        //std::cout << "locus = " << locus << "\n";
        assert(locus >= 0 && locus < SSP->T1_nbBits);

        W_T1_WholeIndividual *= 
            fits[
                THREE * locus + 
                haplo0.getT1_Allele(byte_index,bit_index) +
                haplo1.getT1_Allele(byte_index,bit_index)
            ];
    }
    //std::cout << "W_T1_WholeIndividual = " << W_T1_WholeIndividual << "\n";
    assert(W_T1_WholeIndividual >= 0.0 && W_T1_WholeIndividual <= 1.0);
    

    return W_T1_WholeIndividual;
}

double Individual::CalculateT1EpistaticFitnessOnSubsetOfLoci(const int& Habitat, const std::vector<int>& LociSet)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'CalculateT1EpistaticFitnessOnSubsetOfLoci'\n";
#endif
    double r = 1.0;

    assert(SSP->T1_Epistasis_FitnessEffects.size() > Habitat);
    for (int groupOfLociIndex = 0 ; groupOfLociIndex < SSP->T1_Epistasis_FitnessEffects[Habitat].size() ; groupOfLociIndex++)
    {
        auto& fits = SSP->T1_Epistasis_FitnessEffects[Habitat][groupOfLociIndex];
        auto& loci = SSP->T1_Epistasis_LociIndices[Habitat][groupOfLociIndex];

        // check if all are in LociSet, if some are but not all or if none are. This process could be done only once in SSP for better perfromance

        class comparatorsForSetIntersection
        {
        public:
            bool operator()(int lhs, T1_locusDescription rhs) {return lhs + rhs.locus;}
            bool operator()(T1_locusDescription lhs, int rhs) {return lhs.locus + rhs;}
        };
        comparatorsForSetIntersection comps;

        std::vector<T1_locusDescription> intersection;
        std::set_intersection(
            loci.begin(), loci.end(),
            LociSet.begin(), LociSet.end(),
            std::back_inserter(intersection),
            comps
        );

        if (intersection.size() == LociSet.size() && LociSet.size() > 0)
        {
            // Then all loci of 'LociSet' are in 'loci'
            // -> compute fitness
            int FitnessValueIndex = 0;
            int nbLociUnderEpistasis = loci.size();
            for (int i = 0 ; i < nbLociUnderEpistasis ; i++)
            {
                /*
                Organization of SSP->T1_Epistasis_FitnessEffects[Habitat] example with 3 loci where the first line is the locus SSP->T1_Epistasis_LociIndices[Habitat][0].
                locus 1:    00 00 00 00 00 00 00 00 00 01 01 01 01 01 01 01 01 01 11 11 11 11 11 11 11 11 11
                locus 2:    00 00 00 01 01 01 11 11 11 00 00 00 01 01 01 11 11 11 00 00 00 01 01 01 11 11 11
                locus 3:    00 01 11 00 01 11 00 01 11 00 01 11 00 01 11 00 01 11 00 01 11 00 01 11 00 01 11
                Index:      0  1  2  3  4  5  6  7  8  9  10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 

                Tiny R code to create this kind of sequence
                    cat(paste0(
                        rep(paste0("A",rep(c("0|0", "0|1","1|1"), each=3)), each=3),
                        "_",
                        rep(paste0("B",rep(c("0|0", "0|1","1|1"), 3)), each=3),
                        "_",
                        paste0("B",rep(c("0|0", "0|1","1|1"), 9))
                    ))
                */
                
                T1_locusDescription& T1_locus = loci[i];
                FitnessValueIndex += (haplo0.getT1_Allele(T1_locus.byte_index, T1_locus.bit_index) + haplo1.getT1_Allele(T1_locus.byte_index, T1_locus.bit_index)) * pow(THREE,nbLociUnderEpistasis - i - 1);
            }
            assert(FitnessValueIndex >= 0 && FitnessValueIndex < pow(THREE,nbLociUnderEpistasis));
            assert(fits[FitnessValueIndex] >= 0.0 && fits[FitnessValueIndex] <= 1.0);
            //std::cout << "SSP->T1_Epistasis_FitnessEffects[Habitat][FitnessValueIndex] = " << SSP->T1_Epistasis_FitnessEffects[Habitat][FitnessValueIndex] << std::endl;
            r *= fits[FitnessValueIndex];

        } else if (intersection.size() > 0)
        {
            // Then some but not all loci of 'LociSet' are in 'loci'
            return std::numeric_limits<double>::quiet_NaN();
        } else if (intersection.size() == 0 && LociSet.size() > 0)
        {
            // Then no elements of 'LociSet' are in 'loci'
            // -> Nothing to do (or only 'r *= 1.0;')
        } else
        {
            // Nothing in LociSet
            // -> Nothing to do (or only 'r *= 1.0;')
            assert(LociSet.size() == 0);
        }
    }
        
    return r;
}

double Individual::CalculateT2FitnessOnSubsetOfLoci(const int& Habitat, const std::vector<int>& LociSet)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'CalculateT2FitnessOnSubsetOfLoci'\n";
#endif          
    assert(SSP->T2_FitnessEffects.size() > Habitat);
    auto& fits = SSP->T2_FitnessEffects[Habitat];
    assert(fits.size() == SSP->T2_nbChars);

    double r = 1.0;
    for (const int locus : LociSet)
    {
        r *= pow(fits[locus], haplo0.getT2_Allele(locus));
        r *= pow(fits[locus], haplo1.getT2_Allele(locus));
    }

    //if (r!=1.0)  std::cout << "T2 fit: = " << r << " ";
    assert(r<=1.0 && r>=0.0);
    return r;
}

void Individual::CalculateT3PhenotypeOnSubsetOfLoci(const int& Habitat, const std::vector<int>& LociSet)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'CalculateT3PhenotypeOnSubsetOfLoci'\n";
#endif      
    // Directly write the phenotype in the stat 'T3_IndPhenotype'. 
    // Habitat is used here to allow plasticity. One must not confound plasticity on phenotypes and local selection!

    // set T3_IndPhenotype zero for all dimensions
    T3_IndPhenotype.resize(SSP->T3_PhenoNbDimensions);
    for (auto& DimPheno : T3_IndPhenotype)
    {
        DimPheno = 0.0;
    }

    // Calculate phenotype
    for (const int locus : LociSet)
    {
        double allele = (double) haplo0.getT3_Allele(locus);
        for (int dim = 0 ; dim < SSP->T3_PhenoNbDimensions; dim++)
        {
            std::normal_distribution<double> dist(0.0,SSP->T3_DevelopmentalNoiseStandardDeviation[Habitat][dim]);
            T3_IndPhenotype[dim] += SSP->T3_PhenotypicEffects[Habitat][locus * SSP->T3_PhenoNbDimensions + dim] * allele + dist(GP->mt);
        }
    }
}



double Individual::CalculateT5FitnessMultiplicity(
    const int& Habitat,
    int fitnessMapIndex,
    int& haplo0From,
    int& haplo0To,
    int& haplo1From,
    int& haplo1To
)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'CalculateT5FitnessMultiplicity'\n";
#endif   

    // Get the right fitness array
    auto& fits = SSP->T5_FitnessEffects[Habitat];
    assert(fits.size() == SSP->T5_nbBits);


    if (haplo0.getW_T5(fitnessMapIndex) == -1.0)
    {  
        //std::cout << "Computing W_T5haplo: mutations at ";

        double W_T5haplo0 = 1.0;
        for ( ; haplo0From < haplo0To ; haplo0From++)
        {
            W_T5haplo0 *= fits[haplo0From];
            //std::cout << *it << " ";
        }
        haplo0.setW_T5(W_T5haplo0, fitnessMapIndex);
        //std::cout << "\n";
    }

    if (haplo1.getW_T5(fitnessMapIndex) == -1.0)
    {

        //std::cout << "Computing W_T5hapl1: mutations at ";
        double W_T5haplo1 = 1.0;
        for ( ; haplo1From < haplo1To ; haplo1From++)
        {
            W_T5haplo1 *= fits[haplo1From];
            //std::cout << *it << " ";
        }
        //std::cout << "\n";
        haplo1.setW_T5(W_T5haplo1, fitnessMapIndex);
    }
        
    double r = haplo0.getW_T5(fitnessMapIndex) * haplo1.getW_T5(fitnessMapIndex);  

    assert(r<=1.0 && r>=0.0);

    //std::cout << "r = " << r << "\n";

    return r;
}


double Individual::CalculateT5FitnessNoMultiplicity(const int& Habitat)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'CalculateT5FitnessNoMultiplicity'\n";
#endif   
    auto& fits = SSP->T5_FitnessEffects[Habitat];
    assert(fits.size() == SSP->T5_nbBits * 2);

    double W_T5_WholeIndividual = 1.0;

    auto itHaplo0 = haplo0.T5_AllelesCBegin();
    auto itHaplo1 = haplo1.T5_AllelesCBegin();

    auto itEndHaplo0  = haplo0.T5_AllelesCEnd();
    auto itEndHaplo1  = haplo1.T5_AllelesCEnd();

    /*
    if (itHaplo0 != itEndHaplo0 || itHaplo1 != itEndHaplo1)
    {
        std::cout << "haplo0: ";
        for (auto it = haplo0.T5_AllelesCBegin() ; it != itEndHaplo0 ; it++)
            std::cout << *it << " ";
        std::cout << "\nhaplo1: ";
        for (auto it = haplo1.T5_AllelesCBegin() ; it != itEndHaplo1 ; it++)
            std::cout << *it << " ";
        std::cout << "\n";
    }
    */
    

    bool haplo0ReachedEnd = itHaplo0 == itEndHaplo0;
    bool haplo1ReachedEnd = itHaplo1 == itEndHaplo1;

    while (!haplo0ReachedEnd || !haplo1ReachedEnd )
    {
        //std::cout << "haplo0ReachedEnd = " << std::boolalpha <<  haplo0ReachedEnd << " haplo1ReachedEnd = " << std::boolalpha << haplo1ReachedEnd << "\n";
        if (!haplo0ReachedEnd && !haplo1ReachedEnd)
        {
            //std::cout << "*itHaplo0 = " << *itHaplo0 << " *itHaplo1 = " << *itHaplo1 << "\n";
            // No haplo has reached end of mutations
            if (*itHaplo0 == *itHaplo1)
            {
                // homozygote mutant
                W_T5_WholeIndividual *= fits[*itHaplo0 + 1]; // 1-s
                itHaplo0++;
                itHaplo1++;
                if (itHaplo0 == itEndHaplo0)
                    haplo1ReachedEnd = true;
                if (itHaplo1 == itEndHaplo1)
                    haplo1ReachedEnd = true;
            } else if (*itHaplo0 > *itHaplo1)
            {
                // must wait for haplo1 to eventually find a homozygous
                while (*itHaplo0 > *itHaplo1)
                {
                    W_T5_WholeIndividual *= fits[*itHaplo1]; // 1-hs
                    itHaplo1++;
                    if (itHaplo1 == itEndHaplo1)
                    {
                        haplo1ReachedEnd = true;
                        break;
                    }
                }
            } else // if (*itHaplo0 < *itHaplo1)
            {
                // must wait for haplo0 to eventually find a homozygous
                while (*itHaplo0 < *itHaplo1)
                {
                    W_T5_WholeIndividual *= fits[*itHaplo0]; // 1-hs
                    itHaplo0++;
                    if (itHaplo0 == itEndHaplo0)
                    {
                        haplo0ReachedEnd = true;
                        break;
                    }
                }
            }
        } else
        {
            // One of the two haplos has reached end of mutations
            if (haplo0ReachedEnd)
            {
                //std::cout << "out *itHaplo1 = " << *itHaplo1 << "\n";
                while (itHaplo1 != itEndHaplo1)
                {
                    //std::cout << "out *itHaplo1 = " << *itHaplo1 << "\n";
                    W_T5_WholeIndividual *= fits[*itHaplo1]; // 1-hs
                    itHaplo1++;
                }
                break;
            } else
            {
                //std::cout << "out *itHaplo0 = " << *itHaplo0 << "\n";
                while (itHaplo0 != itEndHaplo0)
                {
                    //std::cout << "*itHaplo0 = " << *itHaplo0 << "\n";
                    W_T5_WholeIndividual *= fits[*itHaplo0]; // 1-hs
                    itHaplo0++;
                }
                break;
            }
        }     
    }


    return W_T5_WholeIndividual;
}

double Individual::CalculateT5FitnessMultiplicityOnSubsetOfLoci(const int& Habitat, const std::vector<int>& LociSet)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'CalculateT5FitnessMultiplicityOnSubsetOfLoci'\n";
#endif   

    assert(SSP->T5_isMultiplicitySelection);

    // Get the right fitness array
    auto& fits = SSP->T5_FitnessEffects[Habitat];
    assert(fits.size() == SSP->T5_nbBits);

    auto itHaplo0 = haplo0.T5_AllelesCBegin();
    auto itHaplo1 = haplo1.T5_AllelesCBegin();

    double r = 1.0;


    for (auto& locus : LociSet)
    {
        itHaplo0 = haplo0.T5_AllelesCiterator(locus, itHaplo0);
        itHaplo1 = haplo1.T5_AllelesCiterator(locus, itHaplo1);

        bool isHaplo0Mutant = itHaplo0 == haplo0.T5_AllelesCEnd() ? false : *itHaplo0 == locus ;
        bool isHaplo1Mutant = itHaplo1 == haplo1.T5_AllelesCEnd() ? false : *itHaplo1 == locus ;

        if (isHaplo0Mutant)
        {
            r *= fits[locus];
        }
        if (isHaplo1Mutant)
        {
            r *= fits[locus];
        }
    }

    return r;
}


double Individual::CalculateT5FitnessNoMultiplicityOnSubsetOfLoci(const int& Habitat, const std::vector<int>& LociSet)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'CalculateT5FitnessNoMultiplicity'\n";
#endif   
    auto& fits = SSP->T5_FitnessEffects[Habitat];
    assert(fits.size() == SSP->T5_nbBits * 2);

    double r = 1.0;

    auto itHaplo0 = haplo0.T5_AllelesCBegin();
    auto itHaplo1 = haplo1.T5_AllelesCBegin();

    for (auto& locus : LociSet)
    {
        itHaplo0 = haplo0.T5_AllelesCiterator(locus, itHaplo0);
        itHaplo1 = haplo1.T5_AllelesCiterator(locus, itHaplo1);

        bool isHaplo0Mutant = itHaplo0 == haplo0.T5_AllelesCEnd() ? false : *itHaplo0 == locus ;
        bool isHaplo1Mutant = itHaplo1 == haplo1.T5_AllelesCEnd() ? false : *itHaplo1 == locus ;

        if (isHaplo0Mutant && isHaplo1Mutant)
        {
            // homozygote mutant
            r *= fits[*itHaplo0 + 1]; // 1-s
            itHaplo0++;
            itHaplo1++;
        } else
        {
            
            if (isHaplo0Mutant)
            {
                r *= fits[locus]; // 1-hs    
                itHaplo0++;
            }

            if (isHaplo1Mutant)
            {
                r *= fits[locus]; // 1-hs
                itHaplo1++;
            }
        }
    }

    return r;
}



