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




double Individual::CalculateT1FitnessNoMultiplicity(const int& Habitat)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'CalculateT1FitnessNoMultiplicity'\n";
#endif   
    auto& fits = SSP->T1_FitnessEffects[Habitat];
    assert(fits.size() == SSP->T1_nbLoci * THREE);


    double W_T1_WholeIndividual = 1.0;

    // All but the last byte
    int fitPosStart = 0;
    for (int byte = 0 ; byte < (SSP->T1_nbChars - 1) ; ++byte)
    {    
        for (auto bit = 0 ; bit < 8 ; ++bit)       
        {
            W_T1_WholeIndividual *= 
                fits[
                    fitPosStart + 
                    haplo0.getT1_Allele(byte, bit) +
                    haplo1.getT1_Allele(byte, bit)
                ];
            fitPosStart += 3;
        }
    }

    // Last byte
    auto lastByte = SSP->T1_nbChars - 1;
    for (auto bit = 0 ; bit < SSP->T1_nbLociLastByte ; ++bit)       
    {
        /*
        assert(fitPosStart == (lastByte * 8 + bit) * 3);
        std::cout << "fitPosStart = " << fitPosStart << "\n";
        std::cout << "index = " << fitPosStart + haplo0.getT1_Allele(lastByte, bit) + haplo1.getT1_Allele(lastByte, bit) << "\n";
        std::cout << "fit = " << fits[fitPosStart + haplo0.getT1_Allele(lastByte, bit) +haplo1.getT1_Allele(lastByte, bit)] << "\n";
        */
        W_T1_WholeIndividual *= 
            fits[
                fitPosStart + 
                haplo0.getT1_Allele(lastByte, bit) +
                haplo1.getT1_Allele(lastByte, bit)
            ];
        fitPosStart += 3;
    }
    assert(fitPosStart == SSP->T1_nbLoci * 3);

    /*for (auto& polymorphicLocus : SSP->simTracker.T1SitesForFitnessNoMultiplicityCalculation)
    {
        W_T1_WholeIndividual *= 
            fits[
                THREE * polymorphicLocus.locus + 
                haplo0.getT1_Allele(polymorphicLocus.byte_index,polymorphicLocus.bit_index) +
                haplo1.getT1_Allele(polymorphicLocus.byte_index,polymorphicLocus.bit_index)
            ];
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
    for (int byte_index = 0 ; byte_index < SSP->T3_nbLoci; byte_index++)
    {
        double allele = (double) haplo0.getT3_Allele(byte_index);
        for (int dim = 0 ; dim < SSP->T3_PhenoNbDimensions; dim++)
        {
            T3_IndPhenotype[dim] += SSP->T3_PhenotypicEffects[Habitat][byte_index * SSP->T3_PhenoNbDimensions + dim] * allele;
            if (SSP->T3_DevelopmentalNoiseStandardDeviation[Habitat][dim] != 0.0)
            {
                std::normal_distribution<double> dist(0.0,SSP->T3_DevelopmentalNoiseStandardDeviation[Habitat][dim]);
                T3_IndPhenotype[dim] += dist(GP->rngw.getRNG());
            }
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


const std::vector<double>& Individual::CalculateFitnessComponents(const int& Habitat)
{
    //std::cout << "Enters in Individual::CalculateFitnessComponents\n";

    if (!SSP->isAnySelection)
    {
        return fitnessComponents; // Should be {1,1,1,1} by default
    }

    double& rT1 = fitnessComponents[0];
    double& rT2 = fitnessComponents[1];
    double& rT1epistasis = fitnessComponents[2];
    double& rT3 = fitnessComponents[3];
    double& rT56 = fitnessComponents[4];

    rT1 = 1.0;
    rT2 = 1.0;
    rT1epistasis = 1.0;
    rT3 = 1.0;
    rT56 = 1.0;
      
    // Trait 1 (excludes epistasis)
    if (SSP->T1_isSelection)
    {
        if (SSP->T1_isMultiplicitySelection)
        {
            rT1 =  haplo0.CalculateT1FitnessMultiplicity(Habitat);
            rT1 *= haplo1.CalculateT1FitnessMultiplicity(Habitat);
        } else
        {
            rT1 = this->CalculateT1FitnessNoMultiplicity(Habitat);
        }
    }
        

    // Trait 1 epistasis
    if (SSP->T1_isEpistasis)
    {
        rT1epistasis = this->CalculateT1EpistaticFitness(Habitat);
    }


    // Trait 2
    if (SSP->T2_isSelection)
    {
        rT2 = haplo0.CalculateT2Fitness(Habitat);
        rT2 *= haplo1.CalculateT2Fitness(Habitat);
    }

    // Trait 3
    if (SSP->T3_isSelection)
    {
        this->CalculateT3Phenotype(Habitat); // save phenotype in the static attribute 'T3_IndPhenotype'
        rT3 = Individual::CalculateT3Fitness(Habitat); // calculate fitness associated to the phenotype saved in the static attribute 'T3_IndPhenotype'
    }

    // Trait 5 
    if (SSP->T56_isSelection)
    {
        if (SSP->T56_isMultiplicitySelection)
        {
            rT56  = haplo0.CalculateT56FitnessMultiplicity(Habitat);
            rT56 *= haplo1.CalculateT56FitnessMultiplicity(Habitat);
        } else
        {
            if (SSP->T5sel_nbLoci)
            {
                //std::cout << "haplo0.T5sel_Alleles.size() = " << haplo0.T5sel_Alleles.size() << "\n";
                //std::cout << "haplo1.T5sel_Alleles.size() = " << haplo1.T5sel_Alleles.size() << "\n";
                rT56 = this->CalculateT56FitnessNoMultiplicity(Habitat, haplo0.T5sel_AllelesBegin(), haplo1.T5sel_AllelesBegin(), haplo0.T5sel_AllelesEnd(), haplo1.T5sel_AllelesEnd());
            }
            else if (SSP->T6sel_nbLoci)
            {
                //std::cout << "haplo0.T6sel_Alleles.size() = " << haplo0.T6sel_Alleles.size() << "\n";
                //std::cout << "haplo1.T6sel_Alleles.size() = " << haplo1.T6sel_Alleles.size() << "\n";
                rT56 = this->CalculateT56FitnessNoMultiplicity(Habitat, haplo0.T6sel_AllelesBegin(), haplo1.T6sel_AllelesBegin(), haplo0.T6sel_AllelesEnd(), haplo1.T6sel_AllelesEnd());
            } else
            {
                std::cout << "Internal error in Individual::CalculateFitnessComponents(const int& Habitat). It appears that individuals has both T5sel and T6sel loci.\n";
                abort();
            }
        }
    }

    //std::cout << "has " << haplo0.T5_howManyMutations() + haplo1.T5_howManyMutations() << " mutations\n";
    //std::cout << "rT56 = " <<rT56 << "\n";
    //std::cout << "Exits in Individual::CalculateFitnessComponents\n";
    return fitnessComponents;
}

double Individual::CalculateFitness(const int& patch_index)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'CalculateFitness'\n";
#endif  

    // Figure out what Habitat it corresponds to
    const int& Habitat = SSP->Habitats[patch_index];

    // Get fitness components
    (void) CalculateFitnessComponents(Habitat); // will set the static fitnessComponents

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
    haplo0 = Haplotype(ShouldReadPopFromBinary);
    haplo1 = Haplotype(ShouldReadPopFromBinary);
}

Individual::Individual(const int patch_index, char Abiogenesis, int ind_index)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'Individual::Individual(const int patch_index, char Abiogenesis)'\n";
#endif
    haplo0 = Haplotype(patch_index,Abiogenesis, ind_index*2);
    haplo1 = Haplotype(patch_index,Abiogenesis, ind_index*2+1);
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

Individual::Individual(Haplotype& h0, Haplotype& h1)
{
    haplo0 = h0;
    haplo1 = h1;
}

Individual::Individual(const Haplotype& knownHaplotype)
{
    haplo0 = knownHaplotype;
    haplo1 = knownHaplotype;   
}


Individual::Individual(){}

Individual::Individual(const Individual& I)
{
    haplo0 = I.haplo0;
    haplo1 = I.haplo1;
}


Individual::Individual(const Individual&& I)
{
    haplo0 = I.haplo0;
    haplo1 = I.haplo1;
}

Individual Individual::operator=(const Individual& I)
{
    haplo0 = I.haplo0;
    haplo1 = I.haplo1;
    return *this;
}


Individual Individual::operator=(const Individual&& I)
{
    haplo0 = I.haplo0;
    haplo1 = I.haplo1;
    return *this;
}






///////////// Fitness On Subset Of Loci

bool Individual::isLocusIsInSet(const int locus, const std::vector<int>& LociSet)
{
    return std::find(LociSet.begin(), LociSet.end(), locus) != LociSet.end();
}







double Individual::CalculateT1FitnessNoMultiplicityOnSubsetOfLoci(const int& Habitat, const std::vector<int>& LociSet)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'CalculateT1FitnessNoMultiplicityOnSubsetOfLoci'\n";
#endif   
    auto& fits = SSP->T1_FitnessEffects[Habitat];
    assert(fits.size() == SSP->T1_nbLoci * THREE);


    double W_T1_WholeIndividual = 1.0;

    for (const int locus : LociSet)
    {
        int byte_index = locus / 8;
        int bit_index = locus % 8;
        //std::cout << "locus = " << locus << "\n";
        assert(locus >= 0 && locus < SSP->T1_nbLoci);

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
            T3_IndPhenotype[dim] += SSP->T3_PhenotypicEffects[Habitat][locus * SSP->T3_PhenoNbDimensions + dim] * allele + dist(GP->rngw.getRNG());
        }
    }
}


template<typename Iterator>
double Individual::CalculateT56FitnessNoMultiplicity(const int& Habitat, Iterator itHaplo0, Iterator itHaplo1, Iterator itHaplo0End, Iterator itHaplo1End)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'CalculateT56FitnessNoMultiplicity'\n";
#endif   
    auto& fits = SSP->T56_FitnessEffects[Habitat];
    assert(fits.size() == SSP->T56sel_nbLoci * 2);

    double w = 1.0;


    while (itHaplo0 != itHaplo0End || itHaplo1 != itHaplo1End)
    {
        if (itHaplo0 != itHaplo0End && itHaplo1 != itHaplo1End)
        {
            if (*itHaplo0 == *itHaplo1)
            {
                w *= fits[*itHaplo0+*itHaplo0 + 1]; // 1-s
                ++itHaplo0;
                ++itHaplo1;
            } else if (*itHaplo0 < *itHaplo1)
            {
                w *= fits[*itHaplo0+*itHaplo0]; // 1-hs
                ++itHaplo0;
            } else // if (*itHaplo0 > *itHaplo1)
            {
                w *= fits[*itHaplo1+*itHaplo1]; // 1-hs
                ++itHaplo1;
            }
        } else if (itHaplo0 != itHaplo0End)
        {
            w *= fits[*itHaplo0+*itHaplo0]; // 1-hs
            ++itHaplo0;
        } else // if (itHaplo1 != itHaplo1End)
        {
            w *= fits[*itHaplo1+*itHaplo1]; // 1-hs
            ++itHaplo1;
        }     
    }

    return w;
}


void Individual::toggleT56LociFromHaplotypes(std::vector<int>& T5ntrlLociToToggle, std::vector<int>& T5selLociToToggle, int Habitat)
{
    if (T5ntrlLociToToggle.size())
    {
        haplo0.toggleT56ntrlLoci(T5ntrlLociToToggle);
        haplo1.toggleT56ntrlLoci(T5ntrlLociToToggle);
    }

    assert(T5selLociToToggle.size() == 0); // No more toggling of selected loci in this version of SimBit
    /*if (T5selLociToToggle.size())
    {
        haplo0.toggleT56selLoci(T5selLociToToggle, Habitat);
        haplo1.toggleT56selLoci(T5selLociToToggle, Habitat);
    }*/ 

}



std::vector<double> Individual::CalculateFitnessComponentsOnSubsetOfLoci(const int& Habitat, const int lociSetIndex)
{
    if (!SSP->isAnySelection)
    {
        return {1.0, 1.0, 1.0, 1.0, 1.0};
    }

    double rT1 = 1.0;
    double rT2 = 1.0;
    double rT1epistasis = 1.0;
    double rT3 = 1.0;
    double rT56 = 1.0;

    auto& lociSet = SSP->subsetT5LociForfitnessSubsetLoci_file[lociSetIndex];
      
    // Trait 1 (excludes epistasis)
    if (SSP->T1_isSelection)
    {
        if (SSP->T1_isMultiplicitySelection)
        {
            rT1 =  haplo0.CalculateT1FitnessMultiplicityOnSubsetOfLoci(Habitat, lociSet);
            rT1 *= haplo1.CalculateT1FitnessMultiplicityOnSubsetOfLoci(Habitat, lociSet);
        } else
        {
            rT1 = this->CalculateT1FitnessNoMultiplicityOnSubsetOfLoci(Habitat, lociSet);
        }
    }
        

    // Trait 1 epistasis
    if (SSP->T1_isEpistasis)
    {
        rT1epistasis = this->CalculateT1EpistaticFitnessOnSubsetOfLoci(Habitat, lociSet);
    }


    // Trait 2
    if (SSP->T2_isSelection)
    {
        rT2 = haplo0.CalculateT2FitnessOnSubsetOfLoci(Habitat, lociSet);
        rT2 *= haplo1.CalculateT2FitnessOnSubsetOfLoci(Habitat, lociSet);
    }

    // Trait 3
    if (SSP->T3_isSelection)
    {
        this->CalculateT3PhenotypeOnSubsetOfLoci(Habitat, lociSet); // save phenotype in the static attribute 'T3_IndPhenotype'
        rT3 = Individual::CalculateT3Fitness(Habitat); // calculate fitness associated to the phenotype saved in the static attribute 'T3_IndPhenotype'
    }

    // Trait 5 
    if (SSP->T56_isSelection)
    {
        if (SSP->T56_isMultiplicitySelection)
        {
            rT56  = haplo0.CalculateT56FitnessMultiplicityOnSubsetOfLoci(Habitat, lociSet);
            rT56 *= haplo1.CalculateT56FitnessMultiplicityOnSubsetOfLoci(Habitat, lociSet);
        } else
        {
            if (SSP->T5sel_nbLoci)
            {
                //std::cout << "haplo0.T5sel_Alleles.size() = " << haplo0.T5sel_Alleles.size() << "\n";
                //std::cout << "haplo1.T5sel_Alleles.size() = " << haplo1.T5sel_Alleles.size() << "\n";
                rT56 = this->CalculateT56FitnessNoMultiplicityOnSubsetOfLoci(Habitat, lociSet, haplo0.T5sel_AllelesBegin(), haplo1.T5sel_AllelesBegin(), haplo0.T5sel_AllelesEnd(), haplo1.T5sel_AllelesEnd());
            }
            else if (SSP->T6sel_nbLoci)
            {
                //std::cout << "haplo0.T6sel_Alleles.size() = " << haplo0.T6sel_Alleles.size() << "\n";
                //std::cout << "haplo1.T6sel_Alleles.size() = " << haplo1.T6sel_Alleles.size() << "\n";
                rT56 = this->CalculateT56FitnessNoMultiplicityOnSubsetOfLoci(Habitat, lociSet, haplo0.T6sel_AllelesBegin(), haplo1.T6sel_AllelesBegin(), haplo0.T6sel_AllelesEnd(), haplo1.T6sel_AllelesEnd());
            } else
            {
                std::cout << "Internal error in Individual::CalculateFitnessComponents(const int& Habitat). It appears that individuals has both T5sel and T6sel loci.\n";
                abort();
            }
        }
    }

    //std::cout << "has " << haplo0.T5_howManyMutations() + haplo1.T5_howManyMutations() << " mutations\n";
    //std::cout << "rT56 = " <<rT56 << "\n";
    //std::cout << "Exits in Individual::CalculateFitnessComponents\n";
    return {rT1, rT2, rT1epistasis, rT3, rT56};
}


template<typename Iterator>
double Individual::CalculateT56FitnessNoMultiplicityOnSubsetOfLoci(const int& Habitat, const std::vector<int>& LociSet, Iterator itHaplo0, Iterator itHaplo1, Iterator itHaplo0End, Iterator itHaplo1End)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'CalculateT56FitnessNoMultiplicityOnSubsetOfLoci'\n";
#endif   
    auto& fits = SSP->T56_FitnessEffects[Habitat];
    assert(fits.size() == SSP->T56_nbLoci * 2);

    double r = 1.0;

    for (auto& locus : LociSet)
    {
        unsigned int impossibleValue = SSP->T56_nbLoci;
        assert(locus < impossibleValue);
        unsigned int ValueFromH0It = impossibleValue;
        size_t value = *itHaplo0;
        while (itHaplo0 != itHaplo0End && locus > value)
        {
            ++itHaplo0;
            ValueFromH0It = value;
        }

        unsigned int ValueFromH1It = impossibleValue;
        value = *itHaplo0;
        while (itHaplo1 != itHaplo1End && locus > value)
        {
            ++itHaplo0;
            ValueFromH1It = value;
        }


        bool isHaplo0Mutant = false;
        bool isHaplo1Mutant = false;
        if (ValueFromH0It != impossibleValue && ValueFromH0It == locus)
        {
            isHaplo0Mutant = true;
        }
        if (ValueFromH1It != impossibleValue && ValueFromH1It == locus)
        {
            isHaplo1Mutant = true;
        }
        

        if (isHaplo0Mutant && isHaplo1Mutant)
        {
            // homozygote mutant
            r *= fits[2*locus + 1]; // 1-s
        } else
        {
            if (isHaplo0Mutant)
            {
                r *= fits[2*(locus)]; // 1-hs    
            }

            if (isHaplo1Mutant)
            {
                r *= fits[2*(locus)]; // 1-hs
            }
        }
    }

    return r;
}
