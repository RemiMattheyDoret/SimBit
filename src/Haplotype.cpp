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



bool Haplotype::isFreeFromMutations(int T1_locusFrom, int T1_locusTo)
{
    assert(SSP->T1_nbChars == T1_Alleles.size());
    
    assert(T1_locusTo <= (T1_Alleles.size() * 8));
    for (int T1_locus = T1_locusFrom ; T1_locus < T1_locusTo ; T1_locus++)
    {
        if (this->getT1_Allele(T1_locus) == 1)
        {
            return false;
        }
    }
     
    return true;   
}

bool Haplotype::isFreeFromMutations()
{
    assert(SSP->T1_nbChars == T1_Alleles.size());
    for (int byte_index = 0 ; byte_index < SSP->T1_nbChars ; byte_index++)
    {
        for (int bit_index = 0 ; bit_index < 8 ; bit_index++)   
        {
            if (this->getT1_Allele(byte_index, bit_index) == 1)
            {
                return false;
            }
        }
    }
    return true;
}

double Haplotype::getW_T1(int fitnessMapIndex)
{
    assert(fitnessMapIndex < W_T1.size());
    return W_T1[fitnessMapIndex];
}
double Haplotype::getW_T2(int fitnessMapIndex)
{
    assert(fitnessMapIndex < W_T2.size());
    return W_T2[fitnessMapIndex];
}

void Haplotype::setAllW_T1(double w)
{
    assert(W_T1.size() == SSP->NbElementsInFitnessMap);
    for (int fitnessMapIndex = 0 ; fitnessMapIndex < SSP->NbElementsInFitnessMap ; fitnessMapIndex++)
    {
        this->setW_T1(w, fitnessMapIndex);
    }
}
void Haplotype::setAllW_T2(double w)
{
    for (int fitnessMapIndex = 0 ; fitnessMapIndex < SSP->NbElementsInFitnessMap ; fitnessMapIndex++)
    {
        this->setW_T2(w, fitnessMapIndex);
    }
}

void Haplotype::setW_T1(double w, int fitnessMapIndex)
{
    assert(W_T1.size() > fitnessMapIndex && fitnessMapIndex >= 0);
    W_T1[fitnessMapIndex] = w;
}

void Haplotype::setW_T2(double w, int fitnessMapIndex)
{
    assert(W_T2.size() > fitnessMapIndex && fitnessMapIndex >= 0);
    W_T2[fitnessMapIndex] = w;
}

template<typename INT>
unsigned char Haplotype::getT1_char(INT char_index)
{
    //assert(char_index >= 0);
    //assert(char_index < SSP->T1_nbChars);
    return T1_Alleles[char_index];
}

bool Haplotype::getT1_Allele(const int T1Locus)
{
    const int char_index = T1Locus / 8;
    const int bit_index  = T1Locus % 8;
    return (this->getT1_Allele(char_index, bit_index));
}

bool Haplotype::getT1_Allele(const int char_index, const int bit_index)
{
    /*
    assert(bit_index < 8);
    assert(bit_index >= 0);
    assert(char_index >= 0);
    assert(char_index < SSP->T1_nbChars);
    */
    /*if ((T1_Alleles[char_index] >> bit_index) & 1)
    {
        std::cout << "got " << char_index << " " << bit_index << "("<<char_index * 8 + bit_index << "). Mut\n";    
    } else
    {
        std::cout << "got " << char_index << " " << bit_index << "("<<char_index * 8 + bit_index << "). No mut\n";
    }*/
    
    return (T1_Alleles[char_index] >> bit_index) & 1;
}

unsigned char Haplotype::getT2_Allele(const int char_index)
{
    //assert(char_index >= 0);
    //assert(char_index < SSP->T2_nbLoci);
    return T2_Alleles[char_index];
}

char Haplotype::getT3_Allele(const int char_index)
{
    //assert(char_index >= 0);
    //assert(char_index < SSP->T2_nbLoci);
    return T3_Alleles[char_index];
}

void Haplotype::setT1_char(int& char_index, unsigned char& c)
{
    //assert(char_index >= 0);
    //assert(char_index < SSP->T1_nbChars);
    T1_Alleles[char_index] = c;
}

void Haplotype::setT1_char(int& char_index, unsigned char&& c)
{
    //assert(char_index >= 0);
    //assert(char_index < SSP->T1_nbChars);
    T1_Alleles[char_index] = c;
}

template<typename valueType>
void Haplotype::setT1_Allele(const int& char_index, const int& bit_index, const valueType& value)
{
    //assert(bit_index < 8);
    //assert(bit_index >= 0);
    //assert(char_index >= 0);
    //assert(char_index < SSP->T1_nbChars);
    T1_Alleles[char_index] = (-value ^ T1_Alleles[char_index]) & (1 << bit_index);
}


void Haplotype::setT1_AlleleToOne(int& char_index, int& bit_index)
{
    //assert(bit_index < 8);
    //assert(bit_index >= 0);
    //assert(char_index >= 0);
    //assert(char_index < SSP->T1_nbChars);
    T1_Alleles[char_index] |= 1 << bit_index;
}

void Haplotype::setT1_AlleleToZero(int& char_index, int& bit_index)
{;
    //assert(bit_index < 8);
    //assert(bit_index >= 0);
    //assert(char_index >= 0);
    //assert(char_index < SSP->T1_nbChars);
    T1_Alleles[char_index] &= ~(1 << bit_index);
}

void Haplotype::setT2_Allele(const int char_index, const unsigned char value)
{
    //assert(char_index >= 0);
    //assert(char_index < SSP->T2_nbLoci);
    T2_Alleles[char_index] = value;
}

void Haplotype::setT3_Allele(const int char_index, const char value)
{
    //assert(char_index >= 0);
    //assert(char_index < SSP->T2_nbLoci);
    T3_Alleles[char_index] = value;
}

void Haplotype::toggleT1_Allele(int& byte_index, int& bit_index)
{
    T1_Alleles[byte_index] ^= 1 << bit_index;
}


void Haplotype::mutateT1_Allele(int MutPosition, int& Habitat)
{
    #ifdef DEBUG
    NBT1MUTATIONS++;
    #endif
    //std::cout << "mutation at " << MutPosition << "\n";
    int byte_index = MutPosition / 8;
    int bit_index  = MutPosition % 8;
    //assert(bit_index < 8);
    //assert(bit_index >= 0);
    //assert(byte_index >= 0);
    //assert(byte_index < SSP->T1_nbChars);


    toggleT1_Allele(byte_index, bit_index);
    //SSP->simTracker.addMutation(byte_index, bit_index, MutPosition);

    if (SSP->T1_isMultiplicitySelection)
    {
        int fitnessMapIndex = SSP->FromLocusToFitnessMapIndex[SSP->FromT1LocusToLocus[MutPosition]];
        double w = this->getW_T1(fitnessMapIndex);
        //if (GP->CurrentGeneration >=50) std::cout << "w=" << w << " ";
        if (w != -1.0)
        {
            // Note that Toggle already happened. If it is mutant, it is because it just happened!
            if ( this->getT1_Allele(byte_index, bit_index) ) // mutation added
            {
                w *= SSP->T1_FitnessEffects[Habitat][byte_index * 8 + bit_index];
                assert(w >= 0.0 && w <= 1.0);
            } else // mutation removed
            {
                if (SSP->T1_FitnessEffects[Habitat][byte_index * 8 + bit_index] == 0)
                { 
                    w = -1.0;
                } else
                {
                    w /= SSP->T1_FitnessEffects[Habitat][byte_index * 8 + bit_index]; 
                    assert(w >= 0.0 && w <= 1.0);
                }
            }
            this->setW_T1( w,   fitnessMapIndex );
        }
    }

   /* if (SSP->T1mutsDirectional)
    {
        if (!getT1_Allele(byte_index, bit_index))
        {
            setT1_AlleleToOne(byte_index, bit_index);
            SSP->simTracker.addMutation(byte_index, bit_index, MutPosition);
            if (SSP->T1_isMultiplicitySelection)
            {
                int fitnessMapIndex = SSP->FromLocusToFitnessMapIndex[SSP->FromT1LocusToLocus[MutPosition]];
                double w = this->getW_T1(fitnessMapIndex);
                if (w != -1.0)
                {
                    w *= SSP->T1_FitnessEffects[Habitat][byte_index * 8 + bit_index];
                    assert(w >= 0.0 && w <= 1.0);
                    this->setW_T1( w,   fitnessMapIndex );
                }
            }
        }            
    } else
    {
        toggleT1_Allele(byte_index, bit_index);
        SSP->simTracker.addMutation(byte_index, bit_index, MutPosition);

        if (SSP->T1_isMultiplicitySelection)
        {
            int fitnessMapIndex = SSP->FromLocusToFitnessMapIndex[SSP->FromT1LocusToLocus[MutPosition]];
            double w = this->getW_T1(fitnessMapIndex);
            if (w != -1.0)
            {
                // Note that Toggle already happened. If it is mutant, it is because it just happened!
                if ( this->getT1_Allele(byte_index, bit_index) ) // mutation added
                {
                    w *= SSP->T1_FitnessEffects[Habitat][byte_index * 8 + bit_index];
                    assert(w >= 0.0 && w <= 1.0);
                } else // mutation removed
                {
                    if (SSP->T1_FitnessEffects[Habitat][byte_index * 8 + bit_index] == 0)
                    {
                        w = -1.0;
                    } else
                    {
                        w /= SSP->T1_FitnessEffects[Habitat][byte_index * 8 + bit_index]; 
                        assert(w >= 0.0 && w <= 1.0);
                    }
                }
                this->setW_T1( w,   fitnessMapIndex );
            }
        }
    }*/
}


void Haplotype::AddMutT2_Allele(int char_index)
{
    //assert(char_index >= 0);
    //assert(char_index < SSP->T2_nbLoci);
    // Add one to number of mutations
    if (T2_Alleles[char_index] >= 250)
    {
        assert(SSP != nullptr);
        if (GP->nbWarningsSentFrom_Haplotype_AddMutT2_Allele <= 100)
        {
            GP->nbWarningsSentFrom_Haplotype_AddMutT2_Allele++;
            std::cout << "WARNING: The T2_char index " << char_index << " of species " << SSP->speciesName << " contains 250 mutations (at generation"<< GP->CurrentGeneration <<")! Almost no more mutation can be added! SimBit will reset the value of all individuals at this locus by substracting the number of mutations by the individuals that has the lowest number of mutations at this locus (slow process that should not happen too often). If there is too much variance on the block for correcting, then SimBit will just abort.\n";
            if (GP->nbWarningsSentFrom_Haplotype_AddMutT2_Allele == 100)
            {
                std::cout << "Further warnings in Haplotype::AddMutT2_Allele for overpassing the maximum number of generation will be silenced.\n";
            }
        }
        Pop::addT2LocusToCorrect(char_index);

        if (T2_Alleles[char_index] == 250)
        {
            T2_Alleles[char_index]++;
        }
    } else
    {
        T2_Alleles[char_index]++;
    }
    if (SSP->T2_isSelection)
    {
        int fitnessMapIndex = SSP->FromLocusToFitnessMapIndex[SSP->FromT2LocusToLocus[char_index]];
        this->setW_T2(-1.0,fitnessMapIndex);
    }
}

void Haplotype::AddMutT2_Allele(int char_index, int Habitat)
{
    //assert(char_index >= 0);
    //assert(char_index < SSP->T2_nbLoci);
    // Add one to number of mutations
    if (T2_Alleles[char_index] >= 254)
    {
        assert(SSP != nullptr);
        if (GP->nbWarningsSentFrom_Haplotype_AddMutT2_Allele <= 100)
        {
            GP->nbWarningsSentFrom_Haplotype_AddMutT2_Allele++;
            std::cout << "WARNING: The T2_char index " << char_index << " of species " << SSP->speciesName << " contains 254 mutations (at generation"<< GP->CurrentGeneration <<")! Almost no more mutation can be added! SimBit will reset the value of all individuals at this locus by substracting the number of mutations by the individuals that has the lowest number of mutations at this locus (slow process that should not happen too often). If there is too much variance on the block for correcting, then SimBit will just abort.\n";
            if (GP->nbWarningsSentFrom_Haplotype_AddMutT2_Allele == 100)
            {
                std::cout << "Further warnings in Haplotype::AddMutT2_Allele for overpassing the maximum number of generation will be silenced.\n";
            }
        }
        Pop::addT2LocusToCorrect(char_index);

        if (T2_Alleles[char_index] == 254)
        {
            T2_Alleles[char_index]++;
        }
    } else
    {
        T2_Alleles[char_index]++;
    }

    // Change Fitness
    if (SSP->T2_isSelection)
    {
        int fitnessMapIndex = SSP->FromLocusToFitnessMapIndex[SSP->FromT2LocusToLocus[char_index]];
        if (this->getW_T2(fitnessMapIndex) != -1.0)
        {
            //std::cout << "SSP->T2_FitnessEffects[Habitat].size() = " << SSP->T2_FitnessEffects[Habitat].size() << "    char_index = " <<  char_index << "\n";
            //assert(SSP->T2_FitnessEffects.size() > Habitat);
            //assert(SSP->T2_FitnessEffects[Habitat].size() > char_index);
            //std::cout << "this->getW_T2() = " << this->getW_T2() << "     SSP->T2_FitnessEffects[Habitat][char_index] = " << SSP->T2_FitnessEffects[Habitat][char_index] << "\n";
            this->setW_T2(this->getW_T2(fitnessMapIndex) * SSP->T2_FitnessEffects[Habitat][char_index], fitnessMapIndex);
        }
    }
}

void Haplotype::AddMutT3_Allele(int char_index)
{
    // Add or substract while preventing
#ifdef DEBUG
std::cout << "T3_Alleles.size() = " << T3_Alleles.size() << "\n";
std::cout << "SSP->T3_nbLoci = " << SSP->T3_nbLoci << "\n";
assert(SSP->T3_nbLoci == T3_Alleles.size());    
assert(char_index < T3_Alleles.size());
assert(char_index >= 0);
#endif
    if (GP->rngw.get_1b())
    {
        if (T3_Alleles[char_index] != CHAR_MAX)
        {
            T3_Alleles[char_index]++;
        }
    } else
    {
        if (T3_Alleles[char_index] != CHAR_MIN)
        {
            T3_Alleles[char_index]--;
        }
    }
}


void Haplotype::copyIntoT1(int from, int to, Haplotype& SourceChromo)
{
    // copy [first, last)
    //std::cout << "In copyIntoT1 from = "<< from << "  to = "<< to << std::endl;
    assert(from >= 0);
    assert(to > from);
    
    int ByteFrom = from / 8; // integer division. It will floor automatically
    int ByteTo   = to / 8; // integer division. It will floor automatically
    assert(ByteTo <= SourceChromo.T1_Alleles.size() + 1);

    int bitIndexFrom = from % 8;
    int bitIndexTo = to % 8;

    if (ByteTo == SSP->T1_nbChars - 1) // This is to deal with the last extra bits. Note that if it copies a few extra bits, it is not an issue.
    {
        bitIndexTo = std::max(bitIndexTo, SSP->T1_nbLociLastByte);
    }


    if (ByteTo == ByteFrom)
    {
        assert(bitIndexFrom < bitIndexTo);
        for (int bit_index = bitIndexFrom ; bit_index < bitIndexTo ; bit_index++)
        {
            if (SourceChromo.getT1_Allele(ByteFrom, bit_index))
            {
                this->setT1_AlleleToOne(ByteFrom, bit_index);
            } else
            {
                this->setT1_AlleleToZero(ByteFrom, bit_index);
            }
        }
    } else if (ByteTo > ByteFrom)
    {
        // Copy first bits
        for (int bit_index = bitIndexFrom ; bit_index < 8 ; bit_index++)
        {
            if (SourceChromo.getT1_Allele(ByteFrom, bit_index))
            {
                this->setT1_AlleleToOne(ByteFrom, bit_index);
            } else
            {
                this->setT1_AlleleToZero(ByteFrom, bit_index);
            }
        }
        
        if (ByteFrom + 1 < ByteTo)
        {
            // Copy what is in between. Only if there is at least a full byte to be copied
            std::copy(
                      SourceChromo.T1_Alleles.begin() + ByteFrom + 1,
                      SourceChromo.T1_Alleles.begin() + ByteTo,
                      this->T1_Alleles.begin() + ByteFrom + 1
                      ); // std::copy is always in the range [from,to)
        }
        
        // Copy last bits
        for (int bit_index = 0 ; bit_index < bitIndexTo ; bit_index++)
        {
            if (SourceChromo.getT1_Allele(ByteTo,bit_index))
            {
                this->setT1_AlleleToOne(ByteTo, bit_index);
            } else
            {
                this->setT1_AlleleToZero(ByteTo, bit_index);
            }
        }
    } else
    {
        std::cout << "ByteFrom is " << ByteFrom << " ByteTo is " << ByteTo << ". ByteTo must be equal or greater to ByteFrom.\n";
        abort();
    }
}

void Haplotype::copyIntoT2(int from, int to, Haplotype& SourceChromo)
{
    // copy [first, last)
    assert(from >= 0);
    assert(to > from);
    //std::cout << "In copyIntoT2: to = " << to << std::endl;
    //std::cout << "In copyIntoT2: SourceChromo.T2_Alleles.size() = " << SourceChromo.T2_Alleles.size() << std::endl;
    assert(to <= SourceChromo.T2_Alleles.size() + 1);
    
    std::copy(
              SourceChromo.T2_Alleles.begin() + from,
              SourceChromo.T2_Alleles.begin() + to,
              this->T2_Alleles.begin() + from
              );
}

void Haplotype::copyIntoT3(int from, int to, Haplotype& SourceChromo)
{
    // copy [first, last)
    assert(from >= 0);
    assert(to > from);
    //std::cout << "to = " << to << "  from = " << from << "  SourceChromo.T3_Alleles.size() = " << SourceChromo.T3_Alleles.size() << "\n";
    assert(to <= SourceChromo.T3_Alleles.size() + 1);
    
    std::copy(
              SourceChromo.T3_Alleles.begin() + from,
              SourceChromo.T3_Alleles.begin() + to,
              this->T3_Alleles.begin() + from
              );
}


Haplotype::Haplotype(){}

void Haplotype::PrintBinaryFile(OutputFile& file)
{
    // Print T1_alleles
    file.writeBinary(reinterpret_cast<const char*>(&T1_Alleles[0]), T1_Alleles.size()*sizeof(unsigned char));

    // Print T2_alleles
    file.writeBinary(reinterpret_cast<const char*>(&T2_Alleles[0]), T2_Alleles.size()*sizeof(unsigned char));

    // Print T3_alleles
    file.writeBinary(reinterpret_cast<const char*>(&T3_Alleles[0]), T3_Alleles.size()*sizeof(char));

    // Print T4_alleles
    if (SSP->T4_nbBits)
        std::cout << "\n\tWARNING: You are printing to a binary file while you have T4 loci. In the current version of SimBit, T4 loci cannot be printed on a binary file (sorry). The binary file will be produced ignoring T4 loci\n";

    // Print T5_alleles (there should not be any flipped meaning normally)
    if (SSP->T56ntrl_nbLoci)
    {
        std::vector<unsigned int> T5ntrlTrueHaplotype = this->getT56ntrlTrueHaplotype();
        T5ntrlTrueHaplotype.insert(T5ntrlTrueHaplotype.begin(), T5ntrlTrueHaplotype.size()); // insert the size of the vector as a start

        file.writeBinary(reinterpret_cast<const char*>(&T5ntrlTrueHaplotype[0]), T5ntrlTrueHaplotype.size()*sizeof(unsigned int));
    }

    if (SSP->T56sel_nbLoci)
    {
        std::vector<unsigned int> T5selTrueHaplotype = this->getT56selTrueHaplotype();
        T5selTrueHaplotype.insert(T5selTrueHaplotype.begin(), T5selTrueHaplotype.size()); // insert the size of the vector as a start
        
        file.writeBinary(reinterpret_cast<const char*>(&T5selTrueHaplotype[0]), T5selTrueHaplotype.size()*sizeof(unsigned int));
    }
}

Haplotype::Haplotype(std::vector<unsigned char> T1_info, std::vector<unsigned char> T2_info, std::vector<char> T3_info, std::vector<unsigned int> T56_info)
{
    assert(SSP != nullptr); 

    /*
    std::cout << "T1_info: ";
    for (auto& e : T1_info) {std::cout << e << " ";}
    std::cout << "\n";
    */

    T1_Alleles = T1_info;
    assert(T1_Alleles.size() == SSP->T1_nbChars);
    T2_Alleles = T2_info;
    assert(T2_Alleles.size() == SSP->T2_nbLoci);
    T3_Alleles = T3_info;
    assert(T3_Alleles.size() == SSP->T3_nbLoci);

    assert(SSP->isT56neutral.size() == SSP->T56_nbLoci);
    assert(T56_info.size() <= SSP->T56_nbLoci);


    // Assign values
    for (unsigned T56_info_index = 0 ; T56_info_index < T56_info.size() ; ++T56_info_index)
    {
        auto& mutPos = T56_info[T56_info_index];

        // Make sure T56_info is strictly increasing
        if (T56_info_index > 0)
        {
            assert(mutPos > T56_info[T56_info_index-1]);
        }

        // assign mutPos
        if (SSP->isT56neutral[mutPos])
        {
            if (SSP->T56ntrl_compress)
            {
                T6ntrl_Alleles.push_back(mutPos);
            } else
            {
                T5ntrl_Alleles.push_back(mutPos);
            }
        } else
        {
            if (SSP->T56sel_compress)
            {
                T6sel_Alleles.push_back(mutPos);
            } else
            {
                T5sel_Alleles.push_back(mutPos);
            }
        }
    }


    if (SSP->T1_isMultiplicitySelection)
    {
        assert(SSP->NbElementsInFitnessMap > 0);
        assert(SSP->T1_nbLoci > 0);
        W_T1.resize(SSP->NbElementsInFitnessMap,-1.0);
    }
    if (SSP->T2_isSelection)
    {
        assert(SSP->NbElementsInFitnessMap > 0);
        assert(SSP->T2_nbLoci > 0);
        W_T2.resize(SSP->NbElementsInFitnessMap,-1.0);
    }
    if (SSP->T56_isMultiplicitySelection)
    {
        assert(SSP->NbElementsInFitnessMap > 0);
        assert(SSP->T56sel_nbLoci > 0);
        W_T56.resize(SSP->NbElementsInFitnessMap,-1.0);
    }
}

Haplotype::Haplotype(bool ShouldReadPopFromBinary)
{
#ifdef CALLENTRANCEFUNCTIONS
    std::cout << "Enters in 'Haplotype::Haplotype(bool ShouldReadPopFromBinary)'\n";
#endif     
    (void) ShouldReadPopFromBinary;

    //.read T1_alleles
    T1_Alleles.resize(SSP->T1_nbChars);
    GP->BinaryFileToRead.read(reinterpret_cast<char*>(&T1_Alleles[0]), SSP->T1_nbChars*sizeof(unsigned char));

    //.read T2_alleles
    T2_Alleles.resize(SSP->T2_nbLoci);
    GP->BinaryFileToRead.read(reinterpret_cast<char*>(&T2_Alleles[0]), SSP->T2_nbLoci*sizeof(unsigned char));

    //.read T3_alleles
    T3_Alleles.resize(SSP->T2_nbLoci);
    GP->BinaryFileToRead.read(reinterpret_cast<char*>(&T3_Alleles[0]), SSP->T3_nbLoci*sizeof(char));

    //.read T4_alleles
    if (SSP->T4_nbBits)
        std::cout << "\n\n\tWARNING: You asked for T4 alleles and asked for data to be read from binary file. As per the current version of SimBit, binary file cannot store T4 alleles data. Those loci will be initiated as perfectly unmutated\n\n";

    //.read T56_alleles
    if (SSP->T56ntrl_nbLoci)
    {
        unsigned int nbElements = 99999;
        GP->BinaryFileToRead.read(reinterpret_cast<char*>(&nbElements), sizeof(unsigned int));
        assert(nbElements >= 0);
        assert(nbElements < SSP->T56ntrl_nbLoci);
        
        if (SSP->T56ntrl_compress)
        {
            T6ntrl_Alleles = CompressedSortedDeque(SSP->T56ntrl_nbLoci);
            std::vector<unsigned int> tmp(nbElements);
            GP->BinaryFileToRead.read(reinterpret_cast<char*>(&tmp[0]), nbElements*sizeof(unsigned int));
            assert(tmp.size() == nbElements);
            CompressedSortedDeque tmp2(tmp, SSP->T6ntrl_nbLoci);
            T6ntrl_Alleles = tmp2;
        } else
        {
            T5ntrl_Alleles.resize(nbElements);
            GP->BinaryFileToRead.read(reinterpret_cast<char*>(&T5ntrl_Alleles[0]), nbElements*sizeof(unsigned int));
        }
    }
        

    if (SSP->T56sel_nbLoci)
    {
        int nbElements = -1;
        GP->BinaryFileToRead.read(reinterpret_cast<char*>(&nbElements), sizeof(unsigned int));
        assert(nbElements >= 0);
        assert(nbElements < SSP->T56sel_nbLoci);
        if (SSP->T56sel_compress)
        {
            T6sel_Alleles = CompressedSortedDeque(SSP->T56sel_nbLoci);
            std::vector<unsigned int> tmp(nbElements);
            GP->BinaryFileToRead.read(reinterpret_cast<char*>(&tmp[0]), nbElements*sizeof(unsigned int));
            assert(tmp.size() == nbElements);
            CompressedSortedDeque tmp2(tmp, SSP->T6sel_nbLoci);
            T6sel_Alleles = tmp2;
        } else
        {
            T5sel_Alleles.resize(nbElements);
            GP->BinaryFileToRead.read(reinterpret_cast<char*>(&T5sel_Alleles[0]), nbElements*sizeof(unsigned int));
        }
    }
        
    


    // Initiate W_T1, W_T2 and W_T5
    if (SSP->T1_isMultiplicitySelection)
        this->W_T1.resize(SSP->NbElementsInFitnessMap, -1.0);

    if (SSP->T2_isSelection)
        this->W_T2.resize(SSP->NbElementsInFitnessMap, -1.0);

    if (SSP->T56_isMultiplicitySelection)
        this->W_T56.resize(SSP->NbElementsInFitnessMap, -1.0);

    /*
    for (auto& elem : T1_Alleles)
    {
        printf("%u ", elem);
    }
    std::cout << "\n";
    */
}

Haplotype::Haplotype(const int patch_index, char Abiogenesis, int indHaplo_index) // indHaplo_index is a double for faster comparison
: T1_Alleles(SSP->T1_nbChars), T2_Alleles(SSP->T2_nbLoci), T3_Alleles(SSP->T3_nbLoci) // No initilization for T5_Alleles and T6_Alleles
{
#ifdef CALLENTRANCEFUNCTIONS
    std::cout << "Enters in 'Haplotype::Haplotype(const int patch_index, char Abiogenesis)'\n";
#endif      
    (void)Abiogenesis; // Does nothing but silence the warning that the char Abiogenesis is not used.
    assert(SSP!=nullptr);
    assert(SSP->T1_nbChars + SSP->T2_nbLoci + SSP->T3_nbLoci + SSP->T4_nbBits + SSP->T56_nbLoci > 0);
    assert(indHaplo_index < 2 * SSP->patchSize[patch_index]);

    // Initiate W_T1 and W_T2
    // Initiate W_T1, W_T2 and W_T5
    if (SSP->T1_isMultiplicitySelection)
        this->W_T1.resize(SSP->NbElementsInFitnessMap, -1.0);

    if (SSP->T2_isSelection)
        this->W_T2.resize(SSP->NbElementsInFitnessMap, -1.0);

    if (SSP->T56_isMultiplicitySelection)
        this->W_T56.resize(SSP->NbElementsInFitnessMap, -1.0);
    
    //std::cout << "created at W_T56.size() = " << W_T56.size() << "\n";
    
    // Initiate T1_Alleles
    if (SSP->T1_nbLoci > 0)
    {        
        if (SSP->T1_Initial_AlleleFreqs_AllZeros) // The whole chromosome must be 0s
        {
            unsigned char c = 0;
            for (int char_index = 0 ; char_index < SSP->T1_nbChars ; char_index++)
            {
                this->setT1_char(char_index,c);
            }
        } else if (SSP->T1_Initial_AlleleFreqs_AllOnes) // The whole chromosome must be 1s
        {
            unsigned char c = 255;
            for (int char_index = 0 ; char_index < SSP->T1_nbChars ; char_index++)
            {
                this->setT1_char(char_index,c);
            }
        } else 
        {
            assert(SSP->T1_Initial_AlleleFreqs.size() > patch_index);
            for (int byte_index = 0 ; byte_index < SSP->T1_nbChars ; byte_index++)
            {
                int bit_index_to;
                if (byte_index == SSP->T1_nbChars - 1)
                {
                    bit_index_to = SSP->T1_nbLociLastByte;
                } else
                {
                    bit_index_to = 8;
                }
                for (int bit_index=0 ; bit_index < bit_index_to ; bit_index++ )
                {
                    int locus = byte_index * 8 + bit_index;

                    assert(SSP->T1_Initial_AlleleFreqs[patch_index].size() > locus);
                    if (SSP->T1_Initial_AlleleFreqs[patch_index][locus] == 1.0)
                    {
                        this->setT1_AlleleToOne(byte_index,bit_index);
                    } else if (SSP->T1_Initial_AlleleFreqs[patch_index][locus] == 0.0)
                    {
                        this->setT1_AlleleToZero(byte_index,bit_index);
                    
                    } else
                    {
                        // Find FakeindHaplo_index. FakeindHaplo_index is a function of locus too so as to reduce LD
                        int locus = byte_index * 8 + bit_index;
                        int FakeindHaplo_index = indHaplo_index - (int)(((long)locus*97) % (long)(2 * SSP->patchSize[patch_index]));
                        //                                         ^ 97 is just a large-ish prime number
                        if (FakeindHaplo_index < 0)
                        {
                            FakeindHaplo_index = 2 * SSP->patchSize[patch_index] + FakeindHaplo_index;
                        }
                        assert(FakeindHaplo_index >= 0);
                        assert(2 * SSP->patchSize[patch_index] > FakeindHaplo_index);

                        // Check if mutation must be put
                        int nbMuts = (int)(SSP->T1_Initial_AlleleFreqs[patch_index][locus] * 2.0 * (double) SSP->patchSize[patch_index]);
                        bool setToOne = false;
                        if (FakeindHaplo_index%2 == 0)
                        {
                            if (FakeindHaplo_index <= nbMuts)
                            {
                                setToOne = true;
                            }
                        } else
                        {
                            if (FakeindHaplo_index <= nbMuts - SSP->patchSize[patch_index])
                            {//                                        ^ This is the max number of mutations that have already been added by first (couple of) if statement(s)       
                                setToOne = true;
                            }
                        }


                        // make mutation if needed
                        if (setToOne)
                        {
                            this->setT1_AlleleToOne(byte_index,bit_index);
                        } else
                        {
                            this->setT1_AlleleToZero(byte_index,bit_index);
                        }
                    }
                }
            }
        }
        assert(T1_Alleles.size() == SSP->T1_nbChars);
    }
        

    // Initiate T2_Alleles
    if (SSP->T2_nbLoci > 0)
    {
        for (int char_index=0 ; char_index < SSP->T2_nbLoci ; char_index++)
        {
            unsigned char c = 0;
            this->setT2_Allele(char_index,c);
        }
        assert(T2_Alleles.size() == SSP->T2_nbLoci);
    }

    // Initiate T3_Alleles
    if (SSP->T3_nbLoci > 0)
    {
        for (int char_index=0 ; char_index < SSP->T3_nbLoci ; char_index++)
        {
            char c = 0;
            this->setT3_Allele(char_index,c);
        }
        assert(T3_Alleles.size() == SSP->T3_nbLoci);
    }

    // Initiate T4_Alleles
    // No initialization. Assume all are fixed at 0

    // Initiate T56_Alleles 
    assert(SSP->T5sel_nbLoci == 0 || SSP->T6sel_nbLoci == 0);
    assert(SSP->T5ntrl_nbLoci == 0 || SSP->T6ntrl_nbLoci == 0);

    if (SSP->T56_nbLoci > 0)
    { 
        // Initialize T6
        if (SSP->T6ntrl_nbLoci)
        {
            T6ntrl_Alleles = CompressedSortedDeque(SSP->T6ntrl_nbLoci);
            (void) T6ntrl_AllelesBegin();
        }

        if (SSP->T6sel_nbLoci)
        {
            T6sel_Alleles = CompressedSortedDeque(SSP->T6sel_nbLoci);
            (void) T6sel_AllelesBegin();
        }


        // Initial Allele Frequencies
        if (SSP->T56_Initial_AlleleFreqs_AllZeros) // The whole chromosome must be 0s
        {
            // Nothing to do
            if (SSP->T5ntrl_nbLoci)
                assert(SSP->T5ntrl_flipped.size() == 0);
        } else if (SSP->T56_Initial_AlleleFreqs_AllOnes) // The whole chromosome must be 1s
        {   
            // Nothing to do
            std::cout << "Internal error: It seems that T56_Initial_AlleleFreqs_AllOnes is true but this should not have been allowed when reading the user input.\n";
            abort();
        } else 
        {

            assert(SSP->T56_Initial_AlleleFreqs.size() > patch_index);
            for (int locus = 0 ; locus < SSP->T56_nbLoci ; locus++)
            {
                auto locusInGender = SSP->FromT56LocusToT56genderLocus[locus].second;
                assert(SSP->T56_Initial_AlleleFreqs[patch_index].size() > locus);
                if (SSP->T56_Initial_AlleleFreqs[patch_index][locus] == 1.0)
                {
                    if (SSP->FromT56LocusToT56genderLocus[locus].first)
                    {
                        if (SSP->T56ntrl_compress)
                        {
                            this->setT6ntrl_AlleleToOne_JustPushBack(locusInGender);
                        } else
                        {
                            this->setT5ntrl_AlleleToOne_JustPushBack(locusInGender);
                        }
                    } else
                    {
                        if (SSP->T56sel_compress)
                        {
                            this->setT6sel_AlleleToOne_JustPushBack(locusInGender);
                        } else
                        {
                            this->setT5sel_AlleleToOne_JustPushBack(locusInGender);
                        }
                    }
                    
                } else if (SSP->T56_Initial_AlleleFreqs[patch_index][locus] == 0.0)
                {
                    // Nothing to do
                    
                } else
                {
                    // Find FakeindHaplo_index. FakeindHaplo_index is a function of locus too so as to reduce LD
                    int FakeindHaplo_index = indHaplo_index - (int)(((long)locus*97) % (long)(2 * SSP->patchSize[patch_index]));
                    //                                                                                   ^ 97 is just a large-ish prime number
                    if (FakeindHaplo_index < 0)
                    {
                        FakeindHaplo_index = 2 * SSP->patchSize[patch_index] + FakeindHaplo_index;
                    }
                    
                    assert(FakeindHaplo_index >= 0);
                    assert(2 * SSP->patchSize[patch_index] > FakeindHaplo_index);

                    // Check if mutation must be put
                    int nbMuts = (int)(SSP->T56_Initial_AlleleFreqs[patch_index][locus] * 2.0 * (double) SSP->patchSize[patch_index]);
                    bool setToOne = false;
                    if (FakeindHaplo_index%2 == 0)
                    {
                        if (FakeindHaplo_index <= nbMuts)
                        {
                            setToOne = true;
                        }
                    } else
                    {
                        if (FakeindHaplo_index <= nbMuts - SSP->patchSize[patch_index])
                        {//                                        ^ This is the max number of mutations that have already been added by first (couple of) if statement(s)       
                            setToOne = true;
                        }
                    }


                    // Make mutation if needed
                    if ( setToOne)
                    {
                        if (SSP->FromT56LocusToT56genderLocus[locus].first)
                        {
                            if (SSP->T56ntrl_compress)
                            {
                                this->setT6ntrl_AlleleToOne_JustPushBack(locusInGender);
                            } else
                            {
                                this->setT5ntrl_AlleleToOne_JustPushBack(locusInGender);
                            }
                        } else
                        {
                            if (SSP->T56sel_compress)
                            {
                                this->setT6sel_AlleleToOne_JustPushBack(locusInGender);
                            } else
                            {
                                this->setT5sel_AlleleToOne_JustPushBack(locusInGender);  
                            }
                        }
                    } // nothing to do for else
                }
            }
        }
    }
}

void Haplotype::setEntireT5_Allele(std::vector<unsigned int>& t5a)
{
    assert(t5a.size() == SSP->T5_nbLoci);
    assert(T5ntrl_Alleles.size() == SSP->T5ntrl_nbLoci);
    assert(T5sel_Alleles.size() == SSP->T5sel_nbLoci);

    size_t ntrlLocus = 0;
    size_t selLocus = 0;
    for (size_t locus = 0 ; locus < SSP->T5_nbLoci ; locus++ )
    {
        if (SSP->FromT56ntrlLocusToLocus[ntrlLocus] < SSP->FromT56selLocusToLocus[selLocus])
        {
            this->T5ntrl_Alleles[ntrlLocus] = t5a[locus];
            ntrlLocus++;

        } else
        {
            this->T5sel_Alleles[selLocus] = t5a[locus];
            selLocus++;            
        }
    }
    assert(selLocus == SSP->T5sel_nbLoci);
    assert(ntrlLocus == SSP->T5ntrl_nbLoci);
}


void Haplotype::setEntireT6_Allele(std::vector<unsigned int>& t6a)
{
    assert(t6a.size() == SSP->T6_nbLoci);

    size_t ntrlLocus = 0;
    size_t selLocus = 0;
    for (size_t locus = 0 ; locus < SSP->T6_nbLoci ; locus++ )
    {
        if (SSP->FromT56ntrlLocusToLocus[ntrlLocus] < SSP->FromT56selLocusToLocus[selLocus])
        {
            this->T6ntrl_Alleles.push_back(t6a[locus]);
            ntrlLocus++;

        } else
        {
            this->T6sel_Alleles.push_back(t6a[locus]);
            selLocus++;            
        }
    }
    assert(selLocus == SSP->T6sel_nbLoci);
    assert(ntrlLocus == SSP->T6ntrl_nbLoci);
}

Haplotype::Haplotype(const Haplotype& other)
{
    if (SSP->T1_nbLoci)
    {
        T1_Alleles = other.T1_Alleles;
        if (SSP->T1_isMultiplicitySelection)
        {
            W_T1 = other.W_T1;
        }
    }

    if (SSP->T2_nbLoci)
    {
        T2_Alleles = other.T2_Alleles;
        W_T2 = other.W_T2;
    }

    if (SSP->T3_nbLoci)
    {
        T3_Alleles = other.T3_Alleles;
    }

    if (SSP->T56ntrl_nbLoci)
    {
        if (SSP->T5ntrl_nbLoci)
        {
            T5ntrl_Alleles = other.T5ntrl_Alleles;
        } else
        {
            T6ntrl_Alleles = other.T6ntrl_Alleles;
        }
    }

    if (SSP->T56sel_nbLoci)
    {
        if (SSP->T5sel_nbLoci)
        {
            T5sel_Alleles = other.T5sel_Alleles;
        } else
        {
            T6sel_Alleles = other.T6sel_Alleles;
        }

        if (SSP->T56_isMultiplicitySelection)
        {
            W_T56 = other.W_T56;
        }
    }
}

Haplotype& Haplotype::operator=(const Haplotype& other) // copy assignment operator
{
    //std::cout << "copy assignment haplotype\n";
if (SSP->T1_nbLoci)
    {
        T1_Alleles = other.T1_Alleles;
        if (SSP->T1_isMultiplicitySelection)
        {
            W_T1 = other.W_T1;
        }
    }

    if (SSP->T2_nbLoci)
    {
        T2_Alleles = other.T2_Alleles;
        W_T2 = other.W_T2;
    }

    if (SSP->T3_nbLoci)
    {
        T3_Alleles = other.T3_Alleles;
    }

    if (SSP->T56ntrl_nbLoci)
    {
        if (SSP->T5ntrl_nbLoci)
        {
            T5ntrl_Alleles = other.T5ntrl_Alleles;
        } else
        {
            T6ntrl_Alleles = other.T6ntrl_Alleles;
        }
    }

    if (SSP->T56sel_nbLoci)
    {
        if (SSP->T5sel_nbLoci)
        {
            T5sel_Alleles = other.T5sel_Alleles;
        } else
        {
            T6sel_Alleles = other.T6sel_Alleles;
        }

        if (SSP->T56_isMultiplicitySelection)
        {
            W_T56 = other.W_T56;
        }
    }
    return *this;
}

Haplotype& Haplotype::operator=(Haplotype&& other) // move assignment operator
{
    T1_Alleles = std::move(other.T1_Alleles);
    T2_Alleles = std::move(other.T2_Alleles);
    T3_Alleles = std::move(other.T3_Alleles);
    T5ntrl_Alleles = std::move(other.T5ntrl_Alleles);
    T5sel_Alleles = std::move(other.T5sel_Alleles);
    T6ntrl_Alleles = std::move(other.T6ntrl_Alleles);
    T6sel_Alleles = std::move(other.T6sel_Alleles);
    W_T1 = std::move(other.W_T1);
    W_T2 = std::move(other.W_T2);
    W_T56 = std::move(other.W_T56);
    //assert(W_T1.size() == other.W_T1.size());
    return *this;
}

void Haplotype::swap(Haplotype& other)
{
    T1_Alleles.swap(other.T1_Alleles);
    T2_Alleles.swap(other.T2_Alleles);
    T3_Alleles.swap(other.T3_Alleles);
    T5ntrl_Alleles.swap(other.T5ntrl_Alleles);
    T5sel_Alleles.swap(other.T5sel_Alleles);
    T6ntrl_Alleles.swap(other.T6ntrl_Alleles);
    T6sel_Alleles.swap(other.T6sel_Alleles);
    W_T1.swap(other.W_T1);
    W_T2.swap(other.W_T2);
    W_T56.swap(other.W_T56);
}


void Haplotype::setW_T56(double w, int fitnessMapIndex)
{
    assert(W_T56.size() > fitnessMapIndex && fitnessMapIndex >= 0);
    W_T56[fitnessMapIndex] = w;
}


void Haplotype::setAllW_T56(double w)
{
    assert(W_T56.size() == SSP->NbElementsInFitnessMap);
    for (int fitnessMapIndex = 0 ; fitnessMapIndex < SSP->NbElementsInFitnessMap ; fitnessMapIndex++)
    {
        this->setW_T56(w, fitnessMapIndex);
    }
}


/*bool Haplotype::getT5sel_Allele(const int locus)
{
    // not to be used too often
    return std::binary_search(T5sel_Alleles.begin(), T5sel_Alleles.end(), locus);
}

bool Haplotype::getT5ntrl_Allele(const int locus)
{
    // not to be used too often
    return std::binary_search(T5ntrl_Alleles.begin(), T5ntrl_Alleles.end(), locus);
}*/


void Haplotype::setT5ntrl_Allele(const int& locus, const bool& value)
{
    // I don't like this function but I don't think it is being used in the code anyway!
    auto haploPosition = std::lower_bound(T5ntrl_Alleles.begin(), T5ntrl_Alleles.end(), locus);
    auto flippedPosition = std::lower_bound(SSP->T5ntrl_flipped.begin(), SSP->T5ntrl_flipped.end(), locus);
    bool haploPresent = haploPosition != T5ntrl_Alleles.end() && *haploPosition == locus ? true : false;
    bool flipped = flippedPosition != SSP->T5ntrl_flipped.end() && *flippedPosition == locus ? true : false;

    if (flipped == haploPresent)
    {
        // There is no mut
        if (value)
        {
            T5ntrl_Alleles.insert(haploPosition, locus);
        }
    } else
    {
        // There is a mut
        if (!value)
        {
            T5ntrl_Alleles.erase(haploPosition);   
        }
    }
}

void Haplotype::setT5sel_Allele(const int& locus, const bool& value)
{
    auto position = std::lower_bound(T5sel_Alleles.begin(), T5sel_Alleles.end(), locus);
   
    if (value)
    {
        if (position == T5sel_Alleles.end())
        {
            T5sel_Alleles.push_back(locus);
        } else if (locus != (*position))
        {
            T5sel_Alleles.insert(position, locus);
        }
    } else
    {
        if (position != T5sel_Alleles.end())
        {
            if (locus == (*position))
            {
                T5sel_Alleles.erase(position);
            }
        }
    }
}



void Haplotype::setT6ntrl_Allele(const int& locus, const bool& value)
{
    // I don't like this function but I don't think it is being used in the code anyway!
    auto haploPosition   = T6ntrl_Alleles.lower_bound(locus);
    auto flippedPosition = SSP->T6ntrl_flipped.lower_bound(locus);
    bool haploPresent = haploPosition != T6ntrl_Alleles.end() && *haploPosition == locus ? true : false;
    bool flipped = flippedPosition != SSP->T6ntrl_flipped.end() && *flippedPosition == locus ? true : false;

    if (flipped == haploPresent)
    {
        // There is no mut
        if (value)
        {
            T6ntrl_Alleles.insert(haploPosition, locus);
        }
    } else
    {
        // There is a mut
        if (!value)
        {
            T6ntrl_Alleles.erase(haploPosition);   
        }
    }
}

void Haplotype::setT6sel_Allele(const int& locus, const bool& value)
{
    auto position = std::lower_bound(T5sel_Alleles.begin(), T5sel_Alleles.end(), locus);
   
    if (value)
    {
        if (position == T5sel_Alleles.end())
        {
            T5sel_Alleles.push_back(locus);
        } else if (locus != (*position))
        {
            T5sel_Alleles.insert(position, locus);
        }
    } else
    {
        if (position != T5sel_Alleles.end())
        {
            if (locus == (*position))
            {
                T5sel_Alleles.erase(position);
            }
        }
    }
}


void Haplotype::setT5_Allele(const int& locus, const bool& value)
{
    auto& TXLocusElement = SSP->FromT56LocusToT56genderLocus[locus];
    TXLocusElement.first ? setT5sel_Allele(TXLocusElement.second, value) : setT5ntrl_Allele(TXLocusElement.second, value);
}

void Haplotype::setT6_Allele(const int& locus, const bool& value)
{
    auto& TXLocusElement = SSP->FromT56LocusToT56genderLocus[locus];
    TXLocusElement.first ? setT6sel_Allele(TXLocusElement.second, value) : setT6ntrl_Allele(TXLocusElement.second, value);
}

void Haplotype::setT5sel_AlleleToOne_JustPushBack(unsigned& locus)
{
    T5sel_Alleles.push_back(locus);
}

void Haplotype::setT5sel_AlleleToOne(int& locus)
{
    this->setT5sel_Allele(locus, true);
}

void Haplotype::setT5sel_AlleleToZero(int& locus)
{
    this->setT5sel_Allele(locus, false);
}

void Haplotype::setT5ntrl_AlleleToOne_JustPushBack(unsigned& locus)
{
    T5ntrl_Alleles.push_back(locus);
}

void Haplotype::setT5ntrl_AlleleToOne(int& locus)
{
    this->setT5ntrl_Allele(locus, true);
}


void Haplotype::setT5ntrl_AlleleToZero(int& locus)
{
    this->setT5ntrl_Allele(locus, false);
}



void Haplotype::setT6sel_AlleleToOne_JustPushBack(unsigned& locus)
{
    T6sel_Alleles.push_back(locus);
}

void Haplotype::setT6sel_AlleleToOne(int& locus)
{
    this->setT6sel_Allele(locus, true);
}

void Haplotype::setT6sel_AlleleToZero(int& locus)
{
    this->setT6sel_Allele(locus, false);
}

void Haplotype::setT6ntrl_AlleleToOne_JustPushBack(unsigned& locus)
{
    T6ntrl_Alleles.push_back(locus);
}

void Haplotype::setT6ntrl_AlleleToOne(int& locus)
{
    this->setT6ntrl_Allele(locus, true);
}


void Haplotype::setT6ntrl_AlleleToZero(int& locus)
{
    this->setT6ntrl_Allele(locus, false);
}



double Haplotype::getW_T56(int fitnessMapIndex)
{
    //std::cout << "fitnessMapIndex = " << fitnessMapIndex << "\n";
    //std::cout << "W_T56.size() = " << W_T56.size() << "\n";
    assert(fitnessMapIndex >= 0 && fitnessMapIndex < W_T56.size());
    return W_T56[fitnessMapIndex];
}

template<typename INT>
void Haplotype::toggleT5ntrl_Allele(INT MutPosition)
{
    assert(MutPosition < SSP->T5ntrl_nbLoci);
    auto position = std::lower_bound(T5ntrl_Alleles.begin(), T5ntrl_Alleles.end(), MutPosition);

    if (position == T5ntrl_Alleles.end())
    {
        // not found
        T5ntrl_Alleles.push_back(MutPosition);
    } else
    {
        if ( MutPosition != (*position))
        {
            // not found    
            T5ntrl_Alleles.insert(position, MutPosition);
        } else
        {
            // found
            T5ntrl_Alleles.erase(position);
        }
    }
}

template<typename INT>
void Haplotype::toggleT5sel_Allele(INT MutPosition)
{
    assert(MutPosition < SSP->T5sel_nbLoci);
    auto position = std::lower_bound(T5sel_Alleles.begin(), T5sel_Alleles.end(), MutPosition);

    if (position == T5sel_Alleles.end())
    {
        // not found
        T5sel_Alleles.push_back(MutPosition);
    } else
    {
        if ( MutPosition != (*position))
        {
            // not found    
            T5sel_Alleles.insert(position, MutPosition);
        } else
        {
            // found
            T5sel_Alleles.erase(position);
        }
    }
}



template<typename INT>
void Haplotype::toggleT6ntrl_Allele(INT MutPosition)
{
    assert(MutPosition < SSP->T6ntrl_nbLoci);
    auto position = T6ntrl_Alleles.lower_bound(MutPosition);

    if (position == T6ntrl_Alleles.end())
    {
        // not found
        T6ntrl_Alleles.push_back(MutPosition);
    } else
    {
        if ( MutPosition != (*position))
        {
            // not found    
            T6ntrl_Alleles.insert(position, MutPosition);
        } else
        {
            // found
            T6ntrl_Alleles.erase(position);
        }
    }
}

template<typename INT>
void Haplotype::toggleT6sel_Allele(INT MutPosition)
{
    assert(MutPosition < SSP->T6sel_nbLoci);
    auto position = T6sel_Alleles.lower_bound(MutPosition);

    if (position == T6sel_Alleles.end())
    {
        // not found
        T6sel_Alleles.push_back(MutPosition);
    } else
    {
        if ( MutPosition != (*position))
        {
            // not found    
            T6sel_Alleles.insert(position, MutPosition);
        } else
        {
            // found
            T6sel_Alleles.erase(position);
        }
    }
}

void Haplotype::printT5sel_Alleles()
{ 
    auto it = this->T5sel_AllelesBegin();
    auto itEnd = this->T5sel_AllelesEnd();
    while (it != itEnd)
    {
        std::cout << *it << " ";
    }
    std::cout << "\n";
}

void Haplotype::printT5ntrl_Alleles()
{ 
    auto it = this->T5ntrl_AllelesBegin();
    auto itEnd = this->T5ntrl_AllelesEnd();
    while (it != itEnd)
    {
        std::cout << *it << " ";
    }
    std::cout << "\n";
}

void Haplotype::printT6sel_Alleles()
{ 
    auto it = this->T6sel_AllelesBegin();
    auto itEnd = this->T6sel_AllelesEnd();
    while (it != itEnd)
    {
        std::cout << *it << " ";
    }
    std::cout << "\n";
}

void Haplotype::printT6ntrl_Alleles()
{ 
    auto it = this->T6ntrl_AllelesBegin();
    auto itEnd = this->T6ntrl_AllelesEnd();
    while (it != itEnd)
    {
        std::cout << *it << " ";
    }
    std::cout << "\n";
}


template<typename INT>
void Haplotype::mutateT5sel_Allele(INT MutPosition, int Habitat)
{   
    ///////////////////
    // Find location //
    ///////////////////
    auto haploP = std::lower_bound(T5sel_Alleles.begin(), T5sel_Alleles.end(), MutPosition);

    ///////////////////////////////
    // Toggle from T5sel_Alleles //
    ///////////////////////////////
    bool isNowFoundinT5sel_Alleles = T56finishMutation(haploP, T5sel_Alleles, MutPosition);
    
    /////////////////////
    // Update fitness //
    ////////////////////

    if (SSP->T56_isMultiplicitySelection)
    {
        updateFitnessAfterT56Mutation(MutPosition, isNowFoundinT5sel_Alleles, Habitat);
    }
}


template<typename INT>
void Haplotype::mutateT6sel_Allele(INT MutPosition, int Habitat)
{
    ///////////////////
    // Find location //
    ///////////////////
    auto haploP = T6sel_Alleles.lower_bound(MutPosition);

    ///////////////////////////////
    // Toggle from T6sel_Alleles //
    ///////////////////////////////
    bool isNowFoundinT6sel_Alleles = T56finishMutation(haploP, T6sel_Alleles, MutPosition);
    
    /////////////////////
    // Update fitness //
    ////////////////////
    if (SSP->T56_isMultiplicitySelection)
    {
        updateFitnessAfterT56Mutation(MutPosition, isNowFoundinT6sel_Alleles, Habitat);
    }
}


template<typename INT>
void Haplotype::mutateT56ntrl_Allele(INT MutPosition)
{
    if (SSP->T56ntrl_compress)
    {
        mutateT6ntrl_Allele(MutPosition);
    } else
    {
        mutateT5ntrl_Allele(MutPosition);
    }
}


template<typename INT>
void Haplotype::mutateT56sel_Allele(INT MutPosition, int Habitat)
{
    if (SSP->T56sel_compress)
    {
        mutateT6sel_Allele(MutPosition, Habitat);
    } else
    {
        mutateT5sel_Allele(MutPosition, Habitat);
    }
}


template<typename INT>
void Haplotype::mutateT5ntrl_Allele(INT MutPosition)
{
    auto position = std::lower_bound(T5ntrl_Alleles.begin(), T5ntrl_Alleles.end(), MutPosition);
    
    if (position == T5ntrl_Alleles.end())
    {
        // not found
        T5ntrl_Alleles.push_back(MutPosition);  
    } else
    {
        if ( MutPosition == (*position))
        {
            // found
            T5ntrl_Alleles.erase(position);
        } else
        {
            // not found
            T5ntrl_Alleles.insert(position, MutPosition);
        }
    }
}

template<typename INT>
void Haplotype::mutateT5ntrl_Allele(std::vector<INT>& MutPositions)
{    
    // MutPositions is already sorted and without duplicates
    size_t prevSize = T5ntrl_Alleles.size();
    T5ntrl_Alleles.insert(T5ntrl_Alleles.end(), MutPositions.begin(), MutPositions.end());
    std::inplace_merge(T5ntrl_Alleles.begin(), T5ntrl_Alleles.begin() + prevSize, T5ntrl_Alleles.end());
    T5ntrl_Alleles.erase(std::unique(T5ntrl_Alleles.begin(), T5ntrl_Alleles.end()), T5ntrl_Alleles.end());
}

template<typename INT>
void Haplotype::mutateT6ntrl_Allele(INT MutPosition)
{
    //T6ntrl_Alleles.assertOrdering("From Haplotype::mutateT6ntrl_Allele, beginning");
    auto position = T6ntrl_Alleles.lower_bound(MutPosition);
    
    if (position == T6ntrl_Alleles.end())
    {
        // not found
        /*if (T6ntrl_Alleles.size())
        {
            if (T6ntrl_Alleles.back() >= MutPosition)
            {
                std::cout << "T6ntrl_Alleles.back() = " << T6ntrl_Alleles.back() << "\n";
                std::cout << "MutPosition = " << MutPosition << "\n";
                std::cout << "T6ntrl_Alleles contains: ";
                for (auto it = T6ntrl_Alleles.begin() ; it < T6ntrl_Alleles.end() ; ++it)
                    std::cout << *it <<  " ";
                std::cout << std::endl;
                abort();
            }
        }*/
        T6ntrl_Alleles.push_back(MutPosition);  
        //T6ntrl_Alleles.assertOrdering("From Haplotype::mutateT6ntrl_Allele. pushed_back");
    } else
    {
        if ( MutPosition == (*position))
        {
            // found
            T6ntrl_Alleles.erase(position);
            //T6ntrl_Alleles.assertOrdering("From Haplotype::mutateT6ntrl_Allele. Erased");
        } else
        {
            // not found
            T6ntrl_Alleles.insert(position, MutPosition);
            //T6ntrl_Alleles.assertOrdering("From Haplotype::mutateT6ntrl_Allele. Inserted");
        }
    }

}

void Haplotype::toggleFromT5ntrl_Allele(int& MutPosition, unsigned int& from)
{
    auto haploP = std::lower_bound(T5ntrl_Alleles.begin() + from, T5ntrl_Alleles.end(), MutPosition);
    from = haploP - T5ntrl_Alleles.begin();
    
    (void) T56finishMutation(haploP, T5ntrl_Alleles, MutPosition);
}

void Haplotype::toggleFromT6ntrl_Allele(int& MutPosition, unsigned int& blockIndexFrom, unsigned int& fromInBlock)
{
    assert(MutPosition >= 0 && MutPosition < SSP->T6ntrl_nbLoci);
    CompressedSortedDeque::iterator haploP = T6ntrl_Alleles.lower_bound_from(MutPosition, blockIndexFrom, fromInBlock); // blockIndexFrom and fromInBlock updated in here
    (void) T56finishMutation(haploP, T6ntrl_Alleles, MutPosition);
}


void Haplotype::toggleFromT5sel_Allele(int& MutPosition, unsigned int& from, int Habitat)
{
    auto haploP = std::lower_bound(T5sel_Alleles.begin() + from, T5sel_Alleles.end(), MutPosition);
    from = haploP - T5sel_Alleles.cbegin();

    if (SSP->T56_isMultiplicitySelection)
        updateFitnessAfterT56Mutation(MutPosition, T56finishMutation(haploP, T5sel_Alleles, MutPosition), Habitat);

}

void Haplotype::toggleFromT6sel_Allele(int& MutPosition, unsigned int& blockIndexFrom, unsigned int& fromInBlock, int Habitat)
{
    auto haploP = T6sel_Alleles.lower_bound_from(MutPosition, blockIndexFrom, fromInBlock); // blockIndexFrom and fromInBlock updated in lower_bound
    if (SSP->T56_isMultiplicitySelection)
        updateFitnessAfterT56Mutation(MutPosition, T56finishMutation(haploP, T6sel_Alleles, MutPosition), Habitat);
}



void Haplotype::clearT56Alleles()
{
    if (SSP->T56ntrl_nbLoci)
    {
        if (SSP->T56ntrl_compress)
        {
            T6ntrl_Alleles.clear();
        } else
        {
            T5ntrl_Alleles.clear();
        }
    }
        
    if (SSP->T56sel_nbLoci)
    {
        if (SSP->T56sel_compress)
        {
            T6sel_Alleles.clear();
        } else
        {
            T5sel_Alleles.clear();
        }
    }
}

void Haplotype::copyIntoT56sel(int from, int to, Haplotype& SourceChromo)
{
    // Must have been emptied first

    assert(from >= 0);
    assert(to > from);
    if (SSP->T56sel_compress)
    {
        T6sel_Alleles.extend(from, to, SourceChromo.T6sel_Alleles);
    } else
    {
        std::vector<unsigned int>::iterator itFrom;
        if (from == 0)
        {
            itFrom = SourceChromo.T5sel_Alleles.begin();
        } else
        {
            itFrom = lower_bound(SourceChromo.T5sel_Alleles.begin(), SourceChromo.T5sel_Alleles.end(), from);
        }
        
        std::vector<unsigned int>::iterator itTo;
        if (to == SSP->T5sel_nbLoci)
        {
            itTo   = SourceChromo.T5sel_Alleles.end();
        } else
        {
            itTo   = lower_bound(itFrom, SourceChromo.T5sel_Alleles.end(), to);
        }
        
        
        T5sel_Alleles.insert(T5sel_Alleles.end(), itFrom, itTo);
    }   
}

void Haplotype::copyIntoT56ntrl(int from, int to, Haplotype& SourceChromo)
{
    // Must have been emptied first

    assert(from >= 0);
    assert(to > from);
    if (SSP->T56ntrl_compress)
    {
        T6ntrl_Alleles.extend(from, to, SourceChromo.T6ntrl_Alleles);
        
    } else
    {
        std::vector<unsigned int>::iterator itFrom;
        if (from == 0)
        {
            itFrom = SourceChromo.T5ntrl_Alleles.begin();
        } else
        {
            itFrom = lower_bound(SourceChromo.T5ntrl_Alleles.begin(), SourceChromo.T5ntrl_Alleles.end(), from);
        }
        
        std::vector<unsigned int>::iterator itTo;
        if (to == SSP->T5ntrl_nbLoci)
        {
            itTo   = SourceChromo.T5ntrl_Alleles.end();
        } else
        {
            itTo   = lower_bound(itFrom, SourceChromo.T5ntrl_Alleles.end(), to);
        }
        
        T5ntrl_Alleles.insert(T5ntrl_Alleles.end(), itFrom, itTo);
    }
}

/*
std::vector<unsigned int>::const_iterator Haplotype::T5sel_AllelesCBegin()
{ 
    return T5sel_Alleles.cbegin();
}

std::vector<unsigned int>::const_iterator Haplotype::T5sel_AllelesCEnd()
{
    return T5sel_Alleles.cend();
}

std::vector<unsigned int>::const_iterator Haplotype::T5sel_AllelesCiterator(int locus, std::vector<unsigned int>::const_iterator from)
{
    
    if (from == T5sel_Alleles.cend())
    {
        return from;
    }

    assert(from <= T5sel_Alleles.cend() );
    assert(from >= T5sel_Alleles.cbegin() );

    if (locus == *from)
    {
        return from;
    }

    if (locus == SSP->T5sel_nbLoci)
    {
        return T5sel_Alleles.cend();
    }

    return std::lower_bound(from, T5sel_Alleles.cend(), locus);
}


std::vector<unsigned int>::const_iterator Haplotype::T5sel_AllelesCiterator(int locus)
{
    if (locus == 0)
    {
        return T5sel_Alleles.cbegin();
    } else if (locus == SSP->T5sel_nbLoci)
    {
        return T5sel_Alleles.cend();
    }
    return std::lower_bound(T5sel_Alleles.cbegin(), T5sel_Alleles.cend(), locus);
}





std::vector<unsigned int>::const_iterator Haplotype::T5ntrl_AllelesCBegin()
{ 
    return T5ntrl_Alleles.cbegin();
}

std::vector<unsigned int>::const_iterator Haplotype::T5ntrl_AllelesCEnd()
{
    return T5ntrl_Alleles.cend();
}

std::vector<unsigned int>::const_iterator Haplotype::T5ntrl_AllelesCiterator(int locus, std::vector<unsigned int>::const_iterator from)
{
    assert(from <= T5ntrl_Alleles.cend() );
    assert(from >= T5ntrl_Alleles.cbegin() );

    if (from == T5ntrl_Alleles.cend())
    {
        return from;
    }

    if (locus == *from)
    {
        return from;
    }

    if (locus == SSP->T5ntrl_nbLoci)
    {
        return T5ntrl_Alleles.cend();
    }

    return std::lower_bound(from, T5ntrl_Alleles.cend(), locus);
}


std::vector<unsigned int>::const_iterator Haplotype::T5ntrl_AllelesCiterator(int locus)
{
    if (locus == 0)
    {
        return T5ntrl_Alleles.cbegin();
    } else if (locus == SSP->T5ntrl_nbLoci)
    {
        return T5ntrl_Alleles.cend();
    }
    return std::lower_bound(T5ntrl_Alleles.cbegin(), T5ntrl_Alleles.cend(), locus);
}
*/
/*int Haplotype::T5_AllelesPosition(int locus, int from)
{
    if (from == T5_Alleles.size())
    {
        return from;   
    }

    if (locus == T5_Alleles[from])
    {
        return from;
    }

    if (locus == SSP->T5_nbLoci)
    {
        return T5_Alleles.size();
    }

    return std::lower_bound(T5_Alleles.begin() + from, T5_Alleles.end(), locus) - T5_Alleles.begin();
}*/

/*int Haplotype::T5_AllelesPosition(int locus)
{
    if (locus == 0)
    {
        return 0;
    } else if (locus == SSP->T5_nbLoci)
    {
        return T5_Alleles.size();
    }
    return std::lower_bound(T5_Alleles.begin(), T5_Alleles.end(), locus) - T5_Alleles.begin();
}*/


/*int Haplotype::T5_howManyMutations()
{
    return T5sel_Alleles.size() + T5ntrl_Alleles.size();
}*/


/*int Haplotype::getT5ntrl_nthMutation(const int n)
{
    return T5ntrl_Alleles[n];
}

int Haplotype::getT5sel_nthMutation(const int n)
{
    return T5sel_Alleles[n];
}*/

void Haplotype::toggleT56ntrlLoci(std::vector<int>& lociToToggle)
{
    // lociToToggle must be sorted and all must be neutral

    if (SSP->T56ntrl_compress)
    {
        unsigned blockIndexFrom = 0;
        unsigned fromInBlock = 0;
        for (auto& locus : lociToToggle)
        {
            toggleFromT6ntrl_Allele(locus, blockIndexFrom,  fromInBlock); // fromInBlock and blockIndexFrom updated in 'toggleT5sel_Allele'
        }    
    } else
    {
        unsigned from = 0;
        for (auto& locus : lociToToggle)
        {
            toggleFromT5ntrl_Allele(locus, from); // from updated in 'toggleT5sel_Allele'
        }    
    }
        
}

/*void Haplotype::toggleT56selLoci(std::vector<int>& lociToToggle, int Habitat)
{
    if (SSP->T56sel_compress)
    {
        unsigned blockIndexFrom = 0;
        unsigned fromInBlock = 0;
        for (auto& locus : lociToToggle)
        {
            toggleFromT6sel_Allele(locus, blockIndexFrom,  fromInBlock, Habitat); // fromInBlock and blockIndexFrom updated in 'toggleT5sel_Allele'
        }    
    } else
    {
        unsigned from = 0;
        for (auto& locus : lociToToggle)
        {
            toggleFromT5sel_Allele(locus, from, Habitat); // from updated in 'toggleT5sel_Allele'
        }    
    }        
}*/

int Haplotype::getNbT5ntrl()
{
    return T5ntrl_Alleles.size();
}

int Haplotype::getNbT6ntrl()
{
    return T6ntrl_Alleles.size();
}

bool Haplotype::isT5ntrlMutation(int locus)
{
    auto it = this->T5ntrl_AllelesBegin();
    auto itEnd = this->T5ntrl_AllelesEnd();
    it.lower_bound(locus);
    if (it != itEnd)
    {
        return *it == locus;
    } else
    {
        return false;
    }
}


bool Haplotype::isT5selMutation(int locus)
{
    std::vector<unsigned int>::iterator itEnd = this->T5sel_AllelesEnd();
    auto it = std::lower_bound(this->T5sel_AllelesBegin(), itEnd, locus);
    if (it != itEnd)
    {
        return *it == locus;
    } else
    {
        return false;
    }
}

bool Haplotype::isT6ntrlMutation(int locus)
{
    auto it = this->T6ntrl_AllelesBegin();
    auto itEnd = this->T6ntrl_AllelesEnd();
    it.lower_bound(locus);
    if (it != itEnd)
    {
        return *it == locus;
    } else
    {
        return false;
    }
}


bool Haplotype::isT6selMutation(int locus)
{
    auto it = T6sel_Alleles.lower_bound(locus);
    auto itEnd = this->T6sel_AllelesEnd();
    if (it != itEnd)
    {
        return *it == locus;
    } else
    {
        return false;
    }
}


void Haplotype::assertT5orderAndUniqueness()
{
    {
        auto& v = T5ntrl_Alleles;
        size_t originalSize = v.size();
        v.erase( std::unique( v.begin(), v.end() ), v.end() );
        assert(v.size() == originalSize);
        assert(std::is_sorted(v.begin(), v.end()));
    }

    {
        auto& v = T5sel_Alleles;
        size_t originalSize = v.size();
        v.erase( std::unique( v.begin(), v.end() ), v.end() );
        assert(v.size() == originalSize);
        assert(std::is_sorted(v.begin(), v.end()));
    }
}



ZipIterator<std::vector<unsigned int>, std::vector<unsigned int>::iterator> Haplotype::T5ntrl_AllelesBegin()
{
    return ZipIterator<std::vector<unsigned int>, std::vector<unsigned int>::iterator>(
        T5ntrl_Alleles.begin(),
        SSP->T5ntrl_flipped.begin(),
        T5ntrl_Alleles.end(),
        SSP->T5ntrl_flipped.end()
    );
}

ZipIterator<std::vector<unsigned int>, std::vector<unsigned int>::iterator> Haplotype::T5ntrl_AllelesEnd()
{
    return ZipIterator<std::vector<unsigned int>, std::vector<unsigned int>::iterator>(
        T5ntrl_Alleles.end(),
        SSP->T5ntrl_flipped.end(),
        T5ntrl_Alleles.end(),
        SSP->T5ntrl_flipped.end()
    );
}

std::vector<unsigned int>::iterator Haplotype::T5sel_AllelesBegin()
{
    return T5sel_Alleles.begin();
}

std::vector<unsigned int>::iterator Haplotype::T5sel_AllelesEnd()
{
    return T5sel_Alleles.end();
}


ZipIterator<CompressedSortedDeque, CompressedSortedDeque::iterator> Haplotype::T6ntrl_AllelesBegin()
{
    return ZipIterator<CompressedSortedDeque, CompressedSortedDeque::iterator>(
        T6ntrl_Alleles.begin(),
        SSP->T6ntrl_flipped.begin(),
        T6ntrl_Alleles.end(),
        SSP->T6ntrl_flipped.end()
    );
}

ZipIterator<CompressedSortedDeque, CompressedSortedDeque::iterator> Haplotype::T6ntrl_AllelesEnd()
{
    return ZipIterator<CompressedSortedDeque, CompressedSortedDeque::iterator>(
        T6ntrl_Alleles.end(),
        SSP->T6ntrl_flipped.end(),
        T6ntrl_Alleles.end(),
        SSP->T6ntrl_flipped.end()
    );
}

CompressedSortedDeque::iterator Haplotype::T6sel_AllelesBegin()
{
    return T6sel_Alleles.begin();
}

CompressedSortedDeque::iterator Haplotype::T6sel_AllelesEnd()
{
    return T6sel_Alleles.end();
}



ZipIterator<std::vector<unsigned int>, std::vector<unsigned int>::iterator> Haplotype::T56ntrl_AllelesBegin(ZipIterator<std::vector<unsigned int>, std::vector<unsigned int>::iterator> nothing)
{
    (void) nothing;
    return T5ntrl_AllelesBegin();
}

ZipIterator<std::vector<unsigned int>, std::vector<unsigned int>::iterator> Haplotype::T56ntrl_AllelesEnd(ZipIterator<std::vector<unsigned int>, std::vector<unsigned int>::iterator> nothing)
{
    (void) nothing;
    return T5ntrl_AllelesEnd();
}


std::vector<unsigned int>::iterator Haplotype::T56sel_AllelesBegin(std::vector<unsigned int>::iterator& nothing)
{
    (void) nothing;
    return T5sel_AllelesBegin();
}

std::vector<unsigned int>::iterator Haplotype::T56sel_AllelesEnd(std::vector<unsigned int>::iterator& nothing)
{
    (void) nothing;
    return T5sel_AllelesEnd();
}

ZipIterator<CompressedSortedDeque, CompressedSortedDeque::iterator> Haplotype::T56ntrl_AllelesBegin(ZipIterator<CompressedSortedDeque, CompressedSortedDeque::iterator> nothing)
{
    (void) nothing;
    return T6ntrl_AllelesBegin();
}

ZipIterator<CompressedSortedDeque, CompressedSortedDeque::iterator> Haplotype::T56ntrl_AllelesEnd(ZipIterator<CompressedSortedDeque, CompressedSortedDeque::iterator> nothing)
{
    (void) nothing;
    return T6ntrl_AllelesEnd();
}


CompressedSortedDeque::iterator Haplotype::T56sel_AllelesBegin(CompressedSortedDeque::iterator& nothing)
{
    (void) nothing;
    return T6sel_AllelesBegin();
}

CompressedSortedDeque::iterator Haplotype::T56sel_AllelesEnd(CompressedSortedDeque::iterator& nothing)
{
    (void) nothing;
    return T6sel_AllelesEnd();
}


ZipIterator<std::vector<unsigned int>, std::vector<unsigned int>::iterator> Haplotype::T56ntrl_AllelesIterator(ZipIterator<std::vector<unsigned int>, std::vector<unsigned int>::iterator> nothing, unsigned value)
{
    (void) nothing;
    return ZipIterator<std::vector<unsigned int>, std::vector<unsigned int>::iterator>(
        std::lower_bound(T5ntrl_Alleles.begin(), T5ntrl_Alleles.end(), value),
        std::lower_bound(SSP->T5ntrl_flipped.begin(), SSP->T5ntrl_flipped.end(), value),
        T5ntrl_Alleles.end(),
        SSP->T5ntrl_flipped.end()
    );
    
}

std::vector<unsigned int>::iterator Haplotype::T56sel_AllelesIterator(std::vector<unsigned int>::iterator& nothing, unsigned value)
{
    (void) nothing;
    return std::lower_bound(T5sel_Alleles.begin(), T5sel_Alleles.end(), value);
}


ZipIterator<CompressedSortedDeque, CompressedSortedDeque::iterator> Haplotype::T56ntrl_AllelesIterator(ZipIterator<CompressedSortedDeque, CompressedSortedDeque::iterator> nothing, unsigned value)
{
    (void) nothing;

    return ZipIterator<CompressedSortedDeque, CompressedSortedDeque::iterator>(
        T6ntrl_Alleles.lower_bound(value),
        SSP->T6ntrl_flipped.lower_bound(value),
        T6ntrl_Alleles.end(),
        SSP->T6ntrl_flipped.end()
    );
}

CompressedSortedDeque::iterator Haplotype::T56sel_AllelesIterator(CompressedSortedDeque::iterator& nothing, unsigned value)
{
    (void) nothing;
    return T6sel_Alleles.lower_bound(value);
}




std::vector<unsigned int> Haplotype::getT56ntrlTrueHaplotype()
{
    std::vector<unsigned int> r;
    if (SSP->T56ntrl_nbLoci)
    {
        if (SSP->T56ntrl_compress)
        {
            auto it = this->T6ntrl_AllelesBegin();
            auto itEnd = this->T6ntrl_AllelesEnd();
            for (; it != itEnd ; ++it)
            {
                r.push_back(*it);
            }
        } else
        {
            auto it = this->T5ntrl_AllelesBegin();
            auto itEnd = this->T5ntrl_AllelesEnd();
            for (; it != itEnd ; ++it)
            {
                r.push_back(*it);
            }
        }
    }
    
    return r;
}

std::vector<unsigned int> Haplotype::getT56selTrueHaplotype()
{
    std::vector<unsigned int> r;
    if (SSP->T56sel_nbLoci)
    {
        if (SSP->T56sel_compress)
        {
            auto it = this->T6sel_AllelesBegin();
            auto itEnd = this->T6sel_AllelesEnd();
            for (; it != itEnd ; ++it)
            {
                r.push_back(*it);
            } 
        } else
        {
            auto it = this->T5sel_AllelesBegin();
            auto itEnd = this->T5sel_AllelesEnd();
            for (; it != itEnd ; ++it)
            {
                r.push_back(*it);
            }
        }
    }
    
    return r;
}


template<typename ITERATOR, typename CONTAINER, typename INT>
bool Haplotype::T56finishMutation(ITERATOR& haploP, CONTAINER& container, INT MutPosition)
{
    bool isNowFoundinContainer;
    if (haploP != container.end())
    {
        if ( MutPosition == *haploP)
        {
            // found
            container.erase(haploP);
            isNowFoundinContainer=false;
        } else
        {   
            container.insert(haploP, MutPosition);
            isNowFoundinContainer=true;
        }
    } else
    {
        // Not found
        container.push_back(MutPosition);
        isNowFoundinContainer=true;
    }

    return isNowFoundinContainer;
}

template<typename INT>
void Haplotype::updateFitnessAfterT56Mutation(INT MutPosition, bool isNowFoundInAlleles, int Habitat)
{
    assert(SSP->FromT56selLocusToLocus.size()>MutPosition);
    int fitnessMapIndex = SSP->FromLocusToFitnessMapIndex[SSP->FromT56selLocusToLocus[MutPosition]];

    double w = this->getW_T56(fitnessMapIndex);
    
    if (w != -1.0)
    {
        /*// Was mutation added or removed?
        bool wasMutationAdded;
        if (SSP->T56sel_compress)
        {
            auto flippedP = SSP->T6sel_flipped.lower_bound(MutPosition);
            if (flippedP != SSP->T6sel_flipped.end() && MutPosition == *flippedP)
            {
                // Found
                wasMutationAdded = !isNowFoundInAlleles;
            } else
            {
                // Not found
                wasMutationAdded = isNowFoundInAlleles;
            }
        } else
        {
            auto flippedP = std::lower_bound(SSP->T5sel_flipped.begin(), SSP->T5sel_flipped.end(), MutPosition);
            if (flippedP != SSP->T5sel_flipped.end() && MutPosition == *flippedP)
            {
                // Found
                wasMutationAdded = !isNowFoundInAlleles;
            } else
            {
                // Not found
                wasMutationAdded = isNowFoundInAlleles;
            }
        }*/
            

        // Update fitness
        if ( isNowFoundInAlleles )
        {
            w *= SSP->T56_FitnessEffects[Habitat][MutPosition];


            assert(w >= 0.0 && w <= 1.0);
        } else
        {
            if (SSP->T56_FitnessEffects[Habitat][MutPosition] == 0.0)
            {
                w = -1.0;
            } else
            {
                w /= SSP->T56_FitnessEffects[Habitat][MutPosition];
                assert(w >= 0.0 && w <= 1.0);
            }
        }
        this->setW_T56( w,   fitnessMapIndex );
    }
}

size_t Haplotype::nbT56muts()
{
    size_t r;

    if (SSP->T56ntrl_compress)
    {
        r = T6ntrl_Alleles.size();
    } else
    {
        r = T5ntrl_Alleles.size();
    }

    if (SSP->T56sel_compress)
    {
        r += T6sel_Alleles.size();
    } else
    {
        r += T5sel_Alleles.size();
    }

    return r;
}

size_t Haplotype::nbT56muts(int fitnessMapIndex)
{
    size_t r = 0;

    if (SSP->T56ntrl_nbLoci)
    {
        if (SSP->T56ntrl_compress)
        {
            std::cout << "Sorry debugging function 'size_t Haplotype::nbT56muts(int fitnessMapIndex)' has not been defined for T6 loci\n";
            abort();
        } else
        {
            if (T5ntrl_Alleles.size() != 0)
            {
                auto from = std::lower_bound(SSP->FromLocusToFitnessMapIndex.begin(), SSP->FromLocusToFitnessMapIndex.end(), fitnessMapIndex);
                auto to   = std::upper_bound(from, SSP->FromLocusToFitnessMapIndex.end(), fitnessMapIndex);
                size_t locusFrom = from - SSP->FromLocusToFitnessMapIndex.begin();
                size_t locusTo = to - SSP->FromLocusToFitnessMapIndex.begin();
                r += 
                    std::lower_bound(T5ntrl_Alleles.begin(), T5ntrl_Alleles.end(), locusTo)
                    -
                    std::lower_bound(T5ntrl_Alleles.begin(), T5ntrl_Alleles.end(), locusFrom)
                ;
            }
        }
    }


    if (SSP->T56sel_nbLoci)
    {
        if (SSP->T56sel_compress)
        {
            std::cout << "Sorry debugging function 'size_t Haplotype::nbT56muts(int fitnessMapIndex)' has not been defined for T6 loci\n";
            abort();
        } else
        {
            if (T5sel_Alleles.size() != 0)
            {
                auto from = std::lower_bound(SSP->FromLocusToFitnessMapIndex.begin(), SSP->FromLocusToFitnessMapIndex.end(), fitnessMapIndex);
                auto to   = std::upper_bound(from, SSP->FromLocusToFitnessMapIndex.end(), fitnessMapIndex);
                size_t locusFrom = from - SSP->FromLocusToFitnessMapIndex.begin();
                size_t locusTo = to - SSP->FromLocusToFitnessMapIndex.begin();
                //std::cout << "fitnessMapIndex = " << fitnessMapIndex << " | locusFrom = " << locusFrom << " | locusTo = " << locusTo << "\n";
                r += 
                    std::lower_bound(T5sel_Alleles.begin(), T5sel_Alleles.end(), locusTo)
                    -
                    std::lower_bound(T5sel_Alleles.begin(), T5sel_Alleles.end(), locusFrom)
                ;
            }
        }
    }

    return r;
}








double Haplotype::CalculateT1FitnessMultiplicity(const int& Habitat)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'CalculateT1FitnessMultiplicity'\n";
#endif   
    assert(SSP->T1_isMultiplicitySelection);

    // Get the right fitness array
    auto& fits = SSP->T1_FitnessEffects[Habitat];
    assert(fits.size() == SSP->T1_nbLoci);

    double r = 1.0;
    for (size_t fitnessMapIndex = 0 ; fitnessMapIndex < SSP->NbElementsInFitnessMap ; ++fitnessMapIndex)
    {
        if (getW_T1(fitnessMapIndex) == -1.0)
        {
            //std::cout << "must recompute fitnessMapIndex " << fitnessMapIndex << "\n";
            double w = 1.0;

            auto to = SSP->FromFitnessMapIndexToTXLocus[fitnessMapIndex].T1;
            //std::cout << "to " << to << "\n";
            assert(to <= SSP->T1_nbLoci);
            auto from = 0;
            if (fitnessMapIndex != 0)
            {
                from = SSP->FromFitnessMapIndexToTXLocus[fitnessMapIndex - 1].T1;
            }

            auto byteFrom = from / 8;
            auto bitFrom = from % 8;
            auto byteTo = to / 8;
            auto bitTo = to % 8;
            auto locus = from;
            /*
            std::cout << "from " << from << "\n";
            std::cout << "to " << to << "\n";
            std::cout << "byteFrom " << byteFrom << "\n";
            std::cout << "bitFrom " << bitFrom << "\n";
            std::cout << "byteTo " << byteTo << "\n";
            std::cout << "bitTo " << bitTo << "\n";
            */

            if (byteFrom == byteTo)
            {
                // Special case where the fitness map block is so small that it is contain within a single byte!
                for (auto bit = bitFrom ; bit < bitTo ; ++bit)
                {
                    if (getT1_Allele(byteFrom, bit))
                    {
                        w *= fits[locus];
                    }
                    ++locus;
                }
            } else
            {
                // first bits
                if (bitFrom != 0)
                {
                    if (getT1_char(byteFrom))
                    {
                        for (auto bit = bitFrom ; bit < 8 ; ++bit)
                        {
                            if (getT1_Allele(byteFrom, bit))
                            {
                                w *= fits[locus];
                            }
                            ++locus;
                        }
                    } else
                    {
                        locus += 8 - bitFrom; // that might be wrong when byteFrom == byteTo but it does not matter for this case
                    }
                    ++byteFrom;
                }
                //std::cout << "after first bits: locus " << locus << "\n";

                // First Bytes
                for (auto byte = byteFrom ; byte < byteTo ; ++byte)
                {
                    if (getT1_char(byte))
                    {
                        for (auto bit = 0 ; bit < 8 ; ++bit)
                        {
                            if (getT1_Allele(byte, bit))
                            {
                                w *= fits[locus];
                            }
                            ++locus;
                        }    
                    } else
                    {
                        locus += 8;
                    }
                }
                //std::cout << "after first bytes: locus " << locus << "\n";

                // Last byte
                if (bitTo != 0)
                {
                    if (getT1_char(byteTo))
                    {
                        for (auto bit = 0 ; bit < bitTo ; ++bit)
                        {
                            if (getT1_Allele(byteTo, bit))
                            {
                                w *= fits[locus];
                            }
                            ++locus;
                        }
                    } else
                    {
                        locus += bitTo; // only needed for assertions below
                    }
                } 
                //std::cout << "after last byte: locus " << locus << "\n";
            }

            assert(locus == to);


            /*for (auto& polymorphicLocus : SSP->simTracker.T1SitesForFitnessMultiplicityCalculation[fitnessMapIndex])
            {
                if (getT1_Allele(polymorphicLocus.byte_index, polymorphicLocus.bit_index))
                {
                    w *= fits[polymorphicLocus.locus];
                }
                
            }*/
            
            assert(w<=1.0 && w>=0.0);
            setW_T1(w, fitnessMapIndex);
            r *= w;
        } else
        {
            r *= getW_T1(fitnessMapIndex);
        }
    }

    assert(r<=1.0 && r>=0.0);
    return r;
}


double Haplotype::CalculateT2Fitness(const int& Habitat)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'CalculateT2Fitness'\n";
#endif          
    assert(SSP->T2_FitnessEffects.size() > Habitat);
    auto& fits = SSP->T2_FitnessEffects[Habitat];
    assert(fits.size() == SSP->T2_nbLoci);


    double r = 1.0;
    int T2_locus = 0;
    for (size_t fitnessMapIndex = 0 ; fitnessMapIndex < SSP->NbElementsInFitnessMap ; ++fitnessMapIndex)
    {
        int T2_locusTo = SSP->FromFitnessMapIndexToTXLocus[fitnessMapIndex].T2;
        if (getW_T2(fitnessMapIndex) == -1.0)
        {
            double w = 1.0;
            for ( ; T2_locus <= T2_locusTo ; ++T2_locus )
            {
                //assert(SSP->T2_FitnessEffects.size() > Habitat);
                //assert(SSP->T2_FitnessEffects[Habitat].size() > char_index);
                w *= pow(fits[T2_locus], getT2_Allele(T2_locus));
            }
            assert(w<=1.0 && w>=0.0);
            setW_T2(w, fitnessMapIndex);
            r *= w;
            //if (W_T2haplo0!=1.0) std::cout << "W_T2haplo0 = " << W_T2haplo0 << std::endl;

        } else
        {
            r *= getW_T2(fitnessMapIndex);
            T2_locus = T2_locusTo;
        }          
    }
    
    assert(r >= 0.0 && r <= 1.0);
    return r;
}


double Haplotype::CalculateT5FitnessMultiplicity(const int& Habitat)
{

    auto& fits = SSP->T56_FitnessEffects[Habitat];
    assert(fits.size() == SSP->T56sel_nbLoci);

    std::vector<unsigned int>::iterator it    = T5sel_Alleles.begin();
    std::vector<unsigned int>::iterator itEnd = T5sel_Alleles.end();

    bool needToRecomputeIt = false;

    double r = 1.0;
    for (size_t fitnessMapIndex = 0 ; fitnessMapIndex < SSP->NbElementsInFitnessMap ; ++fitnessMapIndex)
    {
        if (getW_T56(fitnessMapIndex) != -1.0)
        {
            r *= getW_T56(fitnessMapIndex);
            needToRecomputeIt = true;
        } else
        {
            if (needToRecomputeIt)
            {
                auto locusFrom = SSP->FromFitnessMapIndexToTXLocus[fitnessMapIndex-1].T56sel;
                it = std::lower_bound(it, itEnd, locusFrom);
            }

            unsigned T56_locusTo = SSP->FromFitnessMapIndexToTXLocus[fitnessMapIndex].T56sel;
            double w = 1.0;

            for (; it != itEnd && *it < T56_locusTo ; ++it)
            {                
                w *= fits[*it];
            }
            needToRecomputeIt = false;
            
            assert(w<=1.0 && w>=0.0);
            setW_T56(w, fitnessMapIndex);
            r *= w;
        }
    }
    return r;
}

double Haplotype::CalculateT6FitnessMultiplicity(const int& Habitat)
{

    auto& fits = SSP->T56_FitnessEffects[Habitat];
    assert(fits.size() == SSP->T56sel_nbLoci);

    CompressedSortedDeque::iterator it    = T6sel_Alleles.begin();
    CompressedSortedDeque::iterator itEnd = T6sel_Alleles.end();

    bool needToRecomputeIt = false;

    double r = 1.0;
    for (size_t fitnessMapIndex = 0 ; fitnessMapIndex < SSP->NbElementsInFitnessMap ; ++fitnessMapIndex)
    {
        if (getW_T56(fitnessMapIndex) != -1.0)
        {
            r *= getW_T56(fitnessMapIndex);
            needToRecomputeIt = true;
        } else
        {
            if (needToRecomputeIt)
            {
                auto locusFrom = SSP->FromFitnessMapIndexToTXLocus[fitnessMapIndex-1].T56sel;
                it = it.lower_bound_FromThis(locusFrom);
            }

            unsigned T56_locusTo = SSP->FromFitnessMapIndexToTXLocus[fitnessMapIndex].T56sel;
            double w = 1.0;

            for (; it != itEnd && *it < T56_locusTo ; ++it)
            {                
                w *= fits[*it];
            }
            needToRecomputeIt = false;
            
            assert(w<=1.0 && w>=0.0);
            setW_T56(w, fitnessMapIndex);
            r *= w;
        }
    }
    return r;
}



double Haplotype::CalculateT56FitnessMultiplicity(const int& Habitat)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'CalculateT56FitnessMultiplicity'\n";
#endif   


    double r;
    if (SSP->T5sel_nbLoci)
    {
        r = CalculateT5FitnessMultiplicity(Habitat);
    } else
    {
        assert(SSP->T6sel_nbLoci);
        r = CalculateT6FitnessMultiplicity(Habitat);
    }


    assert(r >= 0.0 && r <= 1.0);
    return r;
}








double Haplotype::CalculateT1FitnessMultiplicityOnSubsetOfLoci(const int& Habitat, const std::vector<int>& LociSet)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'CalculateT1FitnessMultiplicityOnSubsetOfLoci'\n";
#endif   

    assert(SSP->T1_isMultiplicitySelection);

    // Get the right fitness array
    auto& fits = SSP->T1_FitnessEffects[Habitat];
    assert(fits.size() == SSP->T1_nbLoci);


    double r = 1.0;
    for (const int& locus : LociSet)
    {
        int byte_index = locus / 8;
        int bit_index = locus % 8;

        //assert(polymorphicLocus.locus >= T1_locusFrom && polymorphicLocus.locus < T1_locusTo);
        if (getT1_Allele(byte_index,bit_index))
        {
            r *= fits[locus];
        }
    }

    assert(r<=1.0 && r>=0.0);

    return r;
}

double Haplotype::CalculateT2FitnessOnSubsetOfLoci(const int& Habitat, const std::vector<int>& LociSet)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'CalculateT2FitnessOnSubsetOfLoci'\n";
#endif          
    assert(SSP->T2_FitnessEffects.size() > Habitat);
    auto& fits = SSP->T2_FitnessEffects[Habitat];
    assert(fits.size() == SSP->T2_nbLoci);

    double r = 1.0;
    for (const int locus : LociSet)
    {
        r *= pow(fits[locus], getT2_Allele(locus));
    }

    //if (r!=1.0)  std::cout << "T2 fit: = " << r << " ";
    assert(r<=1.0 && r>=0.0);
    return r;
}


template<typename Iterator>
double Haplotype::CalculateT56FitnessMultiplicityOnSubsetOfLoci(const int& Habitat, const std::vector<int>& LociSet, Iterator it, Iterator itEnd)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'CalculateT56FitnessMultiplicityOnSubsetOfLoci'\n";
#endif   

    assert(SSP->T56_isMultiplicitySelection);

    // Get the right fitness array
    auto& fits = SSP->T56_FitnessEffects[Habitat];
    assert(fits.size() == SSP->T5sel_nbLoci);

    double r = 1.0;

    assert(LociSet.size() > 0);
    size_t locusIndex = 0;
    for (; it != itEnd ; ++it)
    {
        while (LociSet[locusIndex] < *it)
        {
            ++locusIndex;
        }
        if (LociSet[locusIndex] == *it)
        {
            r *= fits[LociSet[locusIndex]];
        }
    }

    return r;
}



double Haplotype::CalculateT56FitnessMultiplicityOnSubsetOfLoci(const int& Habitat, const std::vector<int>& LociSet)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'CalculateT56FitnessMultiplicityOnSubsetOfLoci'\n";
#endif   

    double r;
    if (SSP->T5sel_nbLoci)
    {
        r = CalculateT56FitnessMultiplicityOnSubsetOfLoci(Habitat, LociSet, T5sel_Alleles.begin(), T5sel_Alleles.end());
    } else
    {
        r = CalculateT56FitnessMultiplicityOnSubsetOfLoci(Habitat, LociSet, T6sel_Alleles.begin(), T6sel_Alleles.end());
    }
    
    assert(r>=0.0 && r <= 1.0);
    return r;
}


/*bool Haplotype::getT56_Allele(const int locus) // That's a slow function that should not be used too often!
{
    assert(locus < SSP->T56_nbLoci);
    auto T56ntrlLocus = FromLocusToTXLocus[locus].T56ntrl;
    SSP->isT56neutral[]
    if (SSP->FromT56LocusToT56genderLocus[locus].type == "T56ntrl")
    {
        auto T56ntrlLocus = SSP->FromT56LocusToT56genderLocus[locus].T56ntrl;
        if (SSP->T56ntrl_compress)
        {

        } else
        {

        }
    } else
    {
        assert(SSP->FromT56LocusToT56genderLocus[locus].type == "T56sel");
        auto T56selLocus = SSP->FromT56LocusToT56genderLocus[locus].T56sel;
        if (SSP->T56sel_compress)
        {

        } else
        {

        }

    }
}*/
