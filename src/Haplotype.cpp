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
    assert(SSP->Gmap.T1_nbChars == T1_Alleles.size());
    
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
    assert(SSP->Gmap.T1_nbChars == T1_Alleles.size());
    for (int byte_index = 0 ; byte_index < SSP->Gmap.T1_nbChars ; byte_index++)
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

fitnesstype Haplotype::getW_T8()
{
    return W_T8;
}



fitnesstype Haplotype::getW_T1(int fitnessMapIndex)
{
    assert(fitnessMapIndex < W_T1.size());
    return W_T1[fitnessMapIndex];
}
fitnesstype Haplotype::getW_T2(int fitnessMapIndex)
{
    assert(fitnessMapIndex < W_T2.size());
    return W_T2[fitnessMapIndex];
}

size_t Haplotype::getW_T56_size()
{
    return W_T56.size();
}

void Haplotype::setAllW_T1(fitnesstype w)
{
    assert(W_T1.size() == SSP->NbElementsInFitnessMap);
    for (int fitnessMapIndex = 0 ; fitnessMapIndex < SSP->NbElementsInFitnessMap ; fitnessMapIndex++)
    {
        this->setW_T1(w, fitnessMapIndex);
    }
}
void Haplotype::setAllW_T2(fitnesstype w)
{
    for (int fitnessMapIndex = 0 ; fitnessMapIndex < SSP->NbElementsInFitnessMap ; fitnessMapIndex++)
    {
        this->setW_T2(w, fitnessMapIndex);
    }
}

void Haplotype::setW_T8(fitnesstype w)
{
    W_T8 = w;
}

void Haplotype::setW_T1(fitnesstype w, int fitnessMapIndex)
{
    assert(W_T1.size() > fitnessMapIndex && fitnessMapIndex >= 0);
    W_T1[fitnessMapIndex] = w;
}

void Haplotype::setW_T2(fitnesstype w, int fitnessMapIndex)
{
    assert(W_T2.size() > fitnessMapIndex && fitnessMapIndex >= 0);
    W_T2[fitnessMapIndex] = w;
}

template<typename INT>
unsigned char Haplotype::getT1_char(INT char_index)
{
    //assert(char_index >= 0);
    //assert(char_index < SSP->Gmap.T1_nbChars);
    return T1_Alleles[char_index];
}

size_t Haplotype::nbT7Genes()
{
    return T7_Alleles.size();
}

template<typename INT> void Haplotype::removeT7Gene(INT index)
{
    T7_Alleles.erase(T7_Alleles.begin() + index);
}

template<typename INT> void Haplotype::duplicateT7Gene(INT index)
{
    assert(T7_Alleles.size() < SSP->Gmap.T7_nbLoci);
    assert(T7_Alleles.size());
    
    size_t nbPossibilities = SSP->Gmap.T7_nbLoci - T7_Alleles.size();
    size_t possiblityIndex = GP->rngw.uniform_int_distribution(nbPossibilities);
    
    size_t newlocus = SSP->Gmap.T7_nbLoci;
    for (size_t i = 0 ; i < T7_Alleles.size() ; ++i)
    {
        if (possiblityIndex < T7_Alleles[i].T7Locus - i)
        {
            if (i > 0)
            {
                newlocus = T7_Alleles[i-1].T7Locus + possiblityIndex - i;
            } else
            {
                newlocus = possiblityIndex;
            }
            break;
        }
    }
    if (newlocus == SSP->Gmap.T7_nbLoci)
    {
        newlocus = T7_Alleles.back().T7Locus + possiblityIndex - T7_Alleles.size();
    }
    assert(newlocus < SSP->Gmap.T7_nbLoci);

    T7_Alleles.insert(T7_Alleles.begin() + newlocus, T7_Alleles[index]);
}

template<typename INT>
T7Gene& Haplotype::getT7_Allele(const INT index)
{
    return T7_Alleles[index];
}

void Haplotype::clearT7Genes()
{
    std::vector<T7Gene>().swap(T7_Alleles);
}

template<typename INT>
bool Haplotype::getT1_Allele(const INT T1Locus)
{
    const int char_index = T1Locus / 8;
    const int bit_index  = T1Locus % 8;
    return (this->getT1_Allele(char_index, bit_index));
}

template<typename INT1, typename INT2>
bool Haplotype::getT1_Allele(const INT1 char_index, const INT2 bit_index)
{
    /*
    assert(bit_index < 8);
    assert(bit_index >= 0);
    assert(char_index >= 0);
    assert(char_index < SSP->Gmap.T1_nbChars);
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

template<typename INT>
unsigned char Haplotype::getT2_Allele(const INT char_index)
{
    //assert(char_index >= 0);
    //assert(char_index < SSP->Gmap.T2_nbLoci);
    return T2_Alleles[char_index];
}

template<typename INT>
T3type Haplotype::getT3_Allele(const INT index)
{
    //assert(char_index >= 0);
    //assert(char_index < SSP->Gmap.T2_nbLoci);
    return T3_Alleles[index];
}

template<typename INT>
void Haplotype::setT1_char(INT char_index, unsigned char& c)
{
    //assert(char_index >= 0);
    //assert(char_index < SSP->Gmap.T1_nbChars);
    T1_Alleles[char_index] = c;
}

template<typename INT>
void Haplotype::setT1_char(INT char_index, unsigned char&& c)
{
    //assert(char_index >= 0);
    //assert(char_index < SSP->Gmap.T1_nbChars);
    T1_Alleles[char_index] = c;
}

template<typename valueType, typename INT>
void Haplotype::setT1_Allele(const INT char_index, const INT bit_index, const valueType& value)
{
    //assert(bit_index < 8);
    //assert(bit_index >= 0);
    //assert(char_index >= 0);
    //assert(char_index < SSP->Gmap.T1_nbChars);
    T1_Alleles[char_index] = (-value ^ T1_Alleles[char_index]) & (1 << bit_index);
}

template<typename INT>
void Haplotype::setT1_AlleleToOne(INT char_index, INT bit_index)
{
    //assert(bit_index < 8);
    //assert(bit_index >= 0);
    //assert(char_index >= 0);
    //assert(char_index < SSP->Gmap.T1_nbChars);
    T1_Alleles[char_index] |= 1 << bit_index;
}

template<typename INT>
void Haplotype::setT1_AlleleToZero(INT char_index, INT bit_index)
{
    /*
    assert(bit_index < 8);
    assert(bit_index >= 0);
    assert(char_index >= 0);
    assert(char_index < SSP->Gmap.T1_nbChars);
    std::cout << "SSP->Gmap.T1_nbChars = " << SSP->Gmap.T1_nbChars << "\n";
    */
    T1_Alleles[char_index] &= ~(1 << bit_index);
}

template<typename INT>
void Haplotype::setT2_Allele(const INT char_index, const unsigned char value)
{
    //assert(char_index >= 0);
    //assert(char_index < SSP->Gmap.T2_nbLoci);
    T2_Alleles[char_index] = value;
}

template<typename INT>
void Haplotype::setT3_Allele(const INT index, const T3type value)
{
    //assert(char_index >= 0);
    //assert(char_index < SSP->Gmap.T2_nbLoci);
    T3_Alleles[index] = value;
}

template<typename INT>
void Haplotype::toggleT1_Allele(INT byte_index, INT bit_index)
{
    T1_Alleles[byte_index] ^= 1 << bit_index;
}

template<typename INT>
void Haplotype::mutateT1_Allele(INT MutPosition, int& Habitat)
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
    //assert(byte_index < SSP->Gmap.T1_nbChars);

    bool shouldIMutate = true;
    if (SSP->T1_mutDirection != 2)
    {
        if (getT1_Allele(byte_index, bit_index))
        {
            if (SSP->T1_mutDirection == 1) shouldIMutate = false;
        } else
        {
            if (SSP->T1_mutDirection == 0) shouldIMutate = false;
        }  
    }


    if (shouldIMutate)
    {
        toggleT1_Allele(byte_index, bit_index);

        if (SSP->T1_isMultiplicitySelection)
        {
            int fitnessMapIndex = SSP->FromLocusToFitnessMapIndex[SSP->Gmap.FromT1LocusToLocus(MutPosition)];
            fitnesstype w = this->getW_T1(fitnessMapIndex);
            //if (GP->CurrentGeneration >=50) std::cout << "w=" << w << " ";
            if (w != -1.0)
            {
                // Note that Toggle already happened. If it is mutant, it is because it just happened!
                if ( this->getT1_Allele(byte_index, bit_index) ) // mutation added
                {
                    w *= SSP->T1_FitnessEffects[Habitat][byte_index * 8 + bit_index];
                    assert(w >= 0.0);
                } else // mutation removed
                {
                    if (SSP->T1_FitnessEffects[Habitat][byte_index * 8 + bit_index] == 0)
                    { 
                        w = -1.0;
                    } else
                    {
                        w /= SSP->T1_FitnessEffects[Habitat][byte_index * 8 + bit_index]; 
                        assert(w >= 0.0);
                    }
                }
                this->setW_T1( w,   fitnessMapIndex );
            }
        }
    }
}

template<typename INT>
void Haplotype::AddMutT2_Allele(INT char_index)
{
    //assert(char_index >= 0);
    //assert(char_index < SSP->Gmap.T2_nbLoci);
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
        int fitnessMapIndex = SSP->FromLocusToFitnessMapIndex[SSP->Gmap.FromT2LocusToLocus(char_index)];
        this->setW_T2(-1.0,fitnessMapIndex);
    }
}

template<typename INT>
void Haplotype::AddMutT2_Allele(INT char_index, int Habitat)
{
    //assert(char_index >= 0);
    //assert(char_index < SSP->Gmap.T2_nbLoci);
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
        int fitnessMapIndex = SSP->FromLocusToFitnessMapIndex[SSP->Gmap.FromT2LocusToLocus(char_index)];
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

template<typename INT>
void Haplotype::mutateT3_Allele(INT index)
{
    // Add or substract while preventing
    #ifdef DEBUG
    std::cout << "T3_Alleles.size() = " << T3_Alleles.size() << "\n";
    std::cout << "SSP->Gmap.T3_nbLoci = " << SSP->Gmap.T3_nbLoci << "\n";
    assert(SSP->Gmap.T3_nbLoci == T3_Alleles.size());    
    assert(index < T3_Alleles.size());
    assert(index >= 0);
    #endif

    if (SSP->T3_mutationType == '1')
    {
        if (GP->rngw.get_1b())
        {
            if (T3_Alleles[index] != std::numeric_limits<T3type>::max())
            {
                T3_Alleles[index] += SSP->T3_mutationType_effectSizeInfo;
            }
        } else
        {
            if (T3_Alleles[index] != std::numeric_limits<T3type>::lowest())
            {
                T3_Alleles[index] -= SSP->T3_mutationType_effectSizeInfo;
            }
        }
    } else if (SSP->T3_mutationType_effectSizeInfo != 0.0)
    {
        if (SSP->T3_mutationType == 'g')
        {
            T3_Alleles[index] += GP->rngw.normal((T3type)0.0, SSP->T3_mutationType_effectSizeInfo);
        } else if (SSP->T3_mutationType == 'G')
        {
            T3_Alleles[index] = GP->rngw.normal((T3type)0.0, SSP->T3_mutationType_effectSizeInfo);
        } else abort();
    }
}

template<typename INT>
void Haplotype::copyIntoT1(INT from, INT to, Haplotype& SourceChromo)
{
    // copy [first, last)
    //std::cout << "In copyIntoT1 from = "<< from << "  to = "<< to << std::endl;
    assert(from >= 0);
    assert(to > from);
    
    uint32_t ByteFrom = from / 8; // integer division. It will floor automatically
    uint32_t ByteTo   = to / 8; // integer division. It will floor automatically
    assert(ByteTo <= SourceChromo.T1_Alleles.size() + 1);

    uint32_t bitIndexFrom = from % 8;
    uint32_t bitIndexTo = to % 8;

    if (ByteTo == SSP->Gmap.T1_nbChars - 1) // This is to deal with the last extra bits. Note that if it copies a few extra bits, it is not an issue.
    {
        bitIndexTo = std::max(bitIndexTo, SSP->Gmap.T1_nbLociLastByte);
    }


    if (ByteTo == ByteFrom)
    {
        assert(bitIndexFrom < bitIndexTo);
        for (uint32_t bit_index = bitIndexFrom ; bit_index < bitIndexTo ; bit_index++)
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
        for (uint32_t bit_index = bitIndexFrom ; bit_index < 8 ; bit_index++)
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
        for (uint32_t bit_index = 0 ; bit_index < bitIndexTo ; bit_index++)
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

template<typename INT>
void Haplotype::copyIntoT2(INT from, INT to, Haplotype& SourceChromo)
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

template<typename INT>
void Haplotype::copyIntoT3(INT from, INT to, Haplotype& SourceChromo)
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
    

    /*
    // Just wanted to test testing before copying but that's slower
    for (size_t i = from ; i < to; ++i)
    {
        if (this->T3_Alleles[i] != SourceChromo.T3_Alleles[i])
            this->T3_Alleles[i] = SourceChromo.T3_Alleles[i];
    }
    */
    
}


Haplotype::Haplotype()
{
    if (SSP->T1_isMultiplicitySelection)
    {
        assert(SSP->NbElementsInFitnessMap > 0);
        assert(SSP->Gmap.T1_nbLoci > 0);
        W_T1.resize(SSP->NbElementsInFitnessMap,-1.0);
    }
    if (SSP->T2_isSelection)
    {
        assert(SSP->NbElementsInFitnessMap > 0);
        assert(SSP->Gmap.T2_nbLoci > 0);
        W_T2.resize(SSP->NbElementsInFitnessMap,-1.0);
    }
    if (SSP->T56_isMultiplicitySelection)
    {
        assert(SSP->NbElementsInFitnessMap > 0);
        assert(SSP->Gmap.T56sel_nbLoci > 0);
        W_T56.resize(SSP->NbElementsInFitnessMap,-1.0);
    }


    T1_Alleles.resize(SSP->Gmap.T1_nbChars);
    T2_Alleles.resize(SSP->Gmap.T2_nbLoci);
    T3_Alleles.resize(SSP->Gmap.T3_nbLoci);
    //T5sel_Alleles.resize(SSP->Gmap.T5sel_nbLoci);
    //T5ntrl_Alleles.resize(SSP->Gmap.T5ntrl_nbLoci);
    //T6sel_Alleles.resize(SSP->Gmap.T6sel_nbLoci);
    //T6ntrl_Alleles.resize(SSP->Gmap.T6ntrl_nbLoci);
    assert(SSP->Gmap.T7_nbLoci == 0); // Just because I don't know how to initialize that yet
    //T7_Alleles.resize(SSP->Gmap.T7_nbLoci);
}

void Haplotype::PrintBinaryFile(OutputFile& file)
{
    // Print T1_alleles
    if (SSP->Gmap.T1_nbLoci) file.writeBinary(T1_Alleles);

    // Print T2_alleles
    if (SSP->Gmap.T2_nbLoci) file.writeBinary(T2_Alleles);

    // Print T3_alleles
    if (SSP->Gmap.T3_nbLoci) file.writeBinary(T3_Alleles);

    // Print T4 tree aside
    /*
    if (SSP->Gmap.T4_nbLoci)
        std::cout << "\n\tWARNING: You are printing to a binary file while you have T4 loci. In the current version of SimBit, T4 loci cannot be printed on a binary file (sorry). The binary file will be produced ignoring T4 loci\n";
    */

    // Print T56_alleles (there should not be any flipped meaning normally)
    if (SSP->Gmap.T56ntrl_nbLoci)
    {
        std::vector<uint32_t> T5ntrlTrueHaplotype = this->getT56ntrlTrueHaplotype();
        file.writeBinary(T5ntrlTrueHaplotype);
    }

    if (SSP->Gmap.T56sel_nbLoci)
    {
        std::vector<uint32_t> T5selTrueHaplotype = this->getT56selTrueHaplotype();
        file.writeBinary(T5selTrueHaplotype);
    }
}

Haplotype::Haplotype(std::vector<unsigned char> T1_info, std::vector<unsigned char> T2_info, std::vector<T3type> T3_info, uint32_t t4id, std::vector<uint32_t> T56_info)
{
    assert(SSP != nullptr); 

    /*
    std::cout << "T1_info: ";
    for (auto& e : T1_info) {std::cout << e << " ";}
    std::cout << "\n";
    */

    T1_Alleles = T1_info;
    assert(T1_Alleles.size() == SSP->Gmap.T1_nbChars);
    T2_Alleles = T2_info;
    assert(T2_Alleles.size() == SSP->Gmap.T2_nbLoci);
    T3_Alleles = T3_info;
    assert(T3_Alleles.size() == SSP->Gmap.T3_nbLoci);
    assert(T56_info.size() <= SSP->Gmap.T56_nbLoci);


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
        if (SSP->Gmap.isT56neutral(mutPos))
        {
            if (SSP->Gmap.isT56ntrlCompress)
            {
                T6ntrl_Alleles.push_back(mutPos);
            } else
            {
                T5ntrl_Alleles.push_back(mutPos);
            }
        } else
        {
            if (SSP->Gmap.isT56selCompress)
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
        assert(SSP->Gmap.T1_nbLoci > 0);
        W_T1.resize(SSP->NbElementsInFitnessMap,-1.0);
    }
    if (SSP->T2_isSelection)
    {
        assert(SSP->NbElementsInFitnessMap > 0);
        assert(SSP->Gmap.T2_nbLoci > 0);
        W_T2.resize(SSP->NbElementsInFitnessMap,-1.0);
    }
    if (SSP->T56_isMultiplicitySelection)
    {
        assert(SSP->NbElementsInFitnessMap > 0);
        assert(SSP->Gmap.T56sel_nbLoci > 0);
        W_T56.resize(SSP->NbElementsInFitnessMap,-1.0);
    }

    T4ID = t4id;
}

Haplotype::Haplotype(bool ShouldReadPopFromBinary)
{
#ifdef CALLENTRANCEFUNCTIONS
    std::cout << "Enters in 'Haplotype::Haplotype(bool ShouldReadPopFromBinary)'\n";
#endif     
    (void) ShouldReadPopFromBinary;

    //.read T1_alleles
    if (SSP->Gmap.T1_nbLoci) GP->binaryFileToRead.readVector(T1_Alleles);

    
    //.read T2_alleles
    if (SSP->Gmap.T2_nbLoci) GP->binaryFileToRead.readVector(T2_Alleles);

    
    //.read T3_alleles
    if (SSP->Gmap.T3_nbLoci) GP->binaryFileToRead.readVector(T3_Alleles);

    
    //.read T4_alleles
    // nothiung to do. T4Tree read later

    //.read T56_alleles    
    if (SSP->Gmap.T56ntrl_nbLoci)
    {        
        if (SSP->Gmap.isT56ntrlCompress)
        {
            std::vector<uint32_t> tmp;
            GP->binaryFileToRead.readVector(tmp);
            T6ntrl_Alleles = CompressedSortedDeque(tmp, SSP->Gmap.T56ntrl_nbLoci);
            if (T6ntrl_Alleles.size() > SSP->Gmap.T56ntrl_nbLoci)
            {
                std::cout << "Reading population from a binary file. Read a vector of T6ntrl loci that is longer than the number of T6ntrl loci set in this simulation.\n";
                abort();
            }
        } else
        {
            GP->binaryFileToRead.readVector(T5ntrl_Alleles);
            if (T5ntrl_Alleles.size() > SSP->Gmap.T56ntrl_nbLoci)
            {
                std::cout << "Reading population from a binary file. Read a vector of T5ntrl loci that is longer than the number of T5ntrl loci set in this simulation.\n";
                abort();
            }
        }
    }
        

    if (SSP->Gmap.T56sel_nbLoci)
    {        
        if (SSP->Gmap.isT56selCompress)
        {
            std::vector<uint32_t> tmp;
            GP->binaryFileToRead.readVector(tmp);
            T6sel_Alleles = CompressedSortedDeque(tmp, SSP->Gmap.T56sel_nbLoci);
            if (T6sel_Alleles.size() > SSP->Gmap.T56sel_nbLoci)
            {
                std::cout << "Reading population from a binary file. Read a vector of T6sel loci that is longer than the number of T6sel loci set in this simulation.\n";
                abort();
            }
        } else
        {
            GP->binaryFileToRead.readVector(T5sel_Alleles);
            if (T5sel_Alleles.size() > SSP->Gmap.T56sel_nbLoci)
            {
                std::cout << "Reading population from a binary file. Read a vector of T5sel loci that is longer than the number of T5sel loci set in this simulation.\n";
                abort();
            }
        }
    }


    // Initiate W_T1, W_T2 and W_T56
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
: T1_Alleles(SSP->Gmap.T1_nbChars), T2_Alleles(SSP->Gmap.T2_nbLoci), T3_Alleles(SSP->Gmap.T3_nbLoci) // No initilization for T5_Alleles and T6_Alleles
{
#ifdef CALLENTRANCEFUNCTIONS
    std::cout << "Enters in 'Haplotype::Haplotype(const int patch_index, char Abiogenesis)'\n";
#endif   
    (void)Abiogenesis; // Does nothing but silence the warning that the char Abiogenesis is not used.
    assert(SSP!=nullptr);
    assert(SSP->Gmap.T1_nbChars + SSP->Gmap.T2_nbLoci + SSP->Gmap.T3_nbLoci + SSP->Gmap.T4_nbLoci + SSP->Gmap.T56_nbLoci + SSP->Gmap.T8_nbLoci > 0);
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
    if (SSP->Gmap.T1_nbLoci > 0)
    {        
        if (SSP->T1_Initial_AlleleFreqs_AllZeros) // The whole chromosome must be 0s
        {
            unsigned char c = 0;
            for (int char_index = 0 ; char_index < SSP->Gmap.T1_nbChars ; char_index++)
            {
                this->setT1_char(char_index,c);
            }
        } else if (SSP->T1_Initial_AlleleFreqs_AllOnes) // The whole chromosome must be 1s
        {
            unsigned char c = 255;
            for (int char_index = 0 ; char_index < SSP->Gmap.T1_nbChars ; char_index++)
            {
                this->setT1_char(char_index,c);
            }
        } else 
        {
            assert(SSP->T1_Initial_AlleleFreqs.size() > patch_index);
            for (int byte_index = 0 ; byte_index < SSP->Gmap.T1_nbChars ; byte_index++)
            {
                int bit_index_to;
                if (byte_index == SSP->Gmap.T1_nbChars - 1)
                {
                    bit_index_to = SSP->Gmap.T1_nbLociLastByte;
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
        assert(T1_Alleles.size() == SSP->Gmap.T1_nbChars);
    }
        

    // Initiate T2_Alleles
    if (SSP->Gmap.T2_nbLoci > 0)
    {
        for (int char_index=0 ; char_index < SSP->Gmap.T2_nbLoci ; char_index++)
        {
            unsigned char c = 0;
            this->setT2_Allele(char_index,c);
        }
        assert(T2_Alleles.size() == SSP->Gmap.T2_nbLoci);
    }

    // Initiate T3_Alleles
    if (SSP->Gmap.T3_nbLoci > 0)
    {
        for (int i=0 ; i < SSP->Gmap.T3_nbLoci ; i++)
        {
            T3type d = 0.0;
            this->setT3_Allele(i,d);
        }
        assert(T3_Alleles.size() == SSP->Gmap.T3_nbLoci);
    }

    // Initiate T4_Alleles
    // No initialization. Assume all are fixed at 0

    // Initiate T56_Alleles 
    assert(SSP->Gmap.T5sel_nbLoci == 0 || SSP->Gmap.T6sel_nbLoci == 0);
    assert(SSP->Gmap.T5ntrl_nbLoci == 0 || SSP->Gmap.T6ntrl_nbLoci == 0);

    if (SSP->Gmap.T56_nbLoci > 0)
    { 
        // Initialize T6
        if (SSP->Gmap.T6ntrl_nbLoci)
        {
            T6ntrl_Alleles = CompressedSortedDeque(SSP->Gmap.T6ntrl_nbLoci);
            (void) T6ntrl_AllelesBegin();
        }

        if (SSP->Gmap.T6sel_nbLoci)
        {
            T6sel_Alleles = CompressedSortedDeque(SSP->Gmap.T6sel_nbLoci);
            (void) T6sel_AllelesBegin();
        }


        // Initial Allele Frequencies
        if (SSP->T56_Initial_AlleleFreqs_AllZeros) // The whole chromosome must be 0s
        {
            // Nothing to do
            if (SSP->Gmap.T5ntrl_nbLoci)
                assert(SSP->T5ntrl_flipped.size() == 0);
        } else if (SSP->T56_Initial_AlleleFreqs_AllOnes) // The whole chromosome must be 1s
        {   
            // Nothing to do
            std::cout << "Internal error: It seems that T56_Initial_AlleleFreqs_AllOnes is true but this should not have been allowed when reading the user input.\n";
            abort();
        } else 
        {

            assert(SSP->T56_Initial_AlleleFreqs.size() > patch_index);
            for (int locus = 0 ; locus < SSP->Gmap.T56_nbLoci ; locus++)
            {
                auto genderLocus = SSP->Gmap.FromT56LocusToT56genderLocus(locus);
                auto locusInGender = genderLocus.locusInGender;
                assert(SSP->T56_Initial_AlleleFreqs[patch_index].size() > locus);
                if (SSP->T56_Initial_AlleleFreqs[patch_index][locus] == 1.0)
                {
                    if (genderLocus.isNtrl)
                    {
                        if (SSP->Gmap.isT56ntrlCompress)
                        {
                            this->setT6ntrl_AlleleToOne_JustPushBack(locusInGender);
                        } else
                        {
                            this->setT5ntrl_AlleleToOne_JustPushBack(locusInGender);
                        }
                    } else
                    {
                        if (SSP->Gmap.isT56selCompress)
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
                        if (SSP->Gmap.FromT56LocusToT56genderLocus(locus).isNtrl)
                        {
                            if (SSP->Gmap.isT56ntrlCompress)
                            {
                                this->setT6ntrl_AlleleToOne_JustPushBack(locusInGender);
                            } else
                            {
                                this->setT5ntrl_AlleleToOne_JustPushBack(locusInGender);
                            }
                        } else
                        {
                            if (SSP->Gmap.isT56selCompress)
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

    T1_Alleles.shrink_to_fit();
    T2_Alleles.shrink_to_fit();
    T3_Alleles.shrink_to_fit();
}

template<typename INT>
void Haplotype::setEntireT5_Allele(std::vector<INT>& t5a)
{
    assert(t5a.size() == SSP->Gmap.T5_nbLoci);
    assert(T5ntrl_Alleles.size() == SSP->Gmap.T5ntrl_nbLoci);
    assert(T5sel_Alleles.size() == SSP->Gmap.T5sel_nbLoci);

    uint32_t ntrlLocus = 0;
    uint32_t selLocus = 0;
    for (uint32_t locus = 0 ; locus < SSP->Gmap.T5_nbLoci ; locus++ )
    {
        if (SSP->Gmap.FromT56ntrlLocusToLocus(ntrlLocus) < SSP->Gmap.FromT56selLocusToLocus(selLocus))
        {
            this->T5ntrl_Alleles[ntrlLocus] = t5a[locus];
            ntrlLocus++;

        } else
        {
            this->T5sel_Alleles[selLocus] = t5a[locus];
            selLocus++;            
        }
    }
    assert(selLocus == SSP->Gmap.T5sel_nbLoci);
    assert(ntrlLocus == SSP->Gmap.T5ntrl_nbLoci);
}


template<typename INT>
void Haplotype::setEntireT6_Allele(std::vector<INT>& t6a)
{
    assert(t6a.size() == SSP->Gmap.T6_nbLoci);

    uint32_t ntrlLocus = 0;
    uint32_t selLocus = 0;
    for (uint32_t locus = 0 ; locus < SSP->Gmap.T6_nbLoci ; locus++ )
    {
        if (SSP->Gmap.FromT56ntrlLocusToLocus(ntrlLocus) < SSP->Gmap.FromT56selLocusToLocus(selLocus))
        {
            this->T6ntrl_Alleles.push_back(t6a[locus]);
            ntrlLocus++;

        } else
        {
            this->T6sel_Alleles.push_back(t6a[locus]);
            selLocus++;            
        }
    }
    assert(selLocus == SSP->Gmap.T6sel_nbLoci);
    assert(ntrlLocus == SSP->Gmap.T6ntrl_nbLoci);
}

Haplotype::Haplotype(const Haplotype& other)
{
    if (SSP->Gmap.T1_nbLoci)
    {
        T1_Alleles = other.T1_Alleles;
        if (SSP->T1_isMultiplicitySelection)
        {
            W_T1 = other.W_T1;
        }
    }

    if (SSP->Gmap.T2_nbLoci)
    {
        T2_Alleles = other.T2_Alleles;
        W_T2 = other.W_T2;
    }

    if (SSP->Gmap.T3_nbLoci)
    {
        T3_Alleles = other.T3_Alleles;
    }

    if (SSP->Gmap.T56ntrl_nbLoci)
    {
        if (SSP->Gmap.T5ntrl_nbLoci)
        {
            T5ntrl_Alleles = other.T5ntrl_Alleles;
        } else
        {
            T6ntrl_Alleles = other.T6ntrl_Alleles;
        }
    }

    if (SSP->Gmap.T56sel_nbLoci)
    {
        if (SSP->Gmap.T5sel_nbLoci)
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

    if (SSP->Gmap.T4_nbLoci)
    {
        T4ID = other.T4ID;
    }

    if (SSP->Gmap.T8_nbLoci)
    {
        //std::cout << "Copying T8ID " << other.T8ID << "\n";
        T8ID = other.T8ID;
        W_T8 = other.W_T8;
    }
}

/*
Haplotype::Haplotype( Haplotype&& other)
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
    T4ID = std::move(other.T4ID);
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
    T4ID = other.T4ID;
    T8ID = other.T8ID;
}*/

Haplotype& Haplotype::operator=(const Haplotype& other) // copy assignment operator
{
    //std::cout << "copy assignment haplotype\n";
    if (SSP->Gmap.T1_nbLoci)
    {
        T1_Alleles = other.T1_Alleles;
        if (SSP->T1_isMultiplicitySelection)
        {
            W_T1 = other.W_T1;
        }
    }

    if (SSP->Gmap.T2_nbLoci)
    {
        T2_Alleles = other.T2_Alleles;
        W_T2 = other.W_T2;
    }

    if (SSP->Gmap.T3_nbLoci)
    {
        T3_Alleles = other.T3_Alleles;
    }

    if (SSP->Gmap.T56ntrl_nbLoci)
    {
        if (SSP->Gmap.T5ntrl_nbLoci)
        {
            T5ntrl_Alleles = other.T5ntrl_Alleles;
        } else
        {
            T6ntrl_Alleles = other.T6ntrl_Alleles;
        }
    }

    if (SSP->Gmap.T56sel_nbLoci)
    {
        if (SSP->Gmap.T5sel_nbLoci)
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

    if (SSP->Gmap.T4_nbLoci)
    {
        T4ID = other.T4ID;
    }

    if (SSP->Gmap.T8_nbLoci)
    {
        T8ID = other.T8ID;
        W_T8 = other.W_T8;
    }

    return *this;
}
/*
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
    T4ID = std::move(other.T4ID);
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
    T4ID = other.T4ID;
    T8ID = other.T8ID;
    //assert(W_T1.size() == other.W_T1.size());
    return *this;
}*/

void Haplotype::swap(Haplotype& other)
{
    if (SSP->Gmap.T1_nbLoci)
    {
        T1_Alleles.swap(other.T1_Alleles);
        W_T1.swap(other.W_T1);
    }

    if (SSP->Gmap.T2_nbLoci)
    {
        T2_Alleles.swap(other.T2_Alleles);
        W_T2.swap(other.W_T2);
    }

    if (SSP->Gmap.T3_nbLoci)
        T3_Alleles.swap(other.T3_Alleles);

    if (SSP->Gmap.T56sel_nbLoci)    
        W_T56.swap(other.W_T56);

    if (SSP->Gmap.T5ntrl_nbLoci)
        T5ntrl_Alleles.swap(other.T5ntrl_Alleles);

    if (SSP->Gmap.T5sel_nbLoci)
        T5sel_Alleles.swap(other.T5sel_Alleles);

    if (SSP->Gmap.T6ntrl_nbLoci)
        T6ntrl_Alleles.swap(other.T6ntrl_Alleles);
    
    if (SSP->Gmap.T6sel_nbLoci)
        T6sel_Alleles.swap(other.T6sel_Alleles);

    if (SSP->Gmap.T4_nbLoci)    
        std::swap(T4ID, other.T4ID);
    
    if (SSP->Gmap.T8_nbLoci)
    {
       std::swap(T8ID, other.T8ID);
       std::swap(W_T8, other.W_T8);
    }
}


void Haplotype::setW_T56(fitnesstype w, int fitnessMapIndex)
{
    if (!(W_T56.size() > fitnessMapIndex && fitnessMapIndex >= 0))
    {
        std::cout << "SSP->T56_isMultiplicitySelection = " << SSP->T56_isMultiplicitySelection << "\n";
        std::cout << "fitnessMapIndex = " << fitnessMapIndex << "\n";
        std::cout << "SSP->NbElementsInFitnessMap = " << SSP->NbElementsInFitnessMap << "\n";
        std::cout << "W_T56.size() = " << W_T56.size() << "\n";
    }
        
    assert(W_T56.size() > fitnessMapIndex && fitnessMapIndex >= 0);
    W_T56[fitnessMapIndex] = w;
}


void Haplotype::setAllW_T56(fitnesstype w)
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

template<typename INT>
void Haplotype::setT5ntrl_Allele(const INT locus, const bool& value)
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

template<typename INT>
void Haplotype::setT5sel_Allele(const INT locus, const bool& value)
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


template<typename INT>
void Haplotype::setT6ntrl_Allele(const INT locus, const bool& value)
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

template<typename INT>
void Haplotype::setT6sel_Allele(const INT locus, const bool& value)
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

template<typename INT>
void Haplotype::setT5_Allele(const INT locus, const bool& value)
{
    auto& TXLocusElement = SSP->Gmap.FromT56LocusToT56genderLocus(locus);
    TXLocusElement.isNtrl ? setT5sel_Allele(TXLocusElement.locusInGender, value) : setT5ntrl_Allele(TXLocusElement.locusInGender, value);
}

template<typename INT>
void Haplotype::setT6_Allele(const INT locus, const bool& value)
{
    auto& TXLocusElement = SSP->Gmap.FromT56LocusToT56genderLocus(locus);
    TXLocusElement.isNtrl ? setT6sel_Allele(TXLocusElement.locusInGender, value) : setT6ntrl_Allele(TXLocusElement.locusInGender, value);
}

template<typename INT>
void Haplotype::setT5sel_AlleleToOne_JustPushBack(INT locus)
{
    T5sel_Alleles.push_back(locus);
}

template<typename INT>
void Haplotype::setT5sel_AlleleToOne(INT locus)
{
    this->setT5sel_Allele(locus, true);
}

template<typename INT>
void Haplotype::setT5sel_AlleleToZero(INT locus)
{
    this->setT5sel_Allele(locus, false);
}

template<typename INT>
void Haplotype::setT5ntrl_AlleleToOne_JustPushBack(const INT locus)
{
    T5ntrl_Alleles.push_back(locus);
}

template<typename INT>
void Haplotype::setT5ntrl_AlleleToOne(INT locus)
{
    this->setT5ntrl_Allele(locus, true);
}

template<typename INT>
void Haplotype::setT5ntrl_AlleleToZero(INT locus)
{
    this->setT5ntrl_Allele(locus, false);
}


template<typename INT>
void Haplotype::setT6sel_AlleleToOne_JustPushBack(const INT locus)
{
    T6sel_Alleles.push_back(locus);
}

template<typename INT>
void Haplotype::setT6sel_AlleleToOne(INT locus)
{
    this->setT6sel_Allele(locus, true);
}

template<typename INT>
void Haplotype::setT6sel_AlleleToZero(INT locus)
{
    this->setT6sel_Allele(locus, false);
}

template<typename INT>
void Haplotype::setT6ntrl_AlleleToOne_JustPushBack(INT locus)
{
    T6ntrl_Alleles.push_back(locus);
}

template<typename INT>
void Haplotype::setT6ntrl_AlleleToOne(INT locus)
{
    this->setT6ntrl_Allele(locus, true);
}

template<typename INT>
void Haplotype::setT6ntrl_AlleleToZero(INT locus)
{
    this->setT6ntrl_Allele(locus, false);
}



fitnesstype Haplotype::getW_T56(int fitnessMapIndex)
{
    //std::cout << "fitnessMapIndex = " << fitnessMapIndex << "\n";
    //std::cout << "W_T56.size() = " << W_T56.size() << "\n";
    assert(fitnessMapIndex >= 0 && fitnessMapIndex < W_T56.size());
    return W_T56[fitnessMapIndex];
}

template<typename INT>
void Haplotype::toggleT5ntrl_Allele(INT MutPosition)
{
    assert(MutPosition < SSP->Gmap.T5ntrl_nbLoci);
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
    assert(MutPosition < SSP->Gmap.T5sel_nbLoci);
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
    assert(MutPosition < SSP->Gmap.T6ntrl_nbLoci);
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
    assert(MutPosition < SSP->Gmap.T6sel_nbLoci);
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
    auto info = T56finishMutation(haploP, T5sel_Alleles, MutPosition);
        
    /////////////////////
    // Update fitness //
    ////////////////////

    if (info != -1 && SSP->T56_isMultiplicitySelection)
    {
        updateFitnessAfterT56Mutation(MutPosition, info == 1, Habitat);
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
    auto info = T56finishMutation(haploP, T6sel_Alleles, MutPosition);
    
    /////////////////////
    // Update fitness //
    ////////////////////
    if (info != -1 && SSP->T56_isMultiplicitySelection)
    {
        updateFitnessAfterT56Mutation(MutPosition, info == 1, Habitat);
    }
}


template<typename INTORINTVECTOR>
void Haplotype::mutateT56ntrl_Allele(INTORINTVECTOR MutPosition)
{
    if (SSP->Gmap.isT56ntrlCompress)
    {
        mutateT6ntrl_Allele(MutPosition);
    } else
    {
        mutateT5ntrl_Allele(MutPosition);
    }
}


template<typename INTORINTVECTOR>
void Haplotype::mutateT56sel_Allele(INTORINTVECTOR MutPosition, int Habitat)
{
    if (SSP->Gmap.isT56selCompress)
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
    std::cout << "Function void Haplotype::mutateT5ntrl_Allele(std::vector<INT>& MutPositions) should not be usedd at the moment. Internal bug!\n";
    abort();
    /*
    std::vector<size_t> its;
    for (auto& mut : MutPositions)
    {
        its.push_back(std::lower_bound(T5ntrl_Alleles.begin(), T5ntrl_Alleles.end(), mut));
    }

    insertAtPositions(x, values, positions);

    uint32_t prevSize = T5ntrl_Alleles.size();
    T5ntrl_Alleles.insert(T5ntrl_Alleles.end(), MutPositions.begin(), MutPositions.end());
    std::inplace_merge(T5ntrl_Alleles.begin(), T5ntrl_Alleles.begin() + prevSize, T5ntrl_Alleles.end());
    T5ntrl_Alleles.erase(std::unique(T5ntrl_Alleles.begin(), T5ntrl_Alleles.end()), T5ntrl_Alleles.end());
    */
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

void Haplotype::toggleFromT5ntrl_Allele(int& MutPosition, uint32_t& from)
{
    auto haploP = std::lower_bound(T5ntrl_Alleles.begin() + from, T5ntrl_Alleles.end(), MutPosition);
    from = haploP - T5ntrl_Alleles.begin();
    
    (void) T56finishMutation(haploP, T5ntrl_Alleles, MutPosition);
}

void Haplotype::toggleFromT6ntrl_Allele(int& MutPosition, uint32_t& blockIndexFrom, uint32_t& fromInBlock)
{
    assert(MutPosition >= 0 && MutPosition < SSP->Gmap.T6ntrl_nbLoci);
    CompressedSortedDeque::iterator haploP = T6ntrl_Alleles.lower_bound_from(MutPosition, blockIndexFrom, fromInBlock); // blockIndexFrom and fromInBlock updated in here
    (void) T56finishMutation(haploP, T6ntrl_Alleles, MutPosition);
}


void Haplotype::toggleFromT5sel_Allele(int& MutPosition, uint32_t& from, int Habitat)
{
    auto haploP = std::lower_bound(T5sel_Alleles.begin() + from, T5sel_Alleles.end(), MutPosition);
    from = haploP - T5sel_Alleles.cbegin();

    if (SSP->T56_isMultiplicitySelection)
        updateFitnessAfterT56Mutation(MutPosition, T56finishMutation(haploP, T5sel_Alleles, MutPosition) == 1, Habitat);

}

void Haplotype::toggleFromT6sel_Allele(int& MutPosition, uint32_t& blockIndexFrom, uint32_t& fromInBlock, int Habitat)
{
    auto haploP = T6sel_Alleles.lower_bound_from(MutPosition, blockIndexFrom, fromInBlock); // blockIndexFrom and fromInBlock updated in lower_bound
    if (SSP->T56_isMultiplicitySelection)
        updateFitnessAfterT56Mutation(MutPosition, T56finishMutation(haploP, T6sel_Alleles, MutPosition) == 1, Habitat);
}



void Haplotype::clearT56Alleles()
{
    if (SSP->Gmap.T56ntrl_nbLoci)
    {
        if (SSP->Gmap.isT56ntrlCompress)
        {
            T6ntrl_Alleles.clear();
        } else
        {
            T5ntrl_Alleles.clear();
        }
    }
        
    if (SSP->Gmap.T56sel_nbLoci)
    {
        if (SSP->Gmap.isT56selCompress)
        {
            T6sel_Alleles.clear();
        } else
        {
            T5sel_Alleles.clear();
        }
    }
}

template<typename INT>
void Haplotype::copyIntoT56sel(INT from, INT to, Haplotype& SourceChromo)
{
    // Must have been emptied first

    assert(from >= 0);
    assert(to > from);
    if (SSP->Gmap.isT56selCompress)
    {
        T6sel_Alleles.extend(from, to, SourceChromo.T6sel_Alleles);
    } else
    {
        std::vector<uint32_t>::iterator itFrom;
        if (from == 0)
        {
            itFrom = SourceChromo.T5sel_Alleles.begin();
        } else
        {
            itFrom = lower_bound(SourceChromo.T5sel_Alleles.begin(), SourceChromo.T5sel_Alleles.end(), from);
        }
        
        std::vector<uint32_t>::iterator itTo;
        if (to == SSP->Gmap.T5sel_nbLoci)
        {
            itTo   = SourceChromo.T5sel_Alleles.end();
        } else
        {
            itTo   = lower_bound(itFrom, SourceChromo.T5sel_Alleles.end(), to);
        }
        
        
        T5sel_Alleles.insert(T5sel_Alleles.end(), itFrom, itTo);
    }   
}

template<typename INT>
void Haplotype::copyIntoT56ntrl(INT from, INT to, Haplotype& SourceChromo)
{
    // Must have been emptied first

    assert(from >= 0);
    assert(to > from);
    if (SSP->Gmap.isT56ntrlCompress)
    {
        T6ntrl_Alleles.extend(from, to, SourceChromo.T6ntrl_Alleles);
        
    } else
    {
        std::vector<uint32_t>::iterator itFrom;
        if (from == 0)
        {
            itFrom = SourceChromo.T5ntrl_Alleles.begin();
        } else
        {
            itFrom = lower_bound(SourceChromo.T5ntrl_Alleles.begin(), SourceChromo.T5ntrl_Alleles.end(), from);
        }
        
        std::vector<uint32_t>::iterator itTo;
        if (to == SSP->Gmap.T5ntrl_nbLoci)
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
std::vector<uint32_t>::const_iterator Haplotype::T5sel_AllelesCBegin()
{ 
    return T5sel_Alleles.cbegin();
}

std::vector<uint32_t>::const_iterator Haplotype::T5sel_AllelesCEnd()
{
    return T5sel_Alleles.cend();
}

std::vector<uint32_t>::const_iterator Haplotype::T5sel_AllelesCiterator(int locus, std::vector<uint32_t>::const_iterator from)
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

    if (locus == SSP->Gmap.T5sel_nbLoci)
    {
        return T5sel_Alleles.cend();
    }

    return std::lower_bound(from, T5sel_Alleles.cend(), locus);
}


std::vector<uint32_t>::const_iterator Haplotype::T5sel_AllelesCiterator(int locus)
{
    if (locus == 0)
    {
        return T5sel_Alleles.cbegin();
    } else if (locus == SSP->Gmap.T5sel_nbLoci)
    {
        return T5sel_Alleles.cend();
    }
    return std::lower_bound(T5sel_Alleles.cbegin(), T5sel_Alleles.cend(), locus);
}





std::vector<uint32_t>::const_iterator Haplotype::T5ntrl_AllelesCBegin()
{ 
    return T5ntrl_Alleles.cbegin();
}

std::vector<uint32_t>::const_iterator Haplotype::T5ntrl_AllelesCEnd()
{
    return T5ntrl_Alleles.cend();
}

std::vector<uint32_t>::const_iterator Haplotype::T5ntrl_AllelesCiterator(int locus, std::vector<uint32_t>::const_iterator from)
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

    if (locus == SSP->Gmap.T5ntrl_nbLoci)
    {
        return T5ntrl_Alleles.cend();
    }

    return std::lower_bound(from, T5ntrl_Alleles.cend(), locus);
}


std::vector<uint32_t>::const_iterator Haplotype::T5ntrl_AllelesCiterator(int locus)
{
    if (locus == 0)
    {
        return T5ntrl_Alleles.cbegin();
    } else if (locus == SSP->Gmap.T5ntrl_nbLoci)
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

    if (locus == SSP->Gmap.T5_nbLoci)
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
    } else if (locus == SSP->Gmap.T5_nbLoci)
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

    if (SSP->Gmap.isT56ntrlCompress)
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
    if (SSP->Gmap.isT56selCompress)
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

int Haplotype::getNbT5sel()
{
    return T5sel_Alleles.size();
}

int Haplotype::getNbT6sel()
{
    return T6sel_Alleles.size();
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
    std::vector<uint32_t>::iterator itEnd = this->T5sel_AllelesEnd();
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
        uint32_t originalSize = v.size();
        v.erase( std::unique( v.begin(), v.end() ), v.end() );
        assert(v.size() == originalSize);
        assert(std::is_sorted(v.begin(), v.end()));
    }

    {
        auto& v = T5sel_Alleles;
        uint32_t originalSize = v.size();
        v.erase( std::unique( v.begin(), v.end() ), v.end() );
        assert(v.size() == originalSize);
        assert(std::is_sorted(v.begin(), v.end()));
    }
}



ZipIterator<std::vector<uint32_t>, std::vector<uint32_t>::iterator> Haplotype::T5ntrl_AllelesBegin()
{
    return ZipIterator<std::vector<uint32_t>, std::vector<uint32_t>::iterator>(
        T5ntrl_Alleles.begin(),
        SSP->T5ntrl_flipped.begin(),
        T5ntrl_Alleles.end(),
        SSP->T5ntrl_flipped.end()
    );
}

ZipIterator<std::vector<uint32_t>, std::vector<uint32_t>::iterator> Haplotype::T5ntrl_AllelesEnd()
{
    return ZipIterator<std::vector<uint32_t>, std::vector<uint32_t>::iterator>(
        T5ntrl_Alleles.end(),
        SSP->T5ntrl_flipped.end(),
        T5ntrl_Alleles.end(),
        SSP->T5ntrl_flipped.end()
    );
}

std::vector<uint32_t>::iterator Haplotype::T5sel_AllelesBegin()
{
    return T5sel_Alleles.begin();
}

std::vector<uint32_t>::iterator Haplotype::T5sel_AllelesEnd()
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



ZipIterator<std::vector<uint32_t>, std::vector<uint32_t>::iterator> Haplotype::T56ntrl_AllelesBegin(ZipIterator<std::vector<uint32_t>, std::vector<uint32_t>::iterator> nothing)
{
    (void) nothing;
    return T5ntrl_AllelesBegin();
}

ZipIterator<std::vector<uint32_t>, std::vector<uint32_t>::iterator> Haplotype::T56ntrl_AllelesEnd(ZipIterator<std::vector<uint32_t>, std::vector<uint32_t>::iterator> nothing)
{
    (void) nothing;
    return T5ntrl_AllelesEnd();
}


std::vector<uint32_t>::iterator Haplotype::T56sel_AllelesBegin(std::vector<uint32_t>::iterator& nothing)
{
    (void) nothing;
    return T5sel_AllelesBegin();
}

std::vector<uint32_t>::iterator Haplotype::T56sel_AllelesEnd(std::vector<uint32_t>::iterator& nothing)
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


ZipIterator<std::vector<uint32_t>, std::vector<uint32_t>::iterator> Haplotype::T56ntrl_AllelesIterator(ZipIterator<std::vector<uint32_t>, std::vector<uint32_t>::iterator> nothing, unsigned value)
{
    (void) nothing;
    return ZipIterator<std::vector<uint32_t>, std::vector<uint32_t>::iterator>(
        std::lower_bound(T5ntrl_Alleles.begin(), T5ntrl_Alleles.end(), value),
        std::lower_bound(SSP->T5ntrl_flipped.begin(), SSP->T5ntrl_flipped.end(), value),
        T5ntrl_Alleles.end(),
        SSP->T5ntrl_flipped.end()
    );
    
}

std::vector<uint32_t>::iterator Haplotype::T56sel_AllelesIterator(std::vector<uint32_t>::iterator& nothing, unsigned value)
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




std::vector<uint32_t> Haplotype::getT56ntrlTrueHaplotype()
{
    std::vector<uint32_t> r;
    if (SSP->Gmap.T56ntrl_nbLoci)
    {
        if (SSP->Gmap.isT56ntrlCompress)
        {
            auto it = this->T6ntrl_AllelesBegin();
            auto itEnd = this->T6ntrl_AllelesEnd();
            for (; it != itEnd ; ++it)
            {
                r.push_back(*it);
                assert(r.back() < SSP->Gmap.T6ntrl_nbLoci);
            }
        } else
        {
            auto it = this->T5ntrl_AllelesBegin();
            auto itEnd = this->T5ntrl_AllelesEnd();
            for (; it != itEnd ; ++it)
            {
                r.push_back(*it);
                assert(r.back() < SSP->Gmap.T5ntrl_nbLoci);
            }
        }
    }
    
    return r;
}

std::vector<uint32_t> Haplotype::getT56selTrueHaplotype() 
{
    std::vector<uint32_t> r;
    if (SSP->Gmap.T56sel_nbLoci)
    {
        if (SSP->Gmap.isT56selCompress)
        {
            auto it = this->T6sel_AllelesBegin();
            auto itEnd = this->T6sel_AllelesEnd();
            for (; it != itEnd ; ++it)
            {
                r.push_back(*it);
                assert(r.back() < SSP->Gmap.T6sel_nbLoci);
            } 
        } else
        {
            auto it = this->T5sel_AllelesBegin();
            auto itEnd = this->T5sel_AllelesEnd();
            for (; it != itEnd ; ++it)
            {
                r.push_back(*it);
                assert(r.back() < SSP->Gmap.T5sel_nbLoci);
            }
        }
    }
    
    return r;
}


template<typename ITERATOR, typename CONTAINER, typename INT>
char Haplotype::T56finishMutation(ITERATOR& haploP, CONTAINER& container, INT MutPosition)
{
    /*
    info = 0; // Mutation happened and is now not found in container anymore
    info = 1; // Mutation happened and is now found in container 
    info = -1; // Mutation did not happen
    */

    char info;
    if (haploP != container.end())
    {
        if ( MutPosition == *haploP)
        {
            // found
            if (SSP->T56_mutDirection != 1)
            {
                container.erase(haploP);
                info = 0;
            } else
            {
                info = -1;
            }
        } else
        {
            if (SSP->T56_mutDirection != 0)
            {
                container.insert(haploP, MutPosition);
                info = 1;
            } else
            {
                info = -1;
            }
        }
    } else
    {
        // Not found
        if (SSP->T56_mutDirection != 0)
        {
            container.push_back(MutPosition);
            info = 1;
        } else
        {
            info = -1;
        }
    }

    return info;
}

template<typename INT>
void Haplotype::updateFitnessAfterT56Mutation(INT MutPosition, bool isNowFoundInAlleles, int Habitat)
{
    int fitnessMapIndex = SSP->FromLocusToFitnessMapIndex[SSP->Gmap.FromT56selLocusToLocus(MutPosition)];

    fitnesstype w = this->getW_T56(fitnessMapIndex);
    
    if (w != -1.0)
    {
        /*// Was mutation added or removed?
        bool wasMutationAdded;
        if (SSP->Gmap.isT56selCompress)
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


            assert(w >= 0.0);
        } else
        {
            if (SSP->T56_FitnessEffects[Habitat][MutPosition] == 0.0)
            {
                w = -1.0;
            } else
            {
                w /= SSP->T56_FitnessEffects[Habitat][MutPosition];
                assert(w >= 0.0);
            }
        }
        this->setW_T56( w,   fitnessMapIndex );
    }
}

uint32_t Haplotype::nbT56muts()
{
    uint32_t r;

    if (SSP->Gmap.isT56ntrlCompress)
    {
        r = T6ntrl_Alleles.size();
    } else
    {
        r = T5ntrl_Alleles.size();
    }

    if (SSP->Gmap.isT56selCompress)
    {
        r += T6sel_Alleles.size();
    } else
    {
        r += T5sel_Alleles.size();
    }

    return r;
}

uint32_t Haplotype::nbT56muts(int fitnessMapIndex)
{
    uint32_t r = 0;

    if (SSP->Gmap.T56ntrl_nbLoci)
    {
        if (SSP->Gmap.isT56ntrlCompress)
        {
            std::cout << "Sorry debugging function 'uint32_t Haplotype::nbT56muts(int fitnessMapIndex)' has not been defined for T6 loci\n";
            abort();
        } else
        {
            if (T5ntrl_Alleles.size() != 0)
            {
                auto from = std::lower_bound(SSP->FromLocusToFitnessMapIndex.begin(), SSP->FromLocusToFitnessMapIndex.end(), fitnessMapIndex);
                auto to   = std::upper_bound(from, SSP->FromLocusToFitnessMapIndex.end(), fitnessMapIndex);
                uint32_t locusFrom = from - SSP->FromLocusToFitnessMapIndex.begin();
                uint32_t locusTo = to - SSP->FromLocusToFitnessMapIndex.begin();
                r += 
                    std::lower_bound(T5ntrl_Alleles.begin(), T5ntrl_Alleles.end(), locusTo)
                    -
                    std::lower_bound(T5ntrl_Alleles.begin(), T5ntrl_Alleles.end(), locusFrom)
                ;
            }
        }
    }


    if (SSP->Gmap.T56sel_nbLoci)
    {
        if (SSP->Gmap.isT56selCompress)
        {
            std::cout << "Sorry debugging function 'uint32_t Haplotype::nbT56muts(int fitnessMapIndex)' has not been defined for T6 loci\n";
            abort();
        } else
        {
            if (T5sel_Alleles.size() != 0)
            {
                auto from = std::lower_bound(SSP->FromLocusToFitnessMapIndex.begin(), SSP->FromLocusToFitnessMapIndex.end(), fitnessMapIndex);
                auto to   = std::upper_bound(from, SSP->FromLocusToFitnessMapIndex.end(), fitnessMapIndex);
                uint32_t locusFrom = from - SSP->FromLocusToFitnessMapIndex.begin();
                uint32_t locusTo = to - SSP->FromLocusToFitnessMapIndex.begin();
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








fitnesstype Haplotype::CalculateT1FitnessMultiplicity(const int& Habitat)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'CalculateT1FitnessMultiplicity'\n";
#endif   
    assert(SSP->T1_isMultiplicitySelection);

    // Get the right fitness array
    auto& fits = SSP->T1_FitnessEffects[Habitat];
    assert(fits.size() == SSP->Gmap.T1_nbLoci);

    fitnesstype r = 1.0;
    for (uint32_t fitnessMapIndex = 0 ; fitnessMapIndex < SSP->NbElementsInFitnessMap ; ++fitnessMapIndex)
    {
        if (getW_T1(fitnessMapIndex) == -1.0)
        {
            //std::cout << "must recompute fitnessMapIndex " << fitnessMapIndex << "\n";
            fitnesstype w = 1.0;

            auto to = SSP->FromFitnessMapIndexToTXLocus[fitnessMapIndex].T1;
            //std::cout << "to " << to << "\n";
            assert(to <= SSP->Gmap.T1_nbLoci);
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
            
            assert(w>=0.0);
            setW_T1(w, fitnessMapIndex);
            r *= w;
        } else
        {
            r *= getW_T1(fitnessMapIndex);
        }
    }

    assert(r>=0.0);
    return r;
}


fitnesstype Haplotype::CalculateT2Fitness(const int& Habitat)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'CalculateT2Fitness'\n";
#endif          
    assert(SSP->T2_FitnessEffects.size() > Habitat);
    auto& fits = SSP->T2_FitnessEffects[Habitat];
    assert(fits.size() == SSP->Gmap.T2_nbLoci);


    fitnesstype r = 1.0;
    int T2_locus = 0;
    for (uint32_t fitnessMapIndex = 0 ; fitnessMapIndex < SSP->NbElementsInFitnessMap ; ++fitnessMapIndex)
    {
        int T2_locusTo = SSP->FromFitnessMapIndexToTXLocus[fitnessMapIndex].T2;
        if (getW_T2(fitnessMapIndex) == -1.0)
        {
            fitnesstype w = 1.0;
            for ( ; T2_locus <= T2_locusTo ; ++T2_locus )
            {
                //assert(SSP->T2_FitnessEffects.size() > Habitat);
                //assert(SSP->T2_FitnessEffects[Habitat].size() > char_index);
                w *= pow(fits[T2_locus], getT2_Allele(T2_locus));
            }
            assert(w>=0.0);
            setW_T2(w, fitnessMapIndex);
            r *= w;
            //if (W_T2haplo0!=1.0) std::cout << "W_T2haplo0 = " << W_T2haplo0 << std::endl;

        } else
        {
            r *= getW_T2(fitnessMapIndex);
            T2_locus = T2_locusTo;
        }          
    }
    
    assert(r >= 0.0);
    return r;
}


fitnesstype Haplotype::CalculateT5FitnessMultiplicity(const int& Habitat)
{

    auto& fits = SSP->T56_FitnessEffects[Habitat];
    assert(fits.size() == SSP->Gmap.T56sel_nbLoci);

    std::vector<uint32_t>::iterator it    = T5sel_Alleles.begin();
    std::vector<uint32_t>::iterator itEnd = T5sel_Alleles.end();

    bool needToRecomputeIt = false;

    fitnesstype r = 1.0;
    for (uint32_t fitnessMapIndex = 0 ; fitnessMapIndex < SSP->NbElementsInFitnessMap ; ++fitnessMapIndex)
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
            fitnesstype w = 1.0;

            for (; it != itEnd && *it < T56_locusTo ; ++it)
            {                
                w *= fits[*it];
            }
            needToRecomputeIt = false;
            
            assert(w>=0.0);
            setW_T56(w, fitnessMapIndex);
            r *= w;
        }
    }
    return r;
}

fitnesstype Haplotype::CalculateT6FitnessMultiplicity(const int& Habitat)
{

    auto& fits = SSP->T56_FitnessEffects[Habitat];
    assert(fits.size() == SSP->Gmap.T56sel_nbLoci);

    CompressedSortedDeque::iterator it    = T6sel_Alleles.begin();
    CompressedSortedDeque::iterator itEnd = T6sel_Alleles.end();

    bool needToRecomputeIt = false;

    fitnesstype r = 1.0;
    for (uint32_t fitnessMapIndex = 0 ; fitnessMapIndex < SSP->NbElementsInFitnessMap ; ++fitnessMapIndex)
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
            fitnesstype w = 1.0;

            for (; it != itEnd && *it < T56_locusTo ; ++it)
            {                
                w *= fits[*it];
            }
            needToRecomputeIt = false;
            
            assert(w>=0.0);
            setW_T56(w, fitnessMapIndex);
            r *= w;
        }
    }
    return r;
}



fitnesstype Haplotype::CalculateT56FitnessMultiplicity(const int& Habitat)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'CalculateT56FitnessMultiplicity'\n";
#endif   


    fitnesstype r;
    if (SSP->Gmap.T5sel_nbLoci)
    {
        r = CalculateT5FitnessMultiplicity(Habitat);
    } else
    {
        assert(SSP->Gmap.T6sel_nbLoci);
        r = CalculateT6FitnessMultiplicity(Habitat);
    }


    assert(r >= 0.0);
    return r;
}








fitnesstype Haplotype::CalculateT1FitnessMultiplicityOnSubsetOfLoci(const int& Habitat, const std::vector<int>& LociSet)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'CalculateT1FitnessMultiplicityOnSubsetOfLoci'\n";
#endif   

    assert(SSP->T1_isMultiplicitySelection);

    // Get the right fitness array
    auto& fits = SSP->T1_FitnessEffects[Habitat];
    assert(fits.size() == SSP->Gmap.T1_nbLoci);


    fitnesstype r = 1.0;
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

    assert(r>=0.0);

    return r;
}

fitnesstype Haplotype::CalculateT2FitnessOnSubsetOfLoci(const int& Habitat, const std::vector<int>& LociSet)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'CalculateT2FitnessOnSubsetOfLoci'\n";
#endif          
    assert(SSP->T2_FitnessEffects.size() > Habitat);
    auto& fits = SSP->T2_FitnessEffects[Habitat];
    assert(fits.size() == SSP->Gmap.T2_nbLoci);

    fitnesstype r = 1.0;
    for (const int locus : LociSet)
    {
        r *= pow(fits[locus], getT2_Allele(locus));
    }

    //if (r!=1.0)  std::cout << "T2 fit: = " << r << " ";
    assert(r>=0.0);
    return r;
}


template<typename Iterator>
fitnesstype Haplotype::CalculateT56FitnessMultiplicityOnSubsetOfLoci(const int& Habitat, const std::vector<int>& LociSet, Iterator it, Iterator itEnd)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'CalculateT56FitnessMultiplicityOnSubsetOfLoci'\n";
#endif   

    assert(SSP->T56_isMultiplicitySelection);

    // Get the right fitness array
    auto& fits = SSP->T56_FitnessEffects[Habitat];
    assert(fits.size() == SSP->Gmap.T5sel_nbLoci);

    fitnesstype r = 1.0;

    assert(LociSet.size() > 0);
    uint32_t locusIndex = 0;
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



fitnesstype Haplotype::CalculateT56FitnessMultiplicityOnSubsetOfLoci(const int& Habitat, const std::vector<int>& LociSet)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'CalculateT56FitnessMultiplicityOnSubsetOfLoci'\n";
#endif   

    fitnesstype r;
    if (SSP->Gmap.T5sel_nbLoci)
    {
        r = CalculateT56FitnessMultiplicityOnSubsetOfLoci(Habitat, LociSet, T5sel_Alleles.begin(), T5sel_Alleles.end());
    } else
    {
        r = CalculateT56FitnessMultiplicityOnSubsetOfLoci(Habitat, LociSet, T6sel_Alleles.begin(), T6sel_Alleles.end());
    }
    
    assert(r>=0.0);
    return r;
}


/*bool Haplotype::getT56_Allele(const int locus) // That's a slow function that should not be used too often!
{
    assert(locus < SSP->Gmap.T56_nbLoci);
    auto T56ntrlLocus = FromLocusToTXLocus[locus].T56ntrl;
    SSP->Gmap.isT56neutral[]
    if (SSP->Gmap.FromT56LocusToT56genderLocus(locus).type == "T56ntrl")
    {
        auto T56ntrlLocus = SSP->Gmap.FromT56LocusToT56genderLocus(locus).T56ntrl;
        if (SSP->Gmap.isT56ntrlCompress)
        {

        } else
        {

        }
    } else
    {
        assert(SSP->Gmap.FromT56LocusToT56genderLocus(locus).type == "T56sel");
        auto T56selLocus = SSP->Gmap.FromT56LocusToT56genderLocus(locus).T56sel;
        if (SSP->Gmap.isT56selCompress)
        {

        } else
        {

        }

    }
}*/

void Haplotype::shrink_to_fitT56()
{
    if (SSP->Gmap.T5ntrl_nbLoci)
    {
        //T5ntrl_Alleles.clear();
        T5ntrl_Alleles.shrink_to_fit();
    }

    if (SSP->Gmap.T5sel_nbLoci)
    {
        //T5sel_Alleles.clear();
        T5sel_Alleles.shrink_to_fit();
    }

    /*
    if (SSP->Gmap.T6ntrl_nbLoci)
    {
        T6ntrl_Alleles.clear();
        //T6ntrl_Alleles.shrink_to_fit();
    }
    */

    /*
    if (SSP->Gmap.T6sel_nbLoci)
    {
        T6sel_Alleles.clear();
        //T6sel_Alleles.shrink_to_fit();
    }
    */

    if (SSP->T56_isMultiplicitySelection)
    {
        W_T56.shrink_to_fit();
    }
}

void Haplotype::freeT56Memory()
{
    if (SSP->Gmap.T5ntrl_nbLoci)
    {
        T5ntrl_Alleles.clear();
        T5ntrl_Alleles.shrink_to_fit();
    }

    if (SSP->Gmap.T5sel_nbLoci)
    {
        T5sel_Alleles.clear();
        T5sel_Alleles.shrink_to_fit();
    }

    if (SSP->Gmap.T6ntrl_nbLoci)
    {
        T6ntrl_Alleles.clear();
        //T6ntrl_Alleles.shrink_to_fit();
    }

    if (SSP->Gmap.T6sel_nbLoci)
    {
        T6sel_Alleles.clear();
        //T6sel_Alleles.shrink_to_fit();
    }

    if (SSP->T56_isMultiplicitySelection)
    {
        W_T56.clear();
        W_T56.shrink_to_fit();
    }
}

void Haplotype::dealWithT8Info(std::pair<T8ID_type, fitnesstype> pair)
{
    /*
    std::cout << "Dealing with T8 info\n";
    std::cout << "pair.second = " << pair.second <<"\n";
    */
    this->T8ID = pair.first;
    this->setW_T8(pair.second);
}

/*
template<typename T>
void Haplotype::insertEraseAtPositions(
    std::vector<T>& x,                      
    std::vector<size_t>& toInsert,
    std::vector<size_t>& toErase
)
{
    
}
*/