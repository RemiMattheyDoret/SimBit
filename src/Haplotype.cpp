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

unsigned char Haplotype::getT1_char(int& char_index)
{
    //assert(char_index >= 0);
    //assert(char_index < SSP->T1_nbChars);
    return T1_Alleles[char_index];
}

bool Haplotype::getT1_Allele(const int T1Locus)
{
    const int char_index = T1Locus / EIGHT;
    const int bit_index  = T1Locus % EIGHT;
    return (this->getT1_Allele(char_index, bit_index));
}

bool Haplotype::getT1_Allele(const int char_index, const int bit_index)
{
    /*
    assert(bit_index < EIGHT);
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
    //assert(char_index < SSP->T2_nbChars);
    return T2_Alleles[char_index];
}

char Haplotype::getT3_Allele(const int char_index)
{
    //assert(char_index >= 0);
    //assert(char_index < SSP->T2_nbChars);
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

void Haplotype::setT1_Allele(const int& char_index, const int& bit_index, const int& value)
{
    //assert(bit_index < EIGHT);
    //assert(bit_index >= 0);
    //assert(char_index >= 0);
    //assert(char_index < SSP->T1_nbChars);
    T1_Alleles[char_index] = (-value ^ T1_Alleles[char_index]) & (1 << bit_index);
}


void Haplotype::setT1_AlleleToOne(int& char_index, int& bit_index)
{
    //assert(bit_index < EIGHT);
    //assert(bit_index >= 0);
    //assert(char_index >= 0);
    //assert(char_index < SSP->T1_nbChars);
    T1_Alleles[char_index] |= 1 << bit_index;
}

void Haplotype::setT1_AlleleToZero(int& char_index, int& bit_index)
{;
    //assert(bit_index < EIGHT);
    //assert(bit_index >= 0);
    //assert(char_index >= 0);
    //assert(char_index < SSP->T1_nbChars);
    T1_Alleles[char_index] &= ~(1 << bit_index);
}

void Haplotype::setT2_Allele(const int char_index, const unsigned char value)
{
    //assert(char_index >= 0);
    //assert(char_index < SSP->T2_nbChars);
    T2_Alleles[char_index] = value;
}

void Haplotype::setT3_Allele(const int char_index, const char value)
{
    //assert(char_index >= 0);
    //assert(char_index < SSP->T2_nbChars);
    T3_Alleles[char_index] = value;
}


void Haplotype::toggleT1_Allele(int& MutPosition, int Habitat)
{
    #ifdef DEBUG
    NBT1MUTATIONS++;
    #endif
    //std::cout << "mutation at " << MutPosition << "\n";
    int byte_index = MutPosition / 8;
    int bit_index  = MutPosition % 8;
    //assert(bit_index < EIGHT);
    //assert(bit_index >= 0);
    //assert(byte_index >= 0);
    //assert(byte_index < SSP->T1_nbChars);
    T1_Alleles[byte_index] ^= 1 << bit_index;
    //std::cout << "line 180\n";
    SSP->simTracker.addMutation(byte_index, bit_index, MutPosition);
    //std::cout << "line 182\n";

    if (SSP->T1_isMultiplicitySelection)
    {
        //std::cout << "line 186\n";
        //std::cout << "MutPosition = " << MutPosition << "\n";
        //std::cout << "SSP->FromT1LocusToLocus.size() = " << SSP->FromT1LocusToLocus.size() << "\n";
        //assert(MutPosition < SSP->FromT1LocusToLocus.size());
        //std::cout << " SSP->FromT1LocusToLocus[MutPosition] = " << SSP->FromT1LocusToLocus[MutPosition] << "\n";
        //assert(SSP->FromT1LocusToLocus[MutPosition] < SSP->FromLocusToFitnessMapIndex.size());
        int fitnessMapIndex = SSP->FromLocusToFitnessMapIndex[SSP->FromT1LocusToLocus[MutPosition]];
        //std::cout << "line 188\n";
        double w = this->getW_T1(fitnessMapIndex);
        //std::cout << "line 189\n";
        if (w != -1.0)
        {
            //std::cout << "line 191\n";
            if (SSP->T1_isMultiplicitySelection)
            {
                // Note that Toggle already happened. If it is mutant, it is because it just happened!
                //std::cout << "line 196\n";
                if ( this->getT1_Allele(byte_index, bit_index) )
                {
                    w *= SSP->T1_FitnessEffects[Habitat][byte_index * EIGHT + bit_index];
                    assert(w >= 0.0 && w <= 1.0);
                } else
                {
                    if (SSP->T1_FitnessEffects[Habitat][byte_index * EIGHT + bit_index] == 0)
                    {
                        w = -1.0;
                    } else
                    {
                        w /= SSP->T1_FitnessEffects[Habitat][byte_index * EIGHT + bit_index]; 
                        assert(w >= 0.0 && w <= 1.0);
                    }
                }
                this->setW_T1( w,   fitnessMapIndex );
            } else
            {
                this->setW_T1(-1.0, fitnessMapIndex );
            }
        }
    }
}

void Haplotype::toggleT1_Allele(int& MutPosition)
{
    std::cout << "Internal error: I don't think I should be using this function Haplotype::toggleT1_Allele(int& MutPosition)!\n";
    abort();
    //std::cout << "mutation at " << MutPosition << "\n";
    int byte_index = MutPosition / 8;
    int bit_index  = MutPosition % 8;
    //assert(bit_index < EIGHT);
    //assert(bit_index >= 0);
    //assert(byte_index >= 0);
    //assert(byte_index < SSP->T1_nbChars);
    T1_Alleles[byte_index] ^= 1 << bit_index;
    //std::cout << "line 180\n";
    SSP->simTracker.addMutation(byte_index, bit_index, MutPosition);
    //std::cout << "line 182\n";

    if (SSP->T1_isMultiplicitySelection)
    {
        //std::cout << "line 186\n";
        //std::cout << "MutPosition = " << MutPosition << "\n";
        //std::cout << "SSP->FromT1LocusToLocus.size() = " << SSP->FromT1LocusToLocus.size() << "\n";
        //assert(MutPosition < SSP->FromT1LocusToLocus.size());
        //std::cout << " SSP->FromT1LocusToLocus[MutPosition] = " << SSP->FromT1LocusToLocus[MutPosition] << "\n";
        //assert(SSP->FromT1LocusToLocus[MutPosition] < SSP->FromLocusToFitnessMapIndex.size());
        int fitnessMapIndex = SSP->FromLocusToFitnessMapIndex[SSP->FromT1LocusToLocus[MutPosition]];
        //std::cout << "line 189\n";
        this->setW_T1(-1.0, fitnessMapIndex );
    }
}


void Haplotype::AddMutT2_Allele(int& char_index)
{
    //assert(char_index >= 0);
    //assert(char_index < SSP->T2_nbChars);
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

void Haplotype::AddMutT2_Allele(int& char_index, int Habitat)
{
    //assert(char_index >= 0);
    //assert(char_index < SSP->T2_nbChars);
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

void Haplotype::AddMutT3_Allele(int& char_index)
{
    // Add or substract while preventing
#ifdef DEBUG
std::cout << "T3_Alleles.size() = " << T3_Alleles.size() << "\n";
std::cout << "SSP->T3_nbChars = " << SSP->T3_nbChars << "\n";
assert(SSP->T3_nbChars == T3_Alleles.size());    
assert(char_index < T3_Alleles.size());
assert(char_index > 0);
#endif
    if (GP->random_0or1(GP->mt))
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
    
    int ByteFrom = from / EIGHT; // integer division. It will floor automatically
    int ByteTo   = to / EIGHT; // integer division. It will floor automatically
    assert(ByteTo <= SourceChromo.T1_Alleles.size() + 1);

    int bitIndexFrom = from % EIGHT;
    int bitIndexTo = to % EIGHT;

    if (ByteTo == SSP->T1_nbChars - 1) // This is to deal with the last extra bits. Note that if it copies a few extra bits, it is not an issue.
    {
        bitIndexTo = std::max(bitIndexTo, SSP->T1_nbBitsLastByte);
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
        for (int bit_index = bitIndexFrom ; bit_index < EIGHT ; bit_index++)
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

    // Print T5_alleles
    file.writeBinary(reinterpret_cast<const char*>(&T5_Alleles[0]), T5_Alleles.size()*sizeof(size_t));    
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
    T2_Alleles.resize(SSP->T2_nbChars);
    GP->BinaryFileToRead.read(reinterpret_cast<char*>(&T2_Alleles[0]), SSP->T2_nbChars*sizeof(unsigned char));

    //.read T3_alleles
    T3_Alleles.resize(SSP->T2_nbChars);
    GP->BinaryFileToRead.read(reinterpret_cast<char*>(&T3_Alleles[0]), SSP->T3_nbChars*sizeof(char));

    //.read T4_alleles
    if (SSP->T4_nbBits)
        std::cout << "\n\n\tWARNING: You asked for T4 alleles and asked for data to be read from binary file. As per the current version of SimBit, binary file cannot store T4 alleles data. Those loci will be initiated as perfectly unmutated\n\n";

    //.read T5_alleles
    T5_Alleles.resize(SSP->T5_nbBits);
    GP->BinaryFileToRead.read(reinterpret_cast<char*>(&T5_Alleles[0]), SSP->T5_nbBits*sizeof(size_t));

    // Initiate W_T1, W_T2 and W_T5
    if (SSP->T1_isSelection)
    {
        this->W_T1.resize(SSP->NbElementsInFitnessMap, -1.0);
    } else
    {
        this->W_T1.resize(SSP->NbElementsInFitnessMap,1.0);
    }

    if (SSP->T2_isSelection)
    {
        this->W_T2.resize(SSP->NbElementsInFitnessMap, -1.0);
    } else
    {
        this->W_T2.resize(SSP->NbElementsInFitnessMap, 1.0);
    }

    if (SSP->T5_isSelection)
    {
        this->W_T5.resize(SSP->NbElementsInFitnessMap, -1.0);
    } else
    {
        this->W_T5.resize(SSP->NbElementsInFitnessMap, 1.0);
    }

    /*
    for (auto& elem : T1_Alleles)
    {
        printf("%u ", elem);
    }
    std::cout << "\n";
    */
}

Haplotype::Haplotype(const int patch_index, char Abiogenesis)
: T1_Alleles(SSP->T1_nbChars), T2_Alleles(SSP->T2_nbChars), T3_Alleles(SSP->T3_nbChars)
{
#ifdef CALLENTRANCEFUNCTIONS
    std::cout << "Enters in 'Haplotype::Haplotype(const int patch_index, char Abiogenesis)'\n";
#endif      
    (void)Abiogenesis; // Does nothing but silence the warning that the char Abiogenesis is not used.
    assert(SSP!=nullptr);
    assert(SSP->T1_nbChars + SSP->T2_nbChars + SSP->T3_nbChars + SSP->T4_nbBits + SSP->T5_nbBits > 0);

    // Initiate W_T1 and W_T2
    if (SSP->T1_isSelection)
    {
        this->W_T1.resize(SSP->NbElementsInFitnessMap, -1.0);
    } else
    {
        this->W_T1.resize(SSP->NbElementsInFitnessMap, 1.0);
    }
    if (SSP->T2_isSelection)
    {
        this->W_T2.resize(SSP->NbElementsInFitnessMap, -1.0);
    } else
    {
        this->W_T2.resize(SSP->NbElementsInFitnessMap, 1.0);
    }
    if (SSP->T5_isSelection)
    {
        this->W_T5.resize(SSP->NbElementsInFitnessMap, -1.0);
    } else
    {
        this->W_T5.resize(SSP->NbElementsInFitnessMap, 1.0);
    }
    //std::cout << "created at W_T5.size() = " << W_T5.size() << "\n";
    
    // Initiate T1_Alleles
    
    if (SSP->T1_nbBits > 0)
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
                    bit_index_to = SSP->T1_nbBitsLastByte;
                } else
                {
                    bit_index_to = EIGHT;
                }
                for (int bit_index=0 ; bit_index < bit_index_to ; bit_index++ )
                {
                    assert(SSP->T1_Initial_AlleleFreqs[patch_index].size() > byte_index * EIGHT + bit_index);
                    if (SSP->T1_Initial_AlleleFreqs[patch_index][byte_index * EIGHT + bit_index] == 1.0)
                    {
                        this->setT1_AlleleToOne(byte_index,bit_index);
                    } else if (SSP->T1_Initial_AlleleFreqs[patch_index][byte_index * EIGHT + bit_index] == 0.0)
                    {
                        this->setT1_AlleleToZero(byte_index,bit_index);
                        
                    } else
                    {
                        double rnd = GP->random_0and1(GP->mt);
                        if (rnd < SSP->T1_Initial_AlleleFreqs[patch_index][byte_index * EIGHT + bit_index])
                        {
                            this->setT1_AlleleToOne(byte_index,bit_index);
                        }
                        else
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
    if (SSP->T2_nbChars > 0)
    {
        for (int char_index=0 ; char_index < SSP->T2_nbChars ; char_index++)
        {
            unsigned char c = 0;
            this->setT2_Allele(char_index,c);
        }
        assert(T2_Alleles.size() == SSP->T2_nbChars);
    }

    // Initiate T3_Alleles
    if (SSP->T3_nbChars > 0)
    {
        for (int char_index=0 ; char_index < SSP->T3_nbChars ; char_index++)
        {
            char c = 0;
            this->setT3_Allele(char_index,c);
        }
        assert(T3_Alleles.size() == SSP->T3_nbChars);
    }

    // Initiate T4_Alleles
    // No initialization. Assume all are fixed at 0

    // Initiate T5_Alleles
    // No initialization. Assume all are fixed at 0
}

Haplotype::Haplotype(const Haplotype& other)
{
    //std::cout << "copy constructor haplotype\n";
    T1_Alleles = other.T1_Alleles;
    T2_Alleles = other.T2_Alleles;
    T3_Alleles = other.T3_Alleles;
    T5_Alleles = other.T5_Alleles;
    W_T1 = other.W_T1;
    W_T2 = other.W_T2;
    W_T5 = other.W_T5;
}

Haplotype& Haplotype::operator=(const Haplotype& other) // copy assignment operator
{
    //std::cout << "copy assignment haplotype\n";
    T1_Alleles = other.T1_Alleles;
    T2_Alleles = other.T2_Alleles;
    T3_Alleles = other.T3_Alleles;
    T5_Alleles = other.T5_Alleles;
    W_T1 = other.W_T1;
    W_T2 = other.W_T2;
    W_T5 = other.W_T5;
    return *this;
}

Haplotype& Haplotype::operator=(Haplotype&& other) // move assignment operator
{
    //std::cout << "move assignment haplotype\n";
    T1_Alleles = std::move(other.T1_Alleles);
    T2_Alleles = std::move(other.T2_Alleles);
    T3_Alleles = std::move(other.T3_Alleles);
    T5_Alleles = std::move(other.T5_Alleles);
    W_T1 = std::move(other.W_T1);
    W_T2 = std::move(other.W_T2);
    W_T5 = std::move(other.W_T5);
    //assert(W_T1.size() == other.W_T1.size());
    return *this;
}


void Haplotype::setW_T5(double w, int fitnessMapIndex)
{
    assert(W_T5.size() > fitnessMapIndex && fitnessMapIndex >= 0);
    W_T5[fitnessMapIndex] = w;
}

void Haplotype::setAllW_T5(double w)
{
    for (int fitnessMapIndex = 0 ; fitnessMapIndex < SSP->NbElementsInFitnessMap ; fitnessMapIndex++)
    {
        this->setW_T5(w, fitnessMapIndex);
    }
}


bool Haplotype::getT5_Allele(const int locus)
{
    // not to be used too often
    return std::binary_search(T5_Alleles.begin(), T5_Alleles.end(), locus);
}


void Haplotype::setT5_Allele(const int& locus, const bool& value)
{
    auto position = std::lower_bound(T5_Alleles.begin(), T5_Alleles.end(), locus);

    if (value)
    {
        if (position == T5_Alleles.end())
        {
            T5_Alleles.push_back(locus);
        } else if (locus != (*position))
        {
            T5_Alleles.insert(position, locus);
        }
    } else
    {
        if (locus == (*position))
        {
            T5_Alleles.erase(position);
        }
    }
}


void Haplotype::setT5_AlleleToOne(int& locus)
{
    this->setT5_Allele(locus, true);
}


void Haplotype::setT5_AlleleToZero(int& locus)
{
    this->setT5_Allele(locus, false);
}


void Haplotype::toggleT5_Allele(int& MutPosition)
{
    auto position = std::lower_bound(T5_Alleles.begin(), T5_Alleles.end(), MutPosition);
    if (MutPosition != (*position))
    {
        // not found
        T5_Alleles.insert(position, MutPosition);
    } else
    {
        // found
        T5_Alleles.erase(position);
    }

    if (SSP->T5_isMultiplicitySelection)
    {
        //std::cout << "line 186\n";
        //std::cout << "MutPosition = " << MutPosition << "\n";
        //std::cout << "SSP->FromT1LocusToLocus.size() = " << SSP->FromT1LocusToLocus.size() << "\n";
        //assert(MutPosition < SSP->FromT1LocusToLocus.size());
        //std::cout << " SSP->FromT1LocusToLocus[MutPosition] = " << SSP->FromT1LocusToLocus[MutPosition] << "\n";
        //assert(SSP->FromT1LocusToLocus[MutPosition] < SSP->FromLocusToFitnessMapIndex.size());
        int fitnessMapIndex = SSP->FromLocusToFitnessMapIndex[SSP->FromT5LocusToLocus[MutPosition]];
        //std::cout << "line 189\n";
        this->setW_T5(-1.0, fitnessMapIndex );
    }
}

double Haplotype::getW_T5(int fitnessMapIndex)
{
    //std::cout << "fitnessMapIndex = " << fitnessMapIndex << "\n";
    //std::cout << "W_T5.size() = " << W_T5.size() << "\n";
    assert(fitnessMapIndex >= 0 && fitnessMapIndex < W_T5.size());
    return W_T5[fitnessMapIndex];
}



void Haplotype::toggleT5_Allele(int& MutPosition, int Habitat)
{
    //std::cout << "In Haplotype::toggleT5_Allele, MutPosition = " << MutPosition << "\n";
    auto position = std::lower_bound(T5_Alleles.begin(), T5_Alleles.end(), MutPosition);
    bool newAllele;
    //std::cout << "position is " << position - T5_Alleles.begin() << "\n";
    if (position == T5_Alleles.end())
    {
        // not found
        T5_Alleles.push_back(MutPosition);
        newAllele = true;   
    } else
    {
        if ( MutPosition != (*position))
        {
            // not found
            T5_Alleles.insert(position, MutPosition);
            newAllele = true;
            
        } else
        {
            // found
            T5_Alleles.erase(position);
            newAllele = false;
        }

    }

    if (SSP->T5_isMultiplicitySelection)
    {
        //std::cout << "MutPosition = " << MutPosition << "\n";
        //std::cout << "SSP->FromT1LocusToLocus.size() = " << SSP->FromT1LocusToLocus.size() << "\n";
        //assert(MutPosition < SSP->FromT1LocusToLocus.size());
        //std::cout << " SSP->FromT1LocusToLocus[MutPosition] = " << SSP->FromT1LocusToLocus[MutPosition] << "\n";
        //assert(SSP->FromT1LocusToLocus[MutPosition] < SSP->FromLocusToFitnessMapIndex.size());
        int fitnessMapIndex = SSP->FromLocusToFitnessMapIndex[SSP->FromT5LocusToLocus[MutPosition]];
        //std::cout << "In Haplotype::toggleT5_Allele, fitnessMapIndex = "<<fitnessMapIndex<<"\n";

        double w = this->getW_T5(fitnessMapIndex);
        if (w != -1.0)
        {
            // Note that Toggle already happened. If it is mutant, it is because it just happened!
            //std::cout << "line 196\n";
            if ( newAllele )
            {
                w *= SSP->T5_FitnessEffects[Habitat][MutPosition];
                assert(w >= 0.0 && w <= 1.0);
            } else
            {
                if (SSP->T5_FitnessEffects[Habitat][MutPosition] == 0)
                {
                    w = -1.0;
                } else
                {
                    w /= SSP->T5_FitnessEffects[Habitat][MutPosition];
                    assert(w >= 0.0 && w <= 1.0);
                }
            }
            this->setW_T5( w,   fitnessMapIndex );
        
        }
    }
}

void Haplotype::clearT5Alleles()
{
    T5_Alleles.clear();
}

void Haplotype::copyIntoT5(int from, int to, Haplotype& SourceChromo)
{
    // Must have been emptied first
    
    assert(from >= 0);
    assert(to > from);
    
    auto itFrom = lower_bound(SourceChromo.T5_AllelesCBegin(), SourceChromo.T5_AllelesCEnd(), from);
    auto itTo   = lower_bound(itFrom, SourceChromo.T5_AllelesCEnd(), to);
    
    //std::cout << "inserting " << itTo - itFrom << " elements after a T5_Alleles of " << T5_Alleles.size() << " elements\n";
    T5_Alleles.insert(T5_Alleles.end(), itFrom, itTo);
}

std::vector<size_t>::const_iterator Haplotype::T5_AllelesCBegin()
{
    return T5_Alleles.cbegin();
}

std::vector<size_t>::const_iterator Haplotype::T5_AllelesCEnd()
{
    return T5_Alleles.cend();
}

std::vector<size_t>::const_iterator Haplotype::T5_AllelesCiterator(int locus, std::vector<size_t>::const_iterator from)
{
    if (from == T5_Alleles.cend())
    {
        return from;   
    }

    if (locus == *from)
    {
        return from;
    }

    if (locus == SSP->T5_nbBits)
    {
        return T5_Alleles.cend();
    }

    return std::lower_bound(from, T5_Alleles.cend(), locus);
}

std::vector<size_t>::const_iterator Haplotype::T5_AllelesCiterator(int locus)
{
    if (locus == 0)
    {
        return T5_Alleles.cbegin();
    } else if (locus == SSP->T5_nbBits)
    {
        return T5_Alleles.cend();
    }
    return std::lower_bound(T5_Alleles.cbegin(), T5_Alleles.cend(), locus);
}

bool Haplotype::isT5Mutation(size_t locus)
{
    return std::binary_search(T5_Alleles.cbegin(), T5_Alleles.cend(), locus);
}
