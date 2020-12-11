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

class LifeCycle
{
public:

///////////////////
// Inner classes //
///////////////////

    class HaplotypeData;
    class CoupleData;
    class ParentsData;


private:
/////////////////////
// private methods //
/////////////////////

    static void reproduceThroughSwap(Individual& parent, Haplotype& offspringHaplotype, HaplotypeData& parentHaploData, int& patch_index);
    static void reproduceThroughCopy(Individual& parent, Haplotype& TransmittedChrom, HaplotypeData& parentData, int& patch_index);
    //static int  recombination_nbRecs();
    static void recombination_RecPositions(int& nbRecs, Individual& parent);
    static void copyOver(Individual& parent,Haplotype& TransmittedChrom, int segregationIndex);
    static void Mutate_T1(Haplotype& TransmittedChrom, int Habitat);
    static void Mutate_T2(Haplotype& TransmittedChrom, int Habitat);
    static void Mutate_T3(Haplotype& TransmittedChrom);
    static void Mutate_T56(Haplotype& TransmittedChrom, int Habitat);
    static void Mutate_T7(Haplotype& TransmittedChrom);
    static void Mutate(Haplotype& TransmittedChrom, int Habitat);
    static void findAllParents(Pop& pop, const std::vector<int>& patchSizeNextGeneration);
    static CoupleData findCouple(Pop& pop, int& patch_index, int& offspring_index, ParentsData& PD, const int& nextGenPatchSize);

    static std::vector<uint32_t> recombination_breakpoints; // for performance reason. It makes things a little messy though
    static ParentsData PD;


public:
//////////////////////////////
// methods called from main //
//////////////////////////////
    static void BREEDING_SELECTION_DISPERSAL(Pop& Offspring_pop, Pop& Parent_pop); // because called from main
    
};

