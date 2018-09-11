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

class LifeCycle
{
private:
    static void clone(Individual& parent, Individual& offspring, int Habitat, bool HasChangedHabitat);
    static void recombination(Individual& parent,Haplotype& TransmittedChrom);
    static int  recombination_nbRecs();
    static void recombination_RecPositions(int& nbRecs, std::vector<int>& breakpoints, Individual& parent);
    static void recombination_copyOver(Individual& parent,Haplotype& TransmittedChrom, std::vector<int>& breakpoints, int segregationInfo);
    static void Mutate_T1(Haplotype& TransmittedChrom, int Habitat);
    static void Mutate_T2(Haplotype& TransmittedChrom, int Habitat);
    static void Mutate_T3(Haplotype& TransmittedChrom);
    static void Mutate(Haplotype& TransmittedChrom, int Habitat);

public:
    static void BREEDING_SELECTION_DISPERSAL(Pop& Offspring_pop, Pop& Parent_pop); // because called from main
    
};

