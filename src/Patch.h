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


class Patch
{
private:
    std::vector<Individual> inds;
public:
    Individual& getInd(const int& ind_index);
    void setInd(Individual& ind, const int& ind_index);
    Individual& getInd(const int& ind_index, const int& patch_index);
    void removeLastIndividual();
    void AddIndividual(Individual& newInd);
    int getpatchCapacity();
    Patch();
    Patch(const int patch_index, char Abiogenesis);
    Patch(const int patch_index, bool ShouldReadPopFromBinary);
    void PrintBinaryFile(OutputFile& file, int patch_index);
    Patch(const Patch& p);
    Patch operator=(const Patch& p);
    Patch(const Patch&& p);
    Patch operator=(const Patch&& p);

    void toggleT56LociFromEveryone(std::vector<int>& T5ntrlLociToToggle, std::vector<int>& T5selLociToToggle, int Habitat);
    //Patch(Patch&&) = default;
    //Patch& operator=(Patch&& ) = default;

    std::vector<double> computeT5ntrlFrequencies();
    std::vector<double> computeT5selFrequencies();
};


