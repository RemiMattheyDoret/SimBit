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



class Pop
{
private:
    std::vector<Patch> patches;
    static std::vector<int>  T2LociToCorrect;
    std::vector<int>         indexFirstMale;
    
public:
    std::vector<std::vector<std::vector<double>>> CumSumFits; // each patch, each gender, each individual fitness

    Patch& getPatch(const int& patch_index);
    int getNbPatches();
    void AddPatch(Patch& newPatch);
    void RemoveLastPatch();
    int SelectionParent(int& patch_from, int sex);
    void CalculateFitnesses();
    int SelectionOriginPatch(int& patch_to);
    int patchSizeNextGeneration(int patch_index);
    Pop(bool ShouldReadPopFromBinary);
    Pop();

    void PrintBinaryFile();

    static void addT2LocusToCorrect(int T2Locus);
    int correctT2Loci();

    static void updatePops(Pop& pop1, Pop& pop2, int speciesIndex, std::vector<int> previousPatchSizes);
};

