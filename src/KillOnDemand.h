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

class KillOnDemand
{
private:
    std::string functionToCall;

    int isT1LocusFixedAfterGeneration_T1locus;
    int isT1LocusFixedAfterGeneration_generation;
    int isT4LocusFixed_T4Locus;
    std::vector<uint32_t> killWhenT1HaplotypeConfigurationIsFixed_T1loci;
    std::vector<size_t> isTotalPatchSizeGreaterThan_patches;
    size_t isTotalPatchSizeGreaterThan_N;

    bool isT1LocusFixedAfterGeneration(Pop& pop);

    bool isT4LocusFixed(Pop& pop);

    bool isT1HaplotypeConfigurationFixed(Pop& pop);

    bool isTotalPatchSizeGreaterThan(Pop& pop);

public:
    static bool justAboutToKill;
    bool mustKill(Pop& pop);
    void readUserInput(InputReader& input);

    
};

