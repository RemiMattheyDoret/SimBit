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

class ResetGeneticsEvent
{
friend class ResetGenetics;
private:
    int generation;
    char mutationType;
    // 0: set T1 to 0
    // 1: set T1 to 1
    // 2: toggle T1
    // T2 and T3 are just reset to zero
    std::vector<int> T1loci;
    std::vector<int> T2loci;
    std::vector<int> T3loci;
    std::vector<int> patches;
    std::vector<int> haplotypes;
    std::vector<std::vector<int>> individuals;

public:
    ResetGeneticsEvent(
        int g,
        char MT,
        std::vector<int> T1l,
        std::vector<int> T2l,
        std::vector<int> T3l,
        std::vector<int> pa,
        std::vector<int> hap,
        std::vector<std::vector<int>> inds
        );

};

/*
######################################################################
######################################################################
######################################################################
*/

class ResetGenetics
{
private:
    std::vector<ResetGeneticsEvent> events;

public:
    void addEvent(ResetGeneticsEvent event);
    void resetPopIfNeeded(Pop& pop);
    bool isGeneration();
};
