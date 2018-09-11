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


class DispersalData
{
public:
    std::vector<std::vector<std::vector<double>>>   __FullFormForwardMigration;

    std::vector<std::vector<double>>                FullFormForwardMigration; //FullFormForwardMigration[from][to]

    //std::vector<std::vector<double>>                BackwardMigrationOriginal;        // DispMatOriginal     [to][fake_from] returns backward migration
    //std::vector<std::vector<int>>                   BackwardMigrationIndexOriginal;   // DispMatIndexOriginal[to][fake_from] returns index

    std::vector<std::vector<double>>                BackwardMigration;                // DispMatOriginal     [to][fake_from] returns backward migration
    std::vector<std::vector<int>>                   BackwardMigrationIndex;           // DispMatIndexOriginal[to][fake_from] returns index

    bool                                            DispWeightByFitness;

    

    std::vector<int> SetBackwardMigrationAndGetNextGenerationPatchSizes(
        std::vector<std::vector<std::vector<double>>>& CumSumFits,
        bool computingOriginal
    );
    double ComputeEffectOfSpeciesEcologicalInteractions(int patch_index);
    double ComputeEffectiveNumberOfMigrants(
        std::vector<std::vector<std::vector<double>>>& CumSumFits,
        bool computingOriginal,
        int patch_to,
        int patch_from
    );

    void FromFullFormForwardSetReducedFormForward(std::vector<std::vector<double>> FullFormForwardMigration);

    std::vector<std::vector<double>> FromProbaLineToFullFormForwardMigration(std::vector<double>& probs,int center,int CurrentPatchNumber);
    void readDispMat(InputReader& input);
    void updateDispersalData();

    
};

