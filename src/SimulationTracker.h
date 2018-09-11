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



class SimulationTracker
{
private:
    void addMutsSortAndRemoveDuplicateOfT1SitesForFitnessNoMultiplicity();
    void addMutsSortAndRemoveDuplicateOfT1SitesForFitnessMultiplicity();
    void mergeSortRemoveDuplicates(std::vector<T1_locusDescription>& v1, std::vector<T1_locusDescription>& v2);

    bool shouldLocusBeKeptForMultiplicity(Pop& pop, int locus, int byte_index, int bit_index);

    void setToMinus1LociThatShouldNotBeTracked(Pop& pop, std::vector<T1_locusDescription>& trackedLoci );
    void setToMinus1LociThatShouldNotBeTracked(Pop& pop, std::vector<T1_locusDescription*>& trackedLoci );
    void RemoveNoFitnessEffectT1SitesForFitnessNoMultiplicity(Pop& pop);
    void RemoveNoFitnessEffectT1SitesForFitnessMultiplicity(Pop& pop);

    void printT1SitesForFitnessNoMultiplicity(); // for debug purpose
    std::vector<std::vector<T1_locusDescription>>   T1SitesForFitnessMultiplicityCalculationMutationsForNextGeneration;
    std::vector<T1_locusDescription>                T1SitesForFitnessNoMultiplicityCalculationMutationsForNextGeneration;


public:
    std::vector<T1_locusDescription>                T1SitesForFitnessNoMultiplicityCalculation;
    std::vector<std::vector<T1_locusDescription>>   T1SitesForFitnessMultiplicityCalculation;
    bool forceToRecomputeNextTime = false;

    //int MinimalSizeForReducingPolymorphicT1Sites;
    int whenDidExtinctionOccur = -1;
    Genealogy  genealogy;


    void addMutation(int byte_index, int bit_index, int MutPosition);
    void prepareT1SitesForFitness(Pop& pop);
    void prepareFitnessMapIndexForT1SitesForFitnessMultiplicity();
    void recomputeT1SitesForFitnessNoMultiplicityCalculation(Pop& pop); // public because called in resetGenetics
    void recomputeT1SitesForFitnessMultiplicityCalculation(Pop& pop); // public because called in resetGenetics
    

    std::vector<T1_locusDescription> listAllPolymorphicT1Sites(Pop& pop);


    void gotExtinct();
    void Initialization(Pop& pop);

    //void addPolymorphicT1Site();
    //void Cleanup();
    //void CheckIfT1SitesAreStillPolymorphic(Pop& pop);
};
