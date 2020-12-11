
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


class MigrationEvent
{
public:
    double nbInds; // must be a double to avoid rounding issues
    unsigned comingFrom;

    //MigrationEvent(unsigned a, unsigned b):nbInds(a), comingFrom(b){}
    MigrationEvent(double& a, unsigned& b):nbInds(a), comingFrom(b){}
    //MigrationEvent(unsigned&& a, unsigned&& b):nbInds(a), comingFrom(b){}
};



class ListMigrationEvents
{
private:
    std::vector<MigrationEvent> migrantEvents;
    double sumOfIncomingMigrants = 0.0;
    
public:
    void addMigrationEvent(double nbMigrants, unsigned patch_from)
    {
        if (nbMigrants > 0.0)
        {
            migrantEvents.push_back({nbMigrants, patch_from});
            sumOfIncomingMigrants += nbMigrants;
        }
    }

    MigrationEvent& getMigrationEvent(unsigned& migrantEventIndex)
    {
        return migrantEvents[migrantEventIndex];
    }

    unsigned size()
    {
        return migrantEvents.size();
    }

    double getSumOfIncomingMigrants()
    {
        return sumOfIncomingMigrants;
    }
};



class DispersalData
{
private:
    void getMigrationEventsForEachDestinationPatch(std::vector<std::vector<std::vector<double>>>* CumSumFits_p, std::vector<ListMigrationEvents>& migrationEventsForEachDestinationPatch);
    double computeNbOffspringsProducedInPatch(const unsigned patch_from, const double n_t, const double rn_t, const double r);
    void computeBackwardMigrationRates_from_migrationEventsForEachDestinationPatch(std::vector<ListMigrationEvents>& migrationEventsForEachDestinationPatch);
    bool isBackwardRatesReceived;

public:
    std::vector<std::vector<std::vector<double>>>   __forwardMigration;
    std::vector<std::vector<std::vector<unsigned>>> __forwardMigrationIndex;

    std::vector<std::vector<double>>                forwardMigration;       // FullFormForwardMigration[from][fake_to] returns forward migration rate
    std::vector<std::vector<unsigned>>              forwardMigrationIndex;  // FullFormForwardMigration[from][fake_to] returns index to

    std::vector<int>                                nextGenerationPatchSizes;

    //std::vector<std::vector<double>>                BackwardMigrationOriginal;        // DispMatOriginal     [to][fake_from] returns backward migration
    //std::vector<std::vector<int>>                   BackwardMigrationIndexOriginal;   // DispMatIndexOriginal[to][fake_from] returns index

    std::vector<std::vector<double>>                BackwardMigration;                // DispMatOriginal     [to][fake_from] returns backward migration
    std::vector<std::vector<int>>                   BackwardMigrationIndex;           // DispMatIndexOriginal[to][fake_from] returns index

    bool                                            DispWeightByFitness;

    
    const std::vector<int>& setBackwardMigrationIfNeededAndGetNextGenerationPatchSizes(std::vector<std::vector<std::vector<double>>>& CumSumFits);
    void setOriginalBackwardMigrationIfNeeded();

    //void FromFullFormForwardSetReducedFormForward(std::vector<std::vector<double>> FullFormForwardMigration);

    std::vector<std::vector<double>> FromProbaLineToFullFormForwardMigration(std::vector<double>& probs,int center,int CurrentPatchNumber);
    void readDispMat(InputReader& input);
    void pushBack__ForwardMigrationRate(const std::vector<std::vector<double>>& FFFM, const int& CurrentPatchNumber);
    //void updateDispersalData();
    void print();

    
};


