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

class InputReader
{
// It is actually a helper to read. It contains functions that are being called in the.readOptions types of functions in class Parameters
private:
    int VIndex;
    int VIndex_previous_habitat;
    int VIndex_previous_generation;

    std::vector<std::string> V;
    std::string ErrorMessage;

    double readDouble(const std::string& s);
    int readInt(const std::string& s, bool ComingFromMarker);

    int FindHabitatFromString(std::string s_habitat);
    int FindGenerationFromString(std::string s_generation);
    int FindSpeciesIndexFromString(std::string s_species);

public:
    InputReader(std::string entry,std::string ForErrorMessage);
    InputReader(InputReader& fullInput, int from, int to, int speciesIndex);
    void removeAlreadyRead();

    int GetNextElementInt();
    int GetNextElementBool();
    std::string GetErrorMessage();
    double GetNextElementDouble();
    std::string GetNextElementString();
    std::string PeakNextElementString();
    void skipElement();
    int GetNextHabitatMarker(const int habitat);
    int GetNextGenerationMarker(const int generation_index);
    std::vector<std::pair<int,int>> GetRangeOfIndicesForEachSpecies();
    void workDone();
    bool IsThereMoreToRead();
    bool isNextRKeyword();
    std::string print();
    std::string toString();
    int currentVIndex();
    int FindVIndexOfNextMatchingString(std::string toFind, bool throwErrorIfNotFound);
    int numberOfWordsInInput();
    void consideredFullyRead();
    void interpretKeywords();
    int getVIndex();
    size_t getSizeOfV();
    void markedAsRead();
};

