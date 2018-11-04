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


class AllParameters
{
private:
    std::vector<SpeciesSpecificParameters> SSPs;
    GeneralParameters GlobalP;

public:
    GeneralParameters* getGPaddress();
    SpeciesSpecificParameters* getSSPaddress(int speciesIndex);

    void PrintHelpToUser(OptionContainer optionContainer);
    void SetParameters (int argc, char* argv[]);
    void UpdatePopsAndParameters();
    void UpdateParameters();
    void Update(bool updatePopsToo);
    void wrapperOverSpecies(InputReader& fullInput, void (SpeciesSpecificParameters::*readFunction)(InputReader&));
    static std::string readNthLine(const std::string& filename, int N);
    //static void split_string(const std::string &s, std::vector<std::string>& elems);
    void setOptionToDefault(std::string& flag);
    void setOptionToUserInput(std::string& flag, InputReader input);
    bool isSpeciesIndexReadFromBinary(int speciesIndex);
};

