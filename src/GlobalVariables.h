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

#ifndef GLOBALVARIABLES_H // header guards
#define GLOBALVARIABLES_H    

/*
 # # # # # # # # # # # # # # # # # # # #
 # ################################### #
 # ###### Some global constants ###### #
 # ################################### #
 # # # # # # # # # # # # # # # # # # # #
 */

// Compiler variable

// global function
std::string getSimBitVersionLogo();

// allSpecies
extern std::vector<std::pair<Pop,Pop>> allSpecies;

// SimBit Version
extern const std::string VERSION;

// Global allParameters
extern AllParameters allParameters;

// global pointers to attributes of allParameters
extern SpeciesSpecificParameters* SSP;
extern GeneralParameters* GP;   // Also contains patchSizes for all species (duplicate info than from SSP)

//outputWriter

extern OutputWriter outputWriter;

// Thing
const int EIGHT = 8;
const int THREE = 3;

#ifdef DEBUG
extern int NBMIGRANTS;
extern int NBT1MUTATIONS;
#endif

#endif
