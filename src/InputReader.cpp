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

void InputReader::removeWhatPrecedesIndex()
{
    V.erase(V.begin(), V.begin() + VIndex);
    VIndex = 0;
    VIndex_previous_habitat = 0;
    VIndex_previous_generation = 0;
}


InputReader::InputReader(InputReader& fullInput, int from, int to)
{
    if (to == -1)
    {
        to = fullInput.V.size();
    }
    auto first = fullInput.V.begin() + from;
    auto last = fullInput.V.begin() + to;
    this->V = std::vector<std::string>(first, last);
    this->VIndex = 0;
    this->VIndex_previous_habitat = 0;
    this->VIndex_previous_generation = 0;
    this->ErrorMessage = fullInput.ErrorMessage;
}


int InputReader::nextUntilPosition(std::vector<std::string> untils)
{
    for (auto VIndex2 = VIndex ; VIndex2 < V.size() ; ++VIndex2)
    {
        if (std::find(untils.begin(), untils.end(), V[VIndex2]) != untils.end())
        {
            return VIndex2;
        }
    }

    return V.size();
}

std::string InputReader::GetErrorMessage()
{
    return ErrorMessage;
}

std::string InputReader::print()
{
    std::string s;
    for (auto& elem : V)
    {
        std::cout << "\t" << elem << "\n";
        s += elem + " ";
    }
    std::cout << "(current index is " + std::to_string(VIndex) + ")";
    s += "(current index is " + std::to_string(VIndex) + ")";
    return this->toString();
}

std::string InputReader::toString()
{
    std::string s;
    for (auto& elem : V)
        s += elem + " ";
    return s;
}

std::vector<std::pair<int,int>> InputReader::GetRangeOfIndicesForEachSpecies()
{
    std::vector<std::pair<int,int>> r;

    // Does it start with a marker or is the marker implict
    int from = 0;
    if (V[0].at(0) == '@' && V[0].at(1) == 'S')
    {
        int SpeciesIndex = FindSpeciesIndexFromString(V[0]);
        if (SpeciesIndex != 0)
        {
            std::cout << "Message from 'InputReader method GetRangeOfIndicesForEachSpecies': " << ErrorMessage << " the first species marker (@S) is  " << V[0] << ". The first species index (or species name) must be 0 (or "<< GP->speciesNames[0] << ")." << std::endl;
            std::cout << "Please note that using species names instead of species index might be confusing as the species must come in order. The reason is that, like with @G and @H, any missing value, lead to recycling the previous value. Note also that if the very first species index (or species name) is missing, SimBit will just assume it! So, if you don't write any species specific info, then it is assume to be the same for all species.\n";
            std::cout << "Please note that if you want more than one species, then '--speciesNames' should be used before any specific specific argument.\n";
        }
        from = 1;
    }
    
    int currentSpeciesIndex = 0;


    // Get all but the last one
    for (int i = from ; i < V.size() ; i++)
    {
        if (V[i].size() >= 2)
        {
            if (V[i].at(0) == '@' && V[i].at(1) == 'S')
            {
                int SpeciesIndexAtMarker = FindSpeciesIndexFromString(V[i]);
                if (SpeciesIndexAtMarker <= currentSpeciesIndex)
                {
                    std::cout << "Message from 'InputReader method GetRangeOfIndicesForEachSpecies': " << ErrorMessage << " species indexes (or species names) don't come in the right order. The issue was discovered when reading species index " << SpeciesIndexAtMarker << "(species name = "<< GP->speciesNames[SpeciesIndexAtMarker] << " from marker " << V[0] << ")." << std::endl;
                    std::cout << "Please note that using species names instead of species index might be confusing as the species must come in order. The reason is that, like with @G and @H, any missing value, lead to recycling the previous value. Note also that if the very first species index (or species name) is missing, SimBit will just assume it! So, if you don't write any species specific info, then it is assume to be the same for all species.\n";
                    std::cout << "Please note that if you want more than one species, then '--speciesNames' should be used before any specific specific argument.\n";
                    abort();
                }
                int to = i; // to is excluded
                if (from == to)
                {
                    std::cout << "Message from 'InputReader method GetRangeOfIndicesForEachSpecies': " << ErrorMessage << " two @S markers seem to follow each other without entry in between.\n";
                    abort();
                }
                assert(from < to);
                
                std::pair<int, int> pa(from, to);

                for (int SI = currentSpeciesIndex ; SI < SpeciesIndexAtMarker ; SI++)
                {
                    r.push_back(pa);
                }
                from = to + 1; // from is included
                currentSpeciesIndex = SpeciesIndexAtMarker;
            }
        }
    }

    // Do the last one
    int to = V.size();
    assert(from < to);
    std::pair<int, int> pa(from, to);
    for (int SI = currentSpeciesIndex ; SI < GP->nbSpecies ; SI++)
    {
        r.push_back(pa);
    }

    assert(r.size() == GP->nbSpecies);
    return r;
}

bool InputReader::isNextRKeyword()
{
    if (!this->IsThereMoreToRead())
    {
        std::cout << "Message from 'InputReader method isNextRKeyword': " << ErrorMessage << " too few arguments received'" << std::endl;
        abort();
    }
    if (V[VIndex].compare("R")==0)
    {
        VIndex++;
        return true;
    }
    return false;
}

bool InputReader::IsThereMoreToRead()
{
    if (VIndex > V.size())
    {
        std::cout << "Error in InputReader::IsThereMoreToRead. It is an internal error but you might want to check your input parameters anyway.\n";
        std::cout << "input is '" + this->print() + "'\n";
        abort();
    }
    return VIndex != V.size();
}


void InputReader::workDone()
{
    if (this->IsThereMoreToRead())
    {
        std::cout << "Message from 'InputReader method workDone': " << ErrorMessage << "there seems to have more values inputted than expected.\n";
        std::cout << "The entire string being read now is '" << this->toString() << "'\n";
        if (VIndex == 0)
        {
            std::cout << "Nothing has been read from this string. This may suggest an internal bug but please check your input parameters, it is possible that in some circumstance, no input was expected." << std::endl;
        } else
        {
            std::cout << "The last token (or word) to have been read from this string was '"<< V[VIndex-1] << "'' and the following token that has not been read yet is '" << V[VIndex] << "'" << std::endl;
        }
        
        abort();
    }
}



int InputReader::GetNextGenerationMarker(const int generation_index)
{
    int r = -1;
    char markerLetter = 'G';

    // Security
    if (generation_index >= GP->__GenerationChange.size())
    {
        std::cout << "Internal Error. Message from 'InputReader method GetNextGenerationMarker': " << ErrorMessage << " received a generation index (" << generation_index << ") that is bigger than (or as big as) GP->__GenerationChange.size() ("<<GP->__GenerationChange.size()<<").\n";
        abort();
    }
    assert(VIndex <= V.size());


    // If we reached the end
    if (VIndex == V.size())
    {
        r = GP->__GenerationChange[generation_index];
        VIndex = VIndex_previous_generation;
        assert(VIndex >= 0);
        return r;
    }

    // Get string
    std::string s = V[VIndex];

    // Is the string a marker
    bool isItAMarker = false;
    if (s.size() >= 2)
    {
        if (s.at(0) == '@' && s.at(1) == markerLetter)
        {
            isItAMarker = true;
        }
    }

    // If the string is not a marker
    if (!isItAMarker)
    {
        r = GP->__GenerationChange[generation_index];
        VIndex = VIndex_previous_generation;
        return r;
    }


    // get value from marker
    r = FindGenerationFromString(s);


    // if value is too small
    if (r < GP->__GenerationChange[generation_index])
    {
        std::cout << "Message from 'InputReader method GetNextGenerationMarker': " << ErrorMessage << " received '" << s << "', but the generation " << r << " has already been received in a previous generation marker.\n";
        std::cout << "The expected generation was " << GP->__GenerationChange[generation_index] << "\n";
        abort();
    }



    // if value is too big
    if (r > GP->__GenerationChange[generation_index])
    {
        r = GP->__GenerationChange[generation_index];

        // Make sure the generation is realistic
        bool ItIsAllGood = false;
        for (int generation_index_tmp = generation_index ; generation_index_tmp < GP->__GenerationChange.size(); generation_index_tmp++)
        {
            //std::cout << "generation_index_tmp = " << generation_index_tmp << "\n";
            //std::cout << "GP->__GenerationChange[generation_index_tmp] = " << GP->__GenerationChange[generation_index_tmp] << "\n";
            if (GP->__GenerationChange[generation_index_tmp] == r)
            {
                ItIsAllGood = true;
            }
        }
        if (!ItIsAllGood)
        {
            std::cout << "Message from 'InputReader method GetNextGenerationMarker': " << ErrorMessage << " received '" << s << "', but the generation " << r << " has either not been received in '--T (--TemporalChanges)' or has already been received in a previous generation marker. Note that if the option '--T (--TemporalChanges)' is not used, then the default it to consider only @G0 (--T '0') (that is no temporal changes)\n";
            abort();
        }


        VIndex = VIndex_previous_generation;
        return r;
    } 

    // Security
    assert(r >= 0);
    assert(r == GP->__GenerationChange[generation_index]);

    // increment and assign VIndex_previous_generation
    VIndex++;
    VIndex_previous_generation = VIndex;

    return r;
}


int InputReader::GetNextHabitatMarker(const int habitat)
{
     int r = -1;
    char markerLetter = 'H';

    // Security
    if (habitat > SSP->MaxEverHabitat)
    {
        std::cout << "Internal Error. Message from 'InputReader method GetNextHabitatMarker': " << ErrorMessage << " received a habitat  (" << habitat << ") that is bigger than SSP->MaxEverHabitat ("<<SSP->MaxEverHabitat<<").\n";
        abort();
    }
    assert(VIndex <= V.size());


    // If we reached the end
    if (VIndex == V.size())
    {
        r = habitat;
        VIndex = VIndex_previous_habitat;
        assert(VIndex >= 0);
        return r;
    }

    // Get string
    std::string s = V[VIndex];

    // Is the string a marker
    bool isItAMarker = false;
    if (s.size() >= 2)
    {
        if (s.at(0) == '@' && s.at(1) == markerLetter)
        {
            isItAMarker = true;
        }
    }

    // If the string is not a marker
    if (!isItAMarker)
    {
        r = habitat;
        VIndex = VIndex_previous_habitat;
        return r;
    }


    // get value from marker
    r = FindHabitatFromString(s);


    // if value is too small
    if (r < habitat)
    {
        std::cout << "Message from 'InputReader method GetNextHabitatMarker': " << ErrorMessage << " received '" << s << "', but the habitat " << r << " has already been received in a previous habitat marker.\n";
        std::cout << "The expected habitat was " << habitat << "\n";
        abort();
    }



    // if value is too big
    if (r > habitat)
    {
        r = habitat;


        // Make sure the generation is realistic
        bool ItIsAllGood = false;
        for (int habitat_index_tmp = habitat ; habitat_index_tmp <= SSP->MaxEverHabitat; habitat_index_tmp++)
        {
    
            if (habitat_index_tmp == r)
            {
                ItIsAllGood = true;
            }
        }
        if (!ItIsAllGood)
        {
            std::cout << "Message from 'InputReader method GetNextHabitatMarker': " << ErrorMessage << " received '" << s << "', but the habitat " << r << " has either not been received in '--H (--Habitats)' or has already been received in a previous habitat marker. Note that if the option '--H (--Habitats)' is not used, then the default it to consider only @H0 (--H unif 0) (that is all patches share the same habitat)\n";
            abort();
        }


        VIndex = VIndex_previous_habitat;
        return r;
    } 

    // Security
    assert(r >= 0);
    assert(r == habitat);

    // increment and assign VIndex_previous_habitat
    VIndex++;
    VIndex_previous_habitat = VIndex;

    return r;

}




std::pair<int,size_t> InputReader::GetNextLocusInfo()
{
    if (pairsInfo.nb > 0)
    {
        if (pairsInfo.isNextFirst)
        {
            pairsInfo.isNextFirst = false;
            return {pairsInfo.firstType,1};
        } else
        {
            pairsInfo.isNextFirst = true;
            --(pairsInfo.nb);
            if (pairsInfo.nb == 0)
            {
                this->skipElement(); // Skip nb elements now only (it was only peeked before) to make sure isThereMoreToRead says yes.
            }
            return {pairsInfo.secondType,1};
        }
    }


    auto s_type = this->GetNextElementString();

    int type;
    int nbElements;
    if (s_type.size() == 2)
    {
        type = (int) std::stod(s_type.substr(1));
        nbElements = this->GetNextElementInt();
    } else if (s_type.size() == 8)
    {
        if (s_type.substr(2,4) != "pair")
        {
            std::cout << "For option --L (--Loci), received unknown type " << s_type << ". Only types accepted are T1, T2, T3, T4, T5 and T7. Note that the 'T' can be lower case (t1, t2, ...) and can be ignored (1, 2, ...). Note also that T5 might be compressed to T6 (see compression options) but you must still call it 'T5' and not 'T6' (not 'T56' either as SimBit does internally). Finally notet that possibility to indicated two loci ate once with a something like 'T5pairT4'.\n";
            abort();
        }
        
        pairsInfo.firstType = (int) std::stod(s_type.substr(1));
        pairsInfo.secondType = (int) std::stod(s_type.substr(7));
        pairsInfo.isNextFirst = false; // because the first is returned now
        pairsInfo.nb = this->PeakNextElementInt();
        type = pairsInfo.firstType;
        nbElements = 1;
    } else
    {
        if (s_type.size() != 1)
        {
            std::cout << "For option --L (--Loci), received unknown type " << s_type << ". Only types accepted are T1, T2, T3, T4, T5 and T7. Note that the 'T' can be lower case (t1, t2, ...) and can be ignored (1, 2, ...). Note also that T5 might be compressed to T6 (see compression options) but you must still call it 'T5' and not 'T6' (not 'T56' either as SimBit does internally). Finally notet that possibility to indicated two loci ate once with a something like 'T5pairT4'.\n";
            abort();
        }
        type = (int) std::stod(s_type);
        nbElements = this->GetNextElementInt();
    }

    return {type, nbElements};
}


bool InputReader::GetNextElementBool()
{
    std::string r = this->PeakNextElementString();
    if (r.at(0)=='@')
    {
        std::cout << "Message from 'InputReader method GetNextElementBool': " << ErrorMessage << " Expected a string that contains some boolean information (such as 'f', 't', 'false', '0', '1', ...) that is not a species, habitat or generation specific marker but received '" << r << "'" << std::endl;
        abort();
    }
    VIndex++;
    if (r == "f" || r == "false" || r == "False" || r == "FALSE" || r == "0" || r == "F")
    {
        return false;
    } else if  (r == "t" || r == "true" || r == "True" || r == "TRUE" || r == "1" || r == "T")
    {
        return true;
    } else
    {
        std::cout << "Message from 'InputReader method GetNextElementBool': " << ErrorMessage << " Expected a string that contains some boolean information (such as 'f', 't', 'false', '0', '1', ...) but received '" << r << "'" << std::endl;
        abort();   
    }
}

long long int InputReader::PeakNextElementInt()
{
    if (!this->IsThereMoreToRead())
    {
        std::cout << "Message from 'InputReader method GetNextElementInt': " << ErrorMessage << " too few arguments received'" << std::endl;
        abort();
    }
    long long int r = readInt(V[VIndex], false);
    return r;
}

long long int InputReader::GetNextElementInt()
{
    if (!this->IsThereMoreToRead())
    {
        std::cout << "Message from 'InputReader method GetNextElementInt': " << ErrorMessage << " too few arguments received'" << std::endl;
        abort();
    }
    long long int r = readInt(V[VIndex], false);
    VIndex++;
    return r;
}

double InputReader::GetNextElementDouble()
{
    if (!this->IsThereMoreToRead())
    {
        std::cout << "Message from 'InputReader method GetNextElementDouble': " << ErrorMessage << " too few arguments received'" << std::endl;
        abort();
    }
    double r = readDouble(V[VIndex]);
    VIndex++;
    return r;
}

std::string InputReader::PeakNextElementString()
{
    if (!this->IsThereMoreToRead() )
    {
        std::cout << "Message from 'InputReader method PeakNextElementString': " << ErrorMessage << " too few arguments received'" << std::endl;
        abort();
    }
    std::string r = V[VIndex];
    
    return r;
}

void InputReader::skipElement()
{
    VIndex++;
}

std::string InputReader::GetNextElementString()
{
    std::string r = this->PeakNextElementString();
    if (r.at(0)=='@')
    {
        std::cout << "Message from 'InputReader method GetNextElementString': " << ErrorMessage << " Expected a 'string' that is not a species, habitat or generation specific marker but received '" << r << "'. In some occasions, this error message might also be caused by starting a series of marker with not the first marker. Something like '--N @G200 unif 100 @G0 unif 50' instead of '--N @G0 unif 50 @G200 unif 100 '" << std::endl;
        abort();
    }
    VIndex++;
    return r;
}

void InputReader::removeAlreadyRead()
{
    V.erase(V.begin(), V.begin() + VIndex);
    VIndex = 0;
    VIndex_previous_habitat = 0;
    VIndex_previous_generation = 0;
}

InputReader::InputReader(InputReader& fullInput, int from, int to, int speciesIndex)
: VIndex(0), VIndex_previous_habitat(0), VIndex_previous_generation(0)
{
    if (speciesIndex != -1)
    {
        ErrorMessage = fullInput.ErrorMessage + " for species " + GP->speciesNames[speciesIndex] + " (species index = " + std::to_string(speciesIndex) + "), ";
    } else
    {
        ErrorMessage = fullInput.ErrorMessage;
    }
        

    this->VIndex = 0;
    this->VIndex_previous_habitat = 0;
    this->VIndex_previous_generation = 0;

    if (from > to)
    {
        std::cout << ErrorMessage << " internal error in 'InputReader::InputReader(InputReader& fullInput, int from, int to, int speciesIndex)'. The variable 'from' is greater than the variable 'to' (from = "<<from<<", to = "<<to<<").\n";
        abort();
    }

    for (int i = from ; i < to ; i++)
    {
        V.push_back(fullInput.V[i]);
        //std::cout << i << " -> " <<fullInput.V[i] << "\n";
        //std::cout << V.size() << "\n";
    }

    if (to - from > V.size())
    {
        std::cout << ErrorMessage << " internal error in 'InputReader::InputReader(InputReader& fullInput, int from, int to, int speciesIndex)'. The variable 'to' is greater than the size of V (the number of words in input; from = "<<from<<", to = "<<to<<", V.size() = "<<V.size()<<").\n";
        abort();
    }
}


InputReader::InputReader(std::string entry, std::string ForErrorMessage)
: VIndex(0), VIndex_previous_habitat(0), VIndex_previous_generation(0), ErrorMessage(ForErrorMessage)
{
    //std::cout << "Building input reader from entry " << entry << "\n";
    if (entry.size()==0)
    {
        std::cout << "Error Message from the 'InputReader constructor': " << ErrorMessage << " no arguments has been received\n";
        abort();
    }

    assert(entry[0] != ' ');

    std::string word;
    for( size_t i = 0 ; i < entry.size() ; )
    {
        if( entry[i] == ' ' )
        {
            if (word.size())
            {
                V.push_back(word);
                word.clear();
            }
            ++i;
        } else if (entry[i] == '\"' || entry[i] == '\'' )
        {
            //std::cout << "entry["<<i<<"] = " << entry[i] << "\n";
            ++i;
            if (i == entry.size())
            {
                std::cout << "Error Message from the 'InputReader constructor': " << ErrorMessage << "received an oppening quote that at the end of the input. Received either \", \'. Note that btw SimBit does not distinguish between single and double quotes.\n";
                abort();
            }
            if (entry[i] == '\"' || entry[i] == '\'')
            {
                std::cout << "Error Message from the 'InputReader constructor': " << ErrorMessage << "received a closing quote, just after an openning one. Received either \"\", \'\', \"\' or \'\". Note btw that SimBit does not distinguish between single and double quotes.\n";
                abort();
            }
            while ( entry[i] != '\"' && entry[i] != '\'' )
            {
                //std::cout << "entry["<<i<<"] = " << entry[i] << "\n";
                word.push_back(entry[i]);
                ++i;
                if (i == entry.size())
                {
                    std::cout << "Error Message from the 'InputReader constructor': " << ErrorMessage << " received an oppening quote (followed by a number of characters) that does not close. The last character read is '" <<  entry.back() << "'. Note btw that SimBit does not distinguish between single and double quotes.\n";
                    abort();
                }
            }
            ++i;
        } else
        {
            word.push_back(entry[i]);
            ++i;
        }
    }
    if (word.size()) V.push_back(word);
    word.clear();

    if (V.size()  == 0)
    {
        std::cout << "Error Message from the 'InputReader constructor': " << ErrorMessage << "it appears that SimBit failed to build the input reader. This is likely due to an internal bug but you might want to check your input anyway.\n";
        abort();
    }


    #ifdef DEBUG
    std::cout << this->print() << "\n";
    #endif
}

double InputReader::readDouble(const std::string& s)
{
    if (s.at(0)=='@')
    {
        std::cout << "Message from 'InputReader method.readDouble': " << ErrorMessage << " Expected a 'double' value but received a time of habitat specific marker (received '" << s << "')" << std::endl;
            abort();
    }
    if (
            std::count(s.begin(), s.end(), '.') > 1
        || 
             (s.at(0) != '-' && !isdigit(s.at(0)))
        ||
            !isdigit(s.at(s.size() - 1))
    )
    {
        std::cout << "Message from 'InputReader method.readDouble': " << ErrorMessage << " Expected a 'double' value but received '" << s << "'' (error caught at the first security gate) " << std::endl;
        abort();
    }
    double r;
    try
    {
        r = std::stod(s);
    }
    catch (...)
    {
        std::cout << "Message from 'InputReader method.readDouble': " << ErrorMessage << " Expected a 'double' value but received '" << s << "'' (error caught at the second security gate) " << std::endl;
        abort();
    }
    return r;
}

long long int InputReader::readInt(const std::string& s, bool ComingFromMarker)
{
    if (s.at(0)=='@')
    {
        if (ComingFromMarker)
        {
            std::cout << "Message from 'InputReader::readInt': " << ErrorMessage << " Expected an integer value after a species, generation or habitat specific marker (@S, @H or @G) but another '@' instead (received '" << s << "')" << std::endl;
        } else
        {
            std::cout << "Message from 'InputReader::readInt': " << ErrorMessage << " Expected an integer value but received a species, generation or habitat specific marker (@S, @H or @G) (received '" << s << "')" << std::endl;
        }
        std::cout << "It is also possible that this error message comes from an input starting with a marker starting with '@' but does not start with the right one. For example, if you input '--N @G20 100 @100 200' and you forget to specify anything before '@20'. Insread you should do '--N @G0 50 @G20 100 @100 200'\n";
        abort();
    }
    if (
        //std::count(s.begin(), s.end(), '.') > 0
        //|| 
            (s.at(0) != '-' && !isdigit(s.at(0)))
        ||
            !isdigit(s.at(s.size() - 1))
    )
    {
        if (ComingFromMarker)
        {
            std::cout << "Message from 'InputReader::readInt': "<< ErrorMessage << " Expected an integer value after a species, generation or habitat specific marker (@S, @H or @G) but received '" << s << "' (error caught at the first security gate)" <<std::endl;
        } else
        {
            std::cout << "Message from 'InputReader::readInt': "<< ErrorMessage << " Expected an integer value but received '" << s << "' (error caught at the first security gate)" <<std::endl;
        }
        abort();
    }
    long long int r;
    try
    {
        double d = std::stod(s);
        double fraction = d - ((long)d);
        if (fraction > 0.00001 || fraction < -0.00001)
        {
            std::cout << "Message from 'InputReader::readInt': "<< ErrorMessage << " Expected an integer value but received '" << s << "' which seems to be a float number to SimBit" <<std::endl;
            abort();
        }
        if (d > std::numeric_limits<int>::max() || d < std::numeric_limits<int>::min())
        {
            std::cout << "Message from 'InputReader::readInt': "<< ErrorMessage << " received the entry "<< s << " when it was expecting an integer value. Sadly "<< s <<" is outside of the range of value that SimBit represents when reading input (this is a security against overflow). It might mean that SimBit should be able to deal with these numbers and then, the current code should be edited. Please let Remi know about it." <<std::endl;
            abort();
        }
        
        r = (long long int) std::stod(s);
    }
    catch(...)
    {
        if (ComingFromMarker)
        {
            std::cout << "Message from 'InputReader::readInt': "<< ErrorMessage << " Expected an integer value after a species, generation or habitat specific marker (@S, @H or @G) but received '" << s << "'' (error caught at the third security gate) " << std::endl;
        } else
        {
            std::cout << "Message from 'InputReader::readInt': "<< ErrorMessage << " Expected an integer value but received '" << s << "'' (error caught at the second security gate) " << std::endl;
        }
        abort();
    }
    
    return r;
}

int InputReader::FindSpeciesIndexFromString(std::string s_species)
{
    std::string s_SpeciesPrefix("@S");

    // Make sure the option contains '@G'
    if (s_species.find(s_SpeciesPrefix) == std::string::npos)
    {
        std::cout << "Message from 'InputReader method FindSpeciesIndexFromString': "<< ErrorMessage << "  expected info about the species index (zero based counting) or name starting with '@S' but instead, it received '" << s_species << "'" << std::endl;
        std::cout << "Please note that if you want more than one species, then '--speciesNames' should be used before any specific specific argument.\n";
        abort();
    }
    s_species.erase(0,2); // Remove the first two characters (that is '@G')

    int speciesIndex;
    if (isdigit(s_species.at(0)))
    {
       speciesIndex = readInt(s_species, true); // read the species Index
       if (speciesIndex < 0 || speciesIndex >= GP->nbSpecies) // Make sure it is a plausible value
        {
            std::cout << "Message from 'InputReader method FindSpeciesIndexFromString': "<< ErrorMessage << " incorrect speciesIndex received when receiving '" << speciesIndex << "'. Note that in --speciesNames, you have indicated " << GP->nbSpecies << " species. Remember that indices of species (as well as for habitat and genetions) are zero based counting." << std::endl;
            
            abort();
        }
    } else
    {
        std::vector<std::string>::iterator it = std::find(GP->speciesNames.begin(),GP->speciesNames.end(),s_species);
        if (it == GP->speciesNames.end())
        {
            std::cout << ErrorMessage << " after @S SimBit was expecting  a species index (going from '0' to 'nbSpecies - 1') or a species name. As the first character was not a digit, SimBit interpreted the input as the name of a species. The name received is '" << s_species << "' (without simple quote). Under the option --speciesNames, this name was not indicated. You might want to check for capitalization and for hidden characters (such as a tab) (issue found in function 'InputReader::FindSpeciesIndexFromString')\n";
            std::cout << "Please note that using species names instead of species index might be confusing as the species must come in order. The reason is that, like with @G and @H, any missing value, lead to recycling the previous value. Note also that if the very first species index (or species name) is missing, SimBit will just assume it! So, if you don't write any species specific info, then it is assume to be the same for all species.\n";
            abort();
        }
        speciesIndex = it - GP->speciesNames.begin();
    }

    return speciesIndex;
}

int InputReader::FindGenerationFromString(std::string s_generation)
{
    std::string s_GenerationPrefix("@G");

    // Make sure the option contains '@G'
    if (s_generation.find(s_GenerationPrefix) == std::string::npos)
    {
        std::cout << "Message from 'InputReader method FindGenerationFromString': "<< ErrorMessage << "  expected info about the generation starting with '@G' but instead, it received '" << s_generation << "'" << std::endl;
        abort();
    }

    s_generation.erase(0,2); // Remove the first two characters (that is '@G')
    assert(isdigit(s_generation.at(0)));
    int Generation = readInt(s_generation,true); // read the generation

    if (Generation < 0 || Generation > GP->nbGenerations) // Make sure it is a plausible value
    {
        std::cout << "Message from 'InputReader method FindGenerationFromString': "<< ErrorMessage << " incorrect generation received when receiving '" << s_generation << "'" << std::endl;
        abort();
    }

    return Generation;
}

int InputReader::FindHabitatFromString(std::string s_habitat)
{
    std::string s_habitatPrefix("@H");

    // Make sure the option contains '@H'
    if (s_habitat.find(s_habitatPrefix) == std::string::npos)
    {
        std::cout << "Message from 'InputReader method FindHabitatFromString': " << ErrorMessage << "', expected info about the Habitat starting with '@H' but instead, it received '" << s_habitat << "'  Please note that the option '--Habitat' must come before options that are habitat-specific" << std::endl;
        abort();
    }
    s_habitat.erase(0,2); // Remove the first two characters (that is '@H')
    assert(isdigit(s_habitat.at(0)));
    int Habitat = readInt(s_habitat, true); // read the generation
    assert(Habitat > -1 );
    return Habitat;
}

int InputReader::currentVIndex()
{
    return VIndex;
}


int InputReader::FindVIndexOfNextMatchingString(std::string toFind, bool throwErrorIfNotFound)
{
    for (int VIndex_tmp = VIndex ; VIndex_tmp < V.size(); VIndex_tmp++)
    {
        if (V[VIndex_tmp] == toFind)
        {
            return VIndex_tmp;
        }
    }
    // If the keyword has not been found
    if (throwErrorIfNotFound)
    {
        std::cout << "Message from 'InputReader method FindVIndexOfNextMatchingString': " << ErrorMessage << "', expected to find the string '" << toFind << "' after current word index (the current word index is "<< std::to_string(VIndex) <<" and the corresponding word is '"<< V[VIndex] <<"')." << std::endl;
        abort();
    } else
    {
        return V.size();
    }
}


int InputReader::numberOfWordsInInput()
{
    return V.size();
}

void InputReader::consideredFullyRead()
{
    VIndex = V.size();
}

void InputReader::interpretKeywords()
{
    //std::cout << "\n\n\n\nENTERING interpretKeywords with ErrorMessage '"<<ErrorMessage<<"'\n";
    //this->print();
    //std::cout << "V.size() = " << V.size() << "\n";
    //std::cout << "VIndex = " << VIndex << "\n";
    for (int vi = 0; vi < V.size(); ++vi) // vi is int because it must be allowed to go to -1
    {
        std::vector<std::string> toInsert;

        std::string currentKeyword = V[vi];
        if (currentKeyword == "seq")
        {
            if (V.size() < vi+3)
            {
                std::cout << ErrorMessage <<"received keyword 'seq' but could not find three elements (from to by) following the keyword\n";
                abort();
            }
            double from = readDouble(V[vi+1]);
            double to   = readDouble(V[vi+2]);
            double by   = readDouble(V[vi+3]);


            for (double i = from ; i <= to; i += by)
            {
                std::stringstream stream;
                stream << std::fixed << std::setprecision(20) << i;
                toInsert.push_back(stream.str());
            }
                
        } else if (currentKeyword == "seqInt")
        {
            if (V.size() <= vi + 3)
            {
                std::cout << ErrorMessage << "in InputReader::interpretKeywords, received the keyword 'seqInt' followed by less than three values. seqInt expects a 'from', a 'to' and a 'by' value\n";
                abort();
            }
            int from = readInt(V[vi+1], false);
            int to   = readInt(V[vi+2], false);
            int by   = readInt(V[vi+3], false);

            for (auto i = from ; i <= to; i += by)
            {
                toInsert.push_back(std::to_string(i));
            }
                
        } else if (currentKeyword == "rep" || currentKeyword == "repeach")
        {
            if (V.size() < vi+3)
            {
                std::cout << ErrorMessage << "received keyword 'rep' or 'repeach' but could not find two elements (whatToRepeat nbRepeats) following the keyword\n";
                abort();
            }
            auto whatToRepeat = V[vi+1];
            auto nbRepeats    = (int) std::stod(V[vi+2]);
            //std::cout << "whatToRepeat_original = " << whatToRepeat_original << "\n";
            //std::cout << "nbRepeats = " << nbRepeats << "\n";

            

            if (currentKeyword == "rep")
            {
                for (uint32_t repeat_i = 0 ; repeat_i < nbRepeats; repeat_i++)
                {
                    std::istringstream iss (whatToRepeat);
                    auto str_toRep = std::string{};

                    while (iss >> str_toRep)
                    {
                        toInsert.push_back(str_toRep);
                    }
                }
            } else if (currentKeyword == "repeach")
            {
                std::istringstream iss (whatToRepeat);
                auto str_toRep = std::string{};

                while (iss >> str_toRep)
                {
                    for (uint32_t repeat_i = 0 ; repeat_i < nbRepeats; repeat_i++)
                    {
                        toInsert.push_back(str_toRep);
                    }   
                }
            } else
            {
                std::cout << ErrorMessage << "internal error with keywords in InputReader::interpretKeywords\n";
                abort();
            }
                
        }

        // replace string
        if (currentKeyword == "rep" || currentKeyword == "repeach")
        {
            //V.erase(V.begin() + vi, V.begin() + vi + 3);
            //V.insert(std::begin(V) + vi, toInsert.begin(), toInsert.end());
            replace(V, V.begin() + vi, V.begin() + vi + 3, toInsert.begin(), toInsert.end());
        } else if (currentKeyword == "seq" || currentKeyword == "seqInt")
        {
            //V.erase(V.begin() + vi, V.begin() + vi + 4);
            //V.insert(std::begin(V) + vi, toInsert.begin(), toInsert.end());
            replace(V, V.begin() + vi, V.begin() + vi + 4, toInsert.begin(), toInsert.end());
        }        
    }
    /*
    std::cout << "\nEXITING interpretKeywords\n";
    this->print();
    std::cout << "\n\n\n";
    */
}


int InputReader::getVIndex()
{
    return VIndex;
}

uint32_t InputReader::getSizeOfV()
{
    return V.size();
}

void InputReader::markedAsRead()
{
    VIndex = V.size();
}
