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

GeneralParameters* AllParameters::getGPaddress()
{
    return &this->GlobalP;
}


SpeciesSpecificParameters* AllParameters::getSSPaddress(int speciesIndex)
{
    assert(GP == &this->GlobalP);
    assert(speciesIndex < SSPs.size());
    assert(GP->nbSpecies == SSPs.size());
    return &SSPs[speciesIndex];
}

void AllParameters::wrapperOverSpecies(InputReader& fullInput, void (SpeciesSpecificParameters::*readFunction)(InputReader& input))
{
    #ifdef DEBUG
    std::cout << "Enter in wrapperOverSpecies with input "<< fullInput.print() << "\n";
    #endif
    assert(SSPs.size() == GP->nbSpecies);

    std::vector<std::pair<int, int>> rangesToSubsetFullInput = fullInput.GetRangeOfIndicesForEachSpecies();

    assert(rangesToSubsetFullInput.size() == GP->nbSpecies);
    for (int speciesIndex = 0 ; speciesIndex < GP->nbSpecies ; speciesIndex++)
    {
        SSP = &SSPs[speciesIndex]; // It is important to reset this because InputReader uses SSP
        assert(SSP != nullptr);

        int from = rangesToSubsetFullInput[speciesIndex].first;
        int to = rangesToSubsetFullInput[speciesIndex].second;
        assert(from < to);
        InputReader input(fullInput, from, to, speciesIndex);

        (SSP->*readFunction)(input); // call method that readFunction points to on object SSPs[speciesIndex] and give it input in argument.
        input.workDone();
    }
    SSP = nullptr; // Set it back to null to avoid mistakes
    fullInput.consideredFullyRead();
}


void AllParameters::PrintHelpToUser(OptionContainer optionContainer)
{
    std::cout << "Please indicate options in the format '--option1 entry1 --option2 entry2'. See the Manual for help" << std::endl;
    optionContainer.listOptions();
}

std::string AllParameters::readNthLine(const std::string& filename, int N)
{
    N--;
    assert(N >= 0);
    std::ifstream in(filename.c_str());
    if (!in.is_open()) // test if file is open
    { 
        std::cout << "Sorry, SimBit failed to open the optionFile (the optionFile is the file where the input parameters are indicated when you do something like './SimBit f optionFile 1') '" << filename << "'" <<std::endl;
        abort();
    }

    std::string s;
    //for performance
    s.reserve(9999);    

    //skip N lines
    for(int i = 0; i < N; ++i)
        std::getline(in, s);

    std::getline(in,s);
    return s; 
}


/*void AllParameters::split_string(const std::string &s, std::vector<std::string>& elems)
{
    int elem_index = 0;
    bool JustGotASpace = false; // This object allows to ensure that if there are several spaces, then it does not keep empty substrings
    for ( unsigned long int c_index = 0 ; c_index < s.size() ; c_index++ )
    {
        // if elems need to be longer, then just add one element
        if (elem_index == elems.size())
        {
            std::string BeginOneSubstring;
            elems.push_back(BeginOneSubstring);
        }
        if (s.size() <= c_index)
        {
            std::cout << "For some unkown reason SimBit failed to split the input string into its constituants. This is likely because several set of arguments follow each other without option name either because a quote ends too early (something like --Loci '2 12 1' 2 --T1_mu 'A 0.000003 0.000002') or because an option name is missing (something like --Loci '2 12 1 2' 'A 0.000003 0.000002') . Please check your input!\n";
            abort();
        }
        if ( s[c_index] == ' ' )
        {
            if (!JustGotASpace)
            {
                elem_index ++;
            }
            JustGotASpace = true;
        } else
        {
            JustGotASpace = false;
            if (s[c_index] == '\"' || s[c_index] == '\'' )
            {
                c_index++;
                assert(s.size() > c_index && "For some unkown reason SimBit failed to split the input string into its constituants. This is likely because several set of arguments follow each other without option name either because a quote ends too early (something like --Loci '2 12 1' 2 --T1_mu 'A 0.000003 0.000002') or because an option name is missing (something like --Loci '2 12 1 2' 'A 0.000003 0.000002') . Please check your input and eventually report an easily reproducible test case to Remi if you can't find what is wrong!");
                while ( s[c_index] != '\"' && s[c_index] != '\'' )
                {
                    assert(s.size() > c_index && "For some unkown reason SimBit failed to split the input string into its constituants. This is likely because several set of arguments follow each other without option name either because a quote ends too early (something like --Loci '2 12 1' 2 --T1_mu 'A 0.000003 0.000002') or because an option name is missing (something like --Loci '2 12 1 2' 'A 0.000003 0.000002') . Please check your input and eventually report an easily reproducible test case to Remi if you can't find what is wrong!");
                    elems[elem_index] += s[c_index];
                    c_index++;
                }
            } else
            {
                elems[elem_index] += s[c_index];
            }
        }
    }
}*/


void AllParameters::SetParameters (int argc, char* argv[])
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'SetParameters'\n";
#endif   

    OptionContainer optionContainer;

    // set LogfileType to default a priori... This is actually not needed anymore
    outputWriter.LogfileType = 1;

    if (argc <= 1)
    {
        AllParameters::PrintHelpToUser(optionContainer);
        exit(1);
    }
    

    // Make sure there was at least one option after the call of the executable
    if (argc<2)
    {
        std::cout << "Received less than two options (executable incl.). Please look at the manual for help" << std::endl;
        PrintHelpToUser(optionContainer);
        exit(1);
    }

    std::string AllInputInLongString;

    std::string FileOrNot(argv[1]);
    if (
        FileOrNot.compare("optionFile") == 0 ||
        FileOrNot.compare("FILE") == 0 ||
        FileOrNot.compare("File") == 0 ||
        FileOrNot.compare("file") == 0 ||
        FileOrNot.compare("F") == 0 ||
        FileOrNot.compare("f") == 0
        )
    {
        #ifdef DEBUG            
        std::cout << "Reading options from file" << std::endl;
        #endif
        if (argc < 4)
        {
            std::cout << "received only " << argc << " options while trying to read from file. First option should be the executable, the second should be 'FILE' (or something similar), the third should be the line (counting based 1) from the file to consider.\n";
            PrintHelpToUser(optionContainer);
            exit(1);
        }

        // Get optionFilePath
        std::string optionFilePath(argv[2]);

        // Get Line Number from which parameters must be read
        std::string s_LineNumber(argv[3]);
        for (int char_index = 0 ; char_index < s_LineNumber.size() ; ++char_index)
        {
            if (!isdigit(s_LineNumber.at(char_index)))
            {
                std::cout << "When submitting parameters through a file, the option after the path must be an integer (indicates the line number from file, counting zero-based). (received " << s_LineNumber << ")\n";
                abort();
            }
        }
        int LineNumberToRead;
        {
            InputReader input(s_LineNumber, "When reading LineNumber of optionFile,");
            LineNumberToRead = input.GetNextElementInt();
            input.workDone();
        }

        #ifdef DEBUG
        std::cout << "LineNumberToRead = " << LineNumberToRead << std::endl;
        #endif

        if (LineNumberToRead <= 0)
        {
            std::cout << "When submitting parameters through a file, the option after the path must be an integer strictly greater than 0. LineNumber received is " << LineNumberToRead << "\n";
            abort();
        }

        AllInputInLongString = AllParameters::readNthLine(optionFilePath, LineNumberToRead);

        #ifdef DEBUG
        std::cout << "Line read = " << AllInputInLongString << std::endl;
        #endif           

        if (AllInputInLongString.size() < 1)
        {
            std::cout << "The " << LineNumberToRead << "th line from " << optionFilePath << " is empty.\n";
            abort();
        }
    } else
    {
        #ifdef DEBUG            
        std::cout << "Reading options from Command Line" << std::endl;
        #endif       
        for (int i = 1 ; i < argc ; i++)
            AllInputInLongString += std::string(argv[i]) + " ";
    }

    #ifdef DEBUG
    std::cout << "input received = " << AllInputInLongString << "\n";
    #endif

    // Start parsing AllInputInLongString into optionNames and entries
    for (int char_index = 0 ; char_index < AllInputInLongString.size() ; char_index++ )
    {
        if (AllInputInLongString.at(char_index) != ' ')
        {
            if (AllInputInLongString.at(char_index) == '-')
            {
                break;
            } else
            {
                std::cout << "The very first character (excl. white spaces) of the input (after the name of the executable) appears to be '" << AllInputInLongString.at(char_index) << "'. The first character should be '-' as the beginning of the first option (--optionName). In other words, the input does not start with a flag. It is possible that you are trying to indicate a file in which the arguments are found but have forgotten to write the keyword 'file' (or other equivalents such as 'f' or 'FILE') right after the name of the executable.\n";
                abort();
            }
        }
    }

    // Build UserEntries and rename the flags
    std::vector<std::pair<std::string, std::string>> UserEntries;
    bool isItVeryFirstFlag = true;
    bool searchForBeginningOfFlag = true;
    long long flag_from = -1;
    long long previous_flag_to = 0;
    std::string flag;
    //std::cout << "AllInputInLongString = '" << AllInputInLongString << "'\n";

    AllInputInLongString.erase(std::remove_if(AllInputInLongString.begin(), AllInputInLongString.end(),  [](char c){if (c=='\'' || c=='\"') return true; return false;}  ), AllInputInLongString.end());
    //AllInputInLongString.erase(std::remove(AllInputInLongString.begin(), AllInputInLongString.end(), '\''), AllInputInLongString.end());
    //AllInputInLongString.erase(std::remove(AllInputInLongString.begin(), AllInputInLongString.end(), '\"'), AllInputInLongString.end());

    for (int char_index = 0 ; char_index < AllInputInLongString.size() ; char_index++ )
    {
        //bool prout = AllInputInLongString.at(char_index) != '\"' && AllInputInLongString.at(char_index) != '\'';
        //std::cout << char_index << " " << AllInputInLongString.at(char_index) << " " << prout << std::endl;
        if (AllInputInLongString.at(char_index) != '\"' && AllInputInLongString.at(char_index) != '\'')
        {
            if (searchForBeginningOfFlag)
            {
                if (AllInputInLongString.at(char_index) == '-')
                {
                    if (char_index == AllInputInLongString.size() - 1)
                    {
                        std::cout << "Last character of input is '-'. This is weird!\n";
                        abort();
                    }
                    if (AllInputInLongString.at(char_index+1) == '-')
                    {
                        char_index++;
                        flag_from = char_index + 1; // index of first character that follows "--"
                        //std::cout << "flag starts at " << flag_from << "\n";
                        searchForBeginningOfFlag = false;

                        if (!isItVeryFirstFlag)
                        {
                            unsigned long long entryLength = flag_from - 3 - previous_flag_to;
                            //std::cout << "entry starts at " << previous_flag_to + 1 << " and is of length " << entryLength << "\n";
                            
                            if (entryLength <= 0)
                            {
                                std::cout << "Oops. Something went wrong when trying to parse input. This could be caused by two options that follow each other without entry in between or maybe something like '----' or '-- --'. Please check your input.\n";
                                abort();
                            }
                            std::string entry = AllInputInLongString.substr(previous_flag_to + 1, entryLength);

                            reduceString(entry);
                            reduceString(flag);

                            assert(entry.size() > 0);
                            assert(flag.size() > 0);
                            
                            UserEntries.push_back(std::pair<std::string, std::string>(flag, entry));
                        } else
                        {
                            isItVeryFirstFlag = false;
                        }                  
                    }
                }
            } else
            {
                if (AllInputInLongString.at(char_index) == ' ')
                {
                    unsigned long long flag_to = char_index; // to is excluded as usual
                    unsigned long long lengthFlag = flag_to - flag_from; // to is excluded as usual
                    flag = AllInputInLongString.substr(flag_from, lengthFlag);
                    previous_flag_to = flag_to;

                    //std::cout << "flag starts at " << flag_from << " and is of length " << lengthFlag << "\n";

                    
                    flag = optionContainer.renameFlag(flag);

                    searchForBeginningOfFlag = true;
                }
            }
        } else
        {
            std::cout << "Internal error. Found the character " << AllInputInLongString.at(char_index) << " while this character was supposed to be removed from the string by SimBit beforehand\n";
            abort();
        }
    }
    
    if (!searchForBeginningOfFlag)
    {
        std::cout << "The last element of the input appears to be an option name (--Something) without any entry behind. The last flag was '" << AllInputInLongString.substr(flag_from -2) << "'\n";
        abort();
    }

    if (isItVeryFirstFlag)
    {
        std::cout << "No entries were found!\n";
        abort();
    }

    // add last entry
    {
        unsigned long long entryLength = AllInputInLongString.size() - previous_flag_to;
        if (entryLength <= 0)
        {
            std::cout << "Oops. Something went wrong when trying to parse input. There is maybe two options that follow each other without entry in between or maybe something like '----' or '-- --'. Please verify your input.\n";
            abort();
        }
        std::string entry = AllInputInLongString.substr(previous_flag_to, entryLength);

        reduceString(entry);
        reduceString(flag);

        assert(entry.size() > 0);
        assert(flag.size() > 0);
        
        UserEntries.push_back(std::pair<std::string, std::string>(flag, entry));
    }

    // Reorder UserEntries
    /*int UserEntryNextIndexToReorder = 0;
    for (int optionIndex = 0 ; optionIndex < optionContainer.getNbOptions() ; optionIndex++)
    {
        // if everything was reordered
        if (UserEntryNextIndexToReorder >= UserEntries.size()) 
        {
            break;
        }

        // look for flag
        std::string flagToFind = optionContainer.getOptionFirstName(optionIndex);
        for (int UserEntryIndex = UserEntryNextIndexToReorder ; UserEntryIndex < UserEntries.size() ; UserEntryIndex++)
        {
            if (UserEntries[UserEntryIndex].first == flagToFind)
            {
                assert(UserEntryIndex < UserEntries.size());
                std::swap(UserEntries[UserEntryIndex], UserEntries[UserEntryNextIndexToReorder]);
                UserEntryNextIndexToReorder++;
                break;
            }
        }
    }
    if (UserEntryNextIndexToReorder != UserEntries.size())
    {
        std::cout << "Internal error when reorder UserEntries. UserEntryNextIndexToReorder = " << UserEntryNextIndexToReorder << " UserEntries.size() = " << UserEntries.size() << ".\n";
        abort();
    }*/

    #ifdef DEBUG
    std::cout << "Present in UserEntries:\n";
    for (auto& p : UserEntries)
        std::cout << "\t" << p.first << ":\t" << p.second << "\n";
    #endif

    // random SSP security
    assert(SSP == nullptr);

    // Check the temporal changes
    {
        std::vector<int> temporalChanges;
        for (auto& p : UserEntries)
        {
            std::string& entry = p.second;
            size_t Gpos = entry.find("@G", 0);
            while (true)
            {
                Gpos = entry.find("@G", Gpos);
                if (Gpos == std::string::npos)
                    break;
                Gpos += 2;
                size_t start = Gpos;
                size_t end   = entry.find(" ", Gpos);
                size_t endTab = entry.find("\t", Gpos);
                if (endTab < end)
                {
                    end = endTab;
                }
                size_t length = end - start;

                std::string s = entry.substr(start, length);
                if (s.size() == 0)
                {
                    std::cout << "Expected zero (0) or a positive integer after generation specific marker (@G) but ti seems to have received a space (or maybe a quote followed by a space).\n";
                    abort();
                }
                //std::cout << "s = '" << s << "'\n";


                if (
                    std::count(s.begin(), s.end(), '.') > 0
                    || 
                        (s.at(0) != '-' && !isdigit(s.at(0)))
                    ||
                        !isdigit(s.at(s.size() - 1))
                )
                {
                    std::cout << "Expected zero (0) or a positive integer value after a generation specific marker (@G) but received '" << s << "'" <<std::endl;
                    abort();
                }
                int r;
                try
                {
                    double d = std::stod(s);
                    double fraction = d - ((long)d);
                    if (fraction > 0.00001 || fraction < -0.00001)
                    {
                        std::cout << "Expected zero (0) or a positive integer after a generation specific marker (@G) but received '" << s << "' which seems to be a float number to SimBit (but it is surprising it did not get caught the the first security gates)! (error caught at the second security gate in setParameters)" <<std::endl;
                    }
                    r = (int) std::stod(s);
                }
                catch(...)
                {
                    std::cout << "Expected zero (0) or a positive integer after a generation specific marker (@G) but received '" << s << "'' (error caught at the third security gate in setParameters) " << std::endl;
                    abort();
                }

                if (r < 0)
                {
                    std::cout << "Received a generation specific marker (@G) followed by a negative number (followed by the number "<< r <<"). Generations must be zero or positive integers.\n";
                    abort();
                }
                
                temporalChanges.push_back(r);
            }
        }
        GP->readTemporalChanges(temporalChanges); // Will order and remove duplicates
    }
        
    ////// Loop through each possible option
    for (Option& option : optionContainer.options)
    {
        std::string flag = option.optionNames[0]; // THe flags have already been renamed to this first flag

        // Find index in UserEntries
        std::vector<int> UserInputIndexForFlags;
        for (int optionReceivedIndex = 0 ; optionReceivedIndex < UserEntries.size() ; optionReceivedIndex++)
        {
            if (flag == UserEntries[optionReceivedIndex].first)
            {
                UserInputIndexForFlags.push_back(optionReceivedIndex);
            }
        }
        
        option.received(optionContainer);
        //std::cout << "Flag " << flag << "\n";

        // if flag not found in UserEntries
        if (UserInputIndexForFlags.size() == 0) 
        {
            //std::cout << "Set to default...\n";
            this->setOptionToDefault(flag); // The line 'option.received(optionContainer);' is present in the function 'setOptionToDefault' only when something is indeed received
        // if flag was found
        } else 
        {
            for (int UserInputIndexForFlag : UserInputIndexForFlags)
            {
                //std::cout << "Set from user...\n";

                // Check for potential issue
                std::string& entry = UserEntries[UserInputIndexForFlag].second;
                bool issueWithEntry = false;
                if (entry.size() == 0)
                {
                    issueWithEntry = true;
                }
                if (entry.size() == 1)
                {
                    if (entry.at(0) == ' ')
                    {
                        issueWithEntry = true;
                    }
                }
                if (issueWithEntry)
                {
                    std::cout << "For flag / option '" << flag << "', the entry is either empty or contain only whitespaces.\n";
                    abort();
                }

                // Make InputReader
                InputReader input(entry, std::string("In '--") + flag + std::string("', "));
                input.interpretKeywords();

                // Set option
                this->setOptionToUserInput(flag, input);
            }
        }
    }

    // random SSP security
    assert(SSP == nullptr);

    // Security to make sure it initiated all options
    /*if (optionContainer.HowManyOptionsWereInitiated() != optionContainer.getNbOptions())
    {
        std::cout << "Internal error! Some options failed to initiate!\n\n";
        optionContainer.listOptions();
        abort();
    }*/

    #ifdef DEBUG
    std::cout << "Start GP->initializeAllSpeciesPatchSizes()\n" << std::endl;
    #endif   
    GP->initializeAllSpeciesPatchSizes();

    #ifdef DEBUG
    std::cout << "Start outputWriter.SetAllTimes()\n" << std::endl;
    #endif    
    // set All Times
    outputWriter.SetAllTimes();

    #ifdef DEBUG
    std::cout << "Start outputWriter.ShouldSimulationBeDone()\n" << std::endl;
    #endif

    // Test if simulation should be run or if overwritemode and existing files prevent it
    if (!outputWriter.ShouldSimulationBeDone(GP->OverwriteMode))
    {
        std::cout << "\t\tLogfile path:\n\t\t\t" << outputWriter.get_OutputFiles(Logfile)[0].getPath() << "\n\n";
        std::cout << "Looking at existing files, files asked for outputs and the OverwriteMode, the simulation will not run (as figured out by 'outputWriter::ShouldSimulationBeDone'). If you want to overwrite files, please use '--Overwrite 2'\n";

        exit(1);
    }
        

    #ifdef DEBUG
    std::cout << "Start outputWriter.ClearAllFileContent()\n" << std::endl;
    #endif 

    // Clear All files content
    assert(SSP == nullptr);
    for (int speciesIndex = 0 ; speciesIndex < GP->nbSpecies; speciesIndex++)
    {
        SSP = this->getSSPaddress(speciesIndex);
        outputWriter.ClearAllFileContent();
    }
    SSP = nullptr;

    #ifdef DEBUG
    std::cout << "Start outputWriter.AreThereAnyOutput()\n" << std::endl;
    #endif 

    // Check that some output is asked for
    if (outputWriter.AreThereAnyOutput())
    {
    #ifdef DEBUG
    std::cout << "Start outputWriter.IsLastGenerationSampled()\n" << std::endl;
    #endif         
        // Check if last output is sampled
        if (!outputWriter.IsLastGenerationSampled())
        {
            std::cout << "It appears that the last generation is not sampled as this last generation is not indicated in any of the output time. This means that part of the simulation will occur but no output will be given for this part. Is this really on purpose? If yes, please modify this security check in the code.\n\nAs a reminder, the outputs are produced at the end of the generation. Outputs asked for generation '0' are just the population as initialized. Outputs asked for generation '1' is for the first generation. Outputs asked for generation 'nbGenerations + 1' make no sense and should raise an error.\n\nIf you want a dry run you can use '-DryRun'. If you set all output times to '-1' (which is the default), then you will get a warning message but not an error message as here.\n";
            abort();
        }
    } else
    {
        std::cout << "\t\tWARNING: No output has been requested!\n\n";        
    }

    for (auto& SSPi : SSPs)
    {
        #ifdef DEBUG
        std::cout << "Start SSPi.setRandomDistributions()\n" << std::endl;
        #endif      
        SSPi.setRandomDistributions();

        #ifdef DEBUG
        std::cout << "Start SSPi.setFromLocusToFitnessMapIndex()\n" << std::endl;
        #endif 
        SSPi.setFromLocusToFitnessMapIndex();

        // Check initial population size and patch capacity for fecundity = -1
        // Check also that DispWeightByFitness is correct.
        if (SSPi.fecundityForFitnessOfOne == -1.0)
        {
            assert(GP->__PatchNumber[0] > 0);
            assert(GP->__PatchNumber[0] == SSPi.__patchCapacity[0].size());
            assert(GP->__PatchNumber[0] == SSPi.patchSize.size());
            for (int patch_index = 0 ; patch_index < GP->__PatchNumber[0] ; patch_index++)
            {
                if (SSPi.__patchCapacity[0][patch_index] != SSPi.patchSize[patch_index])
                {
                    std::cout << "For patch " << patch_index << " (zero based counting), the initial patchSize does not match the patchCapacity at generation 0, while SSPi.fecundityForFitnessOfOne is set to -1.0. When SSPi.fecundityForFitnessOfOne is set to -1.0, the patchSize must always be at carrying capacity.\n";
                    abort();
                }
            }
        } else
        {
            SSPi.DispWeightByFitness = true;
        }

        // Initialize tree for T4
        if (SSPi.T4_nbBits)
        {
            SSP = &SSPi;
            SSPi.T4Tree.initialize();
            SSP=nullptr;
        }

    }
        
    #ifdef DEBUG
    std::cout << "Start outputWriter.PrintLogfile(AllInputInLongString)\n" << std::endl;
    #endif          
    outputWriter.PrintLogfile(AllInputInLongString);

    // random SSP security
    assert(SSP == nullptr);
}



void AllParameters::UpdatePopsAndParameters()
{
    Update(true);
}

void AllParameters::UpdateParameters()
{
    Update(false);
}

void AllParameters::Update(bool updatePopsToo)
{
    assert(this->GlobalP.nbSpecies > 0);

    /*
    for (auto& g : this->GlobalP.__GenerationChange)
    {
        std::cout << "g = " << g << " ";
    }
    std::cout << "\n";
    */

    auto it = std::find(this->GlobalP.__GenerationChange.begin(), this->GlobalP.__GenerationChange.end(),  this->GlobalP.CurrentGeneration);
    if (it != this->GlobalP.__GenerationChange.end()) // Test if update of parameters must be done at this generation.
    {
        #ifdef DEBUG
        std::cout << "Enters in update part of AllParameters::UpdatePopsAndParameters\n";
        #endif
        // Get 'generation_index'
        int generation_index = it - GP->__GenerationChange.begin();
        assert(generation_index >= 0);
        assert(generation_index < GP->__GenerationChange.size());

        /////////////////////////////////
        //// Change GlobalParameters ////
        /////////////////////////////////

        // assert GP pointer
        assert(GP == this->getGPaddress());

        // Update GP patchNumber
        GP->UpdateParametersPatchNumber(generation_index);


        for (int speciesIndex = 0 ; speciesIndex < this->GlobalP.nbSpecies ; speciesIndex++)
        {
                

            //////////////////////////////////////////
            //// Change SpeciesSpecificParameters ////
            //////////////////////////////////////////

            // Get SSP pointer
            SSP = this->getSSPaddress(speciesIndex);

            // Update SSP
            auto previousPatchSizes = SSP->UpdateParameters(generation_index); // will also call 'GP->saveSSPPatchSize_toGP' and will also update Original of dispersalData

            /////////////////////
            //// Change Pops ////
            /////////////////////

            if (updatePopsToo)
            {
                // get Populations
                Pop& pop1   = allSpecies[speciesIndex].second;
                Pop& pop2   = allSpecies[speciesIndex].first;
                assert(pop1.getNbPatches() == pop2.getNbPatches());

                // Update Pops
                Pop::updatePops(pop1, pop2, speciesIndex, previousPatchSizes); // Will also recompute fitness    
            }

            

            ///////////////////////////
            //// Change SimTracker ////
            ///////////////////////////

            if (GP->CurrentGeneration != 0 && SSP->T1_isSelection)
                SSP->simTracker.forceToRecomputeNextTime = true;

        }
        SSP = nullptr;


#ifdef DEBUG
        std::cout << "Exit after update in AllParameters::UpdatePopsAndParameters\n";
#endif    
    }
}



void AllParameters::setOptionToDefault(std::string& flag)
{
    if (flag == "seed")
    {
        InputReader input("default", "In Default value for --seed,");
        GP->readSeed(input);
    }
    else if (flag == "S")
    {
        std::string sN = "sp";
        SpeciesSpecificParameters ssp(sN, 0);
        SSPs.push_back(std::move(ssp));
        GP->speciesNames.push_back(sN);
        GP->nbSpecies++;
        assert(GP->nbSpecies == 1);
    }
    else if (flag == "nbSubGens")
    {
        for (int speciesIndex = 0 ; speciesIndex < GP->nbSpecies ;speciesIndex++)
        {
            SSPs[speciesIndex].nbSubGenerationsPerGeneration = 1;
        }
    }
    /*else if (flag == "T")
    {
        assert(this->GlobalP.__GenerationChange.size() == 0);
        this->GlobalP.__GenerationChange.push_back(0);
    }*/
    else if (flag == "matingSystem")
    {
        assert(SSPs.size() > 0);
        for (auto& SSPi : this->SSPs)
        {
            InputReader input(std::string("H"), "In Default value for --matingSystem,");
            SSPi.readMatingSystem(input);
        }
    }
    else if (flag == "H")
    {
        assert(SSPs.size() > 0);
        for (auto& SSPi : this->SSPs)
        {
            InputReader input(std::string("@G0 unif 0"), "In Default value for --Habitats,");
            SSPi.readHabitats(input);
        }
    }
    else if (flag == "GP")
    {
        OutputFile::GeneralPath = std::string("");
    }
    else if (flag == "T1_vcf_file")
    {
        // Nothing to do
    }
    else if (flag == "T5_vcf_file")
    {
        // Nothing to do
    }
    else if (flag == "T1_LargeOutput_file")
    {
        // Nothing to do
    }
    else if (flag == "T5_LargeOutput_file")
    {
        // Nothing to do
    }
    else if (flag == "T1_AlleleFreq_file")
    {
        // Nothing to do
    }
    else if (flag == "T5_AlleleFreq_file")
    {
        // Nothing to do
    }
    else if (flag == "Log")
    {
        if (outputWriter.LogfileType != 0)
        {
            OutputFile file("Logfile", Logfile);
            outputWriter.insertOutputFile(std::move(file));
            if (!outputWriter.isFile(Logfile))
            {
                std::cout << "Internal error: Default Logfile failed to be inserted (and failed to be found) in outputWriter\n";
                abort();
            }
        }
    }
    else if (flag == "T1_MeanLD_file")
    {
        // Nothing to do
    }
    else if (flag == "T1_LongestRun_file")
    {
        // Nothing to do
    }
    else if (flag == "T1_HybridIndex_file")
    {
        // Nothing to do
    }
    else if (flag == "T1_ExpectiMinRec_file")
    {
        // Nothing to do
    }
    else if (flag == "T2_LargeOutput_file")
    {
        // Nothing to do
    }
    else if (flag == "SaveBinary_file")
    {
        // Nothing to do
    }
    else if (flag == "T3_LargeOutput_file")
    {
        // Nothing to do
    }
    else if (flag == "T3_MeanVar_file")
    {
        // Nothing to do
    }
    else if (flag == "fitnessSubsetLoci_file")
    {
        // Nothing to do
    }
    else if (flag == "fitness_file")
    {
        // Nothing to do
    }
    else if (flag == "fitnessStats_file")
    {
        // Nothing to do
    }
    else if (flag == "T1_FST_file")
    {
        // Nothing to do
    } else if (flag == "T1_FST_info")
    {
        InputReader input("default", std::string("In default value for '--") + flag + std::string("', "));
        GP->readT1_FST_info(input);
    }
    else if (flag == "extraGeneticInfo_file")
    {
        // Nothing to do
    }
    else if (flag == "patchSize_file")
    {
        // Nothing to do
    }
    else if (flag == "patchSize_file")
    {
        // Nothing to do
    }
    else if (flag == "extinction_file")
    {
        // Nothing to do
    }
    else if (flag == "genealogy_file")
    {
        for (auto& SSPi : this->SSPs)
        {
            SSPi.simTracker.genealogy.setGenealogyToNothing();
        }   
    }
    else if (flag == "coalesce")
    {
        for (auto& SSPi : this->SSPs)
        {
            SSPi.simTracker.genealogy.setcoalesceGenealogyFrequency(0);
        }   
    }
    else if (flag == "T4_LargeOutput_file")
    {
        // Nothing to do
    }
    else if (flag == "T4_vcf_file")
    {
        // Nothing to do
    }
    else if (flag == "T4_SFS_file")
    {
        // Nothing to do
    }
    else if (flag == "T1_SFS_file")
    {
        // Nothing to do
    }
    else if (flag == "T5_SFS_file")
    {
        // Nothing to do
    }
    else if (flag == "T4_printTree")
    {
        // Nothing to do
    } 
    else if (flag == "outputSFSbinSize")
    {
        InputReader input("default","In Default value for --outputSFSbinSize,");
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readOutputSFSbinSize);
    } 
    else if (flag == "LogfileType")
    {
        outputWriter.LogfileType = 1;
    }
    else if (flag == "sequencingErrorRate")
    {
        GP->sequencingErrorRate = 0.0;
    }
    else if (flag == "nbGens")
    {
        std::cout << "--nbGens (--nbGenerations) is missing" << std::endl;
        abort();
    } else if (flag == "startAtGeneration")
    {
        GP->startAtGeneration = 0;
    }
    else if (flag == "PN")
    {
        std::cout << "--PN (--PatchNumber) is missing" << std::endl;
        abort();
    }
    else if (flag == "nbThreads")
    {
        this->GlobalP.nbThreads = 1;   
    }
    else if (flag == "N")
    {
        std::cout << "--N (--patchCapacity) is missing" << std::endl;
        abort();
    }
    else if (flag == "InitialpatchSize")
    {
        // assumes starts at carrying capacity
        InputReader input(std::string("@S0 unif -1"), "In Default value for --InitialpatchSize,");
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readInitialpatchSize);

    }
    else if (flag == "cloningRate")
    {
        for (auto& SSPi : this->SSPs)
            SSPi.cloningRate = 0.0; // no cloning
    }
    else if (flag == "selfingRate")
    {
        for (auto& SSPi : this->SSPs)
            SSPi.selfingRate = -1; // -1 means it selves at rate 1/2N
    }
    else if (flag == "m")
    {
        bool IsThereOnlyOnePatch = true;
        assert(this->GlobalP.__PatchNumber.size() == this->GlobalP.__GenerationChange.size());
        for (int& PN : this->GlobalP.__PatchNumber)
        {
            if (PN != 1) {IsThereOnlyOnePatch = false;}
        }
        if (!IsThereOnlyOnePatch)
        {
            std::cout << "--m (--DispMat) is missing. This option is required as PatchNumber > 0 for at least one moment in the simulation" << std::endl;
            abort();
        } else
        {
            InputReader input(std::string("@S0 @G0 OnePatch"), "In Default value for --m (--DispMat),");
            wrapperOverSpecies(input, &SpeciesSpecificParameters::readDispMat);
        }
    }
    else if (flag == "gameteDispersal")
    {
        InputReader input(std::string("@S0 no"), "In Default value for --gameteDispersal,");
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readGameteDispersal);
    }
    else if (flag == "DispWeightByFitness")
    {
        for (auto& SSPi : this->SSPs)
            SSPi.DispWeightByFitness = false;
    }
    else if (flag == "L")
    {
        std::cout << "--L (--Loci) is missing" << std::endl;
        abort();
    }
    else if (flag == "ploidy")
    {
        for (auto& SSPi : this->SSPs)
            SSPi.ploidy = 2;
    }
    else if (flag == "r")
    {
        for (auto& SSPi : this->SSPs)
        {
            if (SSPi.TotalNbLoci != 1)
            {
                std::cout << "--r (--RecombinationRate) is missing" << std::endl;
                abort();
            }
        }
    }
    else if (flag == "recRateOnMismatch")
    {
        for (auto& SSPi : this->SSPs)
        {
            SSPi.recRateOnMismatch_bool = false;
            SSPi.recRateOnMismatch_halfWindow = -1;
            SSPi.recRateOnMismatch_factor = -1.0;
        }
    }
    else if (flag == "FitnessMapInfo")
    {
        for (auto& SSPi : this->SSPs)
        {
            SSPi.FitnessMapMinimNbLoci = 120;
            SSPi.FitnessMapProbOfEvent = 0.008; // This is a default value but the optimal value will vary a lot depending upon what exactly is being simulated! I can definitely improve the decision of the FitnessMap Boundaries
    
            SSPi.FitnessMapT5WeightProbOfEvent = 0.1;
            SSPi.FitnessMapCoefficient = -9.0; 
        }
    }
    else if (flag == "resetTrackedT1Muts")
    {
        InputReader input("default","In default value for --resetTrackedT1Muts,");
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readResetTrackedT1Muts);  
    }
    else if (flag == "fec")
    {
        for (auto& SSPi : this->SSPs)
            SSPi.fecundityForFitnessOfOne = -1.0;
    }
    else if (flag == "T1_mu")
    {
        for (auto& SSPi : this->SSPs)
        {
            if (SSPi.T1_nbChars)
            {
                std::cout << "You asked for T1 loci for species "<< SSPi.speciesName <<" but option '--T1_MutationRate' is missing!\n";
                abort();
            } else
            {
                SSPi.T1_Total_Mutation_rate = 0.0;
            }
        }
    }
    else if (flag == "T5_mu")
    {
        for (auto& SSPi : this->SSPs)
        {
            if (SSPi.T5_nbBits)
            {
                std::cout << "You asked for T5 loci for species "<< SSPi.speciesName <<" but option '--T5_MutationRate' is missing!\n";
                abort();
            } else
            {
                SSPi.T5_Total_Mutation_rate = 0.0;
            }
        }
    }
    else if (flag == "additiveEffectAmongLoci")
    {
        InputReader input(std::string("@S0 no"), "In Default value for --additiveEffectAmongLoci,");
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readadditiveEffectAmongLoci);   
    } 
    else if (flag == "selectionOn")
    {
        InputReader input(std::string("@S0 fertility"), "In Default value for --selectionOn,");
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readSelectionOn);
    }
    else if (flag == "T1_fit")
    {
        InputReader input(std::string("@S0 @H0 MultiplicityUnif 1.0"), "In Default value for --T1_FitnessEffects,");
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readT1_FitnessEffects);
    }
    else if (flag == "T5_fit")
    {
        InputReader input(std::string("@S0 @H0 MultiplicityUnif 1.0"), "In Default value for --T5_FitnessEffects,");
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readT5_FitnessEffects);
    }
    else if (flag == "T1_ini")
    {
        InputReader input(std::string("@S0 AllZeros"), "In Default value for --T1_ini,");
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readT1_Initial_AlleleFreqs);
    }
    else if (flag == "T5_ini")
    {
        InputReader input(std::string("@S0 AllZeros"), "In Default value for --T5_ini,");
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readT5_Initial_AlleleFreqs);
    }
    else if (flag == "T1_epistasis")
    {
        for (auto& SSPi : this->SSPs)
            SSPi.T1_isEpistasis = false; // Is not even necessary to do
    }
    else if (flag == "T2_mu")
    {
        for (auto& SSPi : this->SSPs)
        {
            if (SSPi.T2_nbChars)
            {
                std::cout << "You asked for T2 loci for species "<< SSPi.speciesName <<" but option '--T2_MutationRate' is missing!\n";
                abort();
            } else
            {
                SSPi.T2_Total_Mutation_rate = 0.0;
            }
        }
    }
    else if (flag == "T2_fit")
    {
        InputReader input(std::string("@S0 @H0 unif 1"), "In Default value for --T2_fit,");
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readT2_FitnessEffects);
    }
    else if (flag == "T3_mu")
    {
        for (auto& SSPi : this->SSPs)
        {
            if (SSPi.T3_nbChars)
            {
                std::cout << "You asked for T3 loci for species "<<SSPi.speciesName<<"but option '--T3_MutationRate' is missing!\n";
                abort();
            } else
            {
                SSPi.T3_Total_Mutation_rate = 0.0;
            }
        }
    }
    else if (flag == "T3_pheno")
    {   
        for (auto& SSPi : this->SSPs)
        {
            if (SSPi.T3_nbChars)
            {
                std::cout << "You asked for T3 loci for species "<<SSPi.speciesName<<" but option '--T3_PhenotypicEffects' is missing!\n";
                abort();
            }
        }
    }
    else if (flag == "T3_fit")
    {
        InputReader input(std::string("NA"), "In Default value for --T3_fit,");
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readT3_FitnessLandscape);
    }
    else if (flag == "T3_DN")
    {
        InputReader input(std::string("unif 0"), "In Default value for --T3_DN,");
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readT3_DevelopmentalNoise);
    }
    else if (flag == "T4_mu")
    {
        for (auto& SSPi : this->SSPs)
        {
            if (SSPi.T4_nbBits)
            {
                std::cout << "You asked for T4 loci for species "<<SSPi.speciesName<<"but option '--T4_MutationRate' is missing!\n";
                abort();
            }
        }
    } else if (flag == "T4_maxAverageNbNodesPerHaplotype")
    {
        InputReader input(std::string("default"), "In Default value for --T4_maxAverageNbNodesPerHaplotype,");
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readT4_maxAverageNbNodesPerHaplotype);
    }
    else if (flag == "eco")
    { 
        assert(this->GlobalP.nbSpecies > 0);
        assert(this->GlobalP.nbSpecies == this->SSPs.size());
        
        //std::cout << "this->GlobalP.nbSpecies = " << this->GlobalP.nbSpecies << "\n";
        for (int speciesIndex_to = 0 ; speciesIndex_to < this->GlobalP.nbSpecies; speciesIndex_to++)
        {
            std::vector<double> fromOneSpecies_effect;
            std::vector<char> fromOneSpecies_type;
            fromOneSpecies_effect.resize(this->GlobalP.nbSpecies);
            fromOneSpecies_type.resize(this->GlobalP.nbSpecies);
            for (int speciesIndex_from = 0 ; speciesIndex_from < this->GlobalP.nbSpecies; speciesIndex_from++)
            {
                fromOneSpecies_type[speciesIndex_from]   = '0';
                fromOneSpecies_effect[speciesIndex_from] = 0.0;
            }
            GP->speciesEcoRel_type.push_back(fromOneSpecies_type);
            GP->speciesEcoRel_effect.push_back(fromOneSpecies_effect);
        }
        
        /*
        std::cout << "Default eco Types:\n";
        for (auto& a : GP->speciesEcoRel_type)
        {
            std::cout << "\t";
            for (auto& b : a)   
            {
                std::cout << "'" << b << "' ";
            }
            std::cout << "\n";
        }
        std::cout << "\n";
        */
    } else if (flag == "growthK")
    {
        assert(GP->speciesEcoRel_type.size() == this->GlobalP.nbSpecies);
        assert(GP->speciesEcoRel_effect.size() == this->GlobalP.nbSpecies);

        InputReader input("unif def", std::string("In '--") + flag + std::string("', "));
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readGrowthK);
    } else if (flag == "resetGenetics")
    {
        // nothing to do
    }
    else if (flag == "Overwrite")
    {
        this->GlobalP.OverwriteMode = 1;
    }
    else if (flag == "readPopFromBinary")
    {
        InputReader input(std::string("false"), "In Default value for --readPopFromBinary,");
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readReadPopFromBinary);
    } 
    else if (flag == "DryRun")
    {
        this->GlobalP.DryRun = false;
    }
    else if (flag == "centralT1LocusForExtraGeneticInfo")
    {
        for (auto& SSPi : this->SSPs)
            SSPi.centralT1LocusForExtraGeneticInfo = -1;
    } else
    {
        std::cout << "Internal error in AllParameters::setOptionToDefault. flag " << flag << " could not be found.\n";
        abort();
    }
}



void AllParameters::setOptionToUserInput(std::string& flag, InputReader input)
{

    if (flag == "GeneralPath" || flag == "GP")
    {

        OutputFile::GeneralPath = input.GetNextElementString();
    } else if (flag == "T1_vcf_file" || flag == "T1_VCF_file")
    {
        OutputFile file(input.GetNextElementString(), T1_vcfFile);
        file.interpretTimeAndSubsetInput(input);
        outputWriter.insertOutputFile(std::move(file));

    } else if (flag == "T5_vcf_file" || flag == "T5_VCF_file")
    {
        OutputFile file(input.GetNextElementString(), T5_vcfFile);
        file.interpretTimeAndSubsetInput(input);
        outputWriter.insertOutputFile(std::move(file));

    } else if (flag == "T1_LargeOutput_file")
    {
        OutputFile file(input.GetNextElementString(), T1_LargeOutputFile);
        file.interpretTimeAndSubsetInput(input);
        outputWriter.insertOutputFile(std::move(file));

    } else if (flag == "T5_LargeOutput_file")
    {
        OutputFile file(input.GetNextElementString(), T5_LargeOutputFile);
        file.interpretTimeAndSubsetInput(input);
        outputWriter.insertOutputFile(std::move(file));

    } else if (flag == "T1_AlleleFreq_file")
    {
        OutputFile file(input.GetNextElementString(), T1_AlleleFreqFile);
        file.interpretTimeAndSubsetInput(input);
        outputWriter.insertOutputFile(std::move(file));
    }  else if (flag == "T5_AlleleFreq_file")
    {
        OutputFile file(input.GetNextElementString(), T5_AlleleFreqFile);
        file.interpretTimeAndSubsetInput(input);
        outputWriter.insertOutputFile(std::move(file));
    } else if (flag == "Log" || flag == "Logfile" || flag == "Logfile_file")
    {
        OutputFile file(input.GetNextElementString(), Logfile);
        outputWriter.insertOutputFile(std::move(file));
    }  else if (flag == "T1_MeanLD_file")
    {
        
        OutputFile file(input.GetNextElementString(), MeanLDFile);
        file.interpretTimeAndSubsetInput(input);
        outputWriter.insertOutputFile(std::move(file));
    }  else if (flag == "T1_LongestRun_file")
    {
        
        OutputFile file(input.GetNextElementString(), LongestRunFile);
        file.interpretTimeAndSubsetInput(input);
        outputWriter.insertOutputFile(std::move(file));

    }  else if (flag == "T1_HybridIndex_file" )
    {
        OutputFile file(input.GetNextElementString(), HybridIndexFile);
        file.interpretTimeAndSubsetInput(input);
        outputWriter.insertOutputFile(std::move(file));
    
    }  else if (flag == "T1_ExpectiMinRec_file" )
    {
        OutputFile file(input.GetNextElementString(), ExpectiMinRecFile);
        file.interpretTimeAndSubsetInput(input);
        outputWriter.insertOutputFile(std::move(file));        

    }  else if (flag == "T2_LargeOutput_file" )
    {     
        OutputFile file(input.GetNextElementString(), T2_LargeOutputFile);
        file.interpretTimeAndSubsetInput(input);
        outputWriter.insertOutputFile(std::move(file));

    }  else if (flag == "SaveBinary_file" )
    {
        OutputFile file(input.GetNextElementString(), SaveBinaryFile);
        file.interpretTimeAndSubsetInput(input);
        outputWriter.insertOutputFile(std::move(file));

    }  else if (flag == "T3_LargeOutput_file" )
    {
        OutputFile file(input.GetNextElementString(), T3_LargeOutputFile);
        file.interpretTimeAndSubsetInput(input);
        outputWriter.insertOutputFile(std::move(file));
        
    }  else if (flag == "T3_MeanVar_file" )
    {
        OutputFile file(input.GetNextElementString(), T3_MeanVarFile);
        file.interpretTimeAndSubsetInput(input);
        outputWriter.insertOutputFile(std::move(file));
         
    }  else if (flag == "fitness_file" )
    {
        OutputFile file(input.GetNextElementString(), fitness);
        file.interpretTimeAndSubsetInput(input);
        outputWriter.insertOutputFile(std::move(file));

    }  else if (flag == "fitnessSubsetLoci_file" )
    {

        OutputFile file(input.GetNextElementString(), fitnessSubsetLoci);

        // Find out the separation between time arguments and the loci to subset
        int VIndexFirstWordOfSSP;
        {
            // false means do not throw error message if word cannot be found. 
            // true means throw error message if word cannot be found. 
            int FirstLociSet = input.FindVIndexOfNextMatchingString("LociSet",true); 
            int FirstSpeciesSpecificIndication = input.FindVIndexOfNextMatchingString("@S",false);

            if (FirstSpeciesSpecificIndication != input.numberOfWordsInInput())
            {
                if (FirstSpeciesSpecificIndication != FirstLociSet - 1)
                {
                    std::cout << "For option --fitnessSubsetLoci_file, received the first species specific marker at word index "<<FirstSpeciesSpecificIndication<<" and the first LociSet at word index " << FirstLociSet << ". Note that if a species specific marker is used, then the keyword 'LociSet' should directly follow it, which is not the case here. Hence SimBit fails to understand this input. Sorry!\n";
                    abort();
                }
                VIndexFirstWordOfSSP = FirstSpeciesSpecificIndication;
            } else
            {
                VIndexFirstWordOfSSP = FirstLociSet;
            }
        } 

        // Subset the inputs
        // -1 indicates this is not specific to a single species (yet)
        //std::cout << "VIndexFirstWordOfSSP = " << VIndexFirstWordOfSSP << "\n";
        InputReader inputForTime(input, input.currentVIndex(), VIndexFirstWordOfSSP, -1);
        InputReader inputForSpeciesSpecificLociSets(input, VIndexFirstWordOfSSP, input.numberOfWordsInInput(), -1);
        //std::cout << "inputForTime = " << inputForTime.print() << "\n";
        //std::cout << "inputForSpeciesSpecificLociSets = " << inputForSpeciesSpecificLociSets.print() << "\n";

        // time
        file.interpretTimeAndSubsetInput(inputForTime);

        // LociSet
        wrapperOverSpecies(inputForSpeciesSpecificLociSets, &SpeciesSpecificParameters::readSubsetLociForfitnessSubsetLoci_file);


        outputWriter.insertOutputFile(std::move(file));

        input.markedAsRead();

    } else if (flag == "fitnessStats_file" )
    {
        OutputFile file(input.GetNextElementString(), fitnessStats);
        file.interpretTimeAndSubsetInput(input);
        outputWriter.insertOutputFile(std::move(file));

    } else if (flag == "T1_FST_file" )
    {
        OutputFile file(input.GetNextElementString(), T1_FST);
        file.interpretTimeAndSubsetInput(input);
        outputWriter.insertOutputFile(std::move(file));

    }  else if (flag == "T1_FST_info")
    {
                GP->readT1_FST_info(input);
    } else if (flag == "extraGeneticInfo_file" )
    {        
        OutputFile file(input.GetNextElementString(), extraGeneticInfo);
        outputWriter.insertOutputFile(std::move(file));
        
    }  else if (flag == "patchSize_file" )
    {
        OutputFile file(input.GetNextElementString(), patchSize);
        file.interpretTimeAndSubsetInput(input);
        outputWriter.insertOutputFile(std::move(file));
     

    }  else if (flag == "extinction_file" )
    {
        OutputFile file(input.GetNextElementString(), extinction);
        outputWriter.insertOutputFile(std::move(file));
        
    } else if (flag == "genealogy_file" )
    {
        OutputFile file(input.GetNextElementString(), genealogy);
        file.interpretTimeAndSubsetInput(input);
                        for (auto& SSPi : this->SSPs)
        {
            SSPi.simTracker.genealogy.setGenealogyTimes(file.getTimes());
        }   
        outputWriter.insertOutputFile(std::move(file));
    }  else if (flag == "coalesce" )
    {
                int coalesceGenealogyFrequency = input.GetNextElementInt();
                if (coalesceGenealogyFrequency < 0)
        {
            std::cout << "In option --coalesced, the values received is negative ( is "<<coalesceGenealogyFrequency<<"). Value of zero mean no coalescence. Any value bigger than 0 indicates the frequency at which SimBit will coalesce. This frequency won't change anything to the output but might well affect the computational time (should not affect the RAM though).\n";
            abort();
        }
        for (auto& SSPi : this->SSPs)
        {
            SSPi.simTracker.genealogy.setcoalesceGenealogyFrequency(coalesceGenealogyFrequency);
        }   
    } 
    else if (flag == "T4_LargeOutput_file")
    {
        OutputFile file(input.GetNextElementString(), T4_LargeOutputFile);
        file.interpretTimeAndSubsetInput(input);
        outputWriter.insertOutputFile(std::move(file));

    }
    else if (flag == "T4_vcf_file" || flag == "T4_VCF_file")
    {
        OutputFile file(input.GetNextElementString(), T4_vcfFile);
        file.interpretTimeAndSubsetInput(input);
        outputWriter.insertOutputFile(std::move(file));

    }
    else if (flag == "T4_SFS_file")
    {
        OutputFile file(input.GetNextElementString(), T4_SFS_file);
        file.interpretTimeAndSubsetInput(input);
        outputWriter.insertOutputFile(std::move(file));

    } 
    else if (flag == "T1_SFS_file")
    {
        OutputFile file(input.GetNextElementString(), T1_SFS_file);
        file.interpretTimeAndSubsetInput(input);
        outputWriter.insertOutputFile(std::move(file));

    }
    else if (flag == "T5_SFS_file")
    {
        OutputFile file(input.GetNextElementString(), T5_SFS_file);
        file.interpretTimeAndSubsetInput(input);
        outputWriter.insertOutputFile(std::move(file));

    }
    else if (flag == "T4_printTree")
    {
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readT4_printTree);
    } 
    else if (flag == "outputSFSbinSize")
    {
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readOutputSFSbinSize);
    }
    else if (flag == "LogfileType")
    {
        
        outputWriter.LogfileType = input.GetNextElementInt();
        if (outputWriter.LogfileType > 2)
        {
            if (outputWriter.LogfileType != 102105116) // 102105116 (stands for 'fit' in ascii) outputs the fitness arrays only
            {
                std::cout << "LogfileType received is " << outputWriter.LogfileType << ". Only 0, 1 and 2 are accepted for the moment.\n";
                abort();
            }
        }
    } else if (flag == "sequencingErrorRate")
    {


        GP->sequencingErrorRate = input.GetNextElementDouble();
        if (GP->sequencingErrorRate < 0.0)
        {
            std::cout << "In --" << flag << ", received a negative error rate (received "<<GP->sequencingErrorRate<<").";
            abort();
        }
    } else if (flag == "nbGens" || flag == "nbGenerations")
    {
        
        GP->nbGenerations = input.GetNextElementInt();
        if (GP->nbGenerations < 0)
        {
            std::cout << "'nbGenerations' is " << GP->nbGenerations << ". 'nbGenerations' cannot be lower than 0." << std::endl;
            abort();
        }
        for (auto& t : GP->__GenerationChange)
        {
            if (t > GP->nbGenerations)
            {
                std::cout << "User asked to simulate " << GP->nbGenerations << " generations (info set by using the option --nbGenerations (--nbGens)). However, somewhere in the command, there are temporal changes at generation " << t <<" (as indicated with generation specific marker @G). It does not make sense to ask for a temporal change at a time after the simulation is over.\n";
                abort();
            }
        }
    } else if (flag == "startAtGeneration")
    {
                GP->startAtGeneration = input.GetNextElementInt();
        assert(GP->nbGenerations >= 0);
        if (GP->startAtGeneration < 0 || GP->startAtGeneration > GP->nbGenerations)
        {
            std::cout << "For option " << flag << ", generation to be received is either negative of larger than the number of generations to be computed. GP->startAtGeneration = " << GP->startAtGeneration << "  GP->nbGenerations = " << GP->nbGenerations  << "\n";
            abort();
        }
    } else if (flag == "PN" || flag == "PatchNumber")
    {
        GP->readPatchNumber(input);
    }  else if (flag == "seed" || flag == "random_seed")
    {

        GP->readSeed(input);
    }  else if (flag == "nbThreads")
    {
        
        
        GP->nbThreads = input.GetNextElementInt();
        //omp_set_num_threads(GP->nbThreads);
    }  else if (flag == "N" || flag == "patchCapacity")
    {
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readpatchCapacity);

    } else if (flag == "nbSubGens" || flag == "nbSubGenerations")
    {
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readnbSubGenerations);

    } else if (flag == "InitialpatchSize")
    {
        // just a security
        for (auto& SSPi : SSPs)
        {
            assert(SSPi.__patchCapacity[0].size() == GP->__PatchNumber[0]);
        }
        // set option
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readInitialpatchSize);

    }  else if (flag == "cloningRate")
    {
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readCloningRate);
        
    }  else if (flag == "selfingRate")
    { 
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readSelfingRate);
        
    }  else if (flag == "DispMat" || flag == "m")
    { 
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readDispMat);
        
    }  else if (flag == "DispWeightByFitness")
    {
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readDispWeightByFitness);
        
    }  else if (flag == "Loci" || flag == "L")
    {
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readLoci);

    }  else if (flag == "ploidy")
    {     
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readPloidy);
    
    }  else if (flag == "RecombinationRate" || flag == "r")
    {
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readRecombinationRate);
        
    }  else if (flag == "recRateOnMismatch")
    {
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readRecRateOnMismatch);
        
    }  else if (flag == "FitnessMapInfo")
    {   
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readFitnessMapInfo);

        
    } else if (flag == "resetTrackedT1Muts")
    {   
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readResetTrackedT1Muts);
        
    }  else if (flag == "fecundityForFitnessOfOne" || flag == "fec")
    {   
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readfecundityForFitnessOfOne);
        
    }  else if (flag == "T1_MutationRate" || flag == "T1_mu")
    {
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readT1_MutationRate);

    } else if (flag == "T5_MutationRate" || flag == "T5_mu")
    {
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readT5_MutationRate);

    } else if (flag == "additiveEffectAmongLoci")
    {
        std::cout << "You are using the option --additiveEffectAmongLoci. The option exists but should not be present in the manual as the option can't be used for the moment. Sorry! Fitness effects are only multiplicative among loci. If you want additivity please, let Remi know and he can eventually code it in for you. It would be quite quick to add this feature in.\n";
        abort();
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readadditiveEffectAmongLoci);
    } else if (flag == "selectionOn")
    {
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readSelectionOn);
    } else if (flag == "T1_FitnessEffects" || flag == "T1_fit")
    {
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readT1_FitnessEffects);
        
    } else if (flag == "T5_FitnessEffects" || flag == "T5_fit")
    {
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readT5_FitnessEffects);
        
    } else if (flag == "T1_Initial_AlleleFreqs" || flag == "T1_ini")
    {   
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readT1_Initial_AlleleFreqs);
        
    } else if (flag == "T5_Initial_AlleleFreqs" || flag == "T5_ini")
    {   
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readT5_Initial_AlleleFreqs);
        
    } else if (flag == "T1_EpistaticFitnessEffects" || flag == "T1_epistasis")
    {   
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readT1_EpistaticFitnessEffects);
    
    }  else if (flag == "T2_MutationRate" || flag == "T2_mu")
    {    
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readT2_MutationRate);
        
    }  else if (flag == "T2_FitnessEffects" || flag == "T2_fit")
    {       
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readT2_FitnessEffects);
  
    }  else if (flag == "T3_MutationRate" || flag == "T3_mu")
    {
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readT3_MutationRate);
        
    }  else if (flag == "T3_PhenotypicEffects" || flag == "T3_pheno")
    {    
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readT3_PhenotypicEffects);
        
    }  else if (flag == "T3_FitnessLandscape" || flag == "T3_fit")
    {    
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readT3_FitnessLandscape);
        
    } else if (flag == "T3_DN" || flag == "T3_DevelopmentalNoise")
    {    
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readT3_DevelopmentalNoise);
        
    } else if (flag == "T4_mu" || flag == "T4_MutationRate")
    {    
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readT4_MutationRate);
        
    } else if (flag == "T4_maxAverageNbNodesPerHaplotype")
    {
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readT4_maxAverageNbNodesPerHaplotype);
    } else if (flag == "Habitats" || flag == "H")
    {
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readHabitats);

    } else if (flag == "matingSystem")
    {
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readMatingSystem);

    }
    /*else if (flag == "TemporalChanges" || flag == "T")
    { 
        GP->readTemporalChanges(input);

    }*/  else if (flag == "Overwrite")
    {
        
        GP->OverwriteMode = input.GetNextElementInt();
        if (
            GP->OverwriteMode != 0 && 
            GP->OverwriteMode != 1 &&
            GP->OverwriteMode != 2  
            )
        {
            std::cout << "'OverwriteMode' = " << GP->OverwriteMode << ". Only values 0, 1 and 2 are accepted:\n\t0 = Don't overwrite\n\t1 = Overwrite even if Logfile is present but not if last output files are present (default)\n\t2 = Overwrite in any case" << std::endl;
            abort();
        }
    }  else if (flag == "readPopFromBinary")
    {
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readReadPopFromBinary);

    } else if (flag == "DryRun")
    {       
        GP->DryRun = true;

    }  else if (flag == "centralT1LocusForExtraGeneticInfo")
    {
        wrapperOverSpecies(input,&SpeciesSpecificParameters::readCentralT1LocusForExtraGeneticInfo);
    
    }  else if (flag == "S" || flag == "species" || flag == "SNames" || flag == "SpeciesNames")
    {
        
        if (SSPs.size() != 0)
        {
            std::cout << "Internal error. SSPs (the attribute of AllParameters containing all species specific parameters) had a length different from 0 before reading -speciesNames. \n";
            abort();
        }
        if (GP->nbSpecies != 0)
        {
            std::cout << "Internal error. GP->nbSpecies was different from 0 before reading -speciesNames (GP->nbSpecies = "<<GP->nbSpecies<<"). \n";
            abort();
        }

        if (input.PeakNextElementString() == "defaultNames")
        {
            input.skipElement();
            GP->nbSpecies = input.GetNextElementInt();
            if (GP->nbSpecies < 1)
            {
                std::cout << "In --speciesNames, received flag 'defaultNames' followed by a number that is lower than 1. SimBit needs at least one species to have something to simulate.\n";
                abort();
            }
            for (int speciesIndex = 0 ; speciesIndex < GP->nbSpecies ; speciesIndex++)
            {
                std::string sN = "sp_" + std::to_string(speciesIndex);
                GP->speciesNames.push_back(sN);
                SpeciesSpecificParameters ssp(sN, GP->nbSpecies);
                SSPs.push_back(ssp);
                assert(GP->speciesNames.size() == SSPs.size());
            }
        } else
        {
            while (input.IsThereMoreToRead())
            {
                std::string sN(input.GetNextElementString());
                if (sN == "seed")
                {
                    std::cout << "You choose to call a species 'seed'. This is the only species name that you cannot choose! Bad luck! The reason is that if you were to print the population to a binary file, then there would be conflict between the paths of the file saving the random seed with the file saving the species 'seed'. Here are some ideas for you on how to name your species; 'Seed', 'SEED', 'nut', 'berry' :D. Note that this error message is printed whether or not you asked to print populations to a binary file.\n";
                    abort();
                }
                SpeciesSpecificParameters ssp(sN, GP->nbSpecies);
                SSPs.push_back(ssp);
                GP->speciesNames.push_back(sN);
                GP->nbSpecies++;
                assert(GP->speciesNames.size() == SSPs.size());
                assert(GP->speciesNames.size() == GP->nbSpecies);
            }
        }

        if (!is_unique(GP->speciesNames))
        {
            std::cout << "Species names received are:\n";
            for (std::string& sN : GP->speciesNames) std::cout << "\t-" << sN << "\n";
            std::cout << "There are duplicate names. Please use unique names!\n";
            abort();
        }

        if (GP->nbSpecies == 0)
        {
            std::cout << "While reading --speciesNames, it appears that zero species names have been indicated. You need at least one species to simulate something.\n";
            abort();
        }
        if (GP->nbSpecies > 1)
        {
            std::cout << "WARNING: The current version is not ready for released because the system of species interaction is currently being modified. Please either simulate a single species at a time or consider a previously released version (or wait a little bit that I finish making these changes). Thank you and sorry!\n";
            abort();   
        }
    }  else if (flag == "eco" || flag == "ecoRelation" || flag == "speciesEcologicalRelationships")
    {   
        
        bool isRandomGauss;
        std::string randomGaussType;
        double randomGaussMean;
        double randomGaussSd;
        if (input.PeakNextElementString() == "randomGauss")
        {
            input.skipElement();
            isRandomGauss = true;
            randomGaussType = input.GetNextElementString();
            randomGaussMean = input.GetNextElementDouble();
            randomGaussSd = input.GetNextElementDouble();
            if (randomGaussSd < 0)
            {
                std::cout << "In --eco, received keyword 'randomGauss'. The keyword was followed by the type ("<<randomGaussType<<"), the mean ("<< randomGaussMean << ") and the SD ("<< randomGaussSd <<"). SD received is negative.\n";
                abort();
            }
        } else
        {
            isRandomGauss = false;
            randomGaussMean = 0.0;
            randomGaussSd = 0.0;
        }

        std::vector<std::vector<char>> transposeOfspeciesEcoRel_type;
        std::vector<std::vector<double>> transposeOfspeciesEcoRel_effect;
        for (int speciesIndex_from = 0 ; speciesIndex_from < GP->nbSpecies; speciesIndex_from++)
        {
            std::vector<char> fromOneSpecies_type;
            std::vector<double> fromOneSpecies_effect;
            fromOneSpecies_type.resize(GP->nbSpecies);
            fromOneSpecies_effect.resize(GP->nbSpecies);
            
            for (int speciesIndex_to = 0 ; speciesIndex_to < GP->nbSpecies; speciesIndex_to++)
            {
                if (speciesIndex_to == speciesIndex_from)
                {
                    fromOneSpecies_type[speciesIndex_from]   = '0';
                    fromOneSpecies_effect[speciesIndex_from] = 0.0;
                } else
                {
                    std::string type;
                    double effect;
                    if (isRandomGauss)
                    {   
                        type = randomGaussType;
                        std::normal_distribution<double> dist(randomGaussMean, randomGaussSd);
                        effect = dist(GP->mt);
                    } else
                    {
                        type = input.GetNextElementString();
                        effect = input.GetNextElementDouble();
                    }

                    if (type != "A" && type != "B" && type != "C" && type != "D" && type != "0")
                    {
                        std::cout << "In option --" << flag << ", the type received for the effect of species '" << this->SSPs[speciesIndex_from].speciesName << "' onto species '" << this->SSPs[speciesIndex_to].speciesName << " is '" << type << "'. Sorry only types are 'A', 'B', 'C' and '0' are recognized.\n\t- 'A' means 'effect is multiplied by the patchSize of the causal species'\n\t- 'B' means 'effect is multiplied by the patchSize of the causal species and divided by the patchSize of the recipient species'\n\t- 'C' means 'effect is divided by the patchSize of the recipient species'\n\t- 'D' means 'effect is multiplied by the patchSize of both the recipient species and the causal species'\n\t- '0' means no effect. It must therefore necessarily be followed by an effect (or magnitude) of 0.0.\n";
                        abort();
                    }
                    assert(type.size() == 1);

                    if (type == "0")
                    {
                        if (effect != 0.0)
                        {
                            std::cout << "In option --" << flag << ", the type received for the effect of species '" << this->SSPs[speciesIndex_from].speciesName << "' onto species '" << this->SSPs[speciesIndex_to].speciesName << " is '" << type << "'. For this type the effect (or magnitude) must be 0.0. Effect received = "<< effect <<"\n";
                        }
                    }
                    if (effect == 0.0)
                    {
                        type = "0";
                    }

                    fromOneSpecies_type[speciesIndex_to] = type.at(0);
                    fromOneSpecies_effect[speciesIndex_to] = effect;
                }
            }
            transposeOfspeciesEcoRel_type.push_back(fromOneSpecies_type);
            transposeOfspeciesEcoRel_effect.push_back(fromOneSpecies_effect);
        }
        GP->speciesEcoRel_type   = transposeSquareMatrix(transposeOfspeciesEcoRel_type);
        GP->speciesEcoRel_effect = transposeSquareMatrix(transposeOfspeciesEcoRel_effect);

        /*
        std::cout << "eco Types:\n";
        for (auto& a : GP->speciesEcoRel_type)
        {
            std::cout << "\t";
            for (auto& b : a)   
            {
                std::cout << "'" << b << "' ";
            }
            std::cout << "\n";
        }
        std::cout << "\n";

        std::cout << "eco Effects:\n";
        for (auto& a : GP->speciesEcoRel_effect)
        {
            std::cout << "\t";
            for (auto& b : a)   
            {
                std::cout << "'" << b << "' ";
            }
            std::cout << "\n";
        }
        std::cout << "\n";
        */
        
        
    } else if (flag == "growthK")
    {
        assert(GP->speciesEcoRel_type.size() == this->GlobalP.nbSpecies);
        assert(GP->speciesEcoRel_effect.size() == this->GlobalP.nbSpecies);
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readGrowthK);

    } else if (flag == "resetGenetics")
    {
        wrapperOverSpecies(input, &SpeciesSpecificParameters::readResetGenetics);

    } else
    {
        std::cout << "Received flag/option " << flag << ". This option is not recognized. Sorry. The error was caught after the renaming which highlights the first security gate has a bug and failed to find this mistake.\n";
        abort();
    }

    input.workDone();
}



bool AllParameters::isSpeciesIndexReadFromBinary(int speciesIndex)
{
    assert(speciesIndex < SSPs.size());
    return SSPs[speciesIndex].readPopFromBinary;
}






