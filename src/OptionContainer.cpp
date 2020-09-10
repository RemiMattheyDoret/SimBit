

/*void OptionContainer::received(std::string& optionNameReceived)
{
    bool isOptionFound = false;
    for (auto& option : options.size())
    {
        else if (optionNameReceived == option)
        {
            isOptionFound = true;
            option->received(this);
        }
    }
    if (!isOptionFound)
    {
        std::cout << "Internal error. OptionContainer::received received optionName " << optionNameReceived << " but no such option could be found in the list of options. While the error should have not been detected that late, you can probably fix the error by checking your input.\n";
        abort();
    }
}*/


void OptionContainer::listOptions()
{
    std::cout << "Here is the entire list of options: \n";
    for (auto& option : options)
    {
        std::cout << "\t" << option.getNamesToPrint() << "\n";
    }
}

int OptionContainer::getNbOptions()
{
    return options.size();
}

std::string OptionContainer::getOptionFirstName(int optionIndex)
{
    return options[optionIndex].optionNames[0];
}

std::string OptionContainer::renameFlag(std::string flag)
{
    for (auto& option : options)
    {
        for (std::string& optionName : option.optionNames)
        {
            if (optionName == flag)
            {
                return option.renameFlag();
            }
        }
    }

    // Levenshtein distance
    std::string* bestPossibleFlag;
    double bestLevenshteinDistance = DBL_MAX;
    for (auto& option : options)
    {
        for (std::string& possibleFlag : option.optionNames)
        {
            double mean_LevenD = mean_levenshtein_distance(possibleFlag, flag);
            if (mean_LevenD < bestLevenshteinDistance)
            {
                bestLevenshteinDistance = mean_LevenD;
                bestPossibleFlag = &possibleFlag;
            }
        }
    }
    assert((*bestPossibleFlag).size() > 0);


    this->listOptions();

    std::cout << "The option '--" << flag << "' was not found. ";
    if (bestLevenshteinDistance < 0.9)
    {
        std::cout << "Did you mean '"<< getOption(*bestPossibleFlag).getNamesToPrint() <<"'?";
    }
    std::cout << "\n";
    abort();
}

uint32_t OptionContainer::minVector(const std::vector<uint32_t> v)
{
    uint32_t min = INT_MAX;
    for (auto& elem : v) min = std::min(min,elem);
    return min;
}

double OptionContainer::mean_levenshtein_distance(const std::string& s1, const std::string& s2) 
{
    const std::uint32_t len1 = s1.size(), len2 = s2.size();
    std::vector<uint32_t> col(len2+1), prevCol(len2+1);
    
    for (uint32_t i = 0; i < prevCol.size(); i++)
        prevCol[i] = i;
    for (uint32_t i = 0; i < len1; i++) {
        col[0] = i+1;
        for (uint32_t j = 0; j < len2; j++)
                        // note that std::min({arg1, arg2, arg3}) works only in C++11,
                        // for C++98 use std::min(std::min(arg1, arg2), arg3)
            col[j+1] = minVector({ prevCol[1 + j] + 1, col[j] + 1, prevCol[j] + (tolower(s1[i])==tolower(s2[j]) ? 0 : 1) });
        col.swap(prevCol);
    }
    return (double) prevCol[len2] / (double) std::max(s1.size(), s2.size());
}



int OptionContainer::HowManyOptionsWereInitiated()
{
    int r = 0;
    for (auto& option : options)
    {
        if (option.howManyTimesReceivedYet > 0)
        {
            r++;
        }
    }
    return r;
}


bool OptionContainer::wasInitiated(std::string name)
{
    for (auto& option : options)
    {
        if (option == name)
        {
            return option.howManyTimesReceivedYet > 0;
        }
    }
    std::cout << "Internal error (although you might want to check your input too)! In  'OptionContainer::hasReceived' could not find option with name " << name << ".\n";
    abort();
}

Option& OptionContainer::getOption(std::string name)
{
    for (auto& option : options)
    {
        if (option == name)
        {
            return option;
        }
    }
    std::cout << "In  'OptionContainer::hasReceived' could not find option with name " << name << ". The error should have been detected when renaming the flags so there is a small internal error here. You might want to check your input anyway.\n\n";
    this->listOptions();
    abort();
}
