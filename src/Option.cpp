
Option::Option(std::vector<std::string> oN, std::vector<std::string> oTMBRB, int CGSI)
: optionNames(oN), optionsThatMustBeReadBefore(oTMBRB), howManyTimesReceivedYet(0), wasRenamedYet(false), canGetSeveralInputs(CGSI)
{}

bool Option::operator==(std::string other)
{
    assert(optionNames.size() > 0);
    for (int nameIndex = 0 ; nameIndex < optionNames.size() ; nameIndex++)
    {
        if (optionNames[nameIndex] == other)
        {
            return true;
        }
    }
    return false;
}

std::string Option::getNamesToPrint()
{
    std::string s;
    s += "--" + optionNames[0];
    if (optionNames.size() > 1)
    {
        s += " (";
        for (int i = 1 ; i < optionNames.size() ; i++)
        {
            s += "--" + optionNames[i];
            if (i != optionNames.size() - 1)
            {
                s += ", ";
            }
        }
        s += ")";
    }
    //s += " ";
    return s;
}

void Option::received(OptionContainer& optionContainer)
{
    // Check it was not received yet
    if (!canGetSeveralInputs && howManyTimesReceivedYet > 0)
    {
        std::cout << "Option " << this->getNamesToPrint() << " was present at least "<<howManyTimesReceivedYet+1<<" times but can only be received once.\n";
        abort();
    }

    howManyTimesReceivedYet++;

    // check the options that must be set before
    for (std::string& preOptionName : optionsThatMustBeReadBefore)
    {
        bool wasOptionFoundSomewhere = false;
        for (auto& option : optionContainer.options ) // 'options' is from outer class
        {
            if (option == preOptionName)
            {
                wasOptionFoundSomewhere = true;
                if (this->canGetSeveralInputs == 0)
                {
                    std::cout << "Internal error: Read option " << this->getNamesToPrint() << " before option " << option.getNamesToPrint() << ".\n";
                    abort();
                }
            }
        }
        if (!wasOptionFoundSomewhere)
        {
            std::cout << "Internal error in OptionContainer, inner class Option. option " << this->getNamesToPrint() << " required option name preOptionName = " << preOptionName << " to be received before but there is no option which name matches this preOptionName.\n";
            abort();
        }
    }
}


std::string Option::renameFlag()
{
    if (this->wasRenamedYet)
    {
        std::cout << "Option " << this->getNamesToPrint() << " was received more than once!\n";
        abort();
    }
    this->wasRenamedYet = true;

    assert(this->optionNames.size() > 0);
    return this->optionNames[0];
}

