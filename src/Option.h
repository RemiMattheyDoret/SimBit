class Option
{
    friend class OptionContainer;
    public:
    std::vector<std::string> optionNames;
    std::vector<std::string> optionsThatMustBeReadBefore;
    bool wasInitiatedYet;
    bool wasRenamedYet;

    Option(std::vector<std::string> oN, std::vector<std::string> oTMBRB);
    
    bool operator==(std::string other);
    std::string getNamesToPrint();

    void received(OptionContainer& optionContainer);
    std::string renameFlag();
};
