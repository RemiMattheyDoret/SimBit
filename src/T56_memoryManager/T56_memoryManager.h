class T56_memoryManager
{
private:
    int attempt_shrinking_every_N_generation;
    bool shouldRecordInfo;
    bool shouldShrink;
    uint32_t maxNtrlSizeLastGeneration;
    uint32_t maxSelSizeLastGeneration;
    

public:
    void doStuff(Haplotype& haplo);
    void setGenerationInfo();
    template<typename INT>
    void set_attempt_shrinking_every_N_generation(INT i);
};