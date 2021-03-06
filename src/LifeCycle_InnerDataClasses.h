class LifeCycle::HaplotypeData
{
public:
    int patch;
    int ind;
    int segregationIndex;
    int nbRecs;
    
    HaplotypeData();
    HaplotypeData(int a ,int b, int c, int d );

    friend bool operator==(HaplotypeData& lhs, HaplotypeData& rhs);
    friend bool operator!=(HaplotypeData& lhs, HaplotypeData& rhs);
};


class LifeCycle::CoupleData
{
public:
    HaplotypeData mother;
    HaplotypeData father;
    
    
    CoupleData();
    CoupleData(HaplotypeData& a, HaplotypeData& b);
    CoupleData(HaplotypeData&& a, HaplotypeData&& b);
};



class LifeCycle::ParentsData
{
public:
    std::vector<std::vector<bool>> cloneInfo;
    std::vector<std::vector<CoupleData>> couples;
    std::vector<std::vector<std::vector<HaplotypeData>>> lastOffspring;
    

    void resizeForNewGeneration(const std::vector<int>& patchSizeNextGeneration);
    ParentsData(); // Does nothing

    void resizeCloneInfo(const std::vector<int>& patchSizeNextGeneration);

    bool shouldIClone(uint32_t patch_index, uint32_t ind_index);
    bool isLastOffspring(HaplotypeData& parent, uint32_t patch_index, uint32_t ind_index, uint32_t segregationIndex);
};

