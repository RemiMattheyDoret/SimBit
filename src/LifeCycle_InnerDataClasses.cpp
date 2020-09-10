/////////////////////
/// HaplotypeData ///
/////////////////////


LifeCycle::HaplotypeData::HaplotypeData():patch(999),ind(9999999),segregationIndex(-9)
{}

LifeCycle::HaplotypeData::HaplotypeData(int a ,int b, int c, int d )
:patch(a), ind(b), segregationIndex(c), nbRecs(d)
{}


bool operator==(LifeCycle::HaplotypeData& lhs, LifeCycle::HaplotypeData& rhs)
{
    return lhs.ind == rhs.ind && lhs.patch == rhs.patch && lhs.segregationIndex == rhs.segregationIndex;
}

bool operator!=(LifeCycle::HaplotypeData& lhs, LifeCycle::HaplotypeData& rhs)
{
    return lhs.ind != rhs.ind || lhs.patch != rhs.patch || lhs.segregationIndex != rhs.segregationIndex;
}
    

//////////////////
/// CoupleData ///
//////////////////
        
LifeCycle::CoupleData::CoupleData(){}
LifeCycle::CoupleData::CoupleData(HaplotypeData& a, HaplotypeData& b)
{
    mother = a;
    father = b;
}
LifeCycle::CoupleData::CoupleData(HaplotypeData&& a, HaplotypeData&& b)
{
    mother = a;
    father = b;
}

///////////////////
/// ParentsData ///
///////////////////


void LifeCycle::ParentsData::resizeCloneInfo(const std::vector<int>& patchSizeNextGeneration)
{
    if (SSP->cloningRate != 0)
    {
        cloneInfo.resize(GP->PatchNumber);
        for (uint32_t patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
        {
            cloneInfo[patch_index].resize(patchSizeNextGeneration[patch_index]);
        }
    }
}        

LifeCycle::ParentsData::ParentsData(){}

void LifeCycle::ParentsData::resizeForNewGeneration(const std::vector<int>& patchSizeNextGeneration)
{
    cloneInfo.resize(GP->PatchNumber);
    couples.resize(GP->PatchNumber);
    lastOffspring.resize(2);
    lastOffspring[0].resize(GP->PatchNumber);
    lastOffspring[1].resize(GP->PatchNumber);
    for (uint32_t patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
    {
        lastOffspring[0][patch_index].resize(SSP->patchSize[patch_index], HaplotypeData(-1,-1,-1,-1));
        lastOffspring[1][patch_index].resize(SSP->patchSize[patch_index], HaplotypeData(-1,-1,-1,-1));
        couples[patch_index].resize(patchSizeNextGeneration[patch_index]);
        if (SSP->cloningRate != 0)
            cloneInfo[patch_index].resize(patchSizeNextGeneration[patch_index]);
    }
}


bool LifeCycle::ParentsData::shouldIClone(uint32_t patch_index, uint32_t ind_index)
{
    if (SSP->cloningRate == 0)
    { 
        return false;
    } else
    {
        return cloneInfo[patch_index][ind_index];
    }
}




bool LifeCycle::ParentsData::isLastOffspring(HaplotypeData& parent, uint32_t patch_index, uint32_t ind_index, uint32_t segregationIndex)
{
    auto& lastParentoffspring = lastOffspring[parent.segregationIndex][parent.patch][parent.ind];
    auto offspring = HaplotypeData(patch_index, ind_index, segregationIndex, 0);
    if (offspring == lastParentoffspring) // Does not compare nbRecs
    {
        return true;
    } else
    {
        return false;
    }
}


