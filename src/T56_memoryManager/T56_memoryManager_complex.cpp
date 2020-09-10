void T56_memoryManager::T5_attempt_shrinking(std::vector<uint32_t>& ntrl, std::vector<uint32_t>& sel) const
{
    if (GP->CurrentGeneration % attempt_shrinking_every_N_generation)
    {
        if (SSP->Gmap.T5ntrl_nbLoci)
        {
            newCapacity = T5ntrl_size_to_shrink_to();
            if (ntrl.capacity() * 1.2 > newCapacity && ntrl.size() < newCapacity)
            {
                ntrl.shrink_to_fit();
                ntrl.reserve(newCapacity);
            }
        }

        if (SSP->Gmap.T5sel_nbLoci)
        {
            newCapacity = SSP->T56_memManager.T5sel_size_to_shrink_to();
            if (sel.capacity() * 1.2 > newCapacity && sel.size() < newCapacity)
            {
                sel.shrink_to_fit();
                sel.reserve(newCapacity);
            }
        }
    }
}

template<typename INT1, typename INT2>
void T5x_reallocate(std::vector<uint32_t>& v, const INT1 sizeNeeded, const INT2 maxSizeLastGeneration, const size_t expected_size_needed) const
{
    if (sizeNeeded > v.capacity())
    {
        if (sizeNeeded > expected_size_needed)
        {
            v.reserve(sizeNeeded * 1.2);
        } else
        {
            auto newCapacity = 4 * v.capacity() > expected_size_needed ? 4 * v.capacity() : expected_size_needed ;
            v.reserve(newCapacity);
        }
    }
}

template<typename INT>
void T56_memoryManager::T5ntrl_reallocate(std::vector<uint32_t>& v, const INT sizeNeeded) const
{
    T5x_reallocate(v, sizeNeeded, T5ntrl_maxSizeLastGeneration, T5ntrl_expected_size_needed);
}

template<typename INT>
void T56_memoryManager::T5sel_reallocate(std::vector<uint32_t>& v, const INT sizeNeeded) const
{
    T5x_reallocate(v, sizeNeeded, T5sel_maxSizeLastGeneration, T5sel_expected_size_needed);
}


template<typename INT>
void T5ntrl_indicate_size(const INT size)
{
    if (T5ntrl_maxSizeLastGeneration < size)
        T5ntrl_maxSizeLastGeneration = size;
}

template<typename INT>
void T5sel_indicate_size(const INT size)
{
    if (T5sel_maxSizeLastGeneration < size)
        T5sel_maxSizeLastGeneration = size;
}

