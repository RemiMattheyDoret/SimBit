class T56_memoryManager
{
private:
    int attempt_shrinking_every_N_generation;
    
    template<typename INT1, typename INT2>
    void T5x_reallocate(std::vector<uint32_t>& v, const INT1 sizeNeeded, const INT2 maxSizeLastGeneration, const size_t expected_size_needed) const;

public:
    void T5_attempt_shrinking(std::vector<uint32_t>& ntrl, std::vector<uint32_t>& sel) const;
    template<typename INT>
    void T5ntrl_reallocate(std::vector<uint32_t>& v, const INT sizeNeeded) const;
    void T5sel_reallocate(std::vector<uint32_t>& v, const INT sizeNeeded) const;
    size_t T5ntrl_size_to_shrink_to() const;
    size_t T5sel_size_to_shrink_to() const;


    template<typename INT>
    void T5ntrl_indicate_size(const INT size);
    template<typename INT>
    void T5sel_indicate_size(const INT size);


};