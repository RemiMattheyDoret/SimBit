
class Walker
{
private:
  std::vector<std::pair<long double, uint32_t>> alias;
  uint32_t N_m_one;
 
public:
  Walker();

  template<typename T>
  Walker(const std::vector<T>& cumSumProbs);

  Walker& operator=(const Walker&& other);
  Walker& operator=(const Walker& other);
  Walker(const Walker&& other);
  Walker(const Walker& other);
  
  uint32_t operator()(RNG_wrapper& rngw) const;
};
