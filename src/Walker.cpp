
Walker::Walker():N_m_one(0){}

Walker::Walker(const std::vector<double>& cumSumProbs)
: N_m_one(cumSumProbs.size()-1)
{
  const size_t N = cumSumProbs.size();

  // Compute probabilities
  std::vector<double> probs;
  probs.reserve(N);
  probs.push_back( cumSumProbs[0] / cumSumProbs.back() );
  for (size_t ind_index = 1 ; ind_index < N ; ++ind_index)
  {
    auto diff = cumSumProbs[ind_index] - cumSumProbs[ind_index-1];
    assert(diff >= 0.0);
    probs.push_back( diff / cumSumProbs.back() );
  }
  //std::cout << "std::accumulate(probs.begin(), probs.end(), 0.0) = " << std::accumulate(probs.begin(), probs.end(), 0.0) << "\n";
    

  // Prepare alias
  alias = std::vector<std::pair<double, size_t>>(N, {0.0, std::numeric_limits<size_t>::max()});
  std::queue<size_t>  small, large;

  for (size_t i = 0; i != N; ++i)
  {
    alias[i].first = N * probs[i]; // This is the expected number of offspring
    if (alias[i].first < 1.0)
    {
      small.push(i);
    } else
    {
      large.push(i);
    }
  }

  while (!small.empty() && !large.empty())
  {
    auto s = small.front(), l = large.front();
    small.pop(), large.pop();
    alias[s].second = l;
    alias[l].first -= (1.0 - alias[s].first);

    if (alias[l].first < 1.0)
    {
      small.push(l);
    } else
    {
      large.push(l);
    }
  }

  assert(alias.size() == N);
}

size_t Walker::operator()(RNG_wrapper& rngw) const
{
  //std::mt19937_64 mt(1253);
  //std::uniform_int_distribution<size_t> d(0, N_m_one);
  //std::uniform_real_distribution<size_t> d2(0,1.0);
  const size_t idx  = rngw.uniform_int_distribution(N_m_one);
  //const size_t idx  = d(mt);
  /*std::cout << "N_m_one = " << N_m_one << "\n";
  std::cout << "idx = " << idx << "\n";
  std::cout << "alias.size() = " << alias.size() << "\n";
  */
  
  assert(alias.size() > idx);
  if (
    alias[idx].second != std::numeric_limits<size_t>::max()
    &&
    rngw.uniform_real_distribution(N_m_one) >= alias[idx].first
  )
  {
    return alias[idx].second;
  } else
  {
    return idx;
  }
}


Walker& Walker::operator=(const Walker&& other)
{
  this->alias = std::move(other.alias);
  this->N_m_one = std::move(other.N_m_one);
  return *this;
}

Walker& Walker::operator=(const Walker& other)
{
  this->alias = other.alias;
  this->N_m_one = other.N_m_one;
  return *this;
}

Walker::Walker(const Walker&& other)
{
  this->alias = std::move(other.alias);
  this->N_m_one = std::move(other.N_m_one);
}

Walker::Walker(const Walker& other)
{
  this->alias = other.alias;
  this->N_m_one = other.N_m_one;
}


