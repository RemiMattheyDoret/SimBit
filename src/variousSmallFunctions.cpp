/*

 Author: Remi Matthey-Doret

    MIT License

    Copyright (c) 2017 Remi Matthey-Doret

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

 */

template <class T>
bool is_unique(std::vector<T> x)
{
    std::sort( x.begin(), x.end() ); // O(N log N)
    return std::adjacent_find( x.begin(), x.end() ) == x.end();
}

template<typename T>
void printVector(std::vector<T> v)
{
  for (T& e : v)
    std::cout << e << " ";
  std::cout << std::endl;
}


template <typename T>
std::vector<uint32_t> reverse_sort_indexes(const std::vector<T> &v) {

  // initialize original index locations
  std::vector<uint32_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  std::sort(idx.begin(), idx.end(),
       [&v](uint32_t i1, uint32_t i2) {return v[i1] > v[i2];});

  return idx;
}

template <typename T>
std::vector<uint32_t> reverse_stable_sort_indexes(const std::vector<T> &v) {

  // initialize original index locations
  std::vector<uint32_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  std::stable_sort(idx.begin(), idx.end(),
       [&v](uint32_t i1, uint32_t i2) {return v[i1] > v[i2];});

  return idx;
}


template <typename T>
std::vector<uint32_t> sort_indexes(const std::vector<T> &v) {

  // initialize original index locations
  std::vector<uint32_t> idx(v.size());
  std::iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  std::sort(idx.begin(), idx.end(),
       [&v](uint32_t i1, uint32_t i2) {return v[i1] < v[i2];});

  return idx;
}

template <typename T>
void reorderNoAssertions(std::vector<T>& v, std::vector<uint32_t>& order)  
{   
    auto v2 = v;
    for (uint32_t i = 0 ; i < v.size() ; i++)
    {
        v[i] = v2[order[i]];
    }
}


template <typename T>
void reorder(std::vector<T>& v, const std::vector<uint32_t>& order, const unsigned char assertions = 0, const int PatchNumber = -1)  
{   
    auto v2 = v;
    T sum = 0;
    for (uint32_t i = 0 ; i < v.size() ; i++)
    {
        v[i] = v2[order[i]];
        if (assertions == 1)
        {
          assert(v[i] >= 0.0 && v[i] <= 1.0 );
          sum += v[i];
        }
        if (assertions == 2) assert(v[i] >= 0 && v[i] < PatchNumber);
    }
    if (v.size() && assertions == 1) assert(abs(sum - 1.0) < 0.00000001);
}

/*template< class T >
void reorder(std::vector<T> &v, std::vector<uint32_t> const &order )  {   
    for ( int s = 1, d; s < order.size(); ++ s ) {
        for ( d = order[s]; d < s; d = order[d] ) ;
        if ( d == s ) while ( d = order[d], d != s ) swap( v[s], v[d] );
    }
}*/


/*std::vector<int> rev_ordered(std::vector<double> const& values) {
    std::vector<int> indices(values.size());
    std::iota(begin(indices), end(indices), static_cast<int>(0));
    
    std::sort(
              begin(indices), end(indices),
              [&](int a, int b) { return values[a] > values[b]; }
              );
    return indices;
}*/


// String manipulation

// trim from start
/*static inline std::string &ltrimString(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(),
            std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
}
*/
static inline std::string &ltrimString(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), [](int c) {return !std::isspace(c);}));
    return s;
}

// trim from end
static inline std::string &rtrimString(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(),
            [](int c) {return !std::isspace(c);}).base(), s.end());
    return s;
}

// trim from both ends
static inline std::string &trimString(std::string &s) {
    return ltrimString(rtrimString(s));
}


void reduceString(std::string& str)
{
    // trim first
    trimString(str);

    // replace '\t' with ' '
    std::replace( str.begin(), str.end(), '\t', ' ');

    // replace sub ranges
    /*auto beginSpace = result.find_first_of(whitespace);
    while (beginSpace != std::string::npos)
    {
        const auto endSpace = result.find_first_not_of(whitespace, beginSpace);
        const auto range = endSpace - beginSpace;

        result.replace(beginSpace, range, fill);

        const auto newStart = beginSpace + fill.length();
        beginSpace = result.find_first_of(whitespace, newStart);
    }*/
}

// matrix manipulation

template< class T >
std::vector<std::vector<T>> transposeSquareMatrix(std::vector<std::vector<T>>& M) // This is slow but it does not matter!
{
  for (int i = 0 ; i < M.size() ; i++) 
  {
    assert(M[i].size() == M.size());
  }

  std::vector<std::vector<T>> r = M;
  for (int i = 0 ; i < M.size() ; i++) 
  {
    for (int j = 0 ; j < M.size() ; j++)
    {
      r[i][j] = M[j][i];
    }
  }
  return r;
}

template <typename T>
std::string to_string_with_precision(const T a_value, const int n = 6)
{
    std::ostringstream out;
    out << std::setprecision(n) << a_value;
    return out.str();
}


/*
Output code:
  A. All elements of `y` are present in `x`
  B. Some but not all elements of `y` are present in `x`
  C. No element of `y` is present in `x`
  U. Undefined (because `y` is empty)
This is a useless function as std::set_intersection does it better with sorted vectors
*/
template<typename T>
char areElementsOfYinX(std::vector<T> x, std::vector<T> y)
{
  if (y.size() == 0)
  {
    return 'U';
  }

  bool whereAnyElementFound = false;
  bool whereAnyElementNotFound = false;
  for (T& y_element : y)
  {
    bool isElementFound = std::find(x.begin(),x.end(),y_element) != x.end();
    if (isElementFound)
    {
      whereAnyElementFound = true;
      if (whereAnyElementNotFound)
      {
        return 'B';
      }
    } else
    {
      whereAnyElementNotFound = true;
      if (whereAnyElementFound)
      {
        return 'B';
      }
    }
  }
  if (whereAnyElementFound)
  {
    return 'A';
  } else if (whereAnyElementNotFound)
  {
    return 'C';
  }
  abort();
}

template<typename T>
void remove_a_from_b(std::vector<T>& a, std::vector<T>& b) 
{
  // remove from vector b all the elements that also exist in vector a.
  // Both vectors must be sorted
  
  auto ia = std::begin(a);
  auto iter = std::remove_if (
       std::begin(b), std::end(b),
       [&ia,&a](T& bo) -> bool {
                       while  (ia != std::end(a) && *ia < bo) ++ia;
                       return (ia != std::end(a) && *ia == bo);
                     });
  (void) iter;
}

template<typename T>
std::vector<T> intersection(std::vector<T>& a, std::vector<T>& b) 
{
  // Both vectors must be sorted
  std::vector<T> r;
  //std::cout << "\na:";
  //for (auto& ae : a)
  //  std::cout << ae.locus << " ";
  //std::cout << "\nb:";
  //for (auto& be : b)
  //  std::cout << be.locus << " ";
  std::set_intersection(a.begin(), a.end(), b.begin(), b.end(), back_inserter(r));
  //std::cout << "\nr:";
  //for (auto& re : r)
  //  std::cout << re.locus << " ";
  //std::cout << "\n";
  return r;
}


template<typename T>
void sortAndRemoveDuplicates(std::vector<T>& vec)
{
  sort( vec.begin(), vec.end() );
  vec.erase( unique( vec.begin(), vec.end() ), vec.end() );
}

template <typename INT, typename T> // INT could be int, uint32_t, char, uint32_t, etc...
void removeIndicesFromVector(std::vector<T>& v, std::vector<INT>& rm )
{
  uint32_t rm_index = 0;
  v.erase(
    std::remove_if(std::begin(v), std::end(v), [&](T& elem)
    {
      if (rm.size() != rm_index && &elem - &v[0] == rm[rm_index])
      {
        rm_index++;
        return true;
      }
      return false;
    }),
    std::end(v)
  );
}


std::vector<uint32_t> whichUnsignedInt(std::vector<bool> b)
{
  std::vector<uint32_t> r;
  for (unsigned i = 0 ; i < b.size() ; ++i)
  {
    if (b[i])
    {
      r.push_back(i);
    }
  }
  return r;
}


std::vector<T1_locusDescription> whichT1_locusDescription(std::vector<bool> b)
{
  std::vector<T1_locusDescription> r;
  for (uint32_t i = 0 ; i < b.size() ; ++i)
  {
    if (b[i])
    {
      r.push_back(T1_locusDescription(i/8, i%8, i));
    }
  }
  return r;
}

template<typename INT, typename INTT>
std::vector<bool> inverseWhich(std::vector<INT> idx, INTT size )
{
  assert(idx.size() <= size);
  std::vector<bool> r(size, false);

  for (auto& i : idx)
  {
    assert(i >= 0 && i < size);
    assert(!r[i]);
    r[i] = true;
  }

  return r;
}


template<typename INT>
std::vector<std::vector<double>> getFullFormMatrix(std::vector<std::vector<double>>& rates, std::vector<std::vector<INT>>& indices)
{
  assert(rates.size() == indices.size());
  auto nbPatches = rates.size(); // Cannot use GP as we aregoing through different generation index
  std::vector<std::vector<double>> r(nbPatches);

  for (uint32_t from_patch_index = 0 ; from_patch_index < nbPatches ; ++from_patch_index)
  {
    r[from_patch_index].resize(nbPatches, 0.0); // important to initialize to zero
    assert(rates[from_patch_index].size() == indices[from_patch_index].size());
    for (uint32_t fake_to_patch_index = 0 ; fake_to_patch_index < rates[from_patch_index].size() ; ++fake_to_patch_index)
    {
      auto& to_patch_index = indices[from_patch_index][fake_to_patch_index];
      auto& rate = rates[from_patch_index][fake_to_patch_index];
      assert(rate >= 0.0 && rate <= 1.0);
      r[from_patch_index][to_patch_index] = rate;
    }

    double sum = 0.0;
    for (uint32_t to_patch_index = 0 ; to_patch_index < nbPatches ; ++to_patch_index)
    {
      sum += r[from_patch_index][to_patch_index];
    }
    assert(fabs(sum - 1.0) < 0.00000001);
  }

  return r;
}


template<typename INT>
std::vector<std::vector<std::vector<double>>> getFullFormMatrix(std::vector<std::vector<std::vector<double>>>& __rates, std::vector<std::vector<std::vector<INT>>>& __indices)
{
  assert(__rates.size() == __indices.size());
  std::vector<std::vector<std::vector<double>>> __r(__indices.size());
  for (uint32_t generation_index = 0 ; generation_index < __rates.size() ; ++generation_index)
  {
    __r[generation_index] = getFullFormMatrix(__rates[generation_index], __indices[generation_index]);
  }

  return __r;
}


template<typename A, typename B>
A myMin(A& a, B& b)
{
  if (a < (A) b)
    return a;
  else
    return (A) b;
}

template<typename A, typename B>
A myMax(A& a, B& b)
{
  if (a < (A) b)
    return (A) b;
  else
    return a;
}

template<typename A, typename B>
A myAbs(A& a, B& b)
{
  if (a > b)
  {
    return a - b;
  } else
  {
    return b - a;
  }
}

void tokenize(std::vector<std::string>& tokens, const std::string& s, const char sep)
{
  size_t pos = 0;
  while (true)
  {
    size_t nextPos = s.find(pos, sep);
    //std::cout << "pos = " << pos << "\n";
    //std::cout << "nextPos = " << nextPos << "\n";
    assert(nextPos > pos);
    //std::cout << "substr is '" << token.substr(pos, nextPos-pos) << "'\n";
    tokens.push_back(s.substr(pos, nextPos-pos));
    if (nextPos == std::string::npos)
      break;
    pos = nextPos + 1;
  }
}



template<typename T>
void replace(typename std::vector<T>& A, typename std::vector<T>::iterator Afrom, typename std::vector<T>::iterator Ato, typename std::vector<T>::iterator Bfrom, typename std::vector<T>::iterator Bto)
{
  /*
  std::cout << "replacing '";
  for (auto it = Afrom ; it != Ato ; ++it)
    std::cout << *it << " ";
  std::cout << "' at position "<< Afrom - A.begin() << " by '";
  for (auto it = Bfrom ; it != Bto ; ++it)
    std::cout << *it << " ";
  std::cout << "'\n";

  std::cout << "before replacing string is '";
  for (size_t i = 0 ; i < A.size() ; ++i)
    std::cout << A[i] << " ";
  std::cout << "'\n";
  */

  auto nbToRemove = Ato - Afrom;
  auto nbToAdd = Bto - Bfrom;
  auto nbToCopy = nbToRemove > nbToAdd ? nbToAdd : nbToRemove ;
  assert(nbToCopy >= 0);
  auto copyTo = Afrom + nbToCopy;
  while ( Afrom < copyTo )
  {
    *Afrom = *Bfrom;
    ++Afrom;
    ++Bfrom;
  }

  if (Bfrom != Bto)
  {
    assert(Afrom == Ato);
    assert(Bfrom < Bto);
    A.insert(Ato, Bfrom, Bto);
  } else if (Afrom != Ato)
  {
    assert(Afrom < Ato);
    assert(Bfrom == Bto);
    A.erase(Afrom, Ato);
  }

  /*
  std::cout << "after replacing string is '";
  for (size_t i = 0 ; i < A.size() ; ++i)
    std::cout << A[i] << " ";
  std::cout << "'\n";
  */
}

bool isLocusConsideredNeutral(double w, double appro)
{
  if (w == 1.0) return true;

  if (w < 1.0)
  {
    auto threshold = 1.0 / (1.0 + appro);
    if (w > threshold) return true; else return false;
  } else // w  >= 1.0
  {
    auto threshold = 1.0 + appro;
    if (w < threshold) return true; else return false;
  }
}


template<typename T>
void mergeSortAppendOddNumber(std::vector<T>& a, std::vector<T>& b, std::vector<T>& c)
{
  /*
    Appends to a, the sorted values from the sorted vectors b and c
  */

  // Assertions
  if (a.size())
  {
    if (b.size())
    {
      assert(a.back() < b.front());
    }

    if (c.size())
    {
      assert(a.back() < c.front());
    }   
  }

  // isThereAnythingToDo
  if (a.size()==0 && b.size()==0 )
    return;


  // Handy variables
  auto bit = b.begin();
  auto cit = c.begin();
  auto bend = b.end();
  auto cend = c.end();

  // Reserve memory in a
  a.reserve(a.size() + b.size() + c.size());


  // Copy stuff
  while (bit != bend && cit != cend)
  {
    if (*bit < *cit)
    {
      if (a.back() == *bit)
      {
        a.pop_back();
      } else
      {
        a.push_back(*bit);
      }
        
      ++bit;
    } else
    {
      if (a.back() == *cit)
      {
        a.pop_back();
      } else
      {
        a.push_back(*cit);
      }

      ++cit;
    }
  }

  // Copy remaining elements of b
  if (bit != cend)
  {
    // first element remaining
    if (a.back() == *bit)
    {
      a.pop_back();
    } else
    {
      a.push_back(*bit);
    }
    ++bit;

    // other elements
    while (bit != cend)
    {
      a.push_back(*bit);
      ++bit;
    }
  }

  // Copy remaining elements of c
  if (cit != cend)
  {
    // first element remaining
    if (a.back() == *cit)
    {
      a.pop_back();
    } else
    {
      a.push_back(*cit);
    }
    ++cit;

    // other elements
    while (cit != cend)
    {
      a.push_back(*cit);
      ++cit;
    }
  }
}


template<typename T>
void filter(std::vector<T>& v, const std::vector<bool>& b)
{
    assert(v.size() == b.size());
    auto it = b.begin();
    v.erase(
      std::remove_if(
        v.begin(),
        v.end(),
        [&](T) { return *it++; }
      ),
      v.end()
    );
}


template<typename T>
void dealWithOddAndEvenDuplicates(std::vector<T>& x)
{
  auto begin = x.begin();
  auto end = x.end();

  if (end - begin <= 1) return;

  std::vector<bool> toErase(end - begin, false);

  auto it = begin;

  T elem = *it;
  size_t nb = 1;
  ++it;
  while (it != end)
  {
    if (*it == elem)
    {
      ++nb;
    }
    else
    {
      if (nb != 1)
      {
        size_t eraseTo = it - begin;
        size_t eraseFrom = it - begin - nb;
        
        if (nb % 2 == 0)
        {
          // even -> keep none
          // nothing to do
        } else
        {
          // odd -> keep one
          ++eraseFrom;
        }

        for (size_t i = eraseFrom ; i < eraseTo ; ++i)
          toErase[i] = true;

        nb = 1;
      }
      elem = *it;
    }
    ++it;
  }

  filter(x, toErase);
  if (x.capacity() > x.size() * 1.2)
    x.shrink_to_fit();
}

//// Insert multiple elements in sorted vector
template<typename T>
void insertInSortedVector(std::vector<T>& x, std::vector<T>& values)
{
  if (values.size() == 0) return;
  if (x.size() == 0) {values.swap(x); return;}
  std::vector<size_t> positions;

  for (auto& value : values)
  {
    positions.push_back(std::lower_bound(x.begin(), x.end(), value));
  }

  insertAtPositions(x, values, positions);
}


//// Insert multiple elements at specified positions into vector
template<typename T>
void insertAtPositions(std::vector<T>& x, std::vector<T>& values, std::vector<size_t>& positions)
{
  // assert values and positions are the same size
  assert(values.size() == positions.size());

  // Special case - values is empty
  if (values.size() == 0) return;

  // Special case - single value to insert
  if (values.size() == 1)
  {
    x.insert(x.begin()+positions.front(), values.front());
    return;
  }

  // sort the values and the positions where those values should be inserted
  auto indices = sort_indexes(positions);
  reorderNoAssertions(positions, indices);
  reorderNoAssertions(values, indices);

  // Special case - x is empty
  if (x.size() == 0)
  {
      x.swap(values);
      return;
  }
  
  // Allocate memory to x
  x.resize(x.size() + values.size());
  
  // Move things to make room for insertions and insert
  int pos_index = positions.size()-1;
  for (size_t i = x.size()-1 ; pos_index >= 0 ; --i)
  {
    if (i == positions[pos_index] + pos_index)
    {
      // A new value should go at index i
      x[i] = std::move(values[pos_index]);
      --pos_index;
    } else
    {
      // The value from index 'i-pos_index-1' must go to index 'i'
      x[i] = std::move(x[i-pos_index-1]);
    }
  }
}

// Insert multiple elements (non-sorted) in sorted vector.
template<typename T>
void insertMultipleElements(std::vector<T>& x, std::vector<T>& values)
{
  if (values.size() == 0 ) return;
  if (values.size() == 1 )
  {
    x.insert(
      std::lower_bound(x.begin(), x.end(), values.front()),
      values.front()
    );

  } else
  {
    /*
    // Technique 1
    I think 'insertAtPositions' is not working properly
    std::vector<size_t> positions;
    positions.reserve(values.size());
    for (auto& value : values)
    {
      positions.push_back(std::lower_bound(x.begin(), x.end(), value) - x.begin());
    }
    insertAtPositions(x, values, positions);
    */

    
    // Technique 2

    for (auto& value : values)
    {
      x.insert(std::lower_bound(x.begin(), x.end(), value), value);
    }
  }
}


//// Erase multiple elements at specified positions into vector
template<typename T>
void eraseAtPositions(std::vector<T>& x, std::vector<size_t>& positions)
{
  // Special case - values is empty
  if (positions.size() == 0) return;

  // Special case - a single value to remove
  if (positions.size() == 1)
  {
    x.erase(positions.front());
    return;
  }

  // sort the positions
  std::sort(positions.begin(), positions.end());
  
  // Move and erase things
  int pos_index = 0;
  for (size_t i = positions.front() ; i < x.size() - positions.size() ; ++i)
  {
    if (pos_index < positions.size() && i == positions[pos_index])
      ++pos_index;
    x[i] = std::move(x[i + pos_index]);
  }
}


template<typename T>
void insertEraseAtPositions(
    std::vector<T>& x,                      
    std::vector<T>& valuesToInsert,        
    std::vector<size_t>& positionsToInsert,
    std::vector<size_t>& positionsToErase
) {
    std::sort(positionsToErase.begin(), positionsToErase.end());

    static std::vector< std::pair< size_t, T > > toInsert;

    toInsert.clear();
    toInsert.reserve( positionsToInsert.size() );
    for (long i = 0; i < positionsToInsert.size(); i++)
        toInsert.emplace_back( positionsToInsert[i], valuesToInsert[i] );
    std::sort(toInsert.begin(), toInsert.end());

    static std::vector<T> y;

    y.clear();
    y.reserve( x.size() - positionsToErase.size() + toInsert.size() );

    long j = 0, k = 0;
    for (long i = 0; i < x.size(); i++)
    {
        if (i == toInsert[k].first)
            y.push_back( toInsert[k++].second );

        if (( j < positionsToErase.size() ) && ( i == positionsToErase[j] ))
            j++;
        else
            y.push_back( x[i] );
    }

    for ( ; k < toInsert.size(); k++)
        y.push_back( toInsert[k].second );

    std::swap(x, y);
}


template<typename T, typename INT>
void insertionSort(std::vector<T>& v, INT firstUnsortedElement = 1)
{
  assert(firstUnsortedElement > 0);
  int i, key, j;
  for (i = firstUnsortedElement; i < v.size(); i++)
  {
    key = v[i];
    j = i - 1;

    /* Move elements of v[0..i-1], that are
    greater than key, to one position ahead
    of their current position */
    while (j >= 0 && v[j] > key)
    {
      v[j + 1] = v[j];
      j = j - 1;
    }
    v[j + 1] = key;
  }
}



/*
//// Insert multiple elements at specified positions into vector
template<typename T>
void insertEraseAtPositions(std::vector<T>& x, std::vector<T>& insertValues, std::vector<size_t>& insertPositions, std::vector<size_t>& erasePositions)
{
  // assert values and positions are the same size
  assert(erasePositions.size() == insertPositions.size());

  // return if nothing to do
  if (insertPositions.size() == 0 && erasePositions.size() == 0) return;

  // If only need to erase
  if (insertPositions.size() == 0)
  {
    eraseAtPositions(x, erasePositions);
    return;
  }

  // If only need to insert
  if (erasePositions.size() == 0)
  {
    insertAtPositions(x, insertValues, insertPositions);
    return;
  }  

  // sort the values and the positions where those values should be inserted
  {
    auto indices = sort_indexes(insertPositions);
    reorderNoAssertions(insertPositions, indices);
    reorderNoAssertions(insertValues, indices);
  }
  
  // Allocate memory to x
  x.resize(x.size() + insertPositions.size());
  
  // Move things to make room for insertions and insert
  int pos_index = positions.size()-1;
  for (size_t i = x.size()-1 ; pos_index >= 0 ; --i)
  {
    if (i == positions[pos_index] + pos_index)
    {
      // A new value should go at index i
      x[i] = std::move(values[pos_index]);
      --pos_index;
    } else
    {
      // The value from index 'i-pos_index-1' must go to index 'i'
      x[i] = std::move(x[i-pos_index-1]);
    }
  }
}
*/


