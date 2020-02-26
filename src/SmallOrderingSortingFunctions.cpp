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
std::vector<size_t> reverse_sort_indexes(const std::vector<T> &v) {

  // initialize original index locations
  std::vector<size_t> idx(v.size());
  iota(idx.begin(), idx.end(), 0);

  // sort indexes based on comparing values in v
  std::sort(idx.begin(), idx.end(),
       [&v](size_t i1, size_t i2) {return v[i1] > v[i2];});

  return idx;
}

template <typename T>
void reorder(std::vector<T>& v, std::vector<size_t>& order)  
{   
    auto v2 = v;
    for (size_t i = 0 ; i < v.size() ; i++)
    {
        v[i] = v2[order[i]];
    }
}

/*template< class T >
void reorder(std::vector<T> &v, std::vector<size_t> const &order )  {   
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
static inline std::string &ltrimString(std::string &s) {
    s.erase(s.begin(), std::find_if(s.begin(), s.end(),
            std::not1(std::ptr_fun<int, int>(std::isspace))));
    return s;
}

// trim from end
static inline std::string &rtrimString(std::string &s) {
    s.erase(std::find_if(s.rbegin(), s.rend(),
            std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
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
  std::set_intersection(a.begin(), a.end(), b.begin(), b.end(), back_inserter(r));
  return r;
}


template<typename T>
void sortAndRemoveDuplicates(std::vector<T>& vec)
{
  sort( vec.begin(), vec.end() );
  vec.erase( unique( vec.begin(), vec.end() ), vec.end() );
}

template <typename INT, typename T> // INT could be int, unsigned int, char, size_t, etc...
void removeIndicesFromVector(std::vector<T>& v, std::vector<INT>& rm )
{
  size_t rm_index = 0;
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


std::vector<unsigned int> whichUnsignedInt(std::vector<bool> b)
{
  std::vector<unsigned int> r;
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
  for (size_t i = 0 ; i < b.size() ; ++i)
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

