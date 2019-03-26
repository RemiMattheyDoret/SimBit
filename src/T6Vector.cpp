// general iterator
template<typename PrefixT, typename SuffixT>
T6VectorIterator<PrefixT, SuffixT>::T6VectorIterator(const T6VectorIterator& other)
{
    this->pPrefix = other.pPrefix;
    this->pSuffix = other.pSuffix;
}


template<typename PrefixT, typename SuffixT>
~T6VectorIterator<PrefixT, SuffixT>::T6VectorIterator()
{
    pPrefix = nullptr;
    pSuffix = nullptr;    
}


template<typename PrefixT, typename SuffixT>
T6VectorIterator& T6VectorIterator<PrefixT, SuffixT>::operator=(const T6VectorIterator& other)
{
    return T6VectorIterator(other);
}


template<typename PrefixT, typename SuffixT>
T6VectorIterator& T6VectorIterator<PrefixT, SuffixT>::operator++() //prefix increment
{
    pSuffix++;
    if ( pSuffix == pPrefix->end() )
    {
        pPrefix++;
        pSuffix = pPrefix->begin();
    }

    return *this;
}

template<typename PrefixT, typename SuffixT>
T6VectorIterator& T6VectorIterator<PrefixT, SuffixT>::operator++() //postfix increment
{
    auto r = *this;
    pSuffix++;
    if ( pSuffix == pPrefix->end() )
    {
        pPrefix++;
        pSuffix = pPrefix->begin();
    }

    return r;
}


template<typename PrefixT, typename SuffixT>
reference T6VectorIterator<PrefixT, SuffixT>::operator*() const
{
    *pPrefix * this->multiplicatorForPrefix + *pSuffix;
}


template<typename PrefixT, typename SuffixT>
friend void T6VectorIterator<PrefixT, SuffixT>::swap(T6VectorIterator& lhs, T6VectorIterator& rhs)
{
    auto tmp = rhs;
    rhs = lhs;
    lhs = tmp;
} //C++11 I think

// bidrectional
template<typename PrefixT, typename SuffixT>
T6VectorIterator& T6VectorIterator<PrefixT, SuffixT>::operator--() //prefix decrement
{
    if ( pSuffix == pPrefix->begin() )
    {
        pPrefix--;
        pSuffix = pPrefix->end() - 1;
    } else
    {
        pSuffix--;
    }

    return *this;
}


template<typename PrefixT, typename SuffixT>
T6VectorIterator T6VectorIterator<PrefixT, SuffixT>::operator--(int) //postfix decrement
{
    auto r = *this;
    if ( pSuffix == pPrefix->begin() )
    {
        pPrefix--;
        pSuffix = pPrefix->end() - 1;
    } else
    {
        pSuffix--;
    }

    return r;
}


//random_access_T6VectorIterator
template<typename PrefixT, typename SuffixT>
friend bool T6VectorIterator<PrefixT, SuffixT>::operator<(const T6VectorIterator& lhs, const T6VectorIterator& rhs)
{
    if (lhs.pPrefix == rhs.pPrefix)
    {
        return lhs.pSuffix < rhs.pSuffix;
    } else
    {
        return lhs.pPrefix < rhs.pPrefix;
    }
}

template<typename PrefixT, typename SuffixT>
friend bool T6VectorIterator<PrefixT, SuffixT>::operator>(const T6VectorIterator& lhs, const T6VectorIterator& rhs)
{
    if (lhs.pPrefix == rhs.pPrefix)
    {
        return lhs.pSuffix > rhs.pSuffix;
    } else
    {
        return lhs.pPrefix > rhs.pPrefix;
    }
}

template<typename PrefixT, typename SuffixT>
friend bool T6VectorIterator<PrefixT, SuffixT>::operator<=(const T6VectorIterator& lhs, const T6VectorIterator& rhs)
{
    if (lhs.pPrefix == rhs.pPrefix)
    {
        return lhs.pSuffix <= rhs.pSuffix;
    } else
    {
        return lhs.pPrefix <= rhs.pPrefix;
    }
}

template<typename PrefixT, typename SuffixT>
friend bool T6VectorIterator<PrefixT, SuffixT>::operator>=(const T6VectorIterator& lhs, const T6VectorIterator& rhs)
{
    if (lhs.pPrefix == rhs.pPrefix)
    {
        return lhs.pSuffix >= rhs.pSuffix;
    } else
    {
        return lhs.pPrefix >= rhs.pPrefix;
    }
}

template<typename PrefixT, typename SuffixT>
T6VectorIterator& T6VectorIterator<PrefixT, SuffixT>::operator+=(size_t i)
{
/*
    for (size_t ii = 0 ; ii < i; ii++)
    {
        this->++;
    }

*/
    size_t firstDiffToEnd = pPrefix->end() - pSuffix;

    if (i < firstDiffToEnd)
    {
        pSuffix += i;
    } else
    {
        i -= firstDiffToEnd;
        while (true)
        {
            pPrefix++;
            size_t diffToEnd = pPrefix->size();
            if (i < diffToEnd)
            {
                pSuffix += i;
                break;
            }
            i -= diffToEnd;
        }
    }
}

template<typename PrefixT, typename SuffixT>
friend T6VectorIterator T6VectorIterator<PrefixT, SuffixT>::operator+(const T6VectorIterator& iterator, size_t i)
{
    T6VectorIterator r(iterator);
    r+=i;
    return r;
}


template<typename PrefixT, typename SuffixT>
friend T6VectorIterator T6VectorIterator<PrefixT, SuffixT>::operator+(size_type i, const T6VectorIterator& iterator)
{
    T6VectorIterator r(iterator);
    r+=i;
    return r;
}


template<typename PrefixT, typename SuffixT>
T6VectorIterator& T6VectorIterator<PrefixT, SuffixT>::operator-=(size_t i)
{
/*
    for (size_t ii = 0 ; ii < i; ii++)
    {
        this->--;
    }

*/
    size_t firstDiffToBegin = pSuffix - pPrefix->begin();

    if (i < firstDiffToBegin)
    {
        pSuffix -= i;
    } else
    {
        i -= firstDiffToBegin;
        while (true)
        {
            pPrefix--;
            size_t diffToBegin = pPrefix->size();
            if (i < diffToBegin)
            {
                pSuffix -= i;
                break;
            }
            i -= diffToBegin;
        }
    }
} 


template<typename PrefixT, typename SuffixT>
friend T6VectorIterator T6VectorIterator<PrefixT, SuffixT>::operator-(const T6VectorIterator& iterator, size_t i)
{
    T6VectorIterator r(iterator);
    r-=i;
    return r;
}

template<typename PrefixT, typename SuffixT>
friend size_t T6VectorIterator<PrefixT, SuffixT>::operator-(T6VectorIterator to, T6VectorIterator from)
{
    /*
    size_t i = 0;
    while (from != to)
    {
        from++;
        i++;
    }
    return i;
    */ 

    size_t r = 0;

    T6VectorIterator newFrom(from);


    while (newFrom.pPrefix != to.pPrefix)
    {
        r += newFrom.pPrefix.size();
        newFrom.pPrefix++;
    }
    r -= from.pPrefix->end() - from.pSuffix;
    r += to.pSuffix - to.pPrefix->begin();

    return r;
}

template<typename PrefixT, typename SuffixT>
T6VectorIterator& T6VectorIterator<PrefixT, SuffixT>::operator[](size_t i)
{
    SuffixT* suf;
    for (auto pref = data.begin() ; true ; pref++)
    {
        if (pref->size() > i)
        {
            suf = pref->begin() + i;
            break;
        } else
        {
            i -= pref->size();
        }
    }

    return T6VectorIterator(pref, suf);    
}






template<typename PrefixT, typename SuffixT>
void T6VectorIterator<PrefixT, SuffixT>::computeMultiplicatorForPrefix()
{
    multiplicatorForPrefix = pow(2,std::numeric_limits<SuffixT>::max());

    // security
    if (SSP->T6_nbLoci > pow(2,std::numeric_limits<PrefixT>::max()) * multiplicatorForPrefix)
    {
        std::cout << "Error detected in void 'T6VectorIterator<PrefixT, SuffixT>::computeMultiplicatorForPrefix()'. There are " << SSP->T6_nbLoci << " but the max value of the prefixType and suffixType are " << pow(2,std::numeric_limits<PrefixT>::max() << " and "<< pow(2,std::numeric_limits<SuffixT>::max() << ", respectively. Their product is therefore " << pow(2,std::numeric_limits<SuffixT>::max() * pow(2,std::numeric_limits<PrefixT>::max() << " which is lower than the number of T6 loci.\n";
        abort();
    }
    if (std::numeric_limits<SuffixT>::max() > std::numeric_limits<PrefixT>::max())
    {
        std::cout << "Error detected in void 'T6VectorIterator<PrefixT, SuffixT>::computeMultiplicatorForPrefix()'. The suffix type has a higher maxx value than the rpefix type. I don't think anyone would want that. Also, in the current setting, it would mess up with 'multiplicatorForPrefix' although that would be easy to fix.\n";
        abort();   
    }
}


