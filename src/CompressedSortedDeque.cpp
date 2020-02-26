
////////////
// Block ///
////////////


// All CSDBlock editing functions will compress. If you don't want compressing, don't use these functions!

void CompressedSortedDeque::CSDBlock::insert(CompressedSortedDeque::iterator it, unsigned int value)
{
    suffixs.insert(it.smallP, compress(value));
}



void CompressedSortedDeque::CSDBlock::push_back(unsigned int value)
{
    suffixs.push_back(compress(value));
}



void CompressedSortedDeque::CSDBlock::insert(std::vector<unsigned short>::iterator it, unsigned int value)
{
    suffixs.insert(it,compress(value));
}




// Constructors
CompressedSortedDeque::CSDBlock::CSDBlock(unsigned int blockIndex)
:adder((unsigned int) blockIndex * (unsigned int) std::numeric_limits<unsigned short>::max())
{
    suffixs.reserve(1);
}

CompressedSortedDeque::CSDBlock::CSDBlock(const CSDBlock& other)
:adder(other.adder), suffixs(other.suffixs)
{
    suffixs.reserve(1);
}

CompressedSortedDeque::CSDBlock::CSDBlock(CSDBlock&& other)
:adder(other.adder), suffixs(other.suffixs)
{
    suffixs.reserve(1);
}

CompressedSortedDeque::CSDBlock& CompressedSortedDeque::CSDBlock::operator=(const CSDBlock& other)
{ 
    adder=other.adder;
    suffixs = other.suffixs;
    return *this;
}



// Other functions

std::vector<unsigned short>::iterator CompressedSortedDeque::CSDBlock::begin()
{
    return suffixs.begin();
}

std::vector<unsigned short>::iterator CompressedSortedDeque::CSDBlock::end()
{
    return suffixs.end();
}

size_t CompressedSortedDeque::CSDBlock::size()
{
    return suffixs.size();
}

unsigned short CompressedSortedDeque::CSDBlock::compress(unsigned int value)
{
    return (unsigned short)(value - adder);
}









///////////////////
// Constructors ///
///////////////////

CompressedSortedDeque::CompressedSortedDeque(){}

CompressedSortedDeque::CompressedSortedDeque(unsigned int pastMaxValue)
// :endIterator(data.begin(), std::vector<unsigned short>::iterator(), data.end()) // That looks ugly but it is because bigPEnd cannot have a default initializer, so I have to initalize it // data(std::vector<CSDBlock>()), 
{
    unsigned short nbBlocks = getBlockForValue(pastMaxValue) + 1;
       
    for (unsigned short blockIndex = 0 ; blockIndex < nbBlocks ; ++blockIndex)
    {
        data.push_back(CSDBlock(blockIndex));
    }
}

CompressedSortedDeque::CompressedSortedDeque(std::vector<unsigned int> in, unsigned int pastMaxValue)
: CompressedSortedDeque(pastMaxValue) 
{
    for (auto& i : in)
    {
        this->push_back(i);
    }
}

CompressedSortedDeque& CompressedSortedDeque::operator=(const CompressedSortedDeque& other)
{
    data = other.data;
    return *this;
}




//////////////////////////
// Get iterator Methods///
/////////////////////////

CompressedSortedDeque::iterator CompressedSortedDeque::begin()
{
    // begin is the position of the first mutation that might not be in the first block
    std::vector<CSDBlock>::iterator bigP;

    bool isEmpty = true;
    // Search beginning
    for (auto it = data.begin() ; it < data.end() ; ++it)
    {
        if (it->size() != 0)
        {
            bigP = it;
            isEmpty = false;
            break;
        }
    }

    // return
    if (isEmpty)
    {
        return end();
    } else
    {
        return iterator(bigP, bigP->begin(), data.begin(), data.end());
    }
}

CompressedSortedDeque::iterator CompressedSortedDeque::end()
{
    return iterator(data.end(), data.back().end(), data.begin(), data.end());   
}


///////////////////
// Edit Methods///
///////////////////



void CompressedSortedDeque::insert(CompressedSortedDeque::iterator& it, unsigned int value)
{
    it.backwardIfWrongBlock(data.begin(), value);
    it.bigP->insert(it.smallP,value); // compressed in here ('it.bigP->suffixs.insert' would not compress)
}   

CompressedSortedDeque::iterator CompressedSortedDeque::insert(unsigned int value)
{
    CompressedSortedDeque::iterator it = CompressedSortedDeque::lower_bound_noForward(value);
    this->insert(it, value);

    return it;
}

void CompressedSortedDeque::push_back(unsigned int value)
{
    size_t blockIndex = getBlockForValue(value);
    data[blockIndex].push_back(value); // compress in here
}




void CompressedSortedDeque::erase(CompressedSortedDeque::iterator& it)
{
    it.bigP->suffixs.erase(it.smallP);
}

template<typename INT>
void CompressedSortedDeque::extend(INT locusFrom, INT locusTo, CompressedSortedDeque& source)
{
   /* std::cout << "~~~~~~\n";
    std::cout << "locusFrom = "<< locusFrom <<"\n";
    std::cout << "locusTo = "<< locusTo <<"\n";
    if (this->size())
    {
        std::cout << "#\n*(this->begin()) = " << *(this->begin()) << "\n";
        std::cout << "this->back() = " << this->back() << "\n#\n";
    }
    std::cout << "~~~~~~\n";*/
    assert(locusFrom < locusTo);
    
    /*std::cout << "source.size() = "<< source.size() <<"\n";
    source.assertOrdering();
    std::cout << "this->size() = "<< this->size() <<"\n";
    this->assertOrdering();*/
    


    auto itFrom = source.lower_bound(locusFrom);
    if (itFrom != source.end() && *itFrom < locusTo )
    {
        auto itTo = source.lower_bound_from(locusTo, itFrom);

        /*std::cout << "####\n";
        if (itFrom != source.end())
        {
            std::cout << "*itFrom = " << *itFrom << "\n";
        }
        if (itTo != source.end())
        {
            std::cout << "*itTo = " << *itTo << "\n";
        }
        std::cout << "####\n";*/

        if (itFrom < itTo)
            extend(itFrom, itTo, source);
            

        /*if (this->size())
        {
            std::cout << "----\n*(this->begin()) = " << *(this->begin()) << "\n";
            std::cout << "this->back() = " << this->back() << "\n";
            std::cout << "----\n";
        }*/
        //std::cout << "second this->size() = "<< this->size() <<"\n";
        //this->assertOrdering();
    }
        
    
        
    
    //std::cout << "this->data.size() = " << this->data.size() << "\n";
    
}

void CompressedSortedDeque::extend(CompressedSortedDeque::iterator& from, CompressedSortedDeque::iterator& to, CompressedSortedDeque& source)
{
    // Copy mutations in first block if the block is not complete

    

    //std::cout << "aaaaaa\n";
    // this->assertOrdering();
    {
        size_t blockIndex = from.bigP - source.data.begin(); 
        if (to.bigP == from.bigP)
        {
            /*std::cout << "*from = " << *from << "\n";
            if (data[blockIndex].suffixs.size() && from.bigP->suffixs.size())
            {
                auto specialIt = iterator(data.begin()+blockIndex, data[blockIndex].suffixs.end()-1, data.begin(), data.end());
                std::cout << "*specialIt = " << *specialIt << "\n";
                assert(*specialIt < *from);           
            }*/
            data[blockIndex].suffixs.insert(data[blockIndex].end(), from.smallP, to.smallP);
            return; // No need to go any further
        } else
        {
            /*std::cout << "*from = " << *from << "\n";
            if (data[blockIndex].suffixs.size() && from.bigP->suffixs.size())
            {
                auto specialIt = iterator(data.begin()+blockIndex, data[blockIndex].suffixs.end()-1, data.begin(), data.end());
                std::cout << "*specialIt = " << *specialIt << "\n";
                assert(*iterator(data.begin()+blockIndex, data[blockIndex].suffixs.end()-1, data.begin(), data.end()) < *from);
            }*/
            data[blockIndex].suffixs.insert(data[blockIndex].end(), from.smallP, from.bigP->end());
            ++from.bigP;
            if (from.bigP == from.bigPEnd) return;
            from.smallP = from.bigP->begin();
        }
    }
    //std::cout << "bbbbbb\n";
    //this->assertOrdering();

    // Copy whole CSDBlock
    {
        if (from.bigP < to.bigP)
        {
            assert(to.bigP <= to.bigPEnd);
            assert(from.bigP < from.bigPEnd);
            
            //size_t blockIndex = from.bigP - source.data.begin(); 
            //auto BlockIt = data.begin() + blockIndex;
            
            //if (BlockIt->suffixs.size());
            {
                size_t i = from.bigP - from.bigPBegin;
                size_t i_end = to.bigP - to.bigPBegin;
                while (i < i_end)
                {
                    data[i] = source.data[i];
                    ++i;
                }
            }
                
            from.bigP = to.bigP;
            if (from.bigP == from.bigPEnd) return;
            from.smallP = from.bigP->begin();
        }
    }
    
    //std::cout << "cccccc\n";      
    //this->assertOrdering();

    // Copy mutation in last block
    {
        if (from.smallP < to.smallP)
        {
            size_t blockIndex = from.bigP - source.data.begin();
            data[blockIndex].suffixs.insert(data[blockIndex].end(), from.smallP, to.smallP);
        }
    }
    //std::cout << "dddddd\n";
    //this->assertOrdering();

}



void CompressedSortedDeque::clear()
{
    for (auto& block : data)
    {
        block.suffixs.clear(); // that will clear the extra block for nothing but whatever
    }
}





/////////////////////
// Search Methods ///
/////////////////////


size_t CompressedSortedDeque::size()
{
    // begin is the position of the first mutation that might not be in the first block
    size_t r = 0;
    for (auto& da : data)
        r += da.suffixs.size();
    return r;
}

unsigned CompressedSortedDeque::back()
{
    assert(!this->isEmpty());
    return *(end()-1);
    /*for (int blockIndex = data.size()-1 ; blockIndex >= 0 ; --blockIndex)
    {
        if (data[blockIndex].suffixs.size())
        {
            return *(iterator(data.begin() + blockIndex, data[blockIndex].suffixs.end()-1, data.begin(), data.end()));
        }
    }

    std::cout << "Internal error in function unsigned CompressedSortedDeque::back(). This function should not have been called as it appears that the container is empty\n";*/
    abort();
}


bool CompressedSortedDeque::isEmpty()
{
    for ( auto& da : data )
    {
        if(da.size())
        {
            return false;
        }
    }
    return true;
}


void CompressedSortedDeque::swap(CompressedSortedDeque& other)
{
    this->data.swap(other.data);
}

CompressedSortedDeque::iterator CompressedSortedDeque::iterator::lower_bound_FromThis(unsigned int value)
{
    unsigned blockIndex = CompressedSortedDeque::getBlockForValue(value);
    auto newBigP = bigPBegin + blockIndex;
    std::vector<unsigned short>::iterator fromSmallP;
    if (this->bigP == newBigP)
    {
        fromSmallP = this->smallP;
    } else
    {
        fromSmallP = newBigP->begin();
    }
    auto newSmallP = std::lower_bound(fromSmallP, newBigP->end(), newBigP->compress(value));
    iterator r = iterator(
        newBigP,
        newSmallP,
        bigPBegin,
        bigPEnd
    );
 
    r.updateIfReachedEndAndThereIsMore();
    return r;
}

CompressedSortedDeque::iterator CompressedSortedDeque::lower_bound(unsigned int value)
{
    unsigned int blockIndex = getBlockForValue(value);
    //std::cout << "blockIndex for " << value << " = "<< blockIndex << "\n";
    auto bigPtoSearchIn = data.begin() + blockIndex;
    iterator r = iterator(
        bigPtoSearchIn,
        std::lower_bound(bigPtoSearchIn->suffixs.begin(), bigPtoSearchIn->suffixs.end(), bigPtoSearchIn->compress(value)),
        data.begin(),
        data.end()
    );
    /*
    assert(r.bigP == bigPtoSearchIn);
    if (r.smallP == r.bigP->suffixs.end())
    {
        std::cout << "In lower_bound, it reached end of block\n";
    } else
    {
        std::cout << "In lower_bound, it did not reach end of block\n";
    }*/

    /*std::cout << "test at line 403\n";
    if (!this->isEmpty() && r == end())
    {
        std::cout << "this->size() = " << this->size() << "\n";
        for (auto it = begin(); it < end() ; ++it)
        {
            std::cout << *it << " ";
        }
        std::cout << "\n";
        std::cout << "value = " << value << "\n";
        std::cout << "back() = " << back() << "\n";
    }
    std::cout << "tested at line 403\n";*/

    r.updateIfReachedEndAndThereIsMore();
    //bool b = r == end();
    //std::cout << "Does lower_bound return the end iterator " << b << "\n";
    return r;
}

CompressedSortedDeque::iterator CompressedSortedDeque::lower_bound_from(unsigned int value, unsigned& blockIndexFrom, unsigned& fromInBlock)
{
    unsigned int blockIndex = getBlockForValue(value);
    auto bigPtoSearchIn = data.begin() + blockIndex;
    
    std::vector<unsigned short>::iterator fromIt;
    if (blockIndex == blockIndexFrom)
    {
        fromIt = bigPtoSearchIn->suffixs.begin() + fromInBlock;
    } else
    {
        fromIt = bigPtoSearchIn->suffixs.begin();
    }
    iterator r = iterator(
        bigPtoSearchIn,
        std::lower_bound(fromIt, bigPtoSearchIn->suffixs.end(), bigPtoSearchIn->compress(value)),
        data.begin(), 
        data.end()
    );
    r.updateIfReachedEndAndThereIsMore();

    fromInBlock = r.getDistanceToBeginningOfBlock();
    blockIndexFrom = blockIndex;

    return r;   
}

CompressedSortedDeque::iterator CompressedSortedDeque::lower_bound_from(unsigned int value, const CompressedSortedDeque::iterator& itFrom)
{
    unsigned int blockIndex = getBlockForValue(value);
    std::vector<unsigned short>::iterator smallPFrom;
    if (itFrom.bigP == data.begin() + blockIndex)
    {
        smallPFrom = itFrom.smallP;
    } else
    {
        smallPFrom = data[blockIndex].suffixs.begin();
    }
    iterator r = iterator(
        data.begin() + blockIndex,
        std::lower_bound(smallPFrom, data[blockIndex].suffixs.end(), data[blockIndex].compress(value)),
        data.begin(), 
        data.end()
    );
    r.updateIfReachedEndAndThereIsMore();
    return r;
}

CompressedSortedDeque::iterator CompressedSortedDeque::lower_bound_noForward(unsigned int value)
{
    unsigned int blockIndex = getBlockForValue(value);
    auto bigPtoSearchIn = data.begin() + blockIndex;
    iterator r = iterator(
        bigPtoSearchIn,
        std::upper_bound(bigPtoSearchIn->suffixs.begin(), bigPtoSearchIn->suffixs.end(), bigPtoSearchIn->compress(value)),
        data.begin(), 
        data.end()
    );
    return r;
}



////////////////////
// Other Methods ///
////////////////////

unsigned short CompressedSortedDeque::getBlockForValue(const unsigned int value)
{
    unsigned int m = (unsigned int) std::numeric_limits<unsigned short>::max();
    return (value + m) / m - 1;
}

std::vector<unsigned> CompressedSortedDeque::toVector()
{
    std::vector<unsigned> r;
    r.reserve(this->size());
    for (auto it = begin() ; it != end() ; ++it)
    {
        r.push_back(*it);
    }
    return r;
}


void CompressedSortedDeque::assertOrdering(std::string message)
{
    if (this->size() > 1)
    {
        auto it = this->begin();
        auto end = this->end();
        unsigned prev = *it;
        ++it;
        for (; it < end ; ++it)
        {
            unsigned curr = *it;
            if (prev >= curr)
            {
                std::cout << "Internal error found in assertOrdering. " <<message <<'\n';
                abort();
            }
                
            prev = curr;
        }
    }
}



///////////////////////
///// Iterator ///////
//////////////////////

CompressedSortedDeque::iterator::iterator(){}

CompressedSortedDeque::iterator::iterator(std::vector<CSDBlock>::iterator b, std::vector<unsigned short>::iterator s, std::vector<CSDBlock>::iterator f, std::vector<CSDBlock>::iterator l)
: bigP(b), smallP(s), bigPBegin(f), bigPEnd(l)
{}

CompressedSortedDeque::iterator::iterator(const CompressedSortedDeque::iterator& other)
:bigP(other.bigP), smallP(other.smallP), bigPBegin(other.bigPBegin), bigPEnd(other.bigPEnd)
{}


CompressedSortedDeque::iterator CompressedSortedDeque::iterator::operator=(const CompressedSortedDeque::iterator& other)
{
    this->bigP = other.bigP;
    this->smallP = other.smallP;
    this->bigPBegin = other.bigPBegin;
    this->bigPEnd = other.bigPEnd;
    return *this;
}



CompressedSortedDeque::iterator CompressedSortedDeque::iterator::operator++() //prefix increment
{
    smallP++;
    this->updateIfReachedEndAndThereIsMore();

    return *this;
}


CompressedSortedDeque::iterator CompressedSortedDeque::iterator::operator++(int) //postfix increment
{
    auto r = *this;
    smallP++;
    this->updateIfReachedEndAndThereIsMore();        

    return r;
}



unsigned CompressedSortedDeque::iterator::operator*() const
{
    //std::cout << "dereference "<< ((unsigned int) bigP->prefix) + (unsigned int)*smallP <<" with prefix " << bigP->prefix << " and short " << *smallP << "\n";
    //std::cout << "bigP->adder = " << bigP->adder << "\n";
    //std::cout << "*smallP = " << *smallP << "\n" << std::flush;
    return (unsigned int) bigP->adder + (unsigned int) (*smallP);
}



// bidrectional

CompressedSortedDeque::iterator CompressedSortedDeque::iterator::operator--() //prefix decrement
{
    if ( smallP == bigP->suffixs.begin() )
    {
        --bigP;
        while (bigP->suffixs.size() == 0 && bigP != bigPBegin)
        {
            --bigP;
        }
        if (bigP->suffixs.size() == 0)
        {
            smallP = bigP->suffixs.begin();
        } else
        {
            smallP = bigP->suffixs.end() - 1;
        }
            
    } else
    {
        --smallP;
    }

    return *this;
}



CompressedSortedDeque::iterator CompressedSortedDeque::iterator::operator--(int) //postfix decrement
{
    auto r = *this;
    this->operator--();

    return r;
}



CompressedSortedDeque::iterator CompressedSortedDeque::iterator::operator+=(size_t i)
{

//   [ ][ ][ ][ ]   [ ][ ][ ][ ][ ]
//       ^  + 5            ~

    size_t diffToEnd = bigP->suffixs.end() - smallP;

    if (i <= diffToEnd)
    {
        smallP += i;
    } else
    {
        i -= diffToEnd;
        while (true)
        {
            bigP++;
            diffToEnd = bigP->suffixs.size();
            if (i <= diffToEnd)
            {
                smallP = bigP->suffixs.begin() + i;
                break;
            }
            i -= diffToEnd;
        }
    }

    this->updateIfReachedEndAndThereIsMore();
    return *this;
}


CompressedSortedDeque::iterator operator+(const CompressedSortedDeque::iterator& it, size_t i)
{
    CompressedSortedDeque::iterator r(it);
    r+=i;
    return r;
}



CompressedSortedDeque::iterator operator+(size_t i, const CompressedSortedDeque::iterator& it)
{
    CompressedSortedDeque::iterator r(it);
    r+=i;
    return r;
}



CompressedSortedDeque::iterator CompressedSortedDeque::iterator::operator-=(size_t i)
{

//   [ ][ ][ ][ ]  [ ][ ][ ]   [ ][ ][ ][ ][ ]
//       ~                      ^  - 1
    if (i == 0)
        return *this;

    if (bigP == bigPEnd)
    {
        --bigP;
    }
    size_t diffToBegin = smallP - bigP->suffixs.begin(); // diff to object past begin looking backward!

    if (i <= diffToBegin)
    {
        smallP -= i;
    } else
    {
        i -= diffToBegin;
        while (true)
        {
            bigP--;
            diffToBegin = bigP->suffixs.size();
            if (i <= diffToBegin)
            {
                smallP = bigP->suffixs.end() - i;
                break;
            }
            i -= diffToBegin; 
        }
    }
    return *this;
} 



CompressedSortedDeque::iterator operator-(const CompressedSortedDeque::iterator& iterator, size_t i)
{
    CompressedSortedDeque::iterator r(iterator);
    r-=i;
    return r;
}


size_t operator-(CompressedSortedDeque::iterator& to, CompressedSortedDeque::iterator& from)
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

    CompressedSortedDeque::iterator newFrom(from);
    while (newFrom.bigP != to.bigP)
    {
        r += newFrom.bigP->suffixs.size();
        newFrom.bigP++;
    }

    r -= from.bigP->suffixs.end() - from.smallP;
    r += to.smallP - to.bigP->suffixs.begin();

    return r;
}


CompressedSortedDeque::iterator CompressedSortedDeque::iterator::operator[](size_t i)
{
    return *this + i;
}



void swap(CompressedSortedDeque::iterator& lhs, CompressedSortedDeque::iterator& rhs)
{
    auto tmp = rhs;
    rhs = lhs;
    lhs = tmp;
} //C++11 I think

//random_access_CompressedSortedDeque::iterator

bool operator<(const CompressedSortedDeque::iterator& lhs, const CompressedSortedDeque::iterator& rhs)
{
    if (lhs.bigP == rhs.bigP)
    {
        return lhs.smallP < rhs.smallP;
    } else
    {
        return lhs.bigP < rhs.bigP;
    }
}


bool operator>(const CompressedSortedDeque::iterator& lhs, const CompressedSortedDeque::iterator& rhs)
{
    if (lhs.bigP == rhs.bigP)
    {
        return lhs.smallP > rhs.smallP;
    } else
    {
        return lhs.bigP > rhs.bigP;
    }
}


bool operator<=(const CompressedSortedDeque::iterator& lhs, const CompressedSortedDeque::iterator& rhs)
{
    if (lhs.bigP == rhs.bigP)
    {
        return lhs.smallP <= rhs.smallP;
    } else
    {
        return lhs.bigP <= rhs.bigP;
    }
}


bool operator>=(const CompressedSortedDeque::iterator& lhs, const CompressedSortedDeque::iterator& rhs)
{
    if (lhs.bigP == rhs.bigP)
    {
        return lhs.smallP >= rhs.smallP;
    } else
    {
        return lhs.bigP >= rhs.bigP;
    }
}

bool operator==(const CompressedSortedDeque::iterator& lhs, const CompressedSortedDeque::iterator& rhs)
{
    return lhs.smallP == rhs.smallP;
}

bool operator!=(const CompressedSortedDeque::iterator& lhs, const CompressedSortedDeque::iterator& rhs)
{
    return lhs.smallP != rhs.smallP;
}

void CompressedSortedDeque::iterator::updateIfReachedEndAndThereIsMore()
{
    while (smallP == bigP->suffixs.end())
    {
        //std::cout << "In updateIfReachedEndAndThereIsMore, it reached end of block\n";
        ++bigP;
        if (bigP == bigPEnd)
        {
            //std::cout << "In updateIfReachedEndAndThereIsMore, it reached bigPend\n";
            break;
        }
        smallP = bigP->suffixs.begin();
    }
}


void CompressedSortedDeque::iterator::backwardIfWrongBlock(std::vector<CSDBlock>::iterator firstBlockIt, const unsigned int value)
{
    if (bigP->suffixs.begin() == smallP)
    {
        std::vector<CSDBlock>::iterator rightBlockIt = firstBlockIt + CompressedSortedDeque::getBlockForValue(value);
        if (rightBlockIt != bigP)
        {
            bigP = rightBlockIt;
            smallP = bigP->suffixs.end();
        }
    }
}

size_t CompressedSortedDeque::iterator::getDistanceToBeginningOfBlock()
{
    
    if (bigP != bigPEnd)
    {
        return this->smallP - this->bigP->suffixs.begin();
    } else
    {
        return 0; // return whatever in fact!
    }
}    




