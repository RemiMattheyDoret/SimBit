
template<typename PrefixT, typename SuffixT>
class T6VectorIterator
{
private:
	PrefixT* pPrefix;
	SuffixT* pSuffix;
    static size_t multiplicatorForPrefix; // pow(2,std::numeric_limits<SuffixT>::max()) 

public:
    

	// general iterator
	T6VectorIterator(const T6VectorIterator&);
    ~T6VectorIterator();
    T6VectorIterator& operator=(const T6VectorIterator&);
    T6VectorIterator& operator++(); //prefix increment
    reference operator*() const;
    friend void swap(T6VectorIterator& lhs, T6VectorIterator& rhs); //C++11 I think


    // bidrectional
    T6VectorIterator& operator--(); //prefix decrement
    T6VectorIterator operator--(int); //postfix decrement


    //random_access_T6VectorIterator
    friend bool operator<(const T6VectorIterator&, const T6VectorIterator&);
    friend bool operator>(const T6VectorIterator&, const T6VectorIterator&);
    friend bool operator<=(const T6VectorIterator&, const T6VectorIterator&);
    friend bool operator>=(const T6VectorIterator&, const T6VectorIterator&);

    T6VectorIterator& operator+=(size_type);
    friend T6VectorIterator operator+(const T6VectorIterator&, size_type);
    friend T6VectorIterator operator+(size_type, const T6VectorIterator&);
    T6VectorIterator& operator-=(size_type);  
    friend T6VectorIterator operator-(const T6VectorIterator&, size_type);
    friend difference_type operator-(T6VectorIterator, T6VectorIterator);

    static void computeMultiplicatorForPrefix();
};