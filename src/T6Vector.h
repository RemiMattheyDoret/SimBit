
template<typename PrefixT, typename SuffixT>
class T6Vector
{
private:
	std::vector<T6VectorElement> data;
	PrefixT multiplicatorForPrefix;

public:
	class T6VectorIterator{};

	// general iterator
	T6Vector(const T6VectorIterator&);
	T6VectorIterator operator=(const T6VectorIterator&);
	size_t get(size_t i);
	void   set(size_t i, size_t value);
	void   push_back(size_t value);
	T6VectorIterator begin();
	T6VectorIterator end();

	void insert(size_t value);
	
	void push_back(size_t value);
	void 
    
    
    static void computeMultiplicatorForPrefix();
};