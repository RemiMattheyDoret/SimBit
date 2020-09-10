
template<typename T>
class ReversibleDeque
{
private:
	std::deque<T> data;
	bool isEndAtEnd;

public:
	ReversibleDeque():isEndAtEnd(true){}
	ReversibleDeque(const ReversibleDeque& other)
	:data(other.data), isEndAtEnd(other.isEndAtEnd)
	{}

	ReversibleDeque(const ReversibleDeque&& other)
	:data(std::move(other.data)), isEndAtEnd(std::move(other.isEndAtEnd))
	{}

	ReversibleDeque& operator=(const ReversibleDeque& other)
	{
		data = other.data;
		isEndAtEnd = other.isEndAtEnd;
		return *this;
	}

	ReversibleDeque& operator=(const ReversibleDeque&& other)
	{
		data = std::move(other.data);
		isEndAtEnd = std::move(other.isEndAtEnd);
		return *this;
	}

	
	void push_back(const T& elem)
	{
		if (isEndAtEnd)
		{
			data.push_back(elem);
		} else
		{
			data.push_front(elem);
		}
	}

	void push_front(const T& elem)
	{
		if (isEndAtEnd)
		{
			data.push_front(elem);
		} else
		{
			data.push_back(elem);
		}
	}

	template<typename INT>
	T& operator[](const INT i)
	{
		if (isEndAtEnd)
		{
			return data[i];
		} else
		{
			return data[data.size()-i-1];
		}	
	}

	void reverse()
	{
		isEndAtEnd = !isEndAtEnd;
	}

	size_t size(){return data.size();} const

	void swap(ReversibleDeque& other)
	{
		this->data.swap(other.data);
	}

	void clear()
	{
		data.clear();
	}

	T& back() 
	{
		if (isEndAtEnd)
		{
			return data.back();
		} else
		{
			return data.front();
		}	
	}

	T& front() 
	{
		if (isEndAtEnd)
		{
			return data.front();
		} else
		{
			return data.back();
		}	
	}
};

