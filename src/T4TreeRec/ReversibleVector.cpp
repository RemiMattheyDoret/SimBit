

template<typename T> T4TreeRec::ReversibleVector<T>::ReversibleVector():isEndAtEnd(true){}
template<typename T> T4TreeRec::ReversibleVector<T>::ReversibleVector(const ReversibleVector<T>& other)
:data_back(other.data_back), data_front(other.data_front), isEndAtEnd(other.isEndAtEnd)
{}

template<typename T> T4TreeRec::ReversibleVector<T>::ReversibleVector(const ReversibleVector<T>&& other)
:data_back(std::move(other.data_back)), data_front(std::move(other.data_front)), isEndAtEnd(std::move(other.isEndAtEnd))
{}

template<typename T> T4TreeRec::ReversibleVector<T>& T4TreeRec::ReversibleVector<T>::operator=(const ReversibleVector<T>& other)
{
	data_back = other.data_back;
	data_front = other.data_front;
	isEndAtEnd = other.isEndAtEnd;
	return *this;
}

template<typename T> T4TreeRec::ReversibleVector<T>& T4TreeRec::ReversibleVector<T>::operator=(const ReversibleVector<T>&& other)
{
	data_back = std::move(other.data_back);
	data_front = std::move(other.data_front);
	isEndAtEnd = std::move(other.isEndAtEnd);
	return *this;
}

template<typename T> void T4TreeRec::ReversibleVector<T>::reserve(uint32_t n)
{
	if (isEndAtEnd)
	{
		data_back.reserve(n);
	} else
	{
		data_front.reserve(n);
	}
}

template<typename T> void T4TreeRec::ReversibleVector<T>::push_back(const T& elem)
{
	if (isEndAtEnd)
	{
		data_back.push_back(elem);
	} else
	{
		data_front.push_back(elem);
	}
}

template<typename T> void T4TreeRec::ReversibleVector<T>::push_front(const T& elem)
{
	if (isEndAtEnd)
	{
		data_front.push_back(elem);
	} else
	{
		data_back.push_back(elem);
	}
}

template<typename T>
template<typename int_type>
T& T4TreeRec::ReversibleVector<T>::operator[](const int_type i) 
{
	/*
		[.....|..]
		[..|.....]
	*/
	if (isEndAtEnd)
	{
		if (i >= data_front.size())
		{
			return data_back[i - data_front.size()];
		} else
		{
			return data_front[data_front.size() - i - 1];
		}
	} else
	{
		if (i >= data_back.size())
		{
			return data_front[i - data_back.size()];
		} else
		{
			return data_back[data_back.size() - i - 1];
		}
	}	
}

template<typename T>
template<typename int_type>
const T& T4TreeRec::ReversibleVector<T>::operator[](const int_type i) const 
{
	/*
		[.....|..]
		[..|.....]
	*/
	if (isEndAtEnd)
	{
		if (i >= data_front.size())
		{
			return data_back[i - data_front.size()];
		} else
		{
			return data_front[data_front.size() - i - 1];
		}
	} else
	{
		if (i >= data_back.size())
		{
			return data_front[i - data_back.size()];
		} else
		{
			return data_back[data_back.size() - i - 1];
		}
	}	
}

template<typename T> void T4TreeRec::ReversibleVector<T>::reverse()
{
	isEndAtEnd = !isEndAtEnd;
	shrink_to_fit();
}

template<typename T> void T4TreeRec::ReversibleVector<T>::shrink_to_fit()
{
	data_front.shrink_to_fit();
	data_back.shrink_to_fit();
}

template<typename T> uint32_t T4TreeRec::ReversibleVector<T>::size() const
{
	return data_front.size() + data_back.size();
} 

template<typename T> void T4TreeRec::ReversibleVector<T>::swap(ReversibleVector<T>& other)
{
	this->data_back.swap(other.data_back);
	this->data_front.swap(other.data_front);
	std::swap(this->isEndAtEnd, other.isEndAtEnd);
}

template<typename T> void T4TreeRec::ReversibleVector<T>::clear()
{
	data_back.clear();
	data_front.clear();
}

template<typename T> T& T4TreeRec::ReversibleVector<T>::back() 
{
	/*
		[.....|..]
		[..|.....]
	*/
	if (isEndAtEnd)
	{
		return data_back.back();
	} else
	{
		return data_front.back();
	}	
}

template<typename T> T& T4TreeRec::ReversibleVector<T>::front() 
{
	/*
		[.....|..]
		[..|.....]
	*/
	if (isEndAtEnd)
	{
		return data_front.back();
	} else
	{
		return data_back.back();
	}	
}

