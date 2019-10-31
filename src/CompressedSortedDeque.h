

class CompressedSortedDeque
{

public:

	//////////////////
	//// iterator ////
	//////////////////

	class CSDBlock;
	class iterator{
	private:
		friend class CompressedSortedDeque;
		std::vector<CSDBlock>::iterator bigP;   // This must always point to the block in which smallP is contained
		std::vector<unsigned short>::iterator smallP; // end() is when bigP == bigPLast, the very end of data
		std::vector<CSDBlock>::iterator bigPBegin;
		std::vector<CSDBlock>::iterator bigPEnd;
		

		void updateIfReachedEndAndThereIsMore();
		void backwardIfWrongBlock(std::vector<CSDBlock>::iterator firstBlockIt, const unsigned int value);

	public:
		iterator();
		iterator(const iterator&);
		iterator(std::vector<CSDBlock>::iterator b, std::vector<unsigned short>::iterator s, std::vector<CSDBlock>::iterator f, std::vector<CSDBlock>::iterator l);
		
		iterator operator=(const iterator&);
		unsigned operator*() const;
		iterator operator++();
		iterator operator++(int);
		iterator operator--(); //prefix decrement
		iterator operator--(int); //postfix decrement

		iterator operator+=(size_t i);
		iterator operator-=(size_t i);
		iterator operator[](size_t i);
		friend iterator operator+(const CompressedSortedDeque::iterator& iterator, size_t i);
		friend iterator operator+(size_t i, const CompressedSortedDeque::iterator& iterator);
		friend iterator operator-(const CompressedSortedDeque::iterator& iterator, size_t i);
		friend size_t operator-(CompressedSortedDeque::iterator& to, CompressedSortedDeque::iterator& from);
		
		friend void swap(CompressedSortedDeque::iterator& lhs, CompressedSortedDeque::iterator& rhs);
		friend bool operator<(const CompressedSortedDeque::iterator& lhs, const CompressedSortedDeque::iterator& rhs);
		friend bool operator>(const CompressedSortedDeque::iterator& lhs, const CompressedSortedDeque::iterator& rhs);
		friend bool operator<=(const CompressedSortedDeque::iterator& lhs, const CompressedSortedDeque::iterator& rhs);
		friend bool operator>=(const CompressedSortedDeque::iterator& lhs, const CompressedSortedDeque::iterator& rhs);
		friend bool operator==(const CompressedSortedDeque::iterator& lhs, const CompressedSortedDeque::iterator& rhs);
		friend bool operator!=(const CompressedSortedDeque::iterator& lhs, const CompressedSortedDeque::iterator& rhs);

		size_t getDistanceToBeginningOfBlock();
		iterator lower_bound_FromThis(unsigned int value);
	};


	//////////////////
	//// CSDBlock ////
	//////////////////

	class CSDBlock
	{
		friend CompressedSortedDeque;
	private:
		unsigned int adder;
		unsigned short compress(unsigned int value);

	public:
		std::vector<unsigned short> suffixs;

		CSDBlock(unsigned int blockIndex);
		CSDBlock(const CSDBlock& other); // Used for operator=
		CSDBlock(CSDBlock&& other); // Used for operator=
		CSDBlock& operator=(const CSDBlock& other);

		std::vector<unsigned short>::iterator begin();
		std::vector<unsigned short>::iterator end();
		size_t size();

		// All CSDBlock editing functions will compress. If you don't want compressing, don't use these functions!

		void insert(CompressedSortedDeque::iterator it, unsigned int value);
		void push_back(unsigned int value);
		void insert(std::vector<unsigned short>::iterator it, unsigned int value);
		void insert(std::vector<unsigned short>::iterator thisFrom, std::vector<unsigned short>::iterator otherFrom, std::vector<unsigned short>::iterator otherTo);
	};


	

private:

	////////////////////
	//// Attributes ////
	////////////////////
	std::vector<CSDBlock> data;


	/////////////////
	//// Methods ////
	/////////////////

	static unsigned short getBlockForValue(const unsigned int value);
	CompressedSortedDeque::iterator lower_bound_noForward(unsigned int value);
	void UpdateEndIterator(CompressedSortedDeque::iterator newPossibleEnd);

	void extend(CompressedSortedDeque::iterator& from, CompressedSortedDeque::iterator& to,  CompressedSortedDeque& source);

public:
	

	CompressedSortedDeque(unsigned int maxValue);
	CompressedSortedDeque(std::vector<unsigned int> in, unsigned int maxValue);
	CompressedSortedDeque(); // Does nothing

	std::vector<unsigned> toVector();

	CompressedSortedDeque& operator=(const CompressedSortedDeque& other);


	iterator begin();
	iterator end();
	size_t size();
	unsigned back();
	bool isEmpty();
	void swap(CompressedSortedDeque& other);

	// Any editing function will compress. The only exception is extend (who takes already compressed data)
	void insert(iterator& it, unsigned int value);
	void erase(iterator& it);
	iterator insert(unsigned int value);
	void push_back(unsigned int value);



	template<typename INT>
	void extend(INT locusFrom, INT locusTo,  CompressedSortedDeque& source);

	//CompressedSortedDeque::iterator lower_bound(CompressedSortedDeque::iterator from, CompressedSortedDeque::iterator to, unsigned int value);
	//CompressedSortedDeque::iterator upper_bound(CompressedSortedDeque::iterator from, CompressedSortedDeque::iterator to, unsigned int value);
	CompressedSortedDeque::iterator lower_bound(unsigned int value);
	CompressedSortedDeque::iterator lower_bound_from(unsigned int value, unsigned& blockIndexFrom, unsigned& fromInBlock);
	CompressedSortedDeque::iterator lower_bound_from(unsigned int value, const CompressedSortedDeque::iterator& itFrom);

	void clear();
	void assertOrdering(std::string message = "");
	
};

