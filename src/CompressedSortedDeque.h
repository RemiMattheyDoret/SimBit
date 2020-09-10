

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
		void backwardIfWrongBlock(std::vector<CSDBlock>::iterator firstBlockIt, const uint32_t value);

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

		iterator operator+=(uint32_t i);
		iterator operator-=(uint32_t i);
		iterator operator[](uint32_t i);
		friend iterator operator+(const CompressedSortedDeque::iterator& iterator, uint32_t i);
		friend iterator operator+(uint32_t i, const CompressedSortedDeque::iterator& iterator);
		friend iterator operator-(const CompressedSortedDeque::iterator& iterator, uint32_t i);
		friend uint32_t operator-(CompressedSortedDeque::iterator& to, CompressedSortedDeque::iterator& from);
		
		friend void swap(CompressedSortedDeque::iterator& lhs, CompressedSortedDeque::iterator& rhs);
		friend bool operator<(const CompressedSortedDeque::iterator& lhs, const CompressedSortedDeque::iterator& rhs);
		friend bool operator>(const CompressedSortedDeque::iterator& lhs, const CompressedSortedDeque::iterator& rhs);
		friend bool operator<=(const CompressedSortedDeque::iterator& lhs, const CompressedSortedDeque::iterator& rhs);
		friend bool operator>=(const CompressedSortedDeque::iterator& lhs, const CompressedSortedDeque::iterator& rhs);
		friend bool operator==(const CompressedSortedDeque::iterator& lhs, const CompressedSortedDeque::iterator& rhs);
		friend bool operator!=(const CompressedSortedDeque::iterator& lhs, const CompressedSortedDeque::iterator& rhs);

		uint32_t getDistanceToBeginningOfBlock();
		iterator lower_bound_FromThis(uint32_t value);
	};


	//////////////////
	//// CSDBlock ////
	//////////////////

	class CSDBlock
	{
		friend CompressedSortedDeque;
	private:
		uint32_t adder;
		unsigned short compress(uint32_t value);

	public:
		std::vector<unsigned short> suffixs;

		CSDBlock(uint32_t blockIndex);
		CSDBlock(const CSDBlock& other); // Used for operator=
		CSDBlock(CSDBlock&& other); // Used for operator=
		CSDBlock& operator=(const CSDBlock& other);

		std::vector<unsigned short>::iterator begin();
		std::vector<unsigned short>::iterator end();
		uint32_t size();

		// All CSDBlock editing functions will compress. If you don't want compressing, don't use these functions!

		void insert(CompressedSortedDeque::iterator it, uint32_t value);
		void push_back(uint32_t value);
		void insert(std::vector<unsigned short>::iterator it, uint32_t value);
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

	static unsigned short getBlockForValue(const uint32_t value);
	CompressedSortedDeque::iterator lower_bound_noForward(uint32_t value);
	void UpdateEndIterator(CompressedSortedDeque::iterator newPossibleEnd);

	void extend(CompressedSortedDeque::iterator& from, CompressedSortedDeque::iterator& to,  CompressedSortedDeque& source);

public:
	

	CompressedSortedDeque(uint32_t maxValue);
	CompressedSortedDeque(std::vector<uint32_t> in, uint32_t maxValue);
	CompressedSortedDeque(); // Does nothing

	std::vector<unsigned> toVector();

	CompressedSortedDeque& operator=(const CompressedSortedDeque& other);


	iterator begin();
	iterator end();
	uint32_t size();
	unsigned back();
	bool isEmpty();
	void swap(CompressedSortedDeque& other);

	// Any editing function will compress. The only exception is extend (who takes already compressed data)
	void insert(iterator& it, uint32_t value);
	void erase(iterator& it);
	iterator insert(uint32_t value);
	void push_back(uint32_t value);



	template<typename INT>
	void extend(INT locusFrom, INT locusTo,  CompressedSortedDeque& source);

	//CompressedSortedDeque::iterator lower_bound(CompressedSortedDeque::iterator from, CompressedSortedDeque::iterator to, uint32_t value);
	//CompressedSortedDeque::iterator upper_bound(CompressedSortedDeque::iterator from, CompressedSortedDeque::iterator to, uint32_t value);
	CompressedSortedDeque::iterator lower_bound(uint32_t value);
	CompressedSortedDeque::iterator lower_bound_from(uint32_t value, unsigned& blockIndexFrom, unsigned& fromInBlock);
	CompressedSortedDeque::iterator lower_bound_from(uint32_t value, const CompressedSortedDeque::iterator& itFrom);

	void clear();
	void assertOrdering(std::string message = "");
	
};

