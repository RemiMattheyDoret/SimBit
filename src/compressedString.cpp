class compressedString
{
	std::string s;
	std::string sep;
	uint32_t counter;
	std::string countingString;

	void addCountingString(bool exactCapacity = false)
	{
		if (counter != 0)
		{
			if (counter == 1)
			{
				if (exactCapacity)
					s.reserve(s.size() + countingString.size());
				s += countingString;
			} else if (counter * countingString.size() > countingString.size() + sep.size() * 3 + 1) // assuming one digit in counter here
			{
				if (exactCapacity) 
					s.reserve(s.size() + sep.size() * 3 + countingString.size() + std::to_string(counter).size());
				s += "R" + countingString + sep + std::to_string(counter) + sep;
			} else
			{
				s.reserve(s.size() + countingString.size() * counter);
				for (size_t c = 0 ; c < counter ; ++c)
					s += countingString;
			}
			counter = 0;	
			countingString = "";
		}
	}

public:
	template<typename INT>
	void reserve(INT i)
	{
		s.reserve(i);
	}

	compressedString(std::string separator = "_")
	: sep(separator), counter(0), countingString("")
	{}

	void add(std::string& i)
	{
		if (i != "")
		{
			if (i == countingString)
			{
				++counter;
			} else
			{
				addCountingString();
				countingString = i;
				counter = 1;
			}
		}
	}

	void add(std::string i)
	{
		if (i != "")
		{
			if (i == countingString)
			{
				++counter;
			} else
			{
				addCountingString();
				countingString = i;
				counter = 1;
			}
		}
	}

	std::string& toString()
	{
		addCountingString(true);
		return s;
	}

};