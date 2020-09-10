
void GeneticMap::shrink_to_fit()
{
    fromLocusToLastChangeIndex.shrink_to_fit();
    fromT1LocusToLastChangeIndex.shrink_to_fit();
    fromT2LocusToLastChangeIndex.shrink_to_fit();
    fromT3LocusToLastChangeIndex.shrink_to_fit();
    fromT4LocusToLastChangeIndex.shrink_to_fit();
    fromT56ntrlLocusToLastChangeIndex.shrink_to_fit();
    fromT56selLocusToLastChangeIndex.shrink_to_fit();
    fromT56LocusToLastChangeIndex.shrink_to_fit();
    TallLocus_atChange.shrink_to_fit();
    whoWillIncreaseFromChange.shrink_to_fit();
    T1Locus_atChange.shrink_to_fit();
    T2Locus_atChange.shrink_to_fit();
    T3Locus_atChange.shrink_to_fit();
    T4Locus_atChange.shrink_to_fit();
    T56ntrlLocus_atChange.shrink_to_fit();
    T56selLocus_atChange.shrink_to_fit();
    T56Locus_atChange.shrink_to_fit();
}


unsigned char GeneticMap::typeFormat(std::string& Type, bool allowSubTypes) const
{
    /*
    std::cout << "Enters GeneticMap::typeFormat(std::string& Type, bool allowSubTypes) with type " << Type << " and allowSubTypes " << allowSubTypes << "\n";
    std::cout << "FT1Type = " << static_cast<unsigned>(FT1Type )<<"\n";
    std::cout << "FT2Type = " << static_cast<unsigned>(FT2Type )<<"\n";
    std::cout << "FT3Type = " << static_cast<unsigned>(FT3Type )<<"\n";
    std::cout << "FT4Type = " << static_cast<unsigned>(FT4Type )<<"\n";
    std::cout << "FT56Type = " << static_cast<unsigned>(FT56Type) <<"\n";
    std::cout << "FT56ntrlType = " << static_cast<unsigned>(FT56ntrlType)<<"\n";
    std::cout << "FT56selType = " << static_cast<unsigned>(FT56selType) <<"\n";
    */
    if (Type == "1" || Type == "T1" || Type == "t1")
    {
        return FT1Type;
    } else if (Type == "2" || Type == "T2" || Type == "t2")
    {
        return FT2Type;
    } else if (Type == "3" || Type == "T3" || Type == "t3")
    {
        return FT3Type;
    } else if (Type == "4" || Type == "T4" || Type == "t4")
    {
        return FT4Type;
    } else if (Type == "5" || Type == "T5" || Type == "t5")
    {
        return FT56Type;
    } else if (allowSubTypes)
    {
        if (Type == "T5ntrl")
        {
            return FT56ntrlType;
        } else if (Type == "T5sel")
        {
            return FT56selType;
        } else
        {
            std::cout << "Received unknown type of locus " << Type << ".\n";
            abort();
        }
    } else
    {
        std::cout << "Received unknown type of locus " << Type << ".\n";
        abort();
    }
}

std::string GeneticMap::typeFormat(unsigned char& c) const
{
    if (c == FT1Type)
    {
        return "T1";
    } else if (c == FT2Type)
    {
        return "T2";
    } else if (c == FT3Type)
    {
        return "T3";
    } else if (c == FT4Type)
    {
        return "T4";
    } else if (c == FT56Type)
    {
        return "T5";
    } else if (c == FT56ntrlType)
    {
        return "T5ntrl";
    } else if (c == FT56selType)
    {
        return "T5sel";
    } else
    {
        std::cout << "Received unknown type of locus " << c << ". Note this is very likely an internal error but you might still want to check your input!\n";
        abort();
    }
}

void GeneticMap::readLoci_SetNewChangeOfType(unsigned char& type, PreviousLocusType& previousLocus)
{
    previousLocus.setType(type);
    std::cout << "adding type " << (size_t)type << " at change_index " << whoWillIncreaseFromChange.size() << "\n";
    whoWillIncreaseFromChange.push_back(type);
    TallLocus_atChange.push_back(TotalNbLoci+1);
    T1Locus_atChange.push_back(T1_nbLoci);
    T2Locus_atChange.push_back(T2_nbLoci);
    T3Locus_atChange.push_back(T3_nbLoci);
    T4Locus_atChange.push_back(T4_nbLoci);
    T56Locus_atChange.push_back(T56_nbLoci);
    T56ntrlLocus_atChange.push_back(T56ntrl_nbLoci);
    T56selLocus_atChange.push_back(T56sel_nbLoci);

    /*if (type == FT1Type)
    {
        ++T1Locus_atChange.back();
    } else if (type == FT2Type)
    {
        ++T2Locus_atChange.back();
    } else if (type == FT3Type)
    {
        ++T3Locus_atChange.back();
    } else if (type == FT4Type)
    {
        ++T4Locus_atChange.back();
    } else if (type == FT56ntrlType)
    {
        ++T56ntrlLocus_atChange.back();
        ++T56Locus_atChange.back();
    } else if (type == FT56selType)
    {
        ++T56selLocus_atChange.back();
        ++T56Locus_atChange.back();
    } else
    {
        if (type == FT56Type)
        {
            std::cout << "Internal error in GeneticMap. When  calling 'readLoci_SetNewChangeOfType', received a type that is T5 without specification of subtype (ntrl or sel).\n";
            abort();
        } else
        {
            std::cout << "Internal error in GeneticMap. When  calling 'readLoci_SetNewChangeOfType', unknown type "<< type <<" received.\n";
            abort();   
        }
    }*/

    std::cout << "T1Locus_atChange.back() = " << T1Locus_atChange.back() << "\n";
    std::cout << "T2Locus_atChange.back() = " << T2Locus_atChange.back() << "\n";
    std::cout << "T3Locus_atChange.back() = " << T3Locus_atChange.back() << "\n";
    std::cout << "T4Locus_atChange.back() = " << T4Locus_atChange.back() << "\n";
    std::cout << "T56ntrlLocus_atChange.back() = " << T56ntrlLocus_atChange.back() << "\n";
    std::cout << "T56Locus_atChange.back() = " << T56Locus_atChange.back() << "\n";
    std::cout << "T56selLocus_atChange.back() = " << T56selLocus_atChange.back() << "\n";
    std::cout << "T56Locus_atChange.back() = " << T56Locus_atChange.back() << "\n";
    std::cout << "whoWillIncreaseFromChange.back() = " << (size_t)whoWillIncreaseFromChange.back() << "\n";
    

    std::cout << "type added\n";
}


void GeneticMap::readLoci(InputReader& input, std::vector<bool>& isT56LocusUnderSelection, bool isT5ntrlCompress, bool isT5selCompress)
{
    this->T1_nbLoci   = 0;
	this->T1_nbChars  = 0;
    this->T2_nbLoci  = 0;
    this->T3_nbLoci  = 0;
    this->T4_nbLoci   = 0;
    this->T5_nbLoci   = 0;
    this->T5sel_nbLoci = 0;
    this->T5ntrl_nbLoci = 0;
    this->T6_nbLoci   = 0;
    this->T6sel_nbLoci = 0;
    this->T6ntrl_nbLoci = 0;
    this->T56_nbLoci   = 0;
    this->T56sel_nbLoci = 0;
    this->T56ntrl_nbLoci = 0;
    this->TotalNbLoci = 0;

    //////////////////////////////////////////////////////
    // Find change positions and set the nbLoci objects //
    //////////////////////////////////////////////////////
    
    PreviousLocusType previousLocus;
    while( input.IsThereMoreToRead() )
    {
        std::string s = input.GetNextElementString();
        auto type = typeFormat(s, false);
        //std::cout << "For string " << s << " the type '" << (size_t)type << "' was matched\n";
        auto nbElements = input.GetNextElementInt();

        if (nbElements)
        {
            if (type == FT56Type)
            {
                assert(isT56LocusUnderSelection.size() > T56_nbLoci+1);
                type = isT56LocusUnderSelection[T56_nbLoci+1] ? FT56selType : FT56ntrlType;
            }
            
            if (!previousLocus.isSameType(type))
            {
                readLoci_SetNewChangeOfType(type, previousLocus);
            }

            

            if (type == FT1Type)
            {
                TotalNbLoci += nbElements;
                T1_nbLoci += nbElements;
            } else if (type == FT2Type)
            {
                TotalNbLoci += nbElements;
                T2_nbLoci += nbElements;
            } else if (type == FT3Type)
            {
                TotalNbLoci += nbElements;
                T3_nbLoci += nbElements;
            } else if (type == FT4Type)
            {
                TotalNbLoci += nbElements;
                T4_nbLoci += nbElements;
            } else if (type == FT56selType || type == FT56ntrlType)
            {
                size_t futur_T56_nbLoci = T56_nbLoci + nbElements;
                while (T56_nbLoci < futur_T56_nbLoci)
                {
                    assert(isT56LocusUnderSelection.size() > T56_nbLoci+1);
                    auto T5currentType = isT56LocusUnderSelection[T56_nbLoci+1] ? FT56selType : FT56ntrlType;
                    if (!previousLocus.isSameType(T5currentType))
                    {
                        readLoci_SetNewChangeOfType(T5currentType, previousLocus);
                    }

                    if (T5currentType == FT56selType)
                    {
                        ++T56sel_nbLoci;
                        if (isT5selCompress)
                        {
                            ++T6sel_nbLoci;
                        }
                        else
                        {
                            ++T5sel_nbLoci;
                        }
                    } else
                    {
                        assert(T5currentType == FT56ntrlType);
                        ++T56ntrl_nbLoci;
                        if (isT5ntrlCompress)
                        {
                            ++T6ntrl_nbLoci;
                        }
                        else
                        {
                            ++T5ntrl_nbLoci;
                        }
                    }

                    ++TotalNbLoci;
                    ++T56_nbLoci;
                }
            } else
            {
                std::cout << "Internal error in GenetiMap::readLoci, received unknowned formatter type " << type << "\n";
                abort();
            }
        }
    }
    ++this->T1_nbLoci;
    ++this->T1_nbChars;
    ++this->T2_nbLoci;
    ++this->T3_nbLoci;
    ++this->T4_nbLoci;
    ++this->T5_nbLoci;
    ++this->T5sel_nbLoci;
    ++this->T5ntrl_nbLoci;
    ++this->T6_nbLoci;
    ++this->T6sel_nbLoci;
    ++this->T6ntrl_nbLoci;
    ++this->T56_nbLoci;
    ++this->T56sel_nbLoci;
    ++this->T56ntrl_nbLoci;
    ++this->TotalNbLoci;
    assert(TotalNbLoci == T1_nbLoci + T2_nbLoci + T3_nbLoci + T4_nbLoci + T56_nbLoci);


        std::cout << "T1_nbLoci = " << T1_nbLoci << "\n";
        std::cout << "T2_nbLoci = " << T2_nbLoci << "\n";
        std::cout << "T3_nbLoci = " << T3_nbLoci << "\n";
        std::cout << "T4_nbLoci = " << T4_nbLoci << "\n";
        std::cout << "T56ntrl_nbLoci = " << T56ntrl_nbLoci << "\n";
        std::cout << "T56sel_nbLoci = " << T56sel_nbLoci << "\n";
        std::cout << "T56_nbLoci = " << T56_nbLoci << "\n";

    T1_nbChars = T1_nbLoci / 8;

    //////////////////
    // T1 last byte //
    //////////////////
    //std::cout << "line 242\n";
    this->T1_nbLociLastByte = this->T1_nbLoci % 8;
    if (this->T1_nbLociLastByte == 0) {this->T1_nbLociLastByte = 8;}
    assert(this->T1_nbLociLastByte >= 0 && this->T1_nbLociLastByte <= 8);


    //////////////////////////////////
    // If nothing has been received //
    //////////////////////////////////
    if (this->TotalNbLoci <= 0)
    {
        std::cout << "In option '--Loci' you seem to have indicated zero loci. The platform requires to have at least one locus (of either Mode) to run.\n";
        abort();
    }

    ////////////////////////////////////////
    // Get rid of types that are not used //
    ////////////////////////////////////////
    //std::cout << "line 261\n";
    if (T1_nbLoci == 0) T1Locus_atChange.clear();
    if (T2_nbLoci == 0) T2Locus_atChange.clear();
    if (T3_nbLoci == 0) T3Locus_atChange.clear();
    if (T4_nbLoci == 0) T4Locus_atChange.clear();
    if (T56ntrl_nbLoci == 0) T56ntrlLocus_atChange.clear();
    if (T56sel_nbLoci == 0) T56selLocus_atChange.clear();
    if (T56_nbLoci == 0) T56Locus_atChange.clear();

    this->shrink_to_fit(); // Just to really make sure everything is shrunk to fit.

    ////////////////
    // Assertions //
    ////////////////

    this->print();

    assert(T1_nbLoci + T2_nbLoci + T3_nbLoci + T4_nbLoci + T56_nbLoci == TotalNbLoci);
    assert(T56_nbLoci == T56sel_nbLoci + T56ntrl_nbLoci);
    if (isT5selCompress)
    {
        assert(T56sel_nbLoci == T6sel_nbLoci);
        assert(0 == T5sel_nbLoci);
    } else
    {
        assert(T56sel_nbLoci == T5sel_nbLoci);
        assert(0 == T6sel_nbLoci);
    }
    if (isT5ntrlCompress)
    {
        assert(T56ntrl_nbLoci == T6ntrl_nbLoci);
        assert(0 == T5ntrl_nbLoci);
    } else
    {
        assert(T56ntrl_nbLoci == T5ntrl_nbLoci);
        assert(0 == T6ntrl_nbLoci);
    }

    /////////////////////
    // Number of types //
    /////////////////////
    {
        /*
        std::cout << "T1_nbLoci = " << T1_nbLoci << "\n";
        std::cout << "T2_nbLoci = " << T2_nbLoci << "\n";
        std::cout << "T3_nbLoci = " << T3_nbLoci << "\n";
        std::cout << "T4_nbLoci = " << T4_nbLoci << "\n";
        std::cout << "T56ntrl_nbLoci = " << T56ntrl_nbLoci << "\n";
        std::cout << "T56sel_nbLoci = " << T56sel_nbLoci << "\n";
        std::cout << "T56_nbLoci = " << T56_nbLoci << "\n";
        */

        size_t nbTypes = (size_t) (T1_nbLoci > 0) + (T2_nbLoci > 0) + (T3_nbLoci > 0) + (T4_nbLoci > 0) + (T56ntrl_nbLoci > 0) + (T56sel_nbLoci > 0);

        std::cout << "nbTypes = " << nbTypes << "\n";
        assert(nbTypes >= 1);
        if (nbTypes == 1)
        {
            isOnlyOneType = true;
            if (T1_nbLoci > 0)
            {
                onlyType = FT1Type;
            } else if (T2_nbLoci > 0)
            {
                onlyType = FT2Type;
            } else if (T3_nbLoci > 0)
            {
                onlyType = FT3Type;
            } else if (T4_nbLoci > 0)
            {
                onlyType = FT4Type;
            } else if (T56ntrl_nbLoci > 0)
            {                
                onlyType = FT56ntrlType;
            } else if (T56sel_nbLoci > 0)
            {
                onlyType = FT56selType;
            } else
            {
                std::cout << "Internal error in geneticMap. There should be only one type present but it cannot find which one it is!\n";
                abort();
            }
        } else
        {
            isOnlyOneType = false;
            onlyType = 250; // Just some non-sense number
        }
    }

    //std::cout << "line 347\n";
    if (!isOnlyOneType)
    {
        ///////////////////////////////////
        // Finish map for faster hashing //
        ///////////////////////////////////
        size_t nbChanges = TallLocus_atChange.size();
        assert(nbChanges == whoWillIncreaseFromChange.size());

        fromLocusToLastChangeIndex.reserve(TotalNbLoci);
        fromT1LocusToLastChangeIndex.reserve(T1_nbLoci);
        fromT2LocusToLastChangeIndex.reserve(T2_nbLoci);
        fromT3LocusToLastChangeIndex.reserve(T3_nbLoci);
        fromT4LocusToLastChangeIndex.reserve(T4_nbLoci);
        fromT56ntrlLocusToLastChangeIndex.reserve(T56ntrl_nbLoci);
        fromT56selLocusToLastChangeIndex.reserve(T56sel_nbLoci);
        fromT56LocusToLastChangeIndex.reserve(T56_nbLoci);

        this->print();
        
        for (size_t change_index = 0 ; change_index < nbChanges ; ++change_index)
        {
            std::cout << "change_index = " << change_index << "\n";
            std::cout << "expandMapToChangeIndex for Tall\n";
            if (TotalNbLoci) expandMapToChangeIndex(fromLocusToLastChangeIndex,   TallLocus_atChange,  250,     change_index); // non-sense type for Tall
            std::cout << "expandMapToChangeIndex for T1\n";
            if (T1_nbLoci) expandMapToChangeIndex(fromT1LocusToLastChangeIndex, T1Locus_atChange,    FT1Type, change_index); 
            std::cout << "expandMapToChangeIndex for T2\n";
            if (T2_nbLoci) expandMapToChangeIndex(fromT2LocusToLastChangeIndex, T2Locus_atChange,    FT2Type, change_index); 
            std::cout << "expandMapToChangeIndex for T3\n";
            if (T3_nbLoci) expandMapToChangeIndex(fromT3LocusToLastChangeIndex, T3Locus_atChange,    FT3Type, change_index); 
            std::cout << "expandMapToChangeIndex for T4\n";
            if (T4_nbLoci) expandMapToChangeIndex(fromT4LocusToLastChangeIndex, T4Locus_atChange,    FT4Type, change_index); 
            std::cout << "expandMapToChangeIndex for T56\n";
            if (T56_nbLoci) expandMapToChangeIndex(fromT56LocusToLastChangeIndex, T56Locus_atChange,   FT56Type, change_index); 
            std::cout << "expandMapToChangeIndex for T56ntrl\n";
            if (T56ntrl_nbLoci) expandMapToChangeIndex(fromT56ntrlLocusToLastChangeIndex, T56ntrlLocus_atChange,    FT56ntrlType, change_index); 
            std::cout << "expandMapToChangeIndex for T56sel\n";
            if (T56sel_nbLoci) expandMapToChangeIndex(fromT56selLocusToLastChangeIndex, T56selLocus_atChange,    FT56selType, change_index); 
        }

        ////////////////
        // Assertions //
        ////////////////

        this->print();

        //std::cout << "line 382\n";
        assert(fromLocusToLastChangeIndex.size()            == TotalNbLoci); 
        assert(fromT1LocusToLastChangeIndex.size()          == T1_nbLoci); 
        assert(fromT2LocusToLastChangeIndex.size()          == T2_nbLoci); 
        assert(fromT3LocusToLastChangeIndex.size()          == T3_nbLoci); 
        assert(fromT4LocusToLastChangeIndex.size()          == T4_nbLoci); 
        assert(fromT56ntrlLocusToLastChangeIndex.size()     == T56ntrl_nbLoci); 
        assert(fromT56selLocusToLastChangeIndex.size()      == T56sel_nbLoci); 
        assert(fromT56LocusToLastChangeIndex.size()         == T56_nbLoci); 

        assert(fromLocusToLastChangeIndex.capacity()        == TotalNbLoci); 
        assert(fromT1LocusToLastChangeIndex.capacity()      == T1_nbLoci); 
        assert(fromT2LocusToLastChangeIndex.capacity()      == T2_nbLoci); 
        assert(fromT3LocusToLastChangeIndex.capacity()      == T3_nbLoci); 
        assert(fromT4LocusToLastChangeIndex.capacity()      == T4_nbLoci); 
        assert(fromT56ntrlLocusToLastChangeIndex.capacity() == T56ntrl_nbLoci); 
        assert(fromT56selLocusToLastChangeIndex.capacity()  == T56sel_nbLoci); 
        assert(fromT56LocusToLastChangeIndex.capacity()     == T56_nbLoci); 
    }
    //std::cout << "line 401\n";
}

template<typename T>
void GeneticMap::print(std::string name, std::vector<T>& v)
{
    std::cout << name << ": ";
    for (auto& e : v) std::cout << (size_t)e << " ";
    std::cout << "\n";
}

template<typename T>
void GeneticMap::print(std::string name, T& v)
{
    std::cout << name << ": " << (size_t) v << "\n";
}


void GeneticMap::print()
{
    std::cout << "\n-------------------------\n\n";
    print("TallLocus_atChange",TallLocus_atChange);
    print("T1Locus_atChange",T1Locus_atChange);
    print("T2Locus_atChange",T2Locus_atChange);
    print("T3Locus_atChange",T3Locus_atChange);
    print("T4Locus_atChange",T4Locus_atChange);
    print("T56ntrlLocus_atChange",T56ntrlLocus_atChange);
    print("T56selLocus_atChange",T56selLocus_atChange);
    print("T56Locus_atChange",T56Locus_atChange);
    print("whoWillIncreaseFromChange",whoWillIncreaseFromChange);
    print("fromLocusToLastChangeIndex",fromLocusToLastChangeIndex);
    print("fromT1LocusToLastChangeIndex",fromT1LocusToLastChangeIndex);
    print("fromT2LocusToLastChangeIndex",fromT2LocusToLastChangeIndex);
    print("fromT3LocusToLastChangeIndex",fromT3LocusToLastChangeIndex);
    print("fromT4LocusToLastChangeIndex",fromT4LocusToLastChangeIndex);
    print("fromT56ntrlLocusToLastChangeIndex",fromT56ntrlLocusToLastChangeIndex);
    print("fromT56selLocusToLastChangeIndex",fromT56selLocusToLastChangeIndex);
    print("fromT56LocusToLastChangeIndex",fromT56LocusToLastChangeIndex);

    print("T1_nbLoci", T1_nbLoci);
    print("T1_nbChars", T1_nbChars);
    print("T1_nbLociLastByte", T1_nbLociLastByte);
    print("T2_nbLoci", T2_nbLoci);
    print("T3_nbLoci", T3_nbLoci);
    print("T4_nbLoci", T4_nbLoci);
    print("T5_nbLoci", T5_nbLoci);
    print("T5sel_nbLoci", T5sel_nbLoci);
    print("T5ntrl_nbLoci", T5ntrl_nbLoci);
    print("T6_nbLoci", T6_nbLoci);
    print("T6sel_nbLoci", T6sel_nbLoci);
    print("T6ntrl_nbLoci", T6ntrl_nbLoci);
    print("T56_nbLoci", T56_nbLoci);
    print("T56sel_nbLoci", T56sel_nbLoci);
    print("T56ntrl_nbLoci", T56ntrl_nbLoci);
    print("TotalNbLoci", TotalNbLoci);
    print("isT5selCompress", isT5selCompress);
    print("isT5ntrlCompress", isT5ntrlCompress);
    std::cout << "\n-------------------------\n\n";
}



template<typename INT>
size_t GeneticMap::getLocusFromTxLocus(INT TxLocus, const std::vector<size_t>& fromTxLocusToLastChangeIndex, const std::vector<size_t>& TxLocus_atChange) const
{
	if (isOnlyOneType) return TxLocus;

	auto& CI = fromTxLocusToLastChangeIndex[TxLocus];
	auto& locusAtChange = TallLocus_atChange[CI];
	auto& TxLocusAtChange = TxLocus_atChange[CI];

	return locusAtChange + (TxLocus - TxLocusAtChange);
}


template<typename INT>
size_t GeneticMap::getTxLocusFromLocus(INT Locus, const std::vector<size_t>& TxLocus_atChange, const unsigned char type) const
{
    std::cout << "getTxLocusFromLocus type = " << type << "\n";
    assert(type == FT1Type && type == FT2Type && type == FT3Type && type == FT4Type && type == FT56Type && type == FT56ntrlType && type == FT56selType );
    if (isOnlyOneType)
    {
        if (type == onlyType)
        {
            return Locus;
        } else
        {
            return 0;
        }
    } 

    auto& CI = TallLocus_atChange[Locus];
    auto& TxLocus = TxLocus_atChange[type];
    if (whoWillIncreaseFromChange[CI] == type)
    {
        return TxLocus + Locus - TallLocus_atChange[CI];
    } else
    {
        return TxLocus;
    }
}






template<typename INT>
T56LocusGender GeneticMap::getT56GLocusFromT56Locus(INT T56Locus) const
{
    if (isOnlyOneType)
    {
        return T56LocusGender(
            onlyType == FT56ntrlType ? true : false,
            T56Locus 
        );
    }

	auto& CI = fromT56LocusToLastChangeIndex[T56Locus];

	if (whoWillIncreaseFromChange[CI] == FT56ntrlType)
	{
		// ntrl
		return T56LocusGender(
			true,
			T56ntrlLocus_atChange[CI] + (T56Locus - T56Locus_atChange[CI])
		);
	} else
	{
		assert(whoWillIncreaseFromChange[CI] == FT56selType);
		// sel
		return T56LocusGender(
			false,
			T56selLocus_atChange[CI] + (T56Locus - T56Locus_atChange[CI])
		);
	}
}


template<typename INT>
bool GeneticMap::isT56neutral(INT T56Locus) const
{
	return getT56GLocusFromT56Locus(T56Locus).isNtrl;
}







template<typename INT>
size_t GeneticMap::getT1LocusFromLocus(INT locus) const
{
	return getTxLocusFromLocus(locus, T1Locus_atChange, FT1Type);
}

template<typename INT>
size_t GeneticMap::getT2LocusFromLocus(INT locus) const
{
	return getTxLocusFromLocus(locus, T2Locus_atChange, FT2Type);
}

template<typename INT>
size_t GeneticMap::getT3LocusFromLocus(INT locus) const
{
	return getTxLocusFromLocus(locus, T3Locus_atChange, FT3Type);
}

template<typename INT>
size_t GeneticMap::getT4LocusFromLocus(INT locus) const
{
	return getTxLocusFromLocus(locus, T4Locus_atChange, FT4Type);
}

template<typename INT>
size_t GeneticMap::getT56ntrlLocusFromLocus(INT locus) const
{
	return getTxLocusFromLocus(locus, T56ntrlLocus_atChange, FT56ntrlType);
}

template<typename INT>
size_t GeneticMap::getT56selLocusFromLocus(INT locus) const
{
	return getTxLocusFromLocus(locus, T56selLocus_atChange, FT56selType);
}





template<typename INT>
size_t GeneticMap::getLocusFromT1Locus(INT T1Locus) const
{
	return getLocusFromTxLocus(T1Locus, fromT1LocusToLastChangeIndex, T1Locus_atChange);
}

template<typename INT>
size_t GeneticMap::getLocusFromT2Locus(INT T2Locus) const
{
	return getLocusFromTxLocus(T2Locus, fromT2LocusToLastChangeIndex, T2Locus_atChange);
}


template<typename INT>
size_t GeneticMap::getLocusFromT3Locus(INT T3Locus) const
{
	return getLocusFromTxLocus(T3Locus, fromT3LocusToLastChangeIndex, T3Locus_atChange);
}

template<typename INT>
size_t GeneticMap::getLocusFromT4Locus(INT T4Locus) const
{ 
	return getLocusFromTxLocus(T4Locus, fromT4LocusToLastChangeIndex, T4Locus_atChange);
}

template<typename INT>
size_t GeneticMap::getLocusFromT56ntrlLocus(INT T56ntrlLocus) const
{
	return getLocusFromTxLocus(T56ntrlLocus, fromT56ntrlLocusToLastChangeIndex, T56ntrlLocus_atChange);
}

template<typename INT>
size_t GeneticMap::getLocusFromT56selLocus(INT T56selLocus) const
{
	return getLocusFromTxLocus(T56selLocus, fromT56selLocusToLastChangeIndex, T56selLocus_atChange);
}

template<typename INT>
size_t GeneticMap::getLocusFromT56Locus(INT T5Locus) const
{
    if (isT56neutral(T5Locus))
    {
        return getLocusFromT56ntrlLocus(T5Locus);
    } else
    {
        return getLocusFromT56selLocus(T5Locus);
    }
}

template<typename INT>
unsigned char GeneticMap::getLocusSubtypeFromLocus(INT locus) const
{
    if (isOnlyOneType)
    {
        if (T1_nbLoci > 0)
        {
            return FT1Type;
        }
        if (T2_nbLoci > 0)
        {
            return FT2Type;
        }
        if (T3_nbLoci > 0)
        {
            return FT3Type;
        }
        if (T4_nbLoci > 0)
        {
            return FT4Type;
        }
        if (T56sel_nbLoci > 0)
        {
            return FT56selType;
        }
        if (T56ntrl_nbLoci > 0)
        {
            return FT56ntrlType;
        }
    }


    auto& CI = TallLocus_atChange[locus];
    return whoWillIncreaseFromChange[CI];
}


/*
void GeneticMap::expandMapToChangeIndexWithRepeatedElement(std::vector<size_t>& v, size_t e, size_t n)
{
    for (size_t i = 0 ; i < n ; ++i)
        v.push_back(e);
}

void GeneticMap::expandMapToChangeIndexWithIncreasingElement(std::vector<size_t>& v, size_t from, size_t to)
{
    for ( ; from < to ; ++from)
        v.push_back(from);
}*/

void GeneticMap::push_back_repeated_value(std::vector<size_t>& v, size_t value, size_t nbRepeats)
{
    v.reserve(v.size() + nbRepeats);
    for (size_t i = 0 ; i < nbRepeats; ++i)
        v.push_back(value);
}

void GeneticMap::expandMapToChangeIndex(std::vector<size_t>& v, std::vector<size_t>& Tx_atChange, unsigned char type, size_t change_index)
{
    assert(Tx_atChange.size() > change_index);
    assert(whoWillIncreaseFromChange.size() > change_index);

    if (
        type ==250 ||
        whoWillIncreaseFromChange[change_index] == type ||
        (type == FT56Type && (whoWillIncreaseFromChange[change_index] == FT56ntrlType || whoWillIncreaseFromChange[change_index] == FT56selType))
    )
    {
        auto from = Tx_atChange[change_index];
        size_t to;
        if (Tx_atChange.size() > change_index + 1)
        {
            to = Tx_atChange[change_index + 1];
        } else
        {
            assert(TallLocus_atChange.size() > change_index);
            size_t nbLociLeftToAdd = TotalNbLoci - TallLocus_atChange[change_index];
            to = from + nbLociLeftToAdd;
            assert(to > from);
        }

        auto nbRepeats = to - from + 1;
        assert(nbRepeats > 0);
        push_back_repeated_value(v, change_index, nbRepeats);
    }
/*
    auto from = Tx_atChange[change_index];
    size_t to;
    if (Tx_atChange.size() > change_index + 1)
    {
        to = Tx_atChange[change_index + 1];

        std::cout << "from = " << from << "\n";
        std::cout << "to = " << to << "\n";

        if (from != -1)
        {
            if (from == to)
            {
                if (type != 250 && type != FT56Type) assert(whoWillIncreaseFromChange[change_index] != type);
            } else
            {
                assert(to > from);
                if (type != 250 && type != FT56Type) assert(whoWillIncreaseFromChange[change_index] == type);
            }
        } else
        {
            assert(to == 0);
        }
    } else
    {
        if (whoWillIncreaseFromChange[change_index] == type)
        {
            size_t nbLociLeftToAdd = TotalNbLoci - TallLocus_atChange[change_index];
            to = from + nbLociLeftToAdd;
            assert(to > from);
        } else
        {
            to = from;
        }
    }

    if (from == to)
    {
        expandMapToChangeIndexWithRepeatedElement(v, from, to - from);
    } else
    {
        assert(to > from);
        expandMapToChangeIndexWithIncreasingElement(v, from, to);
    }*/
}