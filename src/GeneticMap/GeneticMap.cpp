
void GeneticMap::shrink_to_fit()
{
    _FromT1LocusToLocus.shrink_to_fit();
    _FromT2LocusToLocus.shrink_to_fit();
    _FromT3LocusToLocus.shrink_to_fit();
    _FromT4LocusToLocus.shrink_to_fit();
    _FromT56LocusToLocus.shrink_to_fit();
    _FromT56ntrlLocusToLocus.shrink_to_fit();
    _FromT56selLocusToLocus.shrink_to_fit();
    _FromLocusToNextT1Locus.shrink_to_fit();
    _FromLocusToNextT2Locus.shrink_to_fit();
    _FromLocusToNextT3Locus.shrink_to_fit();
    _FromLocusToNextT4Locus.shrink_to_fit();
    _FromLocusToNextT56Locus.shrink_to_fit();
    _FromLocusToNextT56ntrlLocus.shrink_to_fit();
    _FromLocusToNextT56selLocus.shrink_to_fit();
    _FromT56LocusToT56genderLocus.shrink_to_fit();
}


void GeneticMap::setT56GenderLoci(std::vector<std::vector<double>>& T56_fit, bool alreadyKnowItIsAllNtrl, bool alreadyKnowItIsAllSel)
{   
    if (true) //(!(alreadyKnowItIsAllNtrl || alreadyKnowItIsAllSel))
    {
        assert(T56_fit.size());
        uint32_t ntrl_index = 0;
        uint32_t sel_index = 0;
        assert(_FromT56LocusToT56genderLocus.size() == 0);
        _FromT56LocusToT56genderLocus.reserve(T56_fit[0].size());

        assert(isT56ntrlUsed == false);
        assert(isT56selUsed == false);
        nbT56LociGotFromFitnessValues = T56_fit[0].size();

        for (uint32_t i = 0 ; i < nbT56LociGotFromFitnessValues; ++i)
        {
            bool isNtrl = true;
            for (uint32_t habitat = 0 ; habitat < T56_fit.size() ; ++habitat)
            {
                if (T56_fit[habitat][i] != 1.0)
                {
                    isNtrl = false;
                    break;
                }
            }

            if (isNtrl)
            {
                isT56ntrlUsed = true;
                _FromT56LocusToT56genderLocus.push_back({true, ntrl_index});
                ++ntrl_index;
            } else
            {
                isT56selUsed = true;
                _FromT56LocusToT56genderLocus.push_back({false, sel_index});
                ++sel_index;
            }
        }   
        assert(_FromT56LocusToT56genderLocus.size() == nbT56LociGotFromFitnessValues);

        if (!(isT56ntrlUsed && isT56selUsed))
        {
            _FromT56LocusToT56genderLocus.clear();
            _FromT56LocusToT56genderLocus.shrink_to_fit();
        }
    } else
    {
        assert(alreadyKnowItIsAllSel != alreadyKnowItIsAllNtrl);
        if (alreadyKnowItIsAllSel)
        {
            isT56ntrlUsed = false;
            isT56selUsed = true;
        } else
        {
            assert(alreadyKnowItIsAllNtrl);
            isT56ntrlUsed = true;
            isT56selUsed = false;
        }
    }
}


int GeneticMap::readType(const std::string& s_type) const
{
    int type;
    if (s_type.size() == 2)
    {
        type = (int) std::stod(s_type.substr(1));
    } else
    {
        if (s_type.size() != 1)
        {
            std::cout << "For option --L (--Loci), received unknown type " << s_type << ". Only types accepted are T1, T2, T3, T4 and T5. Note that the 'T' can be lower case (t1, t2, ...) and can be ignored (1, 2, ...). Note also that T5 might be compressed to T6 (see compression options) but you must still call it 'T5' and not 'T6' (not 'T56' either as SimBit does internally).\n";
            abort();
        }
        type = (int) std::stod(s_type);
    }
    return type;
}

void GeneticMap::readT56Compression(InputReader& input)
{
    // Ntrl
    if (input.PeakNextElementString() == "default")
    {
        input.skipElement();
        // Nothing to do set it has been set previously in SSP::readT56_fitnessEffects. That's a bit weird but readLoci needs T56_compress, so I can't read --Loci before this and I want to set the compression depending on the number of loci!
    } else
    {
        isT56ntrlCompress = input.GetNextElementBool();
    }   

    // Sel
    if (input.PeakNextElementString() == "default")
    {
        input.skipElement();
        // Nothing to set as it has been set previously in SSP::readT56_fitnessEffects. That's a bit weird but readLoci needs T56_compress, so I can't read --Loci before this and I want to set the compression depending on the number of loci!
    } else
    {
        isT56selCompress = input.GetNextElementBool();
    }
}


void GeneticMap::readLoci(InputReader& input)
{
    /////////////////////////
    // Set nbLoci to zero //
    /////////////////////////

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


    /////////////////////////
    // Figure unused types //
    /////////////////////////

    assert(isT1Used == false);
    assert(isT2Used == false);
    assert(isT3Used == false);
    assert(isT4Used == false);
    assert(isT56Used == false);
    {
        uint32_t sumNbElements = 0;
        auto inputCopy = input;
        while( inputCopy.IsThereMoreToRead() )
        {
            auto s_type = inputCopy.GetNextElementString();
            auto type = readType(s_type);
            auto nbElements = inputCopy.GetNextElementInt();
            sumNbElements += nbElements;
            if (nbElements < 0)
            {
                std::cout << "For option --L (--Loci), received a negative number of loci of type " << s_type << ".\n";
                abort();
            }

            if (nbElements)
            {
                switch (type)
                {
                    case 1:
                        isT1Used = true;
                        break;
                    case 2:
                        isT2Used = true;
                        break;
                    case 3:
                        isT3Used = true;
                        break;
                    case 4:
                        isT4Used = true;
                        break;
                    case 5:
                        isT56Used = true;
                        break;
                    default:
                        std::cout << "For option --L (--Loci), received unknown type " << s_type << ". Only types accepted are T1, T2, T3, T4 and T5. Note that the 'T' can be lower case (t1, t2, ...) and can be ignored (1, 2, ...). Note also that T5 might be compressed to T6 (see compression options) but you must still call it 'T5' and not 'T6' (not 'T56' either as SimBit does internally).\n";
                        abort();
                }
            }
        }
        if (sumNbElements == 0)
        {
            std::cout << "In option '--Loci' you seem to have indicated zero loci. The platform requires to have at least one locus (of any type) to run. This might be a silly limitations if you are only interested in demography but I just have not tested any cases with zero loci and it won't cost you much to just add one T1 or one T5 locus even if you're only interested in demography.\n";
            abort();
        }
    }

    if (isT56Used)
    {
        assert(isT56ntrlUsed || isT56selUsed);
    } else
    {
        assert(!isT56ntrlUsed && !isT56selUsed);
    }


    auto nbTypesUsed = isT1Used + isT2Used + isT3Used + isT4Used + isT56ntrlUsed + isT56selUsed;
    assert(nbTypesUsed > 0);
    if (nbTypesUsed == 1)
    {
        if (isT1Used)
            onlyTypeUsed = 1;
        else if (isT2Used)
            onlyTypeUsed = 3;
        else if (isT3Used)
            onlyTypeUsed = 3;
        else if (isT4Used)
            onlyTypeUsed = 4;
        else if (isT56ntrlUsed)
            onlyTypeUsed = 50;
        else if (isT56selUsed)
            onlyTypeUsed = 51;
        else
        {
            std::cout << "Internal error in GeneticMap::readLoci. SimBit thinks only one type of locus is used but it could not figure which one it is\n";
            abort();
        }
    } else
    {
        // Several different types are used (could be T56ntrl and T56sel only)
        onlyTypeUsed = 250; // non-sense value
    }


    //////////////////////////////////////////////////////
    // Find change positions and set the nbLoci objects //
    //////////////////////////////////////////////////////
    
    while( input.IsThereMoreToRead() )
    {
        auto s_type = input.GetNextElementString();
        auto type = readType(s_type);
        auto nbElements = input.GetNextElementInt();

        if (nbElements)
        {
            for (uint32_t element = 0 ; element < nbElements; ++element)
            {
                switch (type)
                {
                    case 1:
                        if (onlyTypeUsed == 250) _FromT1LocusToLocus.push_back(TotalNbLoci);
                        ++T1_nbLoci;
                        break;
                    case 2:
                        if (onlyTypeUsed == 250)_FromT2LocusToLocus.push_back(TotalNbLoci);
                        ++T2_nbLoci;
                        break;
                    case 3:
                        if (onlyTypeUsed == 250) _FromT3LocusToLocus.push_back(TotalNbLoci);
                        ++T3_nbLoci;
                        break;
                    case 4:
                        if (onlyTypeUsed == 250) _FromT4LocusToLocus.push_back(TotalNbLoci);
                        ++T4_nbLoci;
                        break;
                    case 5:
                        if (nbT56LociGotFromFitnessValues <= T56_nbLoci)
                        {
                            std::cout << "Received at least " << T56_nbLoci << " loci (in option --L aka. --Loci) of type T5 (or T6 if compressed) but received only " << nbT56LociGotFromFitnessValues << " fitness values for T5 loci (with option --T5_fit)\n";
                            abort();
                        }

                        if (isT56neutral(T56_nbLoci))
                        {
                            if (onlyTypeUsed == 250) _FromT56ntrlLocusToLocus.push_back(TotalNbLoci);
                            ++T56ntrl_nbLoci;                        
                        } else
                        {
                            if (onlyTypeUsed == 250) _FromT56selLocusToLocus.push_back(TotalNbLoci);
                            ++T56sel_nbLoci;                        
                        }

                        if (onlyTypeUsed == 250) _FromT56LocusToLocus.push_back(TotalNbLoci);
                        ++T56_nbLoci;
                        break;
                    default:
                        std::cout << "line 277\n";
                        std::cout << "For option --L (--Loci), received unknown type " << s_type << ". Only types accepted are T1, T2, T3, T4 and T5. Note that the 'T' can be lower case (t1, t2, ...) and can be ignored (1, 2, ...). Note also that T5 might be compressed to T6 (see compression options) but you must still call it 'T5' and not 'T6' (not 'T56' either as SimBit does internally). Note this error was found at the second security gate within GeneticMap::readLoci. That sounuds like an internal bug!\n";
                        abort();
                }
                ++TotalNbLoci;

                if (onlyTypeUsed == 250)
                {
                    if (isT1Used) _FromLocusToNextT1Locus.push_back(T1_nbLoci);
                    if (isT2Used) _FromLocusToNextT2Locus.push_back(T2_nbLoci);
                    if (isT3Used) _FromLocusToNextT3Locus.push_back(T3_nbLoci);
                    if (isT4Used) _FromLocusToNextT4Locus.push_back(T4_nbLoci);
                    if (isT56Used) _FromLocusToNextT56Locus.push_back(T56_nbLoci);
                    if (isT56ntrlUsed) _FromLocusToNextT56ntrlLocus.push_back(T56ntrl_nbLoci);
                    if (isT56selUsed) _FromLocusToNextT56selLocus.push_back(T56sel_nbLoci);
                }
            }
        }
    }

    //this->print();
    
    assert(TotalNbLoci == T1_nbLoci + T2_nbLoci + T3_nbLoci + T4_nbLoci + T56_nbLoci);
    assert(T56_nbLoci == T56sel_nbLoci + T56ntrl_nbLoci);
    if (nbT56LociGotFromFitnessValues > T56_nbLoci)
    {
        std::cout << "Received " << T56_nbLoci << " loci (in option --L aka. --Loci) of type T5 (or T6 if compressed) but received only " << nbT56LociGotFromFitnessValues << " fitness values for T5 loci (with option --T5_fit)\n";
        abort();
    }
    
    shrink_to_fit();


    ///////////////
    // T56 stuff //
    ///////////////    

    if (isT56ntrlCompress)
    {
        T6ntrl_nbLoci = T56ntrl_nbLoci;
    } else
    {
        T5ntrl_nbLoci = T56ntrl_nbLoci;
    }

    if (isT56selCompress)
    {
        T6sel_nbLoci = T56sel_nbLoci;
    } else
    {
        T5sel_nbLoci = T56sel_nbLoci;
    }

    T5_nbLoci = T5sel_nbLoci + T5ntrl_nbLoci;
    T6_nbLoci = T6sel_nbLoci + T6ntrl_nbLoci;
    assert(T5ntrl_nbLoci == 0 || T6ntrl_nbLoci == 0);
    assert(T5sel_nbLoci == 0 || T6sel_nbLoci == 0);
    assert(T5ntrl_nbLoci + T6ntrl_nbLoci == T56ntrl_nbLoci);
    assert(T5sel_nbLoci + T6sel_nbLoci == T56sel_nbLoci);
    assert(T5_nbLoci + T6_nbLoci == T56_nbLoci);
    assert(T56_nbLoci == T56sel_nbLoci + T56ntrl_nbLoci);

    //////////////
    // T1 stuff //
    //////////////

    T1_nbChars = ceil((double)T1_nbLoci / 8.0);
    T1_nbLociLastByte = T1_nbLoci % 8;
    if (T1_nbLociLastByte == 0) {T1_nbLociLastByte = 8;}
    assert(T1_nbLociLastByte >= 0 && T1_nbLociLastByte <= 8);
    assert((T1_nbChars - 1) * 8 + T1_nbLociLastByte == T1_nbLoci);


    /////////////////////
    // More Assertions //
    /////////////////////

    assert(TotalNbLoci);
    assert(T1_nbLoci + T2_nbLoci + T3_nbLoci + T4_nbLoci + T56_nbLoci == TotalNbLoci);
    if (isT56selCompress)
    {
        assert(T56sel_nbLoci == T6sel_nbLoci);
        assert(0 == T5sel_nbLoci);
    } else
    {
        assert(T56sel_nbLoci == T5sel_nbLoci);
        assert(0 == T6sel_nbLoci);
    }
    if (isT56ntrlCompress)
    {
        assert(T56ntrl_nbLoci == T6ntrl_nbLoci);
        assert(0 == T5ntrl_nbLoci);
    } else
    {
        assert(T56ntrl_nbLoci == T5ntrl_nbLoci);
        assert(0 == T6ntrl_nbLoci);
    }

    

    if (onlyTypeUsed == 250 && isT1Used)
        assert(_FromT1LocusToLocus.size() == T1_nbLoci); 
    else
        assert(_FromT1LocusToLocus.size() == 0);

    if (onlyTypeUsed == 250 && isT2Used)
        assert(_FromT2LocusToLocus.size() == T2_nbLoci); 
    else
        assert(_FromT2LocusToLocus.size() == 0);

    if (onlyTypeUsed == 250 && isT3Used)
        assert(_FromT3LocusToLocus.size() == T3_nbLoci); 
    else
        assert(_FromT3LocusToLocus.size() == 0);

    if (onlyTypeUsed == 250 && isT4Used)
        assert(_FromT4LocusToLocus.size() == T4_nbLoci); 
    else
        assert(_FromT4LocusToLocus.size() == 0);

    if (onlyTypeUsed == 250 && isT56Used)
        assert(_FromT56LocusToLocus.size() == T56_nbLoci); 
    else
        assert(_FromT56LocusToLocus.size() == 0);

    if (onlyTypeUsed == 250 && isT56ntrlUsed)
        assert(_FromT56ntrlLocusToLocus.size() == T56ntrl_nbLoci); 
    else
        assert(_FromT56ntrlLocusToLocus.size() == 0);

    if (onlyTypeUsed == 250 && isT56selUsed)
        assert(_FromT56selLocusToLocus.size() == T56sel_nbLoci); 
    else
        assert(_FromT56selLocusToLocus.size() == 0);

    if (onlyTypeUsed == 250 && isT1Used)
        assert(_FromLocusToNextT1Locus.size() == TotalNbLoci); 
    else
        assert(_FromLocusToNextT1Locus.size() == 0);

    if (onlyTypeUsed == 250 && isT2Used)
        assert(_FromLocusToNextT2Locus.size() == TotalNbLoci); 
    else
        assert(_FromLocusToNextT2Locus.size() == 0);

    if (onlyTypeUsed == 250 && isT3Used)
        assert(_FromLocusToNextT3Locus.size() == TotalNbLoci); 
    else
        assert(_FromLocusToNextT3Locus.size() == 0);

    if (onlyTypeUsed == 250 && isT4Used)
        assert(_FromLocusToNextT4Locus.size() == TotalNbLoci); 
    else
        assert(_FromLocusToNextT4Locus.size() == 0);

    if (onlyTypeUsed == 250 && isT56Used)
        assert(_FromLocusToNextT56Locus.size() == TotalNbLoci); 
    else
        assert(_FromLocusToNextT56Locus.size() == 0);

    if (onlyTypeUsed == 250 && isT56ntrlUsed)
        assert(_FromLocusToNextT56ntrlLocus.size() == TotalNbLoci); 
    else
        assert(_FromLocusToNextT56ntrlLocus.size() == 0);

    if (onlyTypeUsed == 250 && isT56selUsed)
        assert(_FromLocusToNextT56selLocus.size() == TotalNbLoci); 
    else
        assert(_FromLocusToNextT56selLocus.size() == 0);

    if (onlyTypeUsed == 250 && isT56ntrlUsed && isT56selUsed)
        assert(_FromT56LocusToT56genderLocus.size() == T56_nbLoci);
    else
        assert(_FromT56LocusToT56genderLocus.size() == 0);

    //this->print();
}

template<typename T>
void GeneticMap::print(std::string name, std::vector<T>& v)
{
    std::cout << name << ": ";
    for (auto& e : v) std::cout << (uint32_t)e << " ";
    std::cout << "\n";
}

void GeneticMap::print(std::string name, std::vector<T56LocusGender>& v)
{
    std::cout << name << ": ";
    for (auto& e : v) std::cout << "{" << e.isNtrl << ", " << e.locusInGender << "} ";
    std::cout << "\n";
}

template<typename T>
void GeneticMap::print(std::string name, T& v)
{
    std::cout << name << ": " << (uint32_t) v << "\n";
}


void GeneticMap::print()
{
    std::cout << "\n-------------------------\n\n";
    print("_FromT1LocusToLocus", _FromT1LocusToLocus);
    print("_FromT2LocusToLocus", _FromT2LocusToLocus);
    print("_FromT3LocusToLocus", _FromT3LocusToLocus);
    print("_FromT4LocusToLocus", _FromT4LocusToLocus);
    print("_FromT56LocusToLocus", _FromT56LocusToLocus);
    print("_FromT56ntrlLocusToLocus", _FromT56ntrlLocusToLocus);
    print("_FromT56selLocusToLocus", _FromT56selLocusToLocus);
    print("_FromLocusToT1Locus", _FromLocusToNextT1Locus);
    print("_FromLocusToT2Locus", _FromLocusToNextT2Locus);
    print("_FromLocusToT3Locus", _FromLocusToNextT3Locus);
    print("_FromLocusToT4Locus", _FromLocusToNextT4Locus);
    print("_FromLocusToT56Locus", _FromLocusToNextT56Locus);
    print("_FromLocusToT56ntrlLocus", _FromLocusToNextT56ntrlLocus);
    print("_FromLocusToT56selLocus", _FromLocusToNextT56selLocus);
    print("_FromT56LocusToT56genderLocus", _FromT56LocusToT56genderLocus);

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
    print("isT56selCompress", isT56selCompress);
    print("isT56ntrlCompress", isT56ntrlCompress);
    std::cout << "\n-------------------------\n\n";
}


template<typename INT>
uint32_t GeneticMap::FromT1LocusToLocus(const INT i) const
{
    if (onlyTypeUsed == 250)
    {
        assert(_FromT1LocusToLocus.size() > i);
        return _FromT1LocusToLocus[i];
    } else
    {
        assert(onlyTypeUsed == 1);
        return i;
    }
}

template<typename INT>
uint32_t GeneticMap::FromT2LocusToLocus(const INT i) const
{
    if (onlyTypeUsed == 250)
    {
        assert(_FromT2LocusToLocus.size() > i);
        return _FromT2LocusToLocus[i];
    } else
    {
        assert(onlyTypeUsed == 2);
        return i;
    }
}

template<typename INT>
uint32_t GeneticMap::FromT3LocusToLocus(const INT i) const
{
    if (onlyTypeUsed == 250)
    {
        assert(_FromT3LocusToLocus.size() > i);
        return _FromT3LocusToLocus[i];
    } else
    {
        assert(onlyTypeUsed == 3);
        return i;
    }
}

template<typename INT>
uint32_t GeneticMap::FromT4LocusToLocus(const INT i) const
{
    if (onlyTypeUsed == 250)
    {
        assert(_FromT4LocusToLocus.size() > i);
        return _FromT4LocusToLocus[i];
    } else
    {
        assert(onlyTypeUsed == 4);
        return i;
    }
}

template<typename INT>
uint32_t GeneticMap::FromT56LocusToLocus(const INT i) const
{
    if (onlyTypeUsed == 250)
    {
        assert(_FromT56LocusToLocus.size() > i);
        return _FromT56LocusToLocus[i];
    } else
    {
        assert(onlyTypeUsed == 50 || onlyTypeUsed == 51);
        return i;
    }
}

template<typename INT>
uint32_t GeneticMap::FromT56ntrlLocusToLocus(const INT i) const
{
    if (onlyTypeUsed == 250)
    {
        assert(_FromT56ntrlLocusToLocus.size() > i);
        return _FromT56ntrlLocusToLocus[i];
    } else
    {
        assert(onlyTypeUsed == 50);
        return i;
    }
}

template<typename INT>
uint32_t GeneticMap::FromT56selLocusToLocus(const INT i) const
{
    if (onlyTypeUsed == 250)
    {
        assert(_FromT56selLocusToLocus.size() > i);
        return _FromT56selLocusToLocus[i];
    } else
    {
        assert(onlyTypeUsed == 51);
        return i;
    }
}


template<typename INT>
uint32_t GeneticMap::FromLocusToNextT1Locus(const INT i) const
{
    return FromLocusToNextTxLocus(i, isT1Used, _FromLocusToNextT1Locus, 1);
}

template<typename INT>
uint32_t GeneticMap::FromLocusToNextT2Locus(const INT i) const
{
    return FromLocusToNextTxLocus(i, isT2Used, _FromLocusToNextT2Locus, 2);
}

template<typename INT>
uint32_t GeneticMap::FromLocusToNextT3Locus(const INT i) const
{
    return FromLocusToNextTxLocus(i, isT3Used, _FromLocusToNextT3Locus, 3);
}

template<typename INT>
uint32_t GeneticMap::FromLocusToNextT4Locus(const INT i) const
{
    return FromLocusToNextTxLocus(i, isT4Used, _FromLocusToNextT4Locus, 4);
}

template<typename INT>
uint32_t GeneticMap::FromLocusToNextT56Locus(const INT i) const
{
    return FromLocusToNextTxLocus(i, isT56Used, _FromLocusToNextT56Locus, 5);
}

template<typename INT>
uint32_t GeneticMap::FromLocusToNextT56ntrlLocus(const INT i) const
{
    return FromLocusToNextTxLocus(i, isT56ntrlUsed, _FromLocusToNextT56ntrlLocus, 50);
}

template<typename INT>
uint32_t GeneticMap::FromLocusToNextT56selLocus(const INT i) const
{
    return FromLocusToNextTxLocus(i, isT56selUsed, _FromLocusToNextT56selLocus, 51);
}

template<typename INT>
uint32_t GeneticMap::FromLocusToNextTxLocus(const INT i, const bool isTxUsed, const std::vector<uint32_t>& map, const unsigned char typeValue) const
{
    if (isTxUsed)
    {
        if (onlyTypeUsed == 250)
        {
            assert(map.size() > i);
            return map[i];       
        } else
        {
            assert(onlyTypeUsed == typeValue);
            return i+1;
        }
    } else
    {
        return 0;
    }
}


template<typename INT>
T56LocusGender GeneticMap::FromT56LocusToT56genderLocus(const INT i) const
{
    if (isT56selUsed && isT56ntrlUsed)
    {
        assert(_FromT56LocusToT56genderLocus.size() > i);
        return _FromT56LocusToT56genderLocus[i];
    } else
    {
        assert (isT56selUsed != isT56ntrlUsed);
        if (isT56ntrlUsed)
        {
            return {true, (uint32_t) i};
        } else
        {
            assert(isT56selUsed);
            return {false, (uint32_t) i};
        }
    }
}


template<typename INT>
bool GeneticMap::isT56neutral(const INT i) const
{
    return FromT56LocusToT56genderLocus(i).isNtrl;
}


template<typename INT>
unsigned char GeneticMap::getLocusType(const INT i) const
{
    assert(i < TotalNbLoci);
    if (onlyTypeUsed != 250)
    {
        return onlyTypeUsed;
    } else
    {
        if (i == 0)
        {
            if (isT1Used && FromLocusToNextT1Locus(0) != 0) return 1;
            if (isT2Used && FromLocusToNextT2Locus(0) != 0) return 2;
            if (isT3Used && FromLocusToNextT3Locus(0) != 0) return 3;
            if (isT4Used && FromLocusToNextT4Locus(0) != 0) return 4;
            if (isT56ntrlUsed && FromLocusToNextT56ntrlLocus(0) != 0) return 50;
            if (isT56selUsed && FromLocusToNextT56selLocus(0) != 0) return 51;
            std::cout << "internal error in 'unsigned char GeneticMap::getLocusType(const INT i) const', where i = " << i << ". Could not find locus type.\n";
            abort();
        } else
        {
            if (isT1Used && FromLocusToNextT1Locus(i) != FromLocusToNextT1Locus(i-1)) return 1;
            if (isT2Used && FromLocusToNextT2Locus(i) != FromLocusToNextT2Locus(i-1)) return 2;
            if (isT3Used && FromLocusToNextT3Locus(i) != FromLocusToNextT3Locus(i-1)) return 3;
            if (isT4Used && FromLocusToNextT4Locus(i) != FromLocusToNextT4Locus(i-1)) return 4;
            if (isT56ntrlUsed && FromLocusToNextT56ntrlLocus(i) != FromLocusToNextT56ntrlLocus(i-1)) return 50;
            if (isT56selUsed && FromLocusToNextT56selLocus(i) != FromLocusToNextT56selLocus(i-1)) return 51;
            std::cout << "internal error in 'unsigned char GeneticMap::getLocusType(const INT i) const', where i = " << i << ". Could not find locus type.\n";
            abort();
        }
    }
}


template<typename INT>
LocusDescription GeneticMap::getLocusTypeAndItsIndex(const INT i) const
{
    assert(i < TotalNbLoci);
    if (onlyTypeUsed != 250)
    {
        return {onlyTypeUsed, i};
    } else
    {
        if (i == 0)
        {
            if (isT1Used && FromLocusToNextT1Locus(0) != 0) return {1, 0};
            if (isT2Used && FromLocusToNextT2Locus(0) != 0) return {2, 0};
            if (isT3Used && FromLocusToNextT3Locus(0) != 0) return {3, 0};
            if (isT4Used && FromLocusToNextT4Locus(0) != 0) return {4, 0};
            if (isT56ntrlUsed && FromLocusToNextT56ntrlLocus(0) != 0) return {50, 0};
            if (isT56selUsed && FromLocusToNextT56selLocus(0) != 0) return {51, 0};
            std::cout << "internal error in 'LocusDescription GeneticMap::getLocusTypeAndItsIndex(const INT i) const', where i = " << i << ". Could not find locus type.\n";
            abort();
        } else
        {
            {auto l = FromLocusToNextT1Locus(i-1); if (isT1Used && FromLocusToNextT1Locus(i) != l) return {1,l};}
            {auto l = FromLocusToNextT2Locus(i-1); if (isT2Used && FromLocusToNextT2Locus(i) != l) return {2,l};}
            {auto l = FromLocusToNextT3Locus(i-1); if (isT3Used && FromLocusToNextT3Locus(i) != l) return {3,l};}
            {auto l = FromLocusToNextT4Locus(i-1); if (isT4Used && FromLocusToNextT4Locus(i) != l) return {4,l};}
            {auto l = FromLocusToNextT56ntrlLocus(i-1); if (isT56ntrlUsed && FromLocusToNextT56ntrlLocus(i) != l) return {50,l};}
            {auto l = FromLocusToNextT56selLocus(i-1); if (isT56selUsed && FromLocusToNextT56selLocus(i) != l) return {51,l};}
            std::cout << "internal error in 'LocusDescription GeneticMap::getLocusTypeAndItsIndex(const INT i) const', where i = " << i << ". Could not find locus type.\n";
            abort();
        }
    }
}
