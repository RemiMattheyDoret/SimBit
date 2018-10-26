/*

 Author: Remi Matthey-Doret

    MIT License

    Copyright (c) 2017 Remi Matthey-Doret

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRAreadLineNTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

 */

std::vector<T1_locusDescription> OutputFile::getSubsetToConsider(int speciesIndex)
{
    return this->subset[speciesIndex];
}

template<typename T>
std::vector<T> OutputFile::removeSitesWeDontWant(std::vector<T> sites, int speciesIndex)
{
    // sites must be sorted
    return intersection(this->subset[speciesIndex],sites);
}

const std::vector<std::string> OutputFile::OutputFileTypesNames = 
{
    "Oops... failed to find file type name",
    "Logfile",
    "T1_vcfFile",
    "T1_LargeOutputFile",
    "T1_AlleleFreqFile",
    "MeanLDFile",
    "LongestRunFile",
    "HybridIndexFile",
    "ExpectiMinRecFile", 
    "T2_LargeOutputFile",
    "SaveBinaryFile",
    "T3_MeanVarFile",
    "T3_LargeOutputFile",
    "fitness",
    "fitnessStats",
    "T1_FST",
    "extraGeneticInfo",
    "patchSize",
    "extinction",
    "genealogy",
    "fitnessSubsetLoci"
};    

void OutputFile::openAndReadLine(std::string& line, int generation)
{
    std::ifstream ifs(this->getPath(generation));
    std::getline(ifs,line);
}

void OutputFile::mergeFiles(std::vector<std::string> listOfFiles)
{
    // https://stackoverflow.com/questions/19564450/concatenate-two-huge-files-in-c
    std::string path = this->getPathWithoutGenerationDespiteBeingGenerationSpecific();

    remove(path.c_str());
    ofs.open(path, std::ios_base::app);
    if (!ofs.is_open())
    {
        std::cout << "OutputFileType " << getFileTypeName(OutputFileType) << "(index "<< OutputFileType <<")" << " with path '" << path << "' failed to open in mergeFiles!" << "\n";
        abort();
    }

    for (std::string p : listOfFiles)
    {
        std::ifstream ifs(p);
        ofs << ifs.rdbuf();
    }
    ofs << "\n";

    /*for (std::string p : listOfFiles)
    {
        std::cout << "p = " << p << std::endl;
        std::ifstream ifs(p);
        std::string line;
        std::getline(ifs, line);
        std::cout << "line = "<< line << std::endl;
        ofs << line << "\n";
    }*/

}

bool OutputFile::isLocusInSubset(T1_locusDescription L, int speciesIndex)
{
    bool x = std::binary_search( 
                subset[speciesIndex].begin(),
                subset[speciesIndex].end(),
                L,
                [](const T1_locusDescription& right, const T1_locusDescription& left){return right.locus < left.locus;}
            );
    return x;

    // 'operator==' between 'T1_locusDescription' and 'int' compares with the attribute 'locus'
}

bool OutputFile::isLocusInSubset(int locus, int speciesIndex)
{
    T1_locusDescription L(locus); // minimalist instantiation (no division by eight) used only to search in 'subset'
    return isLocusInSubset(L, speciesIndex);
    // 'operator==' between 'T1_locusDescription' and 'int' compares with the attribute 'locus'
}

std::string OutputFile::getPathForSeed()
{
    return GeneralPath + filename + "_seed" + std::string("_G") + std::to_string(GP->CurrentGeneration) + sequencingErrorStringToAddToFilnames + extension;
}

std::string OutputFile::getPath()
{
    if (this->isSpeciesSpecific)
    {
        assert(SSP != nullptr);
        if (this->isGenerationSpecific)
        {
            // This is the complete path!
            return GeneralPath + filename + "_" + SSP->speciesName + std::string("_G") + std::to_string(GP->CurrentGeneration) + sequencingErrorStringToAddToFilnames + extension;
        } else
        {
            return GeneralPath + filename + "_" + SSP->speciesName + sequencingErrorStringToAddToFilnames + extension;
        }
        
    } else
    {
        if (this->isGenerationSpecific)
        {
            return GeneralPath + filename + std::string("_G") + std::to_string(GP->CurrentGeneration) + sequencingErrorStringToAddToFilnames + extension;
        } else
        {
            return GeneralPath + filename + sequencingErrorStringToAddToFilnames + extension;
        }
    }
}

std::string OutputFile::getPathWithoutGenerationDespiteBeingGenerationSpecific()
{
    if (this->isSpeciesSpecific)
    {
        assert(SSP != nullptr);
        
        return GeneralPath + filename + "_" + SSP->speciesName + sequencingErrorStringToAddToFilnames + extension;
        
    } else
    {
        return GeneralPath + filename + sequencingErrorStringToAddToFilnames + extension;
    }
}

std::string OutputFile::getPath(int generation)
{
    assert(this->isGenerationSpecific);
    if (this->isSpeciesSpecific)
    {
        assert(SSP != nullptr);

        return GeneralPath + filename + "_" + SSP->speciesName + std::string("_G") + std::to_string(generation) + sequencingErrorStringToAddToFilnames + extension;
        
    } else
    {
        return GeneralPath + filename + std::string("_G") + std::to_string(generation) + sequencingErrorStringToAddToFilnames + extension;
    }
}


void OutputFile::interpretTimeAndSubsetInput(InputReader& input)
{
    if (!input.IsThereMoreToRead())
    {
        std::cout << input.GetErrorMessage() << "was expecting time information after file name!\n";
        abort();
    }    

    std::vector<int> v;

    while (input.IsThereMoreToRead() || input.PeakNextElementString() != "subset")
    {
        if (input.PeakNextElementString() == "FromToBy" || input.PeakNextElementString() == "fromtoby" || input.PeakNextElementString() == "fromToBy")
        {
            input.skipElement();
            int from = input.GetNextElementInt();
            int to = input.GetNextElementInt();
            int by = input.GetNextElementInt();


            for (int t = from ; t <= to ; t+=by)
            {
                v.push_back(t);
            }
        } else
        {
            v.push_back(input.GetNextElementInt());
        }
    }
    this->setTimes(v);

    if (input.IsThereMoreToRead())
    {
        std::vector<int> listOfOutputFileTypeThatCanTakeASubset = {
            extraGeneticInfo,
            T1_FST,
            // not fitnessStats
            // not patchSize
            // not fitnessSubsetLoci (obviously)
            // not fitness (fitnessSubsetLoci is here for that)
            T1_AlleleFreqFile,
            MeanLDFile,
            LongestRunFile,
            // not T2_LargeOutputFile
            T1_LargeOutputFile,
            HybridIndexFile,
            ExpectiMinRecFile,
            T1_vcfFile,
            // not T3_LargeOutputFile
            // T3_MeanVarFile
            // extinction
            // not genealogy
            // not SaveBinaryFile
            // not Logfile
        };
        assert(input.GetNextElementString() != "subset");

        if (std::find(listOfOutputFileTypeThatCanTakeASubset.begin(), listOfOutputFileTypeThatCanTakeASubset.end(), this->OutputFileType) != listOfOutputFileTypeThatCanTakeASubset.end())
        {
            std::cout << "Received the 'subset' keyword for outputFile of type " << getFileTypeName(OutputFileType) << " but only the types {";
            for (auto& elem : listOfOutputFileTypeThatCanTakeASubset)
                std::cout << getFileTypeName(elem) << " ";
            std::cout << "}\n";
            abort();
        }

        input.removeAlreadyRead();
        std::vector<std::pair<int, int>> rangesToSubsetFullInput = input.GetRangeOfIndicesForEachSpecies();
        assert(rangesToSubsetFullInput.size() == GP->nbSpecies);
        auto SSP_toReset = SSP;
        for (int speciesIndex = 0 ; speciesIndex < GP->nbSpecies ; speciesIndex++)
        {
            SSP = allParameters.getSSPaddress(speciesIndex); // It is important to reset this because InputReader uses SSP
            assert(SSP != nullptr);

            int from = rangesToSubsetFullInput[speciesIndex].first;
            int to = rangesToSubsetFullInput[speciesIndex].second;
            assert(from < to);
            InputReader inputOneSpecies(input, from, to, speciesIndex);


            while (input.IsThereMoreToRead())
            {
                int locus = inputOneSpecies.GetNextElementInt();
                subset.push_back({locus%8, locus/8, locus});
            }            

        }
        SSP = SSP_toReset;

    } else
    {
        // must list all loci in subset
        std::vector<T1_locusDescription> oneSpeciesSubset;
        for (int speciesIndex = 0 ; speciesIndex < GP->nbSpecies; speciesIndex++)
        {
            for (int locus = 0 ; locus < allParameters.getSSPaddress(speciesIndex)->T1_nbBits; locus++)
            {
                oneSpeciesSubset.push_back({locus%8, locus/8, locus}); 
            }
            subset.push_back(oneSpeciesSubset);
        }
    }
}

bool OutputFile::isEmpty(std::ifstream& pFile)
{
    return pFile.peek() == std::ifstream::traits_type::eof();
}

bool OutputFile::containsRightNumberOfLines(std::ifstream& pFile)
{
    int number_of_lines = 0;
    std::string line;
    while (std::getline(pFile, line))
        ++number_of_lines;

    int nbExpectedLines;
    if (this->OutputFileType == extraGeneticInfo)
    {
        nbExpectedLines = 3;
    } else
    {
        nbExpectedLines = this->getTimes().size() + 1;
    }
    
    if (nbExpectedLines < number_of_lines)
    {
        if (this->OutputFileType != T1_vcfFile)
        {
            std::cout << "\tThe file of type " << getFileTypeName(OutputFileType) << " (index "<< OutputFileType <<")" << " seems to contain more lines than the simulation is expected to produce. Simulation should produce " << nbExpectedLines << " lines (incl. header) and the current file contains " << number_of_lines << " lines." << "\n";
        }
    }
    return nbExpectedLines == number_of_lines;
}

bool OutputFile::DoesFileExist()
{
    // check file exists and contain at least one line
    std::string path = this->getPath();
    std::ifstream f(path);
    if (f.good())
    {
        if (this->containsRightNumberOfLines(f))
        {
            return true;
        }
    }
    return false;
}

bool OutputFile::DoesFileExist(int generation)
{
    // check file exists and contain at least one line
    bool ret = false;
    int ForReset = GP->CurrentGeneration;
    GP->CurrentGeneration = generation;
    std::string path = this->getPath();
    std::ifstream f(path);
    //std::cout << "path = " << path << "\n";
    if (f.good())
    {
        if (this->containsRightNumberOfLines(f))
        {
            ret = true;
        }
    }

    GP->CurrentGeneration = ForReset;
    return ret;

}

bool OutputFile::DoAllFilesOfTypeAlreadyExist()
{
    if (isGenerationSpecific)
    {
        std::vector<int>& filetimes = getTimes();
        for ( int& generation : filetimes)
        {
            if (!this->DoesFileExist(generation))
            {
                return false;
            }
        }
        return true;
    } else
    {
        return this->DoesFileExist();
    }
}

bool OutputFile::DoesAtLeastOneFileOfTypeAlreadyExist()
{
    if (isGenerationSpecific)
    {
        std::vector<int>& filetimes = getTimes();
        for ( int& generation : filetimes)
        {
            if (this->DoesFileExist(generation))
            {
                return true;
            }
        }
        return false;
    } else
    {
        return this->DoesFileExist();
    }
}

std::vector<int>& OutputFile::getTimes()
{
    return times;
}





bool OutputFile::isTime()
{
    if (!doesTimeNeedsToBeSet)
    {
        std::cout << "Internal error: Tried to get times for a file typeDoesTimeNeedsToBeSet who does not need time to be set\n";
        abort();
    }
    return std::find(times.begin(), times.end(),  GP->CurrentGeneration) != times.end();
}

OutputFile::OutputFile(OutputFile&& f)
: filename(std::move(f.filename)), extension(f.extension), OutputFileType(f.OutputFileType), times(f.times), isGenerationSpecific(f.isGenerationSpecific), isSpeciesSpecific(f.isSpeciesSpecific), isNbLinesEqualNbOutputTimes(f.isNbLinesEqualNbOutputTimes), doesTimeNeedsToBeSet(f.doesTimeNeedsToBeSet)
{}

OutputFile::OutputFile(std::string f, OutputFileTypes t)
:filename(f), OutputFileType(t)
{
    if (filename == "NFN" || filename == "nfn")
    {
        filename = "";
    }

    if (t == Logfile)
    {
        this->extension = std::string(".log");
        isGenerationSpecific = false;
        isSpeciesSpecific = false;
        isNbLinesEqualNbOutputTimes = false;
        doesTimeNeedsToBeSet = false;
    } else if (t == T1_vcfFile)
    {
        this->extension = std::string(".T1vcf");
        isGenerationSpecific = true;
        isSpeciesSpecific = true;
        isNbLinesEqualNbOutputTimes = false;
        doesTimeNeedsToBeSet = true;
    } else if (t == T1_LargeOutputFile)
    {
        this->extension = std::string(".T1LO");
        isGenerationSpecific = false;
        isSpeciesSpecific = true;
        isNbLinesEqualNbOutputTimes = true;
        doesTimeNeedsToBeSet = true;
    } else if (t == T1_AlleleFreqFile)
    {
        this->extension = std::string(".T1AllFreq");
        isGenerationSpecific = false;
        isSpeciesSpecific = true;
        isNbLinesEqualNbOutputTimes = true;
        doesTimeNeedsToBeSet = true;
    } else if (t == MeanLDFile)
    {
        this->extension = std::string(".T1LD");
        isGenerationSpecific = false;
        isSpeciesSpecific = true;
        isNbLinesEqualNbOutputTimes = true;
        doesTimeNeedsToBeSet = true;
    } else if (t == LongestRunFile)
    {
        this->extension = std::string(".T1LR");
        isGenerationSpecific = false;
        isSpeciesSpecific = true;
        isNbLinesEqualNbOutputTimes = true;
        doesTimeNeedsToBeSet = true;
    } else if (t == HybridIndexFile)
    {
        this->extension = std::string(".T1HI");
        isGenerationSpecific = false;
        isSpeciesSpecific = true;
        isNbLinesEqualNbOutputTimes = true;
        doesTimeNeedsToBeSet = true;
    } else if (t == ExpectiMinRecFile)
    {
        this->extension = std::string(".T1EMR");
        isGenerationSpecific = false;
        isSpeciesSpecific = true;
        isNbLinesEqualNbOutputTimes = true;
        doesTimeNeedsToBeSet = true;
    } else if (t == T2_LargeOutputFile)
    {
        this->extension = std::string(".T2LO");
        isGenerationSpecific = false;
        isSpeciesSpecific = true;
        isNbLinesEqualNbOutputTimes = true;
        doesTimeNeedsToBeSet = true;
    } else if (t == SaveBinaryFile)
    {
        this->extension = std::string(".bin");
        isGenerationSpecific = true;
        isSpeciesSpecific = true;
        isNbLinesEqualNbOutputTimes = false;
        doesTimeNeedsToBeSet = true;
    } else if (t == T3_MeanVarFile)
    {
        this->extension = std::string(".T3MeanVar");
        isGenerationSpecific = false;
        isSpeciesSpecific = true;
        isNbLinesEqualNbOutputTimes = true;
        doesTimeNeedsToBeSet = true;
    } else if (t == T3_LargeOutputFile)
    {
        this->extension = std::string(".T3LO");
        isGenerationSpecific = false;
        isSpeciesSpecific = true;
        isNbLinesEqualNbOutputTimes = true;
        doesTimeNeedsToBeSet = true;
    } else if (t == fitness)
    {
        this->extension = std::string(".fit");
        isGenerationSpecific = false;
        isSpeciesSpecific = true;
        isNbLinesEqualNbOutputTimes = true;
        doesTimeNeedsToBeSet = true;
    } else if (t == fitnessStats)
    {
        this->extension = std::string(".fitStats");
        isGenerationSpecific = false;
        isSpeciesSpecific = true;
        isNbLinesEqualNbOutputTimes = true;
        doesTimeNeedsToBeSet = true;
    } else if (t == T1_FST)
    {
        this->extension = std::string(".T1_FST");
        isGenerationSpecific = false;
        isSpeciesSpecific = true;
        isNbLinesEqualNbOutputTimes = true;
        doesTimeNeedsToBeSet = true;
    } else if (t == extraGeneticInfo)
    {
        this->extension = std::string(".extraGeneticInfo");
        isGenerationSpecific = false;
        isSpeciesSpecific = true;
        isNbLinesEqualNbOutputTimes = false;
        doesTimeNeedsToBeSet = false;
    }  else if (t == patchSize)
    {
        this->extension = std::string(".patchSize");
        isGenerationSpecific = false;
        isSpeciesSpecific = true;
        isNbLinesEqualNbOutputTimes = true;
        doesTimeNeedsToBeSet = true;
    } else if (t == extinction)
    {
        this->extension = std::string(".extinction");
        isGenerationSpecific = false;
        isSpeciesSpecific = false;
        isNbLinesEqualNbOutputTimes = false;
        doesTimeNeedsToBeSet = false;
    } else if (t == genealogy)
    {
        this->extension = std::string(".genealogy");
        isGenerationSpecific = true;   // Actually the end result is not. SimBit uses getPathWithoutGenerationDespiteBeingGenerationSpecific to produce the last file.
        isSpeciesSpecific = true;
        isNbLinesEqualNbOutputTimes = false;
        doesTimeNeedsToBeSet = true;
    } else if (t == fitnessSubsetLoci)
    {
        this->extension = std::string(".fitSubLoci");
        isGenerationSpecific = false;
        isSpeciesSpecific = true;
        isNbLinesEqualNbOutputTimes = true;
        doesTimeNeedsToBeSet = true;
    } else
    {
        std::cout << "Internal Error: In class 'OutputFile' in 'set_path', unknown fileType\n";
        abort();
    }
}

void OutputFile::openForSeed()
{
    assert(OutputFileType == SaveBinaryFile);
    std::string path = this->getPathForSeed();
    ofs.open(path, std::ios::out | std::ios::binary);

    if (!ofs.is_open())
    {
        std::cout << "When trying to print the random seed to a binary file, the path '" << path << "' failed to open!" << "\n";
        abort();
    }    
}

void OutputFile::open()
{   
    std::string path = this->getPath();
    
    if (OutputFileType == SaveBinaryFile)
    {
        ofs.open(path, std::ios::out | std::ios::binary);
    } else
    {
        ofs.open(path, std::ios_base::app);
    }

    if (!ofs.is_open())
    {
        std::cout << "OutputFileType " << getFileTypeName(OutputFileType) << "(index "<< OutputFileType <<")" << " with path '" << path << "' failed to open!" << "\n";
        abort();
    }
}

void OutputFile::openWithoutGenerationDespiteBeingGenerationSpecific()
{
    std::string path = this->getPathWithoutGenerationDespiteBeingGenerationSpecific();
    
    ofs.open(path, std::ios_base::app);
    if (!ofs.is_open())
    {
        std::cout << "OutputFileType " << getFileTypeName(OutputFileType) << "(index "<< OutputFileType <<")" << " with path '" << path << "' failed to open!" << "\n";
        abort();
    }
}

void OutputFile::open(int generation)
{
    std::string path = this->getPath(generation);
    
    ofs.open(path, std::ios_base::app);
    if (!ofs.is_open())
    {
        std::cout << "OutputFileType " << getFileTypeName(OutputFileType) << "(index "<< OutputFileType <<")" << " with path '" << path << "' failed to open!" << "\n";
        abort();
    }
}

void OutputFile::clearContent()
{
    std::string path;
    
    path = this->getPath();
    
    
    ofs.open(path, std::ofstream::out | std::ofstream::trunc);
    if (!ofs.is_open())
    {
        std::cout << "OutputFileType " << getFileTypeName(OutputFileType) << "(index "<< OutputFileType <<")" << " with path '" << path << "' failed to open!" << "\n";
        abort();
    }
    ofs.close();
}

void OutputFile::clearContentAndLeaveOpen(int generation)
{
    std::string path;
    
    path = this->getPath(generation);
    
    ofs.open(path, std::ofstream::out | std::ofstream::trunc);
    if (!ofs.is_open())
    {
        std::cout << "OutputFileType " << getFileTypeName(OutputFileType) << "(index "<< OutputFileType <<")" << " with path '" << path << "' failed to open!" << "\n";
        abort();
    }
}

bool OutputFile::isOpen()
{
    return ofs.is_open();
}

void OutputFile::write(const std::string& s)
{
    assert(this->isOpen());
    ofs << s;
}

void OutputFile::writeBinary(const char* first, int second)
{
    assert(this->isOpen());
    ofs.write(first, second);
}

void OutputFile::writeBinary(const std::mt19937 x)
{
    assert(this->isOpen());
    ofs << x;   
}

void OutputFile::writeBinary(const int x)
{
    assert(this->isOpen());
    ofs << x;   
}

void OutputFile::close()
{
    ofs.close();
    assert(!this->isOpen());
}

OutputFileTypes OutputFile::getFileType()
{
    return OutputFileType;
}

void OutputFile::setTimes(std::vector<int> x)
{
    // assertions
    for (int i = 0 ; i < x.size(); i++)
    {
        if (i!= 0)
        {
            if (x[i-1] >= x[i])
            {
                std::cout << "For file type '" << getFileTypeName(OutputFileType) << "(index "<< OutputFileType <<")" << "', received a time value that is larger or equal to the previous one. (The "<< i << "th value received is " << x[i] << " while the "<< i-1 <<"th value received is "<< x[i-1] <<"). Please indicate strictly increasing values for the '..._time' options.\n";
                abort();       
            }
        }
        if (x[i] < 0)
        {
            std::cout << "For file type '" << getFileTypeName(OutputFileType) << "(index "<< OutputFileType <<")" << "', received a time value that is negative (Received " << x[i] << ") Note that a generation of 0 refers to sampling before any evolutionary process has happened.\n";
            abort();
        }
        if (x[i] > GP->nbGenerations)
        {
            std::cout << "For file type '" << getFileTypeName(OutputFileType) << "(index "<< OutputFileType <<")" << "', received a time value that is greater than the number of generations (Received " << x[i] << " while nbGenerations = "<< GP->nbGenerations <<"). Note that a generation of 0 refers to sampling before any evolutionary process has happened.\n";
            abort();
        }
    }

    if (OutputFileType == 1 || OutputFileType == 16 || OutputFileType == 18)
    {
        std::cout << "Internal error: Attempt to set times for outputFileType = " << getFileTypeName(OutputFileType) << "(index "<< OutputFileType <<")" << "\n";
        abort();
    }

    
    x.erase(std::remove_if(
        x.begin(), x.end(),
        [](const int& generation)
        { 
            return generation < GP->startAtGeneration;
        }
    ), x.end());
    

    // set times
    this->times = x;
}

bool OutputFile::getDoesTimeNeedsToBeSet()
{
    return doesTimeNeedsToBeSet;
}

std::string OutputFile::getFileTypeName(int fileTypeIndex)
{
    if (fileTypeIndex > OutputFileTypesNames.size() || fileTypeIndex == 0)
    {
        return "Oops... failed to find file type name for file type index " + std::to_string(fileTypeIndex);
    } else
    {
        return OutputFileTypesNames[fileTypeIndex];
    }
}
