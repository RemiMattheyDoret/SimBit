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
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.

 */


typedef enum {
    Logfile = 1,
    T1_vcfFile= 2,
    T1_LargeOutputFile = 3,
    T1_AlleleFreqFile = 4,
    MeanLDFile = 5,
    LongestRunFile = 6,
    HybridIndexFile = 7,
    ExpectiMinRecFile = 8, 
    T2_LargeOutputFile = 9,
    SaveBinaryFile = 10,
    T3_MeanVarFile = 11,
    T3_LargeOutputFile = 12,
    fitness = 13,
    fitnessStats = 14,
    T1_FST = 15,
    extraGeneticInfo = 16,
    patchSize = 17,
    extinction = 18,
    genealogy = 19,
    fitnessSubsetLoci = 20,
    T4_LargeOutputFile = 21,
    T4_vcfFile = 22,
    T4_SFS_file = 23,
    T1_SFS_file = 24,
    T4_printTree = 25, // special type that is not going in outputWriter
    T56_vcfFile= 26,
    T56_SFS_file = 27,
    T56_AlleleFreqFile = 28,
    T56_LargeOutputFile = 29,
    T4CoalescenceFst = 30
} OutputFileTypes; 
/* When creating new types I must not forget to 
    1) tell whether the file name should be time specific
    2) give the extension in 'OutputFile::OutputFile'
    3) write the functions to output to this file
    4) complete OutputFileTypesNames
*/


class OutputFile
{
private:
    std::ofstream ofs;
    std::string filename;
    std::string extension;
    OutputFileTypes OutputFileType;
    std::vector<int> times;
    bool isGenerationSpecific;
    bool isSpeciesSpecific;
    bool isNbLinesEqualNbOutputTimes;
    bool doesTimeNeedsToBeSet;
    std::vector<std::vector<T1_locusDescription>> subset;

    static const std::vector<int> listOfOutputFileTypeThatCanTakeASubset;
    static const std::vector<std::string> OutputFileTypesNames; // initialized in .cpp

public:
    static std::string GeneralPath;
    static std::string sequencingErrorStringToAddToFilnames;

    void assertSubsetSize();
    OutputFile(OutputFile&& f); // move constructor
    OutputFile(std::string f, OutputFileTypes t);
    bool containsRightNumberOfLines(std::ifstream& pFile);
    bool getDoesTimeNeedsToBeSet();
    void open();
    void openForSeed();
    void open(int generation);
    void openWithoutGenerationDespiteBeingGenerationSpecific();
    bool isOpen();
    void write(const std::string& s);
    void writeBinary(const RNG_type x);
    void writeBinary(const int x);
    void writeBinary(const char* first, int second);
    //template<typename T> void writeBinary(const T* first, int second);
    void close();
    void setTimes(std::vector<int> x);
    bool isTime();
    void clearContent();
    void clearContentAndLeaveOpen(int generation);
    bool isEmpty(std::ifstream& pFile);
    bool DoesFileExist();
    bool DoesFileExist(int generation);
    bool DoesAtLeastOneFileOfTypeAlreadyExist();
    bool DoAllFilesOfTypeAlreadyExist();
    void interpretTimeAndSubsetInput(InputReader& input);
    void interpretSubsetInput(InputReader& input);
    std::string getPathForSeed();
    std::string getPath();
    std::string getPath(int generation);
    std::string getPathWithoutGenerationDespiteBeingGenerationSpecific();
    std::vector<int>& getTimes();
    OutputFileTypes getFileType();
    void openAndReadLine(std::string& line, int generation);
    void mergeFiles(std::vector<std::string> listOfFiles);
    std::string getFileTypeName(int fileTypeIndex);
    template<typename T> std::vector<T> removeSitesWeDontWant(std::vector<T> sites, int speciesIndex);
    std::vector<T1_locusDescription> getSubsetToConsider(int speciesIndex);
    bool isLocusInSubset(T1_locusDescription L, int speciesIndex);
    bool isLocusInSubset(int locus, int speciesIndex);

};

