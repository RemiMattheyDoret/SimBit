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
OutputFileTypes
 */

class OutputWriter
{
private:
    // attributes
    std::map<const OutputFileTypes, std::vector<OutputFile>> TypesToPrintOn;

    std::vector<int> AllTimes;
    
    std::clock_t clock_BeforeInitialization;
    std::clock_t clock_BeforeSimulation;
    std::clock_t clock_AfterSimulation;

private:
    template<typename T>
    void ExtendStringForAdvancedLogFile(std::string& s, T& entry, std::string name);
    void ExtendStringForAdvancedLogFile(std::string& s, std::string& entry, std::string name);
    template<typename T>
    void ExtendStringForAdvancedLogFile(std::string& s, std::vector<T>& entry, std::string name);
    template<typename T>
    void ExtendStringForAdvancedLogFile(std::string& s, std::vector<std::vector<T>>& entry, std::string name);
    void ExtendStringForAdvancedLogFile(std::string& s, std::vector<std::vector<std::string>>& entry, std::string name);
    template<typename T>
    void ExtendStringForAdvancedLogFile(std::string& s, std::vector<std::vector<std::vector<T>>>& entry, std::string name);
    void ExtendStringForAdvancedLogFile(std::string& s, std::vector<std::vector<std::vector<T1_locusDescription>>>& entry, std::string name);
    void ExtendStringForAdvancedLogFile(std::string& s, std::vector<FromLocusToTXLocusElement>& entry, std::string name);
    void ExtendStringForAdvancedLogFile(std::string& s, std::vector<std::pair<bool, unsigned>>& entry, std::string name);



    /*template<typename T>
    void ExtendStringForAdvancedLogFile(std::string& s, T& entry, int depth, bool willMoreCome);
    void ExtendStringForAdvancedLogFile(std::string& s, T1_locusDescription& entry, int depth, bool willMoreCome);
    template<typename T>
    void ExtendStringForAdvancedLogFile(std::string& s, std::vector<T>& entry, int depth, bool willMoreCome);
    template<typename T>
    void ExtendStringForAdvancedLogFile(std::string& s, T& entry, std::string name);
    void ExtendStringForAdvancedLogFile(std::string& s, std::string entry, int depth, bool willMoreCome);*/

public:


    // Will contain all the T1_loci that must be outputted
    int                                 LogfileType;


    //Print LogFile
    void FinalizeLogfile();
    void PrintLogfile(std::string& AllInputInLongString);
    

    bool ShouldSimulationBeDone(int OverwriteMode);

    // Clear All files content
    void ClearAllFileContent();
    void RemoveAllFiles();

    // For warning and error messages
    bool AreThereAnyOutput();
    bool IsLastGenerationSampled();
    

    // Print time
    void printTime(long double seconds, std::string message);
    void SetBeforeInitializationTime();
    void SetBeforeSimulationTime();
    void PrintInitializationTimeToLogFile();
    void SetAfterSimulationTime();
    void PrintSimulationTimeToLogFile();

    // Should output be printed at this given generation (generation is in Parameter and is therefore global)
    bool isTime();

    void SetAllTimes();
    void insertOutputFile(OutputFile&& file);

    

    // Get OutputFile
    bool isFile(OutputFileTypes type);
    std::vector<OutputFile>& get_OutputFiles(OutputFileTypes type);

    // WriteOutputs
    void WriteOutputs_patchSize_header(OutputFile& file);
    void WriteOutputs_patchSize(OutputFile& file); // No need to give Pop here as patchSize info is in GP.
    void WriteOutputs_fitness_header(OutputFile& file);
    void WriteOutputs_fitness(Pop& pop, OutputFile& file);
    void WriteOutputs_fitnessSubsetLoci_header(OutputFile& file);
    void WriteOutputs_fitnessSubsetLoci(Pop& pop, OutputFile& file);
    void WriteOutputs_fitnessStats_header(OutputFile& file);
    void WriteOutputs_fitnessStats(Pop& pop, OutputFile& file);
    void WriteOutputs_T1_AlleleFreq_header(OutputFile& file);
    void WriteOutputs_T1_AlleleFreq(Pop& pop, OutputFile& file);
    void WriteOutputs_T56_AlleleFreq_header(OutputFile& file);
    void WriteOutputs_T56_AlleleFreq(Pop& pop, OutputFile& file);
    void WriteOutputs_T1_MeanLD_header(OutputFile& file);
    void WriteOutputs_T1_MeanLD(Pop& pop, OutputFile& file);
    void WriteOutputs_T1_LongestRun_header(OutputFile& file);
    void WriteOutputs_T1_LongestRun(Pop& pop, OutputFile& file);
    void WriteOutputs_T2_LargeOutput_header(OutputFile& file);
    void WriteOutputs_T2_LargeOutput(Pop& pop, OutputFile& file);
    void WriteOutputs_T1_LargeOutput_header(OutputFile& file);
    void WriteOutputs_T1_LargeOutput(Pop& pop, OutputFile& file);
    void WriteOutputs_T4_LargeOutput_header(OutputFile& file, size_t mutPlacingIndex);
    void WriteOutputs_T4_LargeOutput(OutputFile& file, std::vector<std::vector<std::vector<uint32_t>>>& data, size_t mutPlacingIndex);
    void WriteOutputs_T56_LargeOutput_header(OutputFile& file);
    void WriteOutputs_T56_LargeOutput(Pop& pop, OutputFile& file);
    void WriteOutputs_T4_paintedHaplo_file(Pop& pop, OutputFile& file);
    void WriteOutputs_T4_paintedHaploSegmentsDiversity_file(Pop& pop, OutputFile& file);
    void WriteOutputs_T4_paintedHaploSegmentsDiversity_file_header(OutputFile& file);
    void WriteOutputs_T4_SNPfreq_file(OutputFile& file, std::vector<T4TreeRec::PatchCountData>& T4freqs, size_t mutPlacingIndex);
    void WriteOutputs_T4_SNPfreq_file_header(OutputFile& file, size_t mutPlacingIndex);
    void WriteOutputs_Tx_SNPfreq_file_header(OutputFile& file, size_t mutPlacingIndex);
    void WriteOutputs_Tx_SNPfreq_file(OutputFile& file, std::vector<T4TreeRec::PatchCountData>& T4freqs, size_t mutPlacingIndex, Pop& pop);
    template<typename ntrlIterator, typename selIterator>
    void WriteOutputs_T56_LargeOutput_writeData(ntrlIterator ntrlIt, selIterator selIt, ntrlIterator ntrlItEnd, selIterator selItEnd, OutputFile& file, bool printNA);
    void WriteOutputs_T1_HybridIndex_header(OutputFile& file);
    void WriteOutputs_T1_HybridIndex(Pop& pop, OutputFile& file);
    void WriteOutputs_T1_AverageHybridIndex_header(OutputFile& file);
    void WriteOutputs_T1_AverageHybridIndex(Pop& pop, OutputFile& file);
    void WriteOutputs_T1_ExpectiMinRec_header(OutputFile& file);
    void WriteOutputs_T1_ExpectiMinRec(Pop& pop, OutputFile& file);
    void WriteOutputs_T1_vcf(Pop& pop, OutputFile& file);
    void WriteOutputs_T4_vcf(OutputFile& file, std::vector<std::vector<std::vector<uint32_t>>>& data, size_t mutPlacingIndex);
    void writeOutputs_burnInLength(int positiveGeneration);
    //void WriteOutputs_T4CoalescenceFst_header(OutputFile& file);
    //void WriteOutputs_T4CoalescenceFst(OutputFile& file);
    void WriteOutputs_T56_vcf(Pop& pop, OutputFile& file);
    //void WriteOutputs_T56_vcf_writeData(Pop& pop, std::vector<unsigned>& obsFreqs, OutputFile& file);
    //template<typename ntrlIteratorType, typename selIteratorType>
    //void write_T56vcf_forGroupOfLoci(std::vector<uint32_t>& groupOfLoci, std::vector<double>& relFreqsForGroupOfLoci, Pop& pop, ntrlIteratorType& ntrlIteratorTypeInfo, selIteratorType& selIteratorTypeInfo, OutputFile& file);
    void WriteOutputs_T3_LargeOutput_header(OutputFile& file);
    void WriteOutputs_T3_LargeOutput(Pop& pop, OutputFile& file);
    void WriteOutputs_T3_MeanVar_header(OutputFile& file);
    void WriteOutputs_T3_MeanVar(Pop& pop, OutputFile& file);
    void WriteOutputs_T1_haplotypeFreqs(Pop& pop, OutputFile& file);
    void WriteOutputs_extraGeneticInfo(OutputFile& file);
    void WriteOutputs_T1_FST_header(OutputFile& file);
    void WriteOutputs_Tx_FST_header(OutputFile& file);
    void WriteOutputs_T1_FST(Pop& pop, OutputFile& file);
    void WriteOutputs_T1_FST_header_complex(OutputFile& file);
    void WriteOutputs_T1_FST_complex(Pop& pop, OutputFile& file);
    void WriteOutputs_T1SFS(Pop& pop, OutputFile& file);
    void WriteOutputs_T1SFS_header(OutputFile& file);
    void WriteOutputs_T4SFS(OutputFile& file, std::vector<std::vector<std::vector<uint32_t>>>& data, size_t mutPlacingIndex);
    void WriteOutputs_T4SFS_header(OutputFile& file, size_t mutPlacingIndex);
    void WriteOutputs_T56SFS(Pop& pop, OutputFile& file);
    void WriteOutputs_T56SFS_header(OutputFile& file);
    void WriteOutputs_T1or4or56SFS_header(OutputFile& file, size_t mutPlacingIndex = 0);
    void WriteOutputs_T1or4or5SFS(std::vector<std::vector<double>>& obsFreqs, OutputFile& file, size_t mutPlacingIndex = 0);
    void WriteOutputs_extinction(OutputFile& file);
    template<typename ntrlit_t, typename selit_t>
    void WriteOutputs_sampleSequences_process(OutputFile& file, std::vector<uint32_t>& T4data, ntrlit_t it_ntrlBegin, selit_t it_selBegin, ntrlit_t it_ntrlEnd, selit_t it_selEnd, SampleSequenceData& SSD, Haplotype& haplo);
    void WriteOutputs_sampleSequences(std::vector<std::vector<std::vector<uint32_t>>>& T4data, Pop& pop, OutputFile& file, size_t mutPlacingIndex);

    void WriteOutputs_forDefinedPop_general(Pop& pop);
    void WriteOutputs_forDefinedPop_T1(Pop& pop);
    void WriteOutputs_forDefinedPop_T2(Pop& pop);
    void WriteOutputs_forDefinedPop_T3(Pop& pop);
    void WriteOutputs_forDefinedPop_T4AndSampledSeq(Pop& pop);
    void WriteOutputs_forDefinedPop_T56(Pop& pop);

    void WriteOutputs(Pop& realPop);
    void WriteOutputs_forDefinedPop(Pop& pop);
    void imitateSequencingError(Pop& pop);
    void imitateSequencingError(Haplotype& TransmittedChrom);
    void PrintGeneration();
    void PrintGenerationEndOfBurnIn();

    bool shouldNABePrinted(int patch_index);
    bool shouldNABePrinted(int patch_index, int ind_index);

    template<typename CountContainerType> void WriteOutputs_SNPfreq_computeFreq(CountContainerType& counts, Pop& pop, size_t patch_index);
    //void mergeT4countsWithOtherCounts(T4TreeRec::PatchCountData& T4freqs, std::vector<uint32_t>& counts);
    //void mergeT4countsWithOtherCounts(T4TreeRec::PatchCountData& T4freqs, std::map<uint32_t,uint32_t>& counts);

};

