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
    void ExtendStringForAdvancedLogFile(std::string& s, T& entry, int depth, bool willMoreCome);
    void ExtendStringForAdvancedLogFile(std::string& s, T1_locusDescription& entry, int depth, bool willMoreCome);
    template<typename T>
    void ExtendStringForAdvancedLogFile(std::string& s, std::vector<T>& entry, int depth, bool willMoreCome);
    template<typename T>
    void ExtendStringForAdvancedLogFile(std::string& s, T& entry, std::string name);
    void ExtendStringForAdvancedLogFile(std::string& s, std::string entry, int depth, bool willMoreCome);

public:
    // Will contain all the T1_loci that must be outputted
    int                                 LogfileType;


    //Print LogFile
    void FinalizeLogfile();
    void PrintLogfile(std::string& AllInputInLongString);
    

    bool ShouldSimulationBeDone(int OverwriteMode);

    // Clear All files content
    void ClearAllFileContent();

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
    void WriteOutputs_T1_MeanLD_header(OutputFile& file);
    void WriteOutputs_T1_MeanLD(Pop& pop, OutputFile& file);
    void WriteOutputs_T1_LongestRun_header(OutputFile& file);
    void WriteOutputs_T1_LongestRun(Pop& pop, OutputFile& file);
    void WriteOutputs_T2_LargeOutput_header(OutputFile& file);
    void WriteOutputs_T2_LargeOutput(Pop& pop, OutputFile& file);
    void WriteOutputs_T1_LargeOutput_header(OutputFile& file);
    void WriteOutputs_T1_LargeOutput(Pop& pop, OutputFile& file);
    void WriteOutputs_T1_HybridIndex_header(OutputFile& file);
    void WriteOutputs_T1_HybridIndex(Pop& pop, OutputFile& file);
    void WriteOutputs_T1_ExpectiMinRec_header(OutputFile& file);
    void WriteOutputs_T1_ExpectiMinRec(Pop& pop, OutputFile& file);
    void WriteOutputs_T1_vcf(Pop& pop, OutputFile& file);
    void WriteOutputs_T3_LargeOutput_header(OutputFile& file);
    void WriteOutputs_T3_LargeOutput(Pop& pop, OutputFile& file);
    void WriteOutputs_T3_MeanVar_header(OutputFile& file);
    void WriteOutputs_T3_MeanVar(Pop& pop, OutputFile& file);
    void WriteOutputs_extraGeneticInfo(OutputFile& file);
    void WriteOutputs_T1_FST_header(OutputFile& file);
    void WriteOutputs_T1_FST(Pop& pop, OutputFile& file);
    void WriteOutputs_extinction(OutputFile& file);
    void WriteOutputs(Pop& realPop);
    void WriteOutputs_forDefinedPop(Pop& pop);
    void imitateSequencingError(Pop& pop);
    void imitateSequencingError(Haplotype& TransmittedChrom);
    void PrintGeneration();

    bool shouldNABePrinted(int patch_index);
    bool shouldNABePrinted(int patch_index, int ind_index);

};

