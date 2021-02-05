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
    assert(speciesIndex < this->subset.size());
    return this->subset[speciesIndex];
}

template<typename T>
std::vector<T> OutputFile::removeSitesWeDontWant(std::vector<T> sites, int speciesIndex)
{
    // sites must be sorted
    //std::cout << "From OutputFile::removeSitesWeDontWant: speciesIndex = " << speciesIndex << "\n";
    //std::cout << "From OutputFile::removeSitesWeDontWant: this->subset.size() = " << this->subset.size() << "\n";
    assert(speciesIndex < this->subset.size());
    return intersection(this->subset[speciesIndex], sites);
}

const std::vector<std::string> OutputFile::OutputFileTypesNames = {
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
    "fitnessSubsetLoci",
    "T4_LargeOutputFile",
    "T4_vcfFile",
    "T4_SFS_file",
    "T1_SFS_file",
    "T4_printTree",
    "T5_vcfFile",
    "T5_SFS_file",
    "T5_AlleleFreqFile",
    "T5_LargeOutputFile",
    "T4CoalescenceFst",
    "T1_AverageHybridIndexFile",
    "T1_haplotypeFreqs_file",
    "sampleSeq_file"
}; 

const std::vector<int> OutputFile::listOfOutputFileTypeThatCanTakeASubset = {
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
    // not T4_LargeOutputFile,
    // not T4_vcfFile
    // not T4_SFS_file
    T1_SFS_file,
    // not T4_printTree
    // not T5_vcfFile
    // not T5_SFS_file
    // not T5_AlleleFreqFile
    // not T5_LargeOutputFile
    // not T4CoalescenceFst
    T1_AverageHybridIndexFile,
    T1_haplotypeFreqs_file,
    // not sampleSeq_file
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

    std::remove(path.c_str());
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
    assert(speciesIndex < this->subset.size());
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

std::string OutputFile::getPath(std::string patchIndexString, size_t mutPlacementIndex)
{
    if (this->isSpeciesSpecific)
    {
        assert(SSP != nullptr);
        if (this->isGenerationSpecific)
        {
            if (GP->nbSpecies == 1 && SSP->speciesName == "sp")
            {
                if (isT4MutationPlacementSpecific)
                {
                    assert(mutPlacementIndex >= 0);
                    return GeneralPath + filename + std::string("_G") + std::to_string(GP->CurrentGeneration) + patchIndexString + sequencingErrorStringToAddToFilnames + std::string("_mutPlac") + std::to_string(mutPlacementIndex) + extension;
                } else
                {
                    assert(mutPlacementIndex == 0);
                    return GeneralPath + filename + std::string("_G") + std::to_string(GP->CurrentGeneration) + patchIndexString + sequencingErrorStringToAddToFilnames  + extension;
                }
            } else
            {
                if (isT4MutationPlacementSpecific)
                {
                    assert(mutPlacementIndex >= 0);
                    return GeneralPath + filename + "_" + SSP->speciesName + std::string("_G") + std::to_string(GP->CurrentGeneration) + patchIndexString + sequencingErrorStringToAddToFilnames + std::string("_mutPlac") + std::to_string(mutPlacementIndex) + extension;    
                } else
                {
                    assert(mutPlacementIndex == 0);
                    return GeneralPath + filename + "_" + SSP->speciesName + std::string("_G") + std::to_string(GP->CurrentGeneration) + patchIndexString + sequencingErrorStringToAddToFilnames  + extension;    
                }
            }
            
        } else
        {
            if (GP->nbSpecies == 1 && SSP->speciesName == "sp")
            {
                if (isT4MutationPlacementSpecific)
                {
                    assert(mutPlacementIndex >= 0);
                    return GeneralPath + filename + patchIndexString + sequencingErrorStringToAddToFilnames + std::string("_mutPlac") + std::to_string(mutPlacementIndex) + extension;
                } else
                {
                    assert(mutPlacementIndex == 0);
                    return GeneralPath + filename + patchIndexString + sequencingErrorStringToAddToFilnames  + extension;
                }
            } else
            {
                if (isT4MutationPlacementSpecific)
                {
                    assert(mutPlacementIndex >= 0);
                    return GeneralPath + filename + "_" + SSP->speciesName + patchIndexString + sequencingErrorStringToAddToFilnames + std::string("_mutPlac") + std::to_string(mutPlacementIndex) + extension;    
                } else
                {
                    assert(mutPlacementIndex == 0);
                    return GeneralPath + filename + "_" + SSP->speciesName + patchIndexString + sequencingErrorStringToAddToFilnames + extension;    
                }
            }
            
        }
        
    } else
    {
        if (this->isGenerationSpecific)
        {
            if (isT4MutationPlacementSpecific)
            {
                assert(mutPlacementIndex >= 0);
                return GeneralPath + filename + std::string("_G") + std::to_string(GP->CurrentGeneration) + patchIndexString + sequencingErrorStringToAddToFilnames + std::string("_mutPlac") + std::to_string(mutPlacementIndex) + extension;
            } else
            {
                assert(mutPlacementIndex == 0);
                return GeneralPath + filename + std::string("_G") + std::to_string(GP->CurrentGeneration) + patchIndexString + sequencingErrorStringToAddToFilnames  + extension;
            }
        } else
        {
            return GeneralPath + filename + patchIndexString + sequencingErrorStringToAddToFilnames  + extension;
        }
    }       
}

std::string OutputFile::getPathWithoutGenerationDespiteBeingGenerationSpecific(size_t mutPlacementIndex)
{
    if (this->isSpeciesSpecific)
    {
        assert(SSP != nullptr);

        if (GP->nbSpecies == 1 && SSP->speciesName == "sp")
        {
            if (isT4MutationPlacementSpecific)
            {
                assert(mutPlacementIndex >= 0);
                return GeneralPath + filename + sequencingErrorStringToAddToFilnames + std::string("_mutPlac") + std::to_string(mutPlacementIndex) + extension;
            } else
            {
                assert(mutPlacementIndex == 0);
                return GeneralPath + filename + sequencingErrorStringToAddToFilnames  + extension;
            }
        } else
        {
            if (isT4MutationPlacementSpecific)
            {
                assert(mutPlacementIndex >= 0);
                return GeneralPath + filename + "_" + SSP->speciesName + sequencingErrorStringToAddToFilnames + std::string("_mutPlac") + std::to_string(mutPlacementIndex) + extension;
            } else
            {
                assert(mutPlacementIndex == 0);
                return GeneralPath + filename + "_" + SSP->speciesName + sequencingErrorStringToAddToFilnames  + extension;    
            }
        }
        
    } else
    {
        if (isT4MutationPlacementSpecific)
        {
            assert(mutPlacementIndex >= 0);
            return GeneralPath + filename + sequencingErrorStringToAddToFilnames + std::string("_mutPlac") + std::to_string(mutPlacementIndex) + extension;
        } else
        {
            assert(mutPlacementIndex == 0);
            return GeneralPath + filename + sequencingErrorStringToAddToFilnames  + extension;
        }
    }
}

std::string OutputFile::getPath(int generation, std::string patchIndexString, size_t mutPlacementIndex)
{
    assert(this->isGenerationSpecific);
    if (this->isSpeciesSpecific)
    {
        assert(SSP != nullptr);

        if (GP->nbSpecies == 1 && SSP->speciesName == "sp")
        {
            if (isT4MutationPlacementSpecific)
            {
                assert(mutPlacementIndex >= 0);
                return GeneralPath + filename + std::string("_G") + std::to_string(generation) + patchIndexString + sequencingErrorStringToAddToFilnames + std::string("_mutPlac") + std::to_string(mutPlacementIndex) + extension;
            } else
            {
                assert(mutPlacementIndex == 0);
                return GeneralPath + filename + std::string("_G") + std::to_string(generation) + patchIndexString + sequencingErrorStringToAddToFilnames  + extension;
            }
        } else
        {
            if (isT4MutationPlacementSpecific)
            {
                assert(mutPlacementIndex >= 0);
                return GeneralPath + filename + "_" + SSP->speciesName + std::string("_G") + std::to_string(generation) + patchIndexString + sequencingErrorStringToAddToFilnames + std::string("_mutPlac") + std::to_string(mutPlacementIndex) + extension;
            } else
            {
                assert(mutPlacementIndex == 0);
                return GeneralPath + filename + "_" + SSP->speciesName + std::string("_G") + std::to_string(generation) + patchIndexString + sequencingErrorStringToAddToFilnames  + extension;
            }
        }
        
    } else
    {
        if (isT4MutationPlacementSpecific)
        {
            assert(mutPlacementIndex >= 0);
            return GeneralPath + filename + std::string("_G") + std::to_string(generation) + patchIndexString + sequencingErrorStringToAddToFilnames + std::string("_mutPlac") + std::to_string(mutPlacementIndex) + extension;
        } else
        {
            assert(mutPlacementIndex == 0);
            return GeneralPath + filename + std::string("_G") + std::to_string(generation) + patchIndexString + sequencingErrorStringToAddToFilnames  + extension;
        }
    }
}

void OutputFile::interpretTimeForPaintedHaplo(InputReader& input)
{
    std::vector<int> observedGenerations;
    while (input.IsThereMoreToRead())
    {
        // Assert I get keyword 'generations'
        {
            auto keywordgenerations = input.GetNextElementString();
            if (keywordgenerations != "generations")
            {
                std::cout << "For outputs concerning painted haplotypes, expected the keyword 'generations' but received " << keywordgenerations << " instead.\n";
                abort();
            }
        }
   

        auto token = input.GetNextElementString();
        auto sepPos = token.find("-");
        if (sepPos == std::string::npos)
        {
            std::cout << "For outputs concerning painted haplotypes, expected times in format 'paintedGeneration->observedGeneration' (e.g. '50-200'). Received token " << token << " that did not include any '-' (dash) sign. Note if you wrote something like '50 - 200' (note the spaces), then the tokens will be read as '50', '-' and '200' which does not work. So please, don't put spaces around '-'.";
            abort();
        }
        std::string paintedGeneration_s  = token.substr(0, sepPos);
        std::string observedGeneration_s = token.substr(sepPos + 1);
   
        int paintedGeneration;
        int observedGeneration;
        if (paintedGeneration_s == "end")
        {
            std::cout << "For outputs concerning painted haplotypes, expected times in format 'paintedGeneration-observedGeneration' (e.g. '50->200'). The paintedGeneration received is 'end' (end = " << GP->CurrentGeneration << "). Sorry, you cannot paint the haplotype at the last generation as you need to leave at least one generation between the generation when you paint the haplotypes and the generation when you observe how the painted piece have segregated.\n";
            abort();
        } else
        {
            if (paintedGeneration_s == "anc")
            {
                paintedGeneration = std::numeric_limits<int>::lowest();
            } else
            {
                paintedGeneration = (int) std::stod(paintedGeneration_s);
            }
        }
   
        if (observedGeneration_s == "end")
        {
            observedGeneration = GP->nbGenerations;
        } else
        {
            observedGeneration = (int) std::stod(observedGeneration_s);
        }

        if (observedGeneration > GP->nbGenerations)
        {
            std::cout << "For outputs concerning painted haplotypes, expected times in format 'paintedGeneration-observedGeneration' (e.g. '50-200'). The paintedGeneration received is "<<paintedGeneration<<" and the observedGeneration received is "<<observedGeneration<<"'. observedGeneration must be lower or equal to the total nuumber of generations simulated (as indicated with option --nbGens (--nbGenerations), the total number of generations simulated is "<< GP->nbGenerations<<")\n";
            abort();
        }
   
        if (paintedGeneration >= observedGeneration)
        {
            std::cout << "For outputs concerning painted haplotypes, expected times in format 'paintedGeneration-observedGeneration' (e.g. '50-200'). The paintedGeneration received is "<<paintedGeneration<<" and the observedGeneration received is "<<observedGeneration<<"'. observedGeneration must be strictly greater than paintedGeneration\n";
            abort();
        }


        auto generation_index = std::upper_bound(GP->__GenerationChange.begin(), GP->__GenerationChange.end(),  observedGeneration) - GP->__GenerationChange.begin() - 1;
        assert(generation_index < GP->__PatchNumber.size());
        
        // Assert I get keyword 'patch'
        {
            auto keywordPatch = input.GetNextElementString();
            if (keywordPatch != "patch")
            {
                std::cout << "For outputs concerning painted haplotypes, expected the keyword 'patch' but received " << keywordPatch << " instead.\n";
                abort();
            }
        }

        std::vector<int> patch_indices;
        while (input.IsThereMoreToRead() && input.PeakNextElementString() != "nbHaplos" && input.PeakNextElementString() != "generations" && input.PeakNextElementString() != "patch")
        {
            auto patch_index = input.GetNextElementInt();
            patch_indices.push_back(patch_index);

            if (patch_index >= GP->__PatchNumber[generation_index])
            {
                std::cout << "For outputs concerning painted haplotypes, received the patch_index " << patch_index << " for observation at generation " << observedGeneration << " of haplotypes painted at generation " << paintedGeneration << ". generation " << observedGeneration << " there are however only " << GP->__PatchNumber[generation_index] << " patches. It is therefore impossible to sample from patch " << patch_index << ". As a reminder, the patch_index (just like any other indices in SimBit) is zero-based counting.\n";
                abort();
            }
        }


        if (patch_indices.size() == 0)
        {
            std::cout << "For outputs concerning painted haplotypes, expected some input (different from keywords 'patch', 'nbHaplos' or 'generations') after the keyword 'patch'\n";
            abort();
        }


        // Assert I got keyword 'nbHaplos'
        {
            auto keywordnbHaplos = input.GetNextElementString();
            if (keywordnbHaplos != "nbHaplos")
            {
                std::cout << "For outputs concerning painted haplotypes, expected the keyword 'nbHaplos' but received " << keywordnbHaplos << " instead. Note this sounds like an internal bug as the code should have been able to run this line only if the keyword 'nbHaplos' had been found. Weird! Even if you end up finding it is caused by a non-sense input of yours, please let Remi know about it.\n";
                abort();
            }
        }
        

        assert(patch_indices.size());
        std::vector<int> nbHaplotypes;
        while (input.IsThereMoreToRead() && input.PeakNextElementString() != "nbHaplos" && input.PeakNextElementString() != "generations" && input.PeakNextElementString() != "patch")
        {
            auto nbH = input.GetNextElementInt();
            nbHaplotypes.push_back(nbH);

            if (nbHaplotypes.size() > patch_indices.size())
            {
                std::cout << "For outputs concerning painted haplotypes, received more values following the keyword 'nbHaplos' than values following the keyword 'patch'\n";
                abort();
            }

            assert(nbHaplotypes.size() - 1 < patch_indices.size());

            if (GP->nbSpecies == 1)
            {
                if (!allParameters.SSPs[0].T4_paintedHaplo_shouldIgnorePatchSizeSecurityChecks)
                {
                    auto patch_index = patch_indices[nbHaplotypes.size() - 1];

                    if (nbH > 2 * allParameters.SSPs[0].__patchCapacity[generation_index][patch_index])
                    {
                        std::cout << "For outputs concerning painted haplotypes, received the patch_index " << patch_index << " for observation at generation " << observedGeneration << " of haplotypes painted at generation " << paintedGeneration << ". You asked for sampling " << nbH << " haplotypes. The carrying capacity for the patch " << patch_index << " at generation " << observedGeneration << " is only " << allParameters.SSPs[0].__patchCapacity[generation_index][patch_index] << " individuals (or " << 2 * allParameters.SSPs[0].__patchCapacity[generation_index][patch_index]<< " haplotypes).\n";
                        abort();
                    }
                }
            }
        }
        
        if (nbHaplotypes.size() == 0)
        {
            std::cout << "For outputs concerning painted haplotypes, expected some input (different from keywords 'patch', 'nbHaplos' or 'generations') after the keyword 'patch'\n";
            abort();
        }
        if (nbHaplotypes.size() != patch_indices.size())
        {
            std::cout << "For outputs concerning painted haplotypes, received " << patch_indices.size() << " patches (after keyword 'patch') but " << nbHaplotypes.size() << " number of haplotypes per patch values (after keyword 'nbHaplos').\n";
            abort();
        }
        
    
        observedGenerations.push_back(observedGeneration);
        T4TreeRec::generationsToKeepInTheTree.push_back(paintedGeneration);
        this->T4_paintedHaplo_information.push_back(
            {
                (int) paintedGeneration,
                (int) observedGeneration,
                patch_indices,
                nbHaplotypes
            }
        );
    }

    std::sort(observedGenerations.begin(), observedGenerations.end());
    observedGenerations.erase( std::unique( observedGenerations.begin(), observedGenerations.end() ), observedGenerations.end() );    

    this->setTimes(observedGenerations);

    std::sort(T4TreeRec::generationsToKeepInTheTree.begin(), T4TreeRec::generationsToKeepInTheTree.end());
    T4TreeRec::generationsToKeepInTheTree.erase( std::unique( T4TreeRec::generationsToKeepInTheTree.begin(), T4TreeRec::generationsToKeepInTheTree.end() ), T4TreeRec::generationsToKeepInTheTree.end() );


    assert(this->T4_paintedHaplo_information.size());
    std::sort(
        this->T4_paintedHaplo_information.begin(),
        this->T4_paintedHaplo_information.end(),
        [](const OutputFile::T4_paintedHaplo_info& lhs, const OutputFile::T4_paintedHaplo_info& rhs)
        {
            if (lhs.observedGeneration == rhs.observedGeneration)
            {
                return lhs.paintedGeneration < rhs.paintedGeneration;
            } else
            {
                return lhs.observedGeneration < rhs.observedGeneration;
            }
        }
    );

    assert(this->T4_paintedHaplo_information.size());
    /*
    for (size_t i = 1 ; i < this->T4_paintedHaplo_information.size() ; ++i)
    {
        auto& curr = this->T4_paintedHaplo_information[i];
        auto& prev = this->T4_paintedHaplo_information[i-1];
        if (curr.observedGeneration == prev.observedGeneration && curr.paintedGeneration == prev.paintedGeneration && curr.patch_indices == prev.patch_indices)
        {
            std::cout << "For outputs concerning painted haplotypes, received twice the information to sample the sames patches at generation " << prev.observedGeneration << " of haplotypes painted at generation " << curr.paintedGeneration << ". Maybe you only meant to add up the number of haplotypes but I am unsure that's what you meant so I prefer to just raise an error message just in case.\n";
            abort();
        }
    }
    */

    assert(this->T4_paintedHaplo_information.size());
}

bool OutputFile::getIsT4MutationPlacementSpecific()
{
    return isT4MutationPlacementSpecific;
}

void OutputFile::interpretTimeAndSubsetInput(InputReader& input)
{
    if (!input.IsThereMoreToRead())
    {
        std::cout << input.GetErrorMessage() << "was expecting time information after file name!\n";
        abort();
    }    

    std::vector<int> v;
    while (input.IsThereMoreToRead())
    {
        if (input.PeakNextElementString() == "subset" || input.PeakNextElementString() == "sequence")
        {
            break;
        }
        if (input.PeakNextElementString() == "fromtoby" || input.PeakNextElementString() == "FromToBy" || input.PeakNextElementString() == "fromToBy")
        {
            input.skipElement();
            int from = input.GetNextElementInt();
            int to;
            if (input.PeakNextElementString() == "end")
            {
                input.skipElement();
                to = GP->nbGenerations;
            } else
            {
                to = input.GetNextElementInt();
            }                
            int by = input.GetNextElementInt();


            for (int t = from ; t <= to ; t+=by)
            {
                v.push_back(t);
            }
            if (v.back() < to) v.push_back(to);
        } else
        {
            if (input.PeakNextElementString() == "end")
            {
                input.skipElement();
                v.push_back(GP->nbGenerations);
            } else
            {
                v.push_back(input.GetNextElementInt());
            }
                
        }
    }
    this->setTimes(v);
    


    if (std::find(listOfOutputFileTypeThatCanTakeASubset.begin(), listOfOutputFileTypeThatCanTakeASubset.end(), this->OutputFileType) != listOfOutputFileTypeThatCanTakeASubset.end()) 
    {
        if (input.IsThereMoreToRead() && input.PeakNextElementString() != "sequence")
        {
            assert(input.GetNextElementString() == "subset");

            if (std::find(listOfOutputFileTypeThatCanTakeASubset.begin(), listOfOutputFileTypeThatCanTakeASubset.end(), this->OutputFileType) == listOfOutputFileTypeThatCanTakeASubset.end())
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
                assert(from <= to);
                InputReader inputOneSpecies(input, from, to, speciesIndex);

                std::vector<T1_locusDescription> oneSpeciesSubset;

                while (inputOneSpecies.IsThereMoreToRead())
                {
                    int locus = inputOneSpecies.GetNextElementInt();
                    oneSpeciesSubset.push_back({locus/8, locus%8, locus});
                }  
            inputOneSpecies.workDone(); 
                this->subset.push_back(oneSpeciesSubset);

            }
            SSP = SSP_toReset;
        } else
        {
            // empty subset means all loci are in subset
            for (int speciesIndex = 0 ; speciesIndex < GP->nbSpecies; speciesIndex++)
            {
                std::vector<T1_locusDescription> oneSpeciesSubset;
                for (int locus = 0 ; locus < allParameters.getSSPaddress(speciesIndex)->Gmap.T1_nbLoci; locus++)
                    oneSpeciesSubset.push_back({locus/8, locus%8, locus}); 
                this->subset.push_back(oneSpeciesSubset);
            }
        }            
        assert(this->subset.size() == GP->nbSpecies);
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
        if (this->OutputFileType != T1_vcfFile && this->OutputFileType != T4_vcfFile && this->OutputFileType != T56_vcfFile && this->OutputFileType != SaveBinaryFile)
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
    /*
    std::cout << "Trying to find time " << GP->CurrentGeneration << "\n";
    std::cout << "times was set at ";
    for (auto& t : times)
        std::cout << t << " ";
    std::cout << "\n";
    if (std::find(times.begin(), times.end(),  GP->CurrentGeneration) != times.end())
        std::cout << "time has been found!";
    */

    if (KillOnDemand::justAboutToKill)
    {
        if (times.back() > GP->CurrentGeneration)
        {
            return true;
        } else
        {
            return false;
        }
    } else
    {
        return std::find(times.begin(), times.end(),  GP->CurrentGeneration) != times.end();
    }
}

OutputFile::OutputFile(OutputFile&& f)
: filename(std::move(f.filename)), extension(f.extension), OutputFileType(f.OutputFileType), times(f.times), isGenerationSpecific(f.isGenerationSpecific), isSpeciesSpecific(f.isSpeciesSpecific), isPatchSpecific(f.isPatchSpecific), isNbLinesEqualNbOutputTimes(f.isNbLinesEqualNbOutputTimes), doesTimeNeedsToBeSet(f.doesTimeNeedsToBeSet), isT4MutationPlacementSpecific(f.isT4MutationPlacementSpecific), subset(f.subset), T4_paintedHaplo_information(std::move(f.T4_paintedHaplo_information))
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
        isPatchSpecific = false;
        isNbLinesEqualNbOutputTimes = false;
        doesTimeNeedsToBeSet = false;
        isT4MutationPlacementSpecific = false;
    } else if (t == T1_vcfFile)
    {
        this->extension = std::string(".T1vcf");
        isGenerationSpecific = true;
        isSpeciesSpecific = true;
        isPatchSpecific = false;
        isNbLinesEqualNbOutputTimes = false;
        doesTimeNeedsToBeSet = true;
        isT4MutationPlacementSpecific = false;
    } else if (t == T1_LargeOutputFile)
    {
        this->extension = std::string(".T1LO");
        isGenerationSpecific = false;
        isSpeciesSpecific = true;
        isPatchSpecific = false;
        isNbLinesEqualNbOutputTimes = true;
        doesTimeNeedsToBeSet = true;
        isT4MutationPlacementSpecific = false;
    } else if (t == T1_AlleleFreqFile)
    {
        this->extension = std::string(".T1AllFreq");
        isGenerationSpecific = false;
        isSpeciesSpecific = true;
        isPatchSpecific = false;
        isNbLinesEqualNbOutputTimes = true;
        doesTimeNeedsToBeSet = true;
        isT4MutationPlacementSpecific = false;
    } else if (t == MeanLDFile)
    {
        this->extension = std::string(".T1LD");
        isGenerationSpecific = false;
        isSpeciesSpecific = true;
        isPatchSpecific = false;
        isNbLinesEqualNbOutputTimes = true;
        doesTimeNeedsToBeSet = true;
    } else if (t == LongestRunFile)
    {
        this->extension = std::string(".T1LR");
        isGenerationSpecific = false;
        isSpeciesSpecific = true;
        isPatchSpecific = false;
        isNbLinesEqualNbOutputTimes = true;
        doesTimeNeedsToBeSet = true;
        isT4MutationPlacementSpecific = false;
    } else if (t == HybridIndexFile)
    {
        this->extension = std::string(".T1HI");
        isGenerationSpecific = false;
        isSpeciesSpecific = true;
        isPatchSpecific = false;
        isNbLinesEqualNbOutputTimes = true;
        doesTimeNeedsToBeSet = true;
        isT4MutationPlacementSpecific = false;
    } else if (t == T1_AverageHybridIndexFile)
    {
        this->extension = std::string(".T1AverageHI");
        isGenerationSpecific = false;
        isSpeciesSpecific = true;
        isPatchSpecific = false;
        isNbLinesEqualNbOutputTimes = true;
        doesTimeNeedsToBeSet = true;
        isT4MutationPlacementSpecific = false;
    } else if (t == T1_haplotypeFreqs_file)
    {
        this->extension = std::string(".T1haploFreqs");
        isGenerationSpecific = true; // Generation specific because the number of columns will vary
        isSpeciesSpecific = true;
        isPatchSpecific = false;
        isNbLinesEqualNbOutputTimes = true;
        doesTimeNeedsToBeSet = true;
        isT4MutationPlacementSpecific = false;
    } else if (t == ExpectiMinRecFile)
    {
        this->extension = std::string(".T1EMR");
        isGenerationSpecific = false;
        isSpeciesSpecific = true;
        isPatchSpecific = false;
        isNbLinesEqualNbOutputTimes = true;
        doesTimeNeedsToBeSet = true;
        isT4MutationPlacementSpecific = false;
    } else if (t == T2_LargeOutputFile)
    {
        this->extension = std::string(".T2LO");
        isGenerationSpecific = false;
        isSpeciesSpecific = true;
        isPatchSpecific = false;
        isNbLinesEqualNbOutputTimes = true;
        doesTimeNeedsToBeSet = true;
        isT4MutationPlacementSpecific = false;
    } else if (t == SaveBinaryFile)
    {
        this->extension = std::string(".bin");
        isGenerationSpecific = true;
        isSpeciesSpecific = true;
        isPatchSpecific = false;
        isNbLinesEqualNbOutputTimes = false;
        doesTimeNeedsToBeSet = true;
        isT4MutationPlacementSpecific = false;
    } else if (t == T3_MeanVarFile)
    {
        this->extension = std::string(".T3MeanVar");
        isGenerationSpecific = false;
        isSpeciesSpecific = true;
        isPatchSpecific = false;
        isNbLinesEqualNbOutputTimes = true;
        doesTimeNeedsToBeSet = true;
        isT4MutationPlacementSpecific = false;
    } else if (t == T3_LargeOutputFile)
    {
        this->extension = std::string(".T3LO");
        isGenerationSpecific = false;
        isSpeciesSpecific = true;
        isPatchSpecific = false;
        isNbLinesEqualNbOutputTimes = true;
        doesTimeNeedsToBeSet = true;
        isT4MutationPlacementSpecific = false;
    } else if (t == fitness)
    {
        this->extension = std::string(".fit");
        isGenerationSpecific = false;
        isSpeciesSpecific = true;
        isPatchSpecific = false;
        isNbLinesEqualNbOutputTimes = true;
        doesTimeNeedsToBeSet = true;
        isT4MutationPlacementSpecific = false;
    } else if (t == fitnessStats)
    {
        this->extension = std::string(".fitStats");
        isGenerationSpecific = false;
        isSpeciesSpecific = true;
        isPatchSpecific = false;
        isNbLinesEqualNbOutputTimes = true;
        doesTimeNeedsToBeSet = true;
        isT4MutationPlacementSpecific = false;
    } else if (t == T1_FST)
    {
        this->extension = std::string(".T1_FST");
        isGenerationSpecific = false;
        isSpeciesSpecific = true;
        isPatchSpecific = false;
        isNbLinesEqualNbOutputTimes = true;
        doesTimeNeedsToBeSet = true;
        isT4MutationPlacementSpecific = false;
    } else if (t == extraGeneticInfo)
    {
        this->extension = std::string(".extraGeneticInfo");
        isGenerationSpecific = false;
        isSpeciesSpecific = true;
        isPatchSpecific = false;
        isNbLinesEqualNbOutputTimes = false;
        doesTimeNeedsToBeSet = false;
        isT4MutationPlacementSpecific = false;
    }  else if (t == patchSize)
    {
        this->extension = std::string(".patchSize");
        isGenerationSpecific = false;
        isSpeciesSpecific = true;
        isPatchSpecific = false;
        isNbLinesEqualNbOutputTimes = true;
        doesTimeNeedsToBeSet = true;
        isT4MutationPlacementSpecific = false;
    } else if (t == extinction)
    {
        this->extension = std::string(".extinction");
        isGenerationSpecific = false;
        isSpeciesSpecific = false;
        isPatchSpecific = false;
        isNbLinesEqualNbOutputTimes = false;
        doesTimeNeedsToBeSet = false;
        isT4MutationPlacementSpecific = false;
    } else if (t == genealogy)
    {
        this->extension = std::string(".genealogy");
        isGenerationSpecific = true;   // Actually the end result is not. SimBit uses getPathWithoutGenerationDespiteBeingGenerationSpecific to produce the last file.
        isSpeciesSpecific = true;
        isPatchSpecific = false;
        isNbLinesEqualNbOutputTimes = false;
        doesTimeNeedsToBeSet = true;
        isT4MutationPlacementSpecific = false;
    } else if (t == fitnessSubsetLoci)
    {
        this->extension = std::string(".fitSubLoci");
        isGenerationSpecific = false;
        isSpeciesSpecific = true;
        isPatchSpecific = false;
        isNbLinesEqualNbOutputTimes = true;
        doesTimeNeedsToBeSet = true;
        isT4MutationPlacementSpecific = false;
    } else if (t == T4_LargeOutputFile)
    {
        this->extension = std::string(".T4LO");
        isGenerationSpecific = false;
        isSpeciesSpecific = true;
        isPatchSpecific = false;
        isNbLinesEqualNbOutputTimes = true;
        doesTimeNeedsToBeSet = true;
        isT4MutationPlacementSpecific = true;
    } else if (t == T4_vcfFile)
    {
        this->extension = std::string(".T4vcf");
        isGenerationSpecific = true;
        isSpeciesSpecific = true;
        isPatchSpecific = false;
        isNbLinesEqualNbOutputTimes = false;
        doesTimeNeedsToBeSet = true;
        isT4MutationPlacementSpecific = true;
    } else if (t == T4_SFS_file)
    {
        this->extension = std::string(".T4SFS");
        isGenerationSpecific = false;
        isSpeciesSpecific = true;
        isPatchSpecific = false;
        isNbLinesEqualNbOutputTimes = true;
        doesTimeNeedsToBeSet = true;
        isT4MutationPlacementSpecific = true;
    } else if (t == T1_SFS_file)
    {
        this->extension = std::string(".T1SFS");
        isGenerationSpecific = false;
        isSpeciesSpecific = true;
        isPatchSpecific = false;
        isNbLinesEqualNbOutputTimes = true;
        doesTimeNeedsToBeSet = true;
        isT4MutationPlacementSpecific = false;
    } else if (t == T4_printTree)
    {
        this->extension = std::string(".T4tree");
        // The following false/true should not matter as this file is not going into the outputWriter
        isGenerationSpecific = false;
        isSpeciesSpecific = true;
        isPatchSpecific = false;
        isNbLinesEqualNbOutputTimes = false;
        isT4MutationPlacementSpecific = false;
    } else if (t == T56_vcfFile)
    {
        this->extension = std::string(".T5vcf");
        isGenerationSpecific = true;
        isSpeciesSpecific = true;
        isPatchSpecific = true;
        isNbLinesEqualNbOutputTimes = false;
        doesTimeNeedsToBeSet = true;
        isT4MutationPlacementSpecific = false;
    } else if (t == T56_SFS_file)
    {
        this->extension = std::string(".T5SFS");
        isGenerationSpecific = false;
        isSpeciesSpecific = true;
        isPatchSpecific = false;
        isNbLinesEqualNbOutputTimes = true;
        doesTimeNeedsToBeSet = true;
        isT4MutationPlacementSpecific = false;
    } else if (t == T56_AlleleFreqFile)
    {
        this->extension = std::string(".T5AllFreq");
        isGenerationSpecific = false;
        isSpeciesSpecific = true;
        isPatchSpecific = false;
        isNbLinesEqualNbOutputTimes = true;
        doesTimeNeedsToBeSet = true;
        isT4MutationPlacementSpecific = false;
    } else if (t == T56_LargeOutputFile)
    {
        this->extension = std::string(".T5LO");
        isGenerationSpecific = false;
        isSpeciesSpecific = true;
        isPatchSpecific = false;
        isNbLinesEqualNbOutputTimes = true;
        doesTimeNeedsToBeSet = true;
        isT4MutationPlacementSpecific = false;
    } else if (t == T4CoalescenceFst)
    {
        this->extension = std::string(".T4CoalescenceFst");
        isGenerationSpecific = false;
        isSpeciesSpecific = true;
        isPatchSpecific = false;
        isNbLinesEqualNbOutputTimes = true;
        doesTimeNeedsToBeSet = true;
        isT4MutationPlacementSpecific = false;
    } else if (t == sampleSeq_file)
    {
        this->extension = std::string(".seq");
        isGenerationSpecific = true;
        isSpeciesSpecific = true;
        isPatchSpecific = false;
        isNbLinesEqualNbOutputTimes = false;
        doesTimeNeedsToBeSet = true;
        isT4MutationPlacementSpecific = true;
    } else if (t == T4_paintedHaplo_file)
    {
        this->extension = std::string(".paintedHaplo");
        isGenerationSpecific = false;
        isSpeciesSpecific = true;
        isPatchSpecific = false;
        isNbLinesEqualNbOutputTimes = false;
        doesTimeNeedsToBeSet = true;
        isT4MutationPlacementSpecific = false;
    } else if (t == T4_paintedHaploSegmentsDiversity_file)
    {
        this->extension = std::string(".paintedHaploDiversity");
        isGenerationSpecific = false;
        isSpeciesSpecific = true;
        isPatchSpecific = false;
        isNbLinesEqualNbOutputTimes = false;
        doesTimeNeedsToBeSet = true;
        isT4MutationPlacementSpecific = false;
    } else if (t == T4_SNPfreq_file)
    {
        this->extension = std::string(".T4_SNPfreq");
        isGenerationSpecific = false;
        isSpeciesSpecific = true;
        isPatchSpecific = false;
        isNbLinesEqualNbOutputTimes = false;
        doesTimeNeedsToBeSet = true;
        isT4MutationPlacementSpecific = true;
    } else if (t == Tx_SNPfreq_file)
    {
        this->extension = std::string(".Tx_SNPfreq");
        isGenerationSpecific = false;
        isSpeciesSpecific = true;
        isPatchSpecific = false;
        isNbLinesEqualNbOutputTimes = false;
        doesTimeNeedsToBeSet = true;
        isT4MutationPlacementSpecific = true;
    }  else if (t == burnInLength_file)
    {
        this->extension = std::string(".burnInLength");
        isGenerationSpecific = false;
        isSpeciesSpecific = false;
        isPatchSpecific = false;
        isNbLinesEqualNbOutputTimes = false;
        doesTimeNeedsToBeSet = false;
        isT4MutationPlacementSpecific = false;
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

void OutputFile::openPatchSpecific(int patch_index)
{
    std::string patchIndexString("");
    if (patch_index != -1)
    {
        patchIndexString = std::string("_P") + std::to_string(patch_index);
    }

    std::string path = this->getPath(patchIndexString);
    
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

void OutputFile::open()
{       
    std::string path = this->getPath();
    
    if (OutputFileType == SaveBinaryFile)
    {
        ofs.open(path, std::ios::out | std::ios::binary | std::ios::trunc);
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


void OutputFile::openT4(size_t mutPlacementIndex)
{
    std::string emptyString("");       
    std::string path = this->getPath(emptyString, mutPlacementIndex);
    
    
    ofs.open(path, std::ios_base::app);

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

void OutputFile::remove()
{
    std::string path;
    
    if (this->isPatchSpecific)
    {
        for (uint32_t patch_index = 0 ; patch_index < GP->maxEverPatchNumber; ++patch_index)
        {
            assert(!isT4MutationPlacementSpecific);
            std::string patchIndexString = std::string("_P") + std::to_string(patch_index);

            path = this->getPath(patchIndexString);

            ofs.open(path, std::ofstream::out | std::ofstream::trunc);
            if (!ofs.is_open())
            {
                std::cout << "OutputFileType " << getFileTypeName(OutputFileType) << "(index "<< OutputFileType <<")" << " with path '" << path << "' failed to open!" << "\n";
                abort();
            }
            ofs.close();
        
            std::remove(path.c_str());
        }
    } else
    {
        if (isT4MutationPlacementSpecific)
        {
            for (size_t mutPlacingIndex = 0 ; mutPlacingIndex < SSP->T4_nbMutationPlacingsPerOutput ; ++mutPlacingIndex)
            {
                std::string emptyString("");
                path = this->getPath(emptyString, mutPlacingIndex);

                ofs.open(path, std::ofstream::out | std::ofstream::trunc);
                if (!ofs.is_open())
                {
                    std::cout << "OutputFileType " << getFileTypeName(OutputFileType) << "(index "<< OutputFileType <<")" << " with path '" << path << "' failed to open!" << "\n";
                    abort();
                }
                ofs.close();
                
                std::remove(path.c_str());
            }
        } else
        {
            path = this->getPath();

            ofs.open(path, std::ofstream::out | std::ofstream::trunc);
            if (!ofs.is_open())
            {
                std::cout << "OutputFileType " << getFileTypeName(OutputFileType) << "(index "<< OutputFileType <<")" << " with path '" << path << "' failed to open!" << "\n";
                abort();
            }
            ofs.close();
            
            std::remove(path.c_str());
        }
    }
}

void OutputFile::clearContent()
{
    std::string path;
    
    if (this->isPatchSpecific)
    {
        for (uint32_t patch_index = 0 ; patch_index < GP->maxEverPatchNumber; ++patch_index)
        {
            std::string patchIndexString = std::string("_P") + std::to_string(patch_index);

            path = this->getPath(patchIndexString);
        
            ofs.open(path, std::ofstream::out | std::ofstream::trunc);
            if (!ofs.is_open())
            {
                std::cout << "OutputFileType " << getFileTypeName(OutputFileType) << "(index "<< OutputFileType <<")" << " with path '" << path << "' failed to open!" << "\n";
                abort();
            }
            ofs.close();
        }
    } else
    {
        path = this->getPath();
        
        ofs.open(path, std::ofstream::out | std::ofstream::trunc);
        if (!ofs.is_open())
        {
            std::cout << "OutputFileType " << getFileTypeName(OutputFileType) << "(index "<< OutputFileType <<")" << " with path '" << path << "' failed to open!" << "\n";
            abort();
        }
        ofs.close();
    }
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
    //assert(this->isOpen());
    ofs << s;
}

void OutputFile::write(const std::vector<std::string>& s)
{
    for (const auto& elem : s) write(elem);
}


template<typename T>
void OutputFile::writeBinary(const T x)
{
    //std::cout << "sizeof(x) = " << sizeof(x) << "\n";
    assert(this->isOpen());
    ofs.write((char*)&x, sizeof(T));
    //std::cout << "write: " << x << " of type " << typeid(T).name() << "\n";
}

template<typename T>
void OutputFile::writeBinary(const std::vector<T>& x)
{
    assert(this->isOpen());

    this->writeBinary((size_t) x.size());

    /*
    std::cout <<"write: ";
    for (size_t i = 0 ; i < x.size() ; ++i)
    {
        std::cout << x[i] << " ";
    }
    std::cout << " of type " << typeid(T).name() << "\n";
    */
    
    ofs.write(
        reinterpret_cast<const char*>(&x[0]),
        x.size()*sizeof(T)
    );
}

void OutputFile::writeBinary(const char* first, int second)
{
    //std::cout << "write first second\n";
    assert(this->isOpen());
    ofs.write(first, second);
}

void OutputFile::writeBinary(const RNG_type x)
{
    //std::cout << "write RNG\n";
    assert(this->isOpen());
    ofs << x; 
    //std::cout << x << "\n";
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
    /*
    std::cout << "\ttimes set at ";
    for (auto& t : times)
        std::cout << t << " ";
    std::cout << "\n";
    */
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

void OutputFile::assertSubsetSize()
{
    if (std::find(listOfOutputFileTypeThatCanTakeASubset.begin(), listOfOutputFileTypeThatCanTakeASubset.end(), this->OutputFileType) != listOfOutputFileTypeThatCanTakeASubset.end())
    {
        assert(this->subset.size() == GP->nbSpecies);
    }
}
