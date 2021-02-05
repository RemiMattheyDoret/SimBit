
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


void OutputWriter::PrintGenerationEndOfBurnIn()
{
    auto positiveGeneration = GP->CurrentGeneration - std::numeric_limits<int>::lowest();
    assert(positiveGeneration >= 0);
    //std::cout << GP->CurrentGeneration << " -> " << positiveGeneration << "\n";
    std::cout << "\r\t\tBurn in generation: " <<  positiveGeneration << " / " << positiveGeneration << "                  \n\n"<< std::flush;
}

void OutputWriter::PrintGeneration()
{
    //std::cout << "GP->CurrentGeneration = " << GP->CurrentGeneration << "\n";
    //std::cout << "GP->nbGenerations = " << GP->nbGenerations << "\n";
    //std::cout << "GP->printProgress = " << GP->printProgress << "\n";
    if ( GP->CurrentGeneration < 0 )
    {
        if (GP->printProgress)
        {
            auto positiveGeneration = GP->CurrentGeneration - std::numeric_limits<int>::lowest();
            assert(positiveGeneration >= 0);
            //std::cout << GP->CurrentGeneration << " -> " << positiveGeneration << "\n";
            if ( positiveGeneration < 100)
            {
                std::cout << "\r\t\tBurn in generation: " <<  positiveGeneration << "           " << std::flush;
            } else if ( positiveGeneration < 1000 && !( positiveGeneration % 10))
            {
                std::cout << "\r\t\tBurn in generation: " <<  positiveGeneration << "           " << std::flush;
            } else if ( positiveGeneration < 100000 && !(positiveGeneration % 100))
            {
                std::cout << "\r\t\tBurn in generation: " <<  positiveGeneration << "           " << std::flush;
            } else if (!(positiveGeneration % 10000))
            {
                std::cout << "\r\t\tBurn in generation: " <<  positiveGeneration << "           " << std::flush;
            }
        }
    } else
    {
        if ( GP->CurrentGeneration == GP->nbGenerations)
        {
            std::cout << "\r\t\tGeneration: " <<  GP->CurrentGeneration << " / " << GP->nbGenerations << "           \n\n" << std::flush;
        } else if (GP->printProgress)
        {
            if ( GP->CurrentGeneration < 100)
            {
                std::cout << "\r\t\tGeneration: " <<  GP->CurrentGeneration << " / " << GP->nbGenerations << "           " << std::flush;
            } else if ( GP->CurrentGeneration < 1000 && !( GP->CurrentGeneration % 10))
            {
                std::cout << "\r\t\tGeneration: " <<  GP->CurrentGeneration << " / " << GP->nbGenerations << "           " << std::flush;
            } else if ( GP->CurrentGeneration < 100000 && !(GP->CurrentGeneration % 100))
            {
                std::cout << "\r\t\tGeneration: " <<  GP->CurrentGeneration << " / " << GP->nbGenerations << "           " << std::flush;
            } else if (!(GP->CurrentGeneration % 10000))
            {
                std::cout << "\r\t\tGeneration: " <<  GP->CurrentGeneration << " / " << GP->nbGenerations << "           " << std::flush;
            }
        }
    }
}

/*
template<typename T>
void OutputWriter::ExtendStringForAdvancedLogFile(std::string& s, T& entry, int depth, bool willMoreCome)
{
    std::string tabs;
    for (int i = 0 ; i < depth ; i++)
        tabs += "\t";
    s += tabs + std::to_string(entry);
    if (!willMoreCome)
    {
      s+= "\n";
    } else
    {
      s+=" ";
    }
}

void OutputWriter::ExtendStringForAdvancedLogFile(std::string& s, std::string entry, int depth, bool willMoreCome)
{
    std::string tabs;
    for (int i = 0 ; i < depth ; i++)
        tabs += "\t";
    s += tabs + entry;
    if (!willMoreCome)
    {
      s+= "\n";
    } else
    {
      s+=" ";
    }
}

template<typename T>
void OutputWriter::ExtendStringForAdvancedLogFile(std::string& s, std::vector<T>& entry, int depth, bool willMoreCome)
{
    std::string tabs;
    for (int i = 0 ; i < depth ; i++)
        tabs += "\t";
    s += tabs + "{\n";
    if (entry.size() > 0)
    {
        for (int i = 0 ; i < (entry.size() - 1) ; i++)
          ExtendStringForAdvancedLogFile(s,entry[i],depth+1,true);
        ExtendStringForAdvancedLogFile(s,entry.back(),depth+1,false);
    }

    s += tabs + "}\n";
}

void OutputWriter::ExtendStringForAdvancedLogFile(std::string& s, T1_locusDescription& entry, int depth, bool willMoreCome)
{
    std::string tabs;
    for (int i = 0 ; i < depth ; i++)
        tabs += "\t";
    s += tabs + entry.toString();
    if (!willMoreCome)
    {
      s+= "\n";
    } else
    {
      s+=" ";
    }
}


template<typename T>
void OutputWriter::ExtendStringForAdvancedLogFile(std::string& s, T& entry, std::string name)
{
    s += name + "\n";
    ExtendStringForAdvancedLogFile(s,entry,0,false);
    s+= "\n\n";
}*/


template<typename T>
void OutputWriter::ExtendStringForAdvancedLogFile(std::string& s, T& entry, std::string name)
{
    s += name + "\n" + std::to_string(entry) + "\n\n";
}

void OutputWriter::ExtendStringForAdvancedLogFile(std::string& s, std::string& entry, std::string name)
{
    s += name + "\n'" + entry + "'\n\n";
}

template<typename T>
void OutputWriter::ExtendStringForAdvancedLogFile(std::string& s, std::vector<T>& entry, std::string name)
{
    s += name + "\nc(";
    for (auto& elem : entry)
        s += std::to_string(elem) + ",";
    if (s.back() == ',') s.pop_back();
    s +=  ")\n\n";
}

template<typename T>
void OutputWriter::ExtendStringForAdvancedLogFile(std::string& s, std::vector<std::vector<T>>& entry, std::string name)
{
    s += name + "\nlist("; 
    for (auto& elemA : entry)
    {
        s += "c(";
        for (auto& elemB : elemA)
            s += std::to_string(elemB) + ",";
        if (s.back() == ',') s.pop_back();
        s += "),";
    }
    if (s.back() == ',') s.pop_back();
    s +=  ")\n\n";
}

void OutputWriter::ExtendStringForAdvancedLogFile(std::string& s, std::vector<std::vector<std::string>>& entry, std::string name)
{
    s += name + "\nlist("; 
    for (auto& elemA : entry)
    {
        s += "c(";
        for (auto& elemB : elemA)
            s += "'" + elemB + "',";
        if (s.back() == ',') s.pop_back();
        s += "),";
    }
    if (s.back() == ',') s.pop_back();
    s +=  ")\n\n";
}


template<typename T>
void OutputWriter::ExtendStringForAdvancedLogFile(std::string& s, std::vector<std::vector<std::vector<T>>>& entry, std::string name)
{
    s += name + "\nlist("; 
    for (auto& elemA : entry)
    {
        s += "list(";
        for (auto& elemB : elemA)
        {
            s += "c(";
            for (auto& elemC : elemB)
            {
                s += std::to_string(elemC) + ",";
            }
            if (s.back() == ',') s.pop_back();
            s += "),";
        }
        if (s.back() == ',') s.pop_back();
        s += "),";
    }
    if (s.back() == ',') s.pop_back();
    s +=  ")\n\n";
}

void OutputWriter::ExtendStringForAdvancedLogFile(std::string& s, std::vector<std::vector<std::vector<T1_locusDescription>>>& entry, std::string name)
{
    s += name + "\nlist("; 
    for (auto& elemA : entry)
    {
        s += "list(";
        for (auto& elemB : elemA)
        {
            s += "list(";
            for (auto& elemC : elemB)
                s += "c(" + std::to_string(elemC.byte_index) + "," + std::to_string(elemC.bit_index) + "," + std::to_string(elemC.locus) + "),";
            if (s.back() == ',') s.pop_back();
            s += "),";
        }
        if (s.back() == ',') s.pop_back();
        s += "),";
    }
    if (s.back() == ',') s.pop_back();
    s +=  ")\n\n";
}

void OutputWriter::ExtendStringForAdvancedLogFile(std::string& s, std::vector<FromLocusToTXLocusElement>& entry, std::string name)
{
    s += name + "\nlist("; 
    for (auto& elem : entry)
        s += "c(" + std::to_string(elem.T1) + "," + std::to_string(elem.T2) + "," + std::to_string(elem.T3) + "," + std::to_string(elem.T4) + "," + std::to_string(elem.T56ntrl) + "," + std::to_string(elem.T56sel) + "," + std::to_string(elem.TType) + "),";
    if (s.back() == ',') s.pop_back();
    s +=  ")\n\n";
}

void OutputWriter::ExtendStringForAdvancedLogFile(std::string& s, std::vector<std::pair<bool, unsigned>>& entry, std::string name)
{
    s += name + "\nlist("; 
    for (auto& elem : entry)
        s += "c(" + std::to_string(elem.first) + "," + std::to_string(elem.second) + "),";
    if (s.back() == ',') s.pop_back();
    s +=  ")\n\n";
}




bool OutputWriter::AreThereAnyOutput()
{
    bool AreThereAnyOutput = false;
    for (auto& pair : TypesToPrintOn)
    {
        OutputFileTypes type = pair.first;
        auto& files = get_OutputFiles(type);
        for (auto& file : files)
        {
            if (file.getTimes().size() != 0)
            {
                AreThereAnyOutput = true;
                break;
            }
        }
    }
    return AreThereAnyOutput;
}

bool OutputWriter::IsLastGenerationSampled()
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'IsLastGenerationSampled'\n";
#endif  
    bool IsLastGenerationSampled = false;
    for (auto& pair : TypesToPrintOn)
    {
        OutputFileTypes type = pair.first;
        auto& files = get_OutputFiles(type);
        for (auto& file : files)
        {
            if (!file.getDoesTimeNeedsToBeSet())
            {
                
                IsLastGenerationSampled = true;
                break;
            } else
            {
                assert(file.getTimes().size() > 0);
                if (file.getTimes().back() == GP->nbGenerations)
                {
                    IsLastGenerationSampled = true;
                    break;
                }
            }
        }
    }
    return IsLastGenerationSampled;
}

bool OutputWriter::ShouldSimulationBeDone(int OverwriteMode)
{
    if (OverwriteMode == 2)
    {
        SSP = nullptr;
        return true;
    }
    else if (OverwriteMode == 1)
    {
        assert(TypesToPrintOn.size() >= 0);
        if (TypesToPrintOn.size() == 1) // only the logfile
        {
            SSP = nullptr;
            return true;
        }

        for (int speciesIndex = 0 ;speciesIndex < GP->nbSpecies; speciesIndex++)
        {
            SSP =  allParameters.getSSPaddress(speciesIndex);
            
            for (auto& pair : TypesToPrintOn)
            {
                OutputFileTypes type = pair.first;
                assert(isFile(type));
                auto& files = get_OutputFiles(type);
                for (auto& file : files)
                {
                    if (file.getDoesTimeNeedsToBeSet())
                    {
                        if (!file.DoAllFilesOfTypeAlreadyExist())
                        {
                            SSP = nullptr;
                            return true;
                        }
                    }
                }
            }
        }
        SSP = nullptr;
        return false;
    }
    else if (OverwriteMode == 0)
    {
        for (int speciesIndex = 0 ;speciesIndex < GP->nbSpecies; speciesIndex++)
        {
            SSP =  allParameters.getSSPaddress(speciesIndex);

            for (auto& pair : TypesToPrintOn)
            {
                OutputFileTypes type = pair.first;
                assert(isFile(type));
                auto& files = get_OutputFiles(type);
                for (auto& file : files)
                {
                    if (file.DoesAtLeastOneFileOfTypeAlreadyExist())
                    {
                        SSP = nullptr;
                        return false;
                    }
                }
            }
        }
        SSP = nullptr;

        return true;
    } else
    {
        std::cout << "Unkown OverwriteMode discovered in 'OutputWriter::ShouldSimulationBeDone'. This is liekly an internal error but you might want to check your input data.\n";
        abort();
    }
}

void OutputWriter::RemoveAllFiles()
{
    int ForReset = GP->CurrentGeneration;
    for (auto& pair : TypesToPrintOn)
    {
        OutputFileTypes type = pair.first;
        assert(isFile(type));
        auto& files = get_OutputFiles(type);
        for (auto& file : files)
        {
            if (file.getDoesTimeNeedsToBeSet())
            {
                auto& tt = file.getTimes();
                if (tt.size() == 0)
                {
                    std::cout << "In 'OutputWriter::ClearAllFileContent', 'times' was not set for the OutputFileType " << type << ". This could well be an internal issue but you might want to check your input parameters.\n";
                    abort();
                }

                for (auto& t : tt)
                {
                    GP->CurrentGeneration = t;
                    file.remove();
                    
                }
            } else if (file.getFileType() == Logfile && this->LogfileType != 0)
            {
                file.remove();
            } else if (file.getFileType() == burnInLength_file) 
            {
                file.remove();
            }
        }
    }
    GP->CurrentGeneration = ForReset;
}

void OutputWriter::ClearAllFileContent()
{
    int ForReset = GP->CurrentGeneration;
    for (auto& pair : TypesToPrintOn)
    {
        OutputFileTypes type = pair.first;
        assert(isFile(type));
        auto& files = get_OutputFiles(type);
        for (auto& file : files)
        {
            if (file.getDoesTimeNeedsToBeSet())
            {
                auto& tt = file.getTimes();
                if (tt.size() == 0)
                {
                    std::cout << "In 'OutputWriter::ClearAllFileContent', 'times' was not set for the OutputFileType " << type << ". This could well be an internal issue but you might want to check your input parameters.\n";
                    abort();
                }
                for (auto& t : tt)
                {
                    GP->CurrentGeneration = t;
                    file.clearContent();    
                }
            } else if (this->LogfileType != 0)
            {
                file.clearContent();
            }  
        }
    }
    GP->CurrentGeneration = ForReset;
}

void OutputWriter::FinalizeLogfile()
{
    if (this->LogfileType!=0)
    {
        OutputFile& file = this->get_OutputFiles(Logfile)[0];
        file.open();
        file.write(std::string("\t\t------ Simulation is over. Good Job! ----\n"));
        file.close();
    }
}

void OutputWriter::PrintLogfile(std::string& AllInputInLongString)
{
#ifdef DEBUG
    std::cout << "Entering 'PrintLogfile'" << std::endl;
#endif    
    if (this->LogfileType!=0)
    {
        assert(this->isFile(Logfile));
        OutputFile& file = get_OutputFiles(Logfile)[0];
        file.open();

        // write logo and version
        file.write(getSimBitVersionLogo());

        // write seed
        file.write(std::string("Initial Random seed:\n\t") + std::to_string(GP->rngw.getInitialSeed()) + std::string("\n\n")); 

        // Write User Input
        file.write(std::string("User Input:\n\t") + AllInputInLongString + std::string("\n\n"));
    
        // Write extra stuff
        if (this->LogfileType==2)
        {
        #ifdef DEBUG
        std::cout << "In 'PrintLogfile', enters in if (this->LogfileType==2)" << std::endl;
        #endif 

            std::string s;
            s.reserve(10000);
        
            s += "################### General Parameters ###################\n\n";

            this->ExtendStringForAdvancedLogFile(s, OutputFile::GeneralPath, "GeneralPath");
            
            this->ExtendStringForAdvancedLogFile(s, GP->__PatchNumber,"__PatchNumber");
            this->ExtendStringForAdvancedLogFile(s, GP->__GenerationChange,"__GenerationChange");

            // nbGenerations
            this->ExtendStringForAdvancedLogFile(s, GP->nbGenerations,"nbGenerations");
            this->ExtendStringForAdvancedLogFile(s, GP->CurrentGeneration,"CurrentGeneration");

            // nbSpecies
            this->ExtendStringForAdvancedLogFile(s, GP->nbSpecies,"nbSpecies");

            assert(SSP == nullptr);

            s += "################### Species Specific Parameters ###################\n\n";

            for (int species_index = 0 ; species_index < GP->nbSpecies ; species_index++)
            {
                SSP = allParameters.getSSPaddress(species_index);

                // Demarcation
                s += "##### species: " + SSP->speciesName + "\n\n";


                // Outputs
                for (auto& pair : TypesToPrintOn)
                {
                    OutputFileTypes type = pair.first;
                    auto& files = get_OutputFiles(type);
                    for (auto& file : files) 
                    {
                        std::string path(file.getPath());
                        this->ExtendStringForAdvancedLogFile(s, path, "Output File Paths");  
                        //this->ExtendStringForAdvancedLogFile(s, get_OutputFiles(type).getTimes(), "T1_AlleleFreq_time");  
                    }
                }
                
                // Basic Demography
                this->ExtendStringForAdvancedLogFile(s, SSP->__patchCapacity,"__patchCapacity");
                this->ExtendStringForAdvancedLogFile(s, SSP->patchSize,"InitialpatchSize");
                this->ExtendStringForAdvancedLogFile(s, SSP->TotalpatchCapacity,"TotalpatchCapacity");
                this->ExtendStringForAdvancedLogFile(s, SSP->nbSubGenerationsPerGeneration,"nbSubGenerationsPerGeneration");
                this->ExtendStringForAdvancedLogFile(s, SSP->cloningRate,"cloningRate");
                this->ExtendStringForAdvancedLogFile(s, SSP->selfingRate,"selfingRate");
                this->ExtendStringForAdvancedLogFile(s, SSP->malesAndFemales,"malesAndFemales");
                this->ExtendStringForAdvancedLogFile(s, SSP->sexRatio,"sexRatio");
                
                
                // Dispersal
                this->ExtendStringForAdvancedLogFile(s, SSP->dispersalData.__forwardMigration,"dispersalData.__forwardMigration");
                this->ExtendStringForAdvancedLogFile(s, SSP->dispersalData.__forwardMigrationIndex,"dispersalData.__forwardMigrationIndex");
                {
                    auto FFM = getFullFormMatrix(SSP->dispersalData.__forwardMigration,SSP->dispersalData.__forwardMigrationIndex);
                    this->ExtendStringForAdvancedLogFile(s, FFM, "fullFormForwardMigrationMatrix");
                }
                
                
                // Genetic Map
                this->ExtendStringForAdvancedLogFile(s, SSP->NbElementsInFitnessMap,"NbElementsInFitnessMap");
                this->ExtendStringForAdvancedLogFile(s, SSP->FitnessMapProbOfEvent,"FitnessMapProbOfEvent");

                // Genetics and selection Both T1 and T2
                this->ExtendStringForAdvancedLogFile(s, SSP->ploidy,"ploidy");
                this->ExtendStringForAdvancedLogFile(s, SSP->RecombinationRate,"RecombinationRate");
                this->ExtendStringForAdvancedLogFile(s, SSP->TotalRecombinationRate,"TotalRecombinationRate [M]");
                this->ExtendStringForAdvancedLogFile(s, SSP->ChromosomeBoundaries,"ChromosomeBoundaries");
                this->ExtendStringForAdvancedLogFile(s, SSP->Gmap.TotalNbLoci,"TotalNbLoci");

                this->ExtendStringForAdvancedLogFile(s, SSP->recRateOnMismatch_bool,"recRateOnMismatch_bool");
                this->ExtendStringForAdvancedLogFile(s, SSP->recRateOnMismatch_halfWindow,"recRateOnMismatch_halfWindow");
                this->ExtendStringForAdvancedLogFile(s, SSP->recRateOnMismatch_factor,"recRateOnMismatch_factor");
                this->ExtendStringForAdvancedLogFile(s, SSP->fecundityForFitnessOfOne,"fecundityForFitnessOfOne");

                this->ExtendStringForAdvancedLogFile(s, SSP->__growthK,"__growthK");

                
                // Genetics and selection T1
                this->ExtendStringForAdvancedLogFile(s, SSP->T1_isMultiplicitySelection,"T1_isMultiplicitySelection");
                this->ExtendStringForAdvancedLogFile(s, SSP->T1_Initial_AlleleFreqs,"T1_Initial_AlleleFreqs");
                
                this->ExtendStringForAdvancedLogFile(s, SSP->Gmap.T1_nbChars,"T1_nbChars");
                this->ExtendStringForAdvancedLogFile(s, SSP->Gmap.T1_nbLoci,"T1_nbLoci");
                this->ExtendStringForAdvancedLogFile(s, SSP->T1_FitnessEffects,"T1_FitnessEffects");
                this->ExtendStringForAdvancedLogFile(s, SSP->T1_isSelection,"T1_isSelection");
                this->ExtendStringForAdvancedLogFile(s, SSP->T1_isEpistasis,"T1_isEpistasis");
                this->ExtendStringForAdvancedLogFile(s, SSP->T1_MutationRate,"T1_MutationRate");
                this->ExtendStringForAdvancedLogFile(s, SSP->T1_Total_Mutation_rate,"T1_Total_Mutation_rate");
                
                //this->ExtendStringForAdvancedLogFile(s, SSP->ProgrammedT1Mutations,"ProgrammedT1Mutations");
                //this->ExtendStringForAdvancedLogFile(s, SSP->ProgrammedT1MutationsIndexToDo,"ProgrammedT1MutationsIndexToDo");

                this->ExtendStringForAdvancedLogFile(s, SSP->T1_Epistasis_LociIndices,"T1_Epistasis_LociIndices");
                this->ExtendStringForAdvancedLogFile(s, SSP->T1_Epistasis_FitnessEffects,"T1_Epistasis_FitnessEffects");
                this->ExtendStringForAdvancedLogFile(s, SSP->Gmap.T1_nbLociLastByte,"T1_nbLociLastByte");
                this->ExtendStringForAdvancedLogFile(s, SSP->T1_Epistasis_FitnessEffects,"T1_Epistasis_FitnessEffects");

                // Genetics and selection T2
                this->ExtendStringForAdvancedLogFile(s, SSP->Gmap.T2_nbLoci,"T2_nbLoci");
                this->ExtendStringForAdvancedLogFile(s, SSP->T2_FitnessEffects,"T2_FitnessEffects");
                this->ExtendStringForAdvancedLogFile(s, SSP->T2_isSelection,"T2_isSelection");
                this->ExtendStringForAdvancedLogFile(s, SSP->T2_MutationRate,"T2_MutationRate");
                this->ExtendStringForAdvancedLogFile(s, SSP->T2_Total_Mutation_rate,"T2_Total_Mutation_rate");


                // Genetics and selection T3
                this->ExtendStringForAdvancedLogFile(s, SSP->Gmap.T3_nbLoci,"T3_nbLoci");
                this->ExtendStringForAdvancedLogFile(s, SSP->T3_PhenotypicEffects,"T3_PhenotypicEffects");
                this->ExtendStringForAdvancedLogFile(s, SSP->T3_fitnessLandscapeOptimum,"T3_fitnessLandscapeOptimum");
                this->ExtendStringForAdvancedLogFile(s, SSP->T3_fitnessLandscapeLinearGradient,"T3_fitnessLandscapeLinearGradient");
                this->ExtendStringForAdvancedLogFile(s, SSP->T3_isSelection,"T3_isSelection");
                this->ExtendStringForAdvancedLogFile(s, SSP->T3_MutationRate,"T3_MutationRate");
                this->ExtendStringForAdvancedLogFile(s, SSP->T3_Total_Mutation_rate,"T3_Total_Mutation_rate");
                this->ExtendStringForAdvancedLogFile(s, SSP->T3_fitnessLandscapeType,"T3_fitnessLandscapeType");
                this->ExtendStringForAdvancedLogFile(s, SSP->T3_PhenoNbDimensions,"T3_PhenoNbDimensions");

                // Genetics T4
                this->ExtendStringForAdvancedLogFile(s, SSP->Gmap.T4_nbLoci,"T4_nbLoci");
                this->ExtendStringForAdvancedLogFile(s, SSP->T4_simplifyEveryNGenerations,"T4_simplifyEveryNGenerations");
                this->ExtendStringForAdvancedLogFile(s, SSP->T4_MutationRate,"T4_MutationRate");


                // Genetics T5
                this->ExtendStringForAdvancedLogFile(s, SSP->Gmap.T5_nbLoci,"T5_nbLoci");
                this->ExtendStringForAdvancedLogFile(s, SSP->Gmap.T5ntrl_nbLoci,"T5ntrl_nbLoci");
                this->ExtendStringForAdvancedLogFile(s, SSP->Gmap.T5sel_nbLoci,"T5sel_nbLoci");
                this->ExtendStringForAdvancedLogFile(s, SSP->Gmap.T6_nbLoci,"T6_nbLoci");
                this->ExtendStringForAdvancedLogFile(s, SSP->Gmap.T6ntrl_nbLoci,"T6ntrl_nbLoci");
                this->ExtendStringForAdvancedLogFile(s, SSP->Gmap.T6sel_nbLoci,"T6sel_nbLoci");
                this->ExtendStringForAdvancedLogFile(s, SSP->Gmap.T56_nbLoci,"T56_nbLoci");
                this->ExtendStringForAdvancedLogFile(s, SSP->Gmap.T56ntrl_nbLoci,"T56ntrl_nbLoci");
                this->ExtendStringForAdvancedLogFile(s, SSP->Gmap.T56sel_nbLoci,"T56sel_nbLoci");
                this->ExtendStringForAdvancedLogFile(s, SSP->T56_FitnessEffects,"T56_FitnessEffects");
                this->ExtendStringForAdvancedLogFile(s, SSP->T56_MutationRate,"T56_MutationRate");
                this->ExtendStringForAdvancedLogFile(s, SSP->T56_isSelection,"T56_isSelection");
                this->ExtendStringForAdvancedLogFile(s, SSP->T56_isMultiplicitySelection,"T56_isMultiplicitySelection");


                // Ecology
                this->ExtendStringForAdvancedLogFile(s, SSP->__Habitats,"__Habitats");
                this->ExtendStringForAdvancedLogFile(s, SSP->MaxHabitat,"MaxHabitat");

                // Other
                this->ExtendStringForAdvancedLogFile(s, SSP->centralT1LocusForExtraGeneticInfo,"centralT1LocusForExtraGeneticInfo");
                this->ExtendStringForAdvancedLogFile(s, SSP->subsetT1LociForfitnessSubsetLoci_file,"subsetT1LociForfitnessSubsetLoci_file");
                this->ExtendStringForAdvancedLogFile(s, SSP->subsetT2LociForfitnessSubsetLoci_file,"subsetT2LociForfitnessSubsetLoci_file");
                this->ExtendStringForAdvancedLogFile(s, SSP->subsetT3LociForfitnessSubsetLoci_file,"subsetT3LociForfitnessSubsetLoci_file");
                this->ExtendStringForAdvancedLogFile(s, SSP->subsetT1epistasisLociForfitnessSubsetLoci_file,"subsetT1epistasisLociForfitnessSubsetLoci_file");

            }

            SSP = nullptr;
            s += "###################\n\n";

            file.write(s);
        }

        if (this->LogfileType==102105116)
        {
            // fitness info only
            std::string s;
            s.reserve(1000);
        
            s += "################### General Parameters ###################\n\n";

            this->ExtendStringForAdvancedLogFile(s, OutputFile::GeneralPath, "GeneralPath");
            
            this->ExtendStringForAdvancedLogFile(s, GP->__PatchNumber,"__PatchNumber");
            this->ExtendStringForAdvancedLogFile(s, GP->__GenerationChange,"__GenerationChange");

            // nbGenerations
            this->ExtendStringForAdvancedLogFile(s, GP->nbGenerations,"nbGenerations");
            this->ExtendStringForAdvancedLogFile(s, GP->CurrentGeneration,"CurrentGeneration");

            // nbSpecies
            this->ExtendStringForAdvancedLogFile(s, GP->nbSpecies,"nbSpecies");

            assert(SSP == nullptr);

            s += "################### Species Specific Parameters ###################\n\n";

            for (int species_index = 0 ; species_index < GP->nbSpecies ; species_index++)
            {
                SSP = allParameters.getSSPaddress(species_index);

                // Demarcation
                s += "##### species: " + SSP->speciesName + "\n\n";

                // Genetics 
                this->ExtendStringForAdvancedLogFile(s, SSP->Gmap.TotalNbLoci,"TotalNbLoci");

                
                // Genetics and selection T1
                this->ExtendStringForAdvancedLogFile(s, SSP->T1_isMultiplicitySelection,"T1_isMultiplicitySelection");
                
                this->ExtendStringForAdvancedLogFile(s, SSP->Gmap.T1_nbChars,"T1_nbChars");
                this->ExtendStringForAdvancedLogFile(s, SSP->Gmap.T1_nbLoci,"T1_nbLoci");
                this->ExtendStringForAdvancedLogFile(s, SSP->T1_FitnessEffects,"T1_FitnessEffects");
                this->ExtendStringForAdvancedLogFile(s, SSP->T1_isSelection,"T1_isSelection");
                this->ExtendStringForAdvancedLogFile(s, SSP->T1_isEpistasis,"T1_isEpistasis");
                
                //this->ExtendStringForAdvancedLogFile(s, SSP->ProgrammedT1Mutations,"ProgrammedT1Mutations");
                //this->ExtendStringForAdvancedLogFile(s, SSP->ProgrammedT1MutationsIndexToDo,"ProgrammedT1MutationsIndexToDo");

                this->ExtendStringForAdvancedLogFile(s, SSP->T1_Epistasis_LociIndices,"T1_Epistasis_LociIndices");
                this->ExtendStringForAdvancedLogFile(s, SSP->T1_Epistasis_FitnessEffects,"T1_Epistasis_FitnessEffects");
                this->ExtendStringForAdvancedLogFile(s, SSP->Gmap.T1_nbLociLastByte,"T1_nbLociLastByte");
                this->ExtendStringForAdvancedLogFile(s, SSP->T1_Epistasis_FitnessEffects,"T1_Epistasis_FitnessEffects");

                // Genetics and selection T2
                this->ExtendStringForAdvancedLogFile(s, SSP->Gmap.T2_nbLoci,"T2_nbLoci");
                this->ExtendStringForAdvancedLogFile(s, SSP->T2_FitnessEffects,"T2_FitnessEffects");
                this->ExtendStringForAdvancedLogFile(s, SSP->T2_isSelection,"T2_isSelection");


                // Genetics and selection T3
                this->ExtendStringForAdvancedLogFile(s, SSP->Gmap.T3_nbLoci,"T3_nbLoci");
                this->ExtendStringForAdvancedLogFile(s, SSP->T3_PhenotypicEffects,"T3_PhenotypicEffects");
                this->ExtendStringForAdvancedLogFile(s, SSP->T3_fitnessLandscapeOptimum,"T3_fitnessLandscapeOptimum");
                this->ExtendStringForAdvancedLogFile(s, SSP->T3_fitnessLandscapeLinearGradient,"T3_fitnessLandscapeLinearGradient");
                this->ExtendStringForAdvancedLogFile(s, SSP->T3_isSelection,"T3_isSelection");
                this->ExtendStringForAdvancedLogFile(s, SSP->T3_fitnessLandscapeType,"T3_fitnessLandscapeType");
                this->ExtendStringForAdvancedLogFile(s, SSP->T3_PhenoNbDimensions,"T3_PhenoNbDimensions");

                // Genetics T4
                this->ExtendStringForAdvancedLogFile(s, SSP->Gmap.T4_nbLoci,"T4_nbLoci");

                // Genetics T5
                this->ExtendStringForAdvancedLogFile(s, SSP->T56_FitnessEffects,"T56_FitnessEffects");
                this->ExtendStringForAdvancedLogFile(s, SSP->T56_isSelection,"T56_isSelection");
                this->ExtendStringForAdvancedLogFile(s, SSP->T56_isMultiplicitySelection,"T56_isMultiplicitySelection");


                // Ecology
                this->ExtendStringForAdvancedLogFile(s, SSP->__Habitats,"__Habitats");
                this->ExtendStringForAdvancedLogFile(s, SSP->MaxHabitat,"MaxHabitat");
            }

            SSP = nullptr;
            s += "###################\n\n";

            file.write(s);
        }
        
    file.close();
    }
}

void OutputWriter::insertOutputFile(OutputFile&& file)
{
    auto iter = TypesToPrintOn.find(file.getFileType());
    
    /*std::cout << "In TypesToPrintOn, there is: ";
    for( typename std::map<const OutputFileTypes, std::vector<OutputFile>>::const_iterator it = TypesToPrintOn.begin(); it != TypesToPrintOn.end(); ++it )
    {
        std::cout << it->first << " ";
    }
    std::cout << std::endl;*/
    //std::cout << "Try to insert outputFile of type " << file.getFileType() << ".\n";
    file.assertSubsetSize();

    if (iter == TypesToPrintOn.end())
    {
        //std::cout << "there is no file of this type of the moment.\n";        
        std::vector<OutputFile> v;
        v.push_back(std::move(file));
        TypesToPrintOn.insert(
            std::pair<const OutputFileTypes, std::vector<OutputFile>>(file.getFileType(), std::move(v))
        );
    } else
    {
        //std::cout << "there is a file of this type already exists.\n";
        (*iter).second.push_back(std::move(file));
    }
    (void) TypesToPrintOn.at(file.getFileType());
}

bool OutputWriter::isTime()
{
    if (KillOnDemand::justAboutToKill)
    {
        if (AllTimes.back() > GP->CurrentGeneration)
        {
            return true;
        } else
        {
            return false;
        }
    } else
    {
        return std::find(AllTimes.begin(), AllTimes.end(),  GP->CurrentGeneration) != AllTimes.end();
    }
}

void OutputWriter::SetAllTimes()
{
    assert(AllTimes.size() == 0);
    AllTimes.push_back(0);
    for (auto& pair : TypesToPrintOn)
    {
        OutputFileTypes type = pair.first;
        auto& files = get_OutputFiles(type);
        for (auto& file : files)
        {
            if (file.getDoesTimeNeedsToBeSet())
            {
                auto& tt = file.getTimes();
                if (tt.size() == 0)
                {
                    std::cout << "In 'OutputWriter::SetAllTimes', 'times' was not set for the OutputFileType " << type << ". This could well be an internal issue but you might want to check your input parameters.\n";
                    abort();
                }
                for (auto& t : tt)
                {
                    AllTimes.push_back(t);
                }
            }
        }
    }

    // sort and remove duplicates
    std::sort( AllTimes.begin(), AllTimes.end() );
    AllTimes.erase( unique( AllTimes.begin(), AllTimes.end() ), AllTimes.end() );

    if (AllTimes.back() != GP->nbGenerations)
    {
        AllTimes.push_back(GP->nbGenerations);
    }

    //shrink to fit
    AllTimes.shrink_to_fit();
}

bool OutputWriter::isFile(OutputFileTypes type)
{
    return TypesToPrintOn.find(type) != TypesToPrintOn.end();
}

std::vector<OutputFile>& OutputWriter::get_OutputFiles(OutputFileTypes type)
{
    /*if (!isFile(type))
    {
        std::cout << "Internal Error: In 'OutputWriter::get_OutputFiles' failed to get the requested type.\n";
        abort();
    }*/
    /*
    if (type == T1_vcfFile)
    {
        std::cout << "in 'OutputWriter::get_OutputFiles' at generation " << GP->CurrentGeneration  << std::endl;
    }
    */

    return TypesToPrintOn.at(type);
}

void OutputWriter::SetBeforeInitializationTime()
{
    clock_BeforeInitialization = clock();
}

void OutputWriter::SetBeforeSimulationTime()
{
    clock_BeforeSimulation = clock();
}

void OutputWriter::SetAfterSimulationTime()
{
    clock_AfterSimulation = clock();
}

void OutputWriter::PrintInitializationTimeToLogFile()
{
    long double nbSeconds = ((long double) clock_BeforeSimulation - (long double) clock_BeforeInitialization) / (long double) CLOCKS_PER_SEC;
    printTime(nbSeconds, std::string("Real time for initialization: "));    
}
void OutputWriter::PrintSimulationTimeToLogFile()
{
    long double nbSeconds = ((long double) clock_AfterSimulation - (long double) clock_BeforeSimulation) / (long double) CLOCKS_PER_SEC;
    printTime(nbSeconds, std::string("Real time for the simulation: "));
}


void OutputWriter::printTime(long double seconds, std::string message)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'printTime'\n";
#endif      

    int intMinutes = (int) floor(seconds / 60);
    int intHours   = floor(intMinutes / 60);
    intMinutes     -= 60 * intHours;
    int intSeconds = (int) floor((seconds - 3600 * intHours - 60 * intMinutes) + 0.5);
    assert(intMinutes >= 0 && intMinutes <= 60);
    assert(intSeconds >= 0 && intSeconds <= 60);
    std::string s("\t\t");
    if (intHours > 0)
    {
        s += message + std::to_string(intHours) + std::string(" hours ") + std::to_string(intMinutes) + std::string(" minutes ")  + std::to_string(intSeconds) + std::string(" seconds ");
    } else if (intMinutes > 0)
    {
        s += message + std::to_string(intMinutes) + std::string(" minutes ")  + std::to_string(intSeconds) + std::string(" seconds ");
    } else
    {
        s += message + std::to_string(intSeconds) + std::string(" seconds ");
    }
    s += "(" + to_string_with_precision(seconds, 4) + " seconds)\n\n";

    //Print on Console
    std::cout << s;

    //Print on Logfile
    if (this->LogfileType != 0)
    {
        assert(this->isFile(Logfile));
        OutputFile& file = this->get_OutputFiles(Logfile)[0];
        file.open();
        file.write(s);
        file.close();
    }
}

void OutputWriter::WriteOutputs_extraGeneticInfo(OutputFile& file)
{

    assert(SSP->MaxHabitat == 0); // Only makes sense in absence of local selection

    file.open();
    std::string s;


    /*
    ##################
    ######## B #######
    ##################
    */
        

    // Copy SSP->RecombinationRate. This is done because SSP->RecombinationRate could be of length 1 even if TotalNbLoci>1 if they all have the same recombination rate (I think)
    std::vector<double> PRecRatCopy;
    PRecRatCopy.resize(SSP->Gmap.TotalNbLoci-1);
    if (SSP->RecombinationRate.size() > 1)
    {
        PRecRatCopy = SSP->RecombinationRate;
    } else
    {
        for (int interlocus = 1 ; interlocus < SSP->RecombinationRate.size() ; interlocus++)
        {
            PRecRatCopy[interlocus] = SSP->RecombinationRate[0];
        }
    }


    std::vector<double> T1MR; // non cumulative mutation rate
    if (SSP->Gmap.T1_nbLoci > 0)
    {
        T1MR.resize(SSP->T1_MutationRate.size());
        T1MR[0] = SSP->T1_MutationRate[0];
        for (int locus = 1 ; locus < SSP->T1_MutationRate.size() ; locus++)
        {
            T1MR[locus] = SSP->T1_MutationRate[locus] - SSP->T1_MutationRate[locus-1];
            assert(T1MR[locus] >= 0);
            /*if (T1MR[locus] < 0.0)
            {
                assert(T1MR[locus] > -0.00000001);
                T1MR[locus] = 0.0;
            }*/
        }
    }
    assert(T1MR.size() == SSP->Gmap.T1_nbLoci);
    assert(SSP->Gmap.FromLocusToNextT1Locus(SSP->Gmap.TotalNbLoci-1) == SSP->Gmap.T1_nbLoci);


    std::vector<double> T2MR; // non cumulative mutation rate
    if (SSP->Gmap.T2_nbLoci > 0)
    {
        T2MR.resize(SSP->T2_MutationRate.size());
        T2MR[0] = SSP->T2_MutationRate[0];
        for (int locus = 1 ; locus < SSP->T2_MutationRate.size()-1 ; locus++)
        {
            T2MR[locus] = SSP->T2_MutationRate[locus] - SSP->T2_MutationRate[locus-1];
            assert(T2MR[locus] >= 0);
            /*if (T2MR[locus] < 0.0)
            {
                assert(T2MR[locus] > -0.00000001);
                T2MR[locus] = 0.0;
            }*/
        }
    }
    assert(T2MR.size() == SSP->Gmap.T2_nbLoci);
    assert(SSP->Gmap.FromLocusToNextT2Locus(SSP->Gmap.TotalNbLoci-1) == SSP->Gmap.T2_nbLoci );

    std::vector<double> T3MR; // non cumulative mutation rate
    if (SSP->Gmap.T3_nbLoci > 0)
    {
        T3MR.resize(SSP->T3_MutationRate.size());
        T3MR[0] = SSP->T3_MutationRate[0];
        for (int locus = 1 ; locus < SSP->T3_MutationRate.size() ; locus++)
        {
            T2MR[locus] = SSP->T3_MutationRate[locus] - SSP->T3_MutationRate[locus-1];
            assert(T3MR[locus] >= 0);
            /*if (T3MR[locus] < 0.0)
            {
                assert(T3MR[locus] > -0.00000001);
                T3MR[locus] = 0.0;
            }*/
        }
    }
    assert(T3MR.size() == SSP->Gmap.T3_nbLoci);
    assert(SSP->Gmap.FromLocusToNextT3Locus(SSP->Gmap.TotalNbLoci-1) == SSP->Gmap.T3_nbLoci);
     /*   
    std::cout << "T1MR[" << 0 << "] = " << T1MR[0] << "\n";
    std::cout << "T1MR[" << 1 << "] = " << T1MR[1] << "\n";
    std::cout << "T1MR[" << 2 << "] = " << T1MR[2] << "\n";
    std::cout << "T1MR[" << 3 << "] = " << T1MR[3] << "\n";
    std::cout << "T1MR[" << T1MR.size()-4 << "] = " << T1MR[T1MR.size()-4] << "\n";
    std::cout << "T1MR[" << T1MR.size()-3 << "] = " << T1MR[T1MR.size()-3] << "\n";
    std::cout << "T1MR[" << T1MR.size()-2 << "] = " << T1MR[T1MR.size()-2] << "\n";
    std::cout << "T1MR[" << T1MR.size()-1 << "] = " << T1MR[T1MR.size()-1] << "\n";

    std::cout << "\n\n\n";*/

    double B_HudsonKaplan_mean;
    double B_NordborgEtAl_mean;

    if (!SSP->T1_isSelection && !SSP->T1_isEpistasis && !SSP->T2_isSelection && !SSP->T3_isSelection )
    {
        B_HudsonKaplan_mean = 1.0;
        B_NordborgEtAl_mean = 1.0;
    } else
    {
        B_HudsonKaplan_mean = 0.0;
        B_NordborgEtAl_mean = 0.0;

        for (auto& centralLocusObject : file.getSubsetToConsider(SSP->speciesIndex))
        {
            int centralLocus = centralLocusObject.locus;

             // Start computing E
            double E_HudsonKaplan = 0.0;
            double E_NordborgEtAl = 0.0;

            int previous_T1Locus = -1;
            int previous_T2Locus = -1;
            int previous_T3Locus = -1;

            double currentDistanceIncM = 0.0;
            

            assert(centralLocus >= 0 && centralLocus <= SSP->Gmap.TotalNbLoci);
            for (int locus = 0 ; locus < SSP->Gmap.TotalNbLoci ; locus++)
            { 
                // compute distance
                currentDistanceIncM = std::abs(PRecRatCopy[locus] - PRecRatCopy[centralLocus]);
                

        /*        if (std::abs(locus - centralLocus) < 8)
                {
                    std::cout << "central is " << locus - centralLocus << " compared 'centralLocus' (negative is left of centralLocus)\n";
                    //if (centralIsOnTheRight) {} else {std::cout << "central is on the left\n";}
                    std::cout << "currentDistanceIncM = " << currentDistanceIncM << "\n";
                }
                    */


                // Find out what locus type it is and compute its appropriate fitness effect
                int T1_locus = SSP->Gmap.FromLocusToNextT1Locus(locus) - 1;
                int T2_locus = SSP->Gmap.FromLocusToNextT2Locus(locus) - 1;
                int T3_locus = SSP->Gmap.FromLocusToNextT3Locus(locus) - 1;


                //std::cout << "T1_locus = " << T1_locus << " T2_locus = " << T2_locus << " T3_locus = " << T3_locus << "\n";

                // Get value of 't' (heterozygote selection coefficient) and 'currentMutationRate'
                double currentMutationRate;
                double t;
                if (T1_locus > previous_T1Locus)
                {
                    assert(T1_locus >= 0);
                    assert(previous_T1Locus + 1 == T1_locus);
                    assert(previous_T2Locus == T2_locus);
                    assert(previous_T3Locus == T3_locus);
                    previous_T1Locus = T1_locus;

                    assert(T1_locus <= T1MR.size());
                    currentMutationRate = T1MR[T1_locus];

                    //std::cout << "SSP->T1_isMultiplicitySelection = " << SSP->T1_isMultiplicitySelection << std::endl;
                    if (SSP->T1_isMultiplicitySelection)
                    {
                        //std::cout << " SSP->T1_FitnessEffects[0][locus] = " << SSP->T1_FitnessEffects[0][locus] << std::endl;
                        t = 1.0 - SSP->T1_FitnessEffects[0][T1_locus];
                    } else
                    {
                        //std::cout << " SSP->T1_FitnessEffects[0][THREE * locus + 1] = " << SSP->T1_FitnessEffects[0][THREE * locus + 1] << std::endl;
                        t = 1.0 - SSP->T1_FitnessEffects[0][THREE * T1_locus + 1];
                    }
                } else if (T2_locus > previous_T2Locus)
                {
                    assert(T2_locus >= 0);
                    assert(previous_T1Locus == T1_locus);
                    assert(previous_T2Locus + 1 == T2_locus);
                    assert(previous_T3Locus == T3_locus);
                    previous_T2Locus = T2_locus;

                    assert(T2_locus <= T2MR.size());
                    currentMutationRate = T2MR[T2_locus];

                    //std::cout << " SSP->T2_FitnessEffects[0][locus] = " << SSP->T2_FitnessEffects[0][T2_locus] << std::endl;
                    t = 1 - SSP->T2_FitnessEffects[0][T2_locus];
                } else if (T3_locus > previous_T3Locus)
                {
                    assert(T3_locus >= 0);
                    assert(previous_T1Locus == T1_locus);
                    assert(previous_T2Locus == T2_locus);
                    assert(previous_T3Locus + 1 == T3_locus);
                    previous_T3Locus = T3_locus;

                    assert(T3_locus <= T3MR.size());
                    currentMutationRate = T3MR[T3_locus];

                    std::cout << "Oops 'extraGeneticInfo', when computing expected effect of background selection is undefined for the moment for loci of type 3\n";
                    abort();
                } else
                {
                    std::cout << "Internal error. Should not enter here code 47954";
                    abort();
                }
                
                // compute the expected effect for this locus
                //std::cout << "t = " << t << " SSP->T1_MutationRate[locus] = " << SSP->T1_MutationRate[locus] << " currentDistanceIncM = " << currentDistanceIncM << std::endl;
                assert(t >= 0.0);
                /*if (t < 0.0)
                {
                    std::cout << "t = " << t << std::endl;
                    assert(t > -0.0000001);
                    t = 0.0;
                }*/

                
                assert(currentMutationRate >= 0.0);
                assert(currentDistanceIncM >= 0.0);
                double currentDistanceInRate = (1 - exp(-2*currentDistanceIncM)) / 2;
                assert(currentDistanceInRate >= 0.0);
                assert(currentDistanceInRate <= 0.5);
                assert(t <= 1.0);
                assert(t >= 0.0);

                if (t > 0.0)
                {
                    double addToE_HudsonKaplan;
                    double addToE_NordborgEtAl;

                    // E_HudsonKaplan
                    if (E_HudsonKaplan > 9900000)
                    {
                        addToE_HudsonKaplan = 0.0;
                    } else
                    {
                        addToE_HudsonKaplan = (currentMutationRate/2 * t) / (2*(pow(t + currentDistanceInRate, 2)));
                    }

                    // E_NordborgEtAl
                    if (E_NordborgEtAl > 9900000)
                    {
                        addToE_NordborgEtAl = 0.0;
                    } else
                    {
                        addToE_NordborgEtAl = (currentMutationRate) / (t * pow(1 + (1-t) * currentDistanceInRate / t, 2));
                    }

                    // Add
                    E_HudsonKaplan += addToE_HudsonKaplan;
                    E_NordborgEtAl += addToE_NordborgEtAl;
                }
            }
            assert(E_HudsonKaplan >= 0.0);
            assert(E_NordborgEtAl >= 0.0);

            // B is the expected nucleotide diversity at a focal netural site j, relative to its value in the absence of selection.
            B_HudsonKaplan_mean               += exp(-E_HudsonKaplan);
            B_NordborgEtAl_mean               += exp(-E_NordborgEtAl);
        }
        B_HudsonKaplan_mean                 /= file.getSubsetToConsider(SSP->speciesIndex).size();
        B_NordborgEtAl_mean                 /= file.getSubsetToConsider(SSP->speciesIndex).size();
    }
    
       
    assert(B_HudsonKaplan_mean <= 1.0 && B_HudsonKaplan_mean >= 0.0);
    assert(B_NordborgEtAl_mean <= 1.0 && B_NordborgEtAl_mean >= 0.0);

    s += std::string("B_HudsonKaplan = ") + std::to_string(B_HudsonKaplan_mean) + std::string("\n");
    s += std::string("B_NordborgEtAl = ") + std::to_string(B_NordborgEtAl_mean) + std::string("\n");


    /*
    ##################
    ####### mu #######
    ##################
    */

    double averageMuOverSubsetToConsider;
    if (SSP->Gmap.T1_nbLoci > 0)
    {
        averageMuOverSubsetToConsider = SSP->T1_MutationRate[file.getSubsetToConsider(SSP->speciesIndex)[file.getSubsetToConsider(SSP->speciesIndex).size()-2].locus] - SSP->T1_MutationRate[file.getSubsetToConsider(SSP->speciesIndex)[0].locus];
        averageMuOverSubsetToConsider /= file.getSubsetToConsider(SSP->speciesIndex).size()-2;
        //std::cout << averageMuOverSubsetToConsider << "\n";
    } else
    {
        averageMuOverSubsetToConsider = -1.0;
    }   
    
    std::ostringstream s_tmpForPrecision;
    s_tmpForPrecision << std::setprecision(50) << averageMuOverSubsetToConsider;
    s += std::string("averageMuOverSubsetToConsider = ") + s_tmpForPrecision.str() + std::string("\n");
    //s += std::string("averageMuOverSubsetToConsider = ") + std::to_string(averageMuOverSubsetToConsider) + std::string("\n");


    file.write(s);
    file.close();
}



void OutputWriter::WriteOutputs_Tx_FST_header(OutputFile& file)
{
    const std::string s_tab("\t");
    std::string s("Generation\tFST_WeirCockerham_ratioOfAverages\tFST_Nei_ratioOfAverages\tFST_WeirCockerham_averageOfRatios\tFST_Nei_averageOfRatios\n");
    file.open();
    file.write(s);
    file.close();
}

void OutputWriter::WriteOutputs_T1_FST_header(OutputFile& file)
{
    WriteOutputs_Tx_FST_header(file);
}
/*
void OutputWriter::WriteOutputs_T5_FST_header(OutputFile& file)
{
    WriteOutputs_Tx_FST_header(OutputFile& file);
}*/

void OutputWriter::writeOutputs_burnInLength(int positiveGeneration)
{
    std::vector<OutputFile>& files = this->get_OutputFiles(burnInLength_file);
    assert(files.size());
    auto& file = files[0];

    file.open();
    std::string s = std::to_string(positiveGeneration) + std::string("\n");
    file.write(s);
    file.close();
}


void OutputWriter::WriteOutputs_T1_FST(Pop& pop, OutputFile& file)
{
    std::string s = std::to_string(GP->CurrentGeneration);
    std::string s_tab("\t");
    auto data = pop.computeT1RelativeFrequenciesPerPatch();

    /*
    for (auto& dataPatch : data)
    {
        for (auto& dataLocus : dataPatch)
        {
            std::cout << dataLocus << " ";
        }   
        std::cout << "\n";
    }
    */

    auto lociToConsider = file.subset[SSP->speciesIndex];


    double sumNumeratorNei = 0.0;
    double sumNumeratorWC = 0.0;
    double sumDenominatorNei = 0.0;
    double sumDenominatorWC = 0.0;
    double sumRatiosNei = 0.0;
    double sumRatiosWC = 0.0;
    assert(lociToConsider.size());
    for (auto& locusDescription : lociToConsider)
    {
        auto locus = locusDescription.locus;

        // Compute Hs, Ht and Hb for the specific locus
        double Hs = 0.0;
        double meanFreq = 0.0;
        double Hb = 0.0;
        uint32_t nbNonEmptyPatches = 0;
        uint32_t HbNbComparisons = 0;
        for (uint32_t patch_index = 0 ; patch_index < GP->PatchNumber ; patch_index++)
        {
            if (SSP->patchSize[patch_index])
            {
                ++nbNonEmptyPatches;
                auto x = data[patch_index][locus];
                assert(x >= 0.0 && x <= 1.0);
                Hs += 2 * x * (1 - x);
                meanFreq += x;

                for (uint32_t other_patch_index = patch_index+1 ; other_patch_index < GP->PatchNumber ; other_patch_index++)
                {
                    if (SSP->patchSize[other_patch_index])
                    {
                        ++HbNbComparisons;
                        auto other_x = data[other_patch_index][locus];
                        assert(other_x >= 0.0 && other_x <= 1.0);
                        Hb += other_x * (1 - x) + x * (1 - other_x);
                    }
                }
            }
        }

        //std::cout << "locus = " << locus << "\n";
        //std::cout << "\tnbNonEmptyPatches = " << nbNonEmptyPatches << "\n";
        //std::cout << "\tHbNbComparisons = " << HbNbComparisons << "\n";
        

        double Ht;
        if (nbNonEmptyPatches == 0)
        {
            Hs = std::numeric_limits<double>::quiet_NaN();
            Ht = std::numeric_limits<double>::quiet_NaN();
        } else
        {
            Hs /= nbNonEmptyPatches;
            assert(Hs >= 0.0 && Hs <= 1.0);

            meanFreq /= nbNonEmptyPatches;
            assert(meanFreq >= 0.0 && meanFreq <= 1.0);
            Ht = 2 * meanFreq * (1-meanFreq);
            assert(Ht >= 0.0 && Ht <= 1.0);
        }

        if (HbNbComparisons == 0)
        {
            Hb = std::numeric_limits<double>::quiet_NaN();
        } else
        {
            Hb /= HbNbComparisons;
            assert(Hb >= 0.0 && Hb <= 1.0);
        }
        

        //std::cout << "\tHs = " << Hs << "\n";
        //std::cout << "\tHt = " << Ht << "\n";
        //std::cout << "\tHb = " << Hb << "\n";

        // Compute stats to compute FST
        sumNumeratorNei += Ht-Hs;
        sumNumeratorWC += Hb-Hs;
        sumDenominatorNei += Ht;
        sumDenominatorWC += Hb;

        if (Ht != 0.0)
            sumRatiosNei += (Ht-Hs)/Ht;
        
        if (Hb != 0.0)
            sumRatiosWC += (Hb-Hs)/Hb;
    }


    double FST_Nei_averageOfRatios = sumRatiosNei / lociToConsider.size();
    double FST_WC_averageOfRatios = sumRatiosWC / lociToConsider.size();
    double FST_Nei_ratioOfAverages = sumNumeratorNei / sumDenominatorNei;
    double FST_WC_ratioOfAverages = sumNumeratorWC / sumDenominatorWC;

    s += s_tab + std::to_string(FST_WC_ratioOfAverages) + s_tab + std::to_string(FST_Nei_ratioOfAverages) + s_tab + std::to_string(FST_WC_averageOfRatios) + s_tab + std::to_string(FST_Nei_averageOfRatios) + "\n"; 

    file.open();
    file.write(s);
    file.close();
}

/*void OutputWriter::WriteOutputs_T5_FST(Pop& pop, OutputFile& file)
{
    std::cout << "Sorry outputs for T5_FST have not been written yet\n";
    abort();
}*/

void OutputWriter::WriteOutputs_T1_FST_header_complex(OutputFile& file)
{
    assert(GP->output_FST_nbPatchesToConsider.size() > 0);

    file.open();

    std::string s;

    
    s += std::string("Generation");

    for (int nbPatchesToConsider : GP->output_FST_nbPatchesToConsider)
    {
        std::string patchesToConsider(nbPatchesToConsider,1);   // nbPatchesToConsider leading 1's. nbPatchesToConsider is K in binomial coef
        patchesToConsider.resize(GP->maxEverPatchNumber,0);     // GP->maxEverPatchNumber-nbPatchesToConsider trailing 0's. GP->maxEverPatchNumber is N in binomial coef

        do
        {
            // Get indices of patches to consider
            std::vector<int> patchesToConsiderIndices;
            for (int patch_index = 0 ; patch_index < patchesToConsider.size() ;patch_index++)
            {
                if (patchesToConsider[patch_index])
                {
                    patchesToConsiderIndices.push_back(patch_index);
                }
            }
            assert(patchesToConsiderIndices.size() == nbPatchesToConsider);

            std::string patchCombinationString("\t");
            for (int patch_index_index =0 ; patch_index_index < (patchesToConsiderIndices.size() - 1); patch_index_index++)
            {
                int patch_index = patchesToConsiderIndices[patch_index_index];
                patchCombinationString += std::string("P") + std::to_string(patch_index);
            }
            patchCombinationString += std::string("P") + std::to_string(patchesToConsiderIndices.back());

            s += patchCombinationString + std::string("_averageOfRatios_GST");
            s += patchCombinationString + std::string("_averageOfRatios_WCST");
            s += patchCombinationString + std::string("_ratioOfAverages_GST");
            s += patchCombinationString + std::string("_ratioOfAverages_WCST");
        } while (std::prev_permutation(patchesToConsider.begin(), patchesToConsider.end()));
    }

    s += std::string("\n");
    
    file.write(s);
    file.close();
}

void OutputWriter::WriteOutputs_T1_FST_complex(Pop& pop, OutputFile& file)
{
    file.open();
    const std::string s_tab("\t");
    std::string s;
    s += std::to_string(GP->CurrentGeneration);

    // compute alleleFreq at each patch
    std::vector<std::vector<double>> alleleFreqs;  // alleleFreqs[patch][locus]
    alleleFreqs.resize(GP->PatchNumber);
    auto polymorphicLoci = file.removeSitesWeDontWant(pop.listT1PolymorphicLoci(), SSP->speciesIndex);
    
    for (int patch_index = 0 ; patch_index < GP->PatchNumber ; patch_index++)
    {
        alleleFreqs[patch_index].resize(polymorphicLoci.size());


        for (uint32_t SNPindex = 0 ; SNPindex < polymorphicLoci.size() ; SNPindex++)
        {
            int byte_index = polymorphicLoci[SNPindex].byte_index;
            int bit_index = polymorphicLoci[SNPindex].bit_index;

            if (SSP->patchSize[patch_index] == 0)
            {
                alleleFreqs[patch_index][SNPindex] += std::numeric_limits<double>::quiet_NaN();
            } else
            {
                alleleFreqs[patch_index][SNPindex] = 0.0;
                for (uint32_t ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ind_index++)
                {
                    alleleFreqs[patch_index][SNPindex] += pop.getPatch(patch_index).getInd(ind_index).getHaplo(0).getT1_Allele(byte_index,bit_index);
                    alleleFreqs[patch_index][SNPindex] += pop.getPatch(patch_index).getInd(ind_index).getHaplo(1).getT1_Allele(byte_index,bit_index);
                }
                alleleFreqs[patch_index][SNPindex] /= 2 * SSP->patchSize[patch_index];
                assert(alleleFreqs[patch_index][SNPindex] >= 0.0 && alleleFreqs[patch_index][SNPindex] <= 1.0);
            }
        }
    }
    assert(alleleFreqs.size() == GP->PatchNumber);


    // Compute all FSTs
    for (int nbPatchesToConsider : GP->output_FST_nbPatchesToConsider)
    {
        std::string patchesToConsider(nbPatchesToConsider,1);   // nbPatchesToConsider leading 1's. nbPatchesToConsider is K in binomial coef
        patchesToConsider.resize(GP->maxEverPatchNumber,0);     // GP->maxEverPatchNumber-nbPatchesToConsider trailing 0's. GP->maxEverPatchNumber is N in binomial coef

        do
        {
            /////////////////
            //// Get indices of patches to consider
            /////////////////
            std::vector<int> patchesToConsiderIndices;
            for (int patch_index = 0 ; patch_index < patchesToConsider.size() ;patch_index++)
            {
                if (patchesToConsider[patch_index])
                {
                    patchesToConsiderIndices.push_back(patch_index);
                }
            }
            assert(patchesToConsiderIndices.size() == nbPatchesToConsider);

            /////////////////
            //// Compute FSTs
            /////////////////
            // Initialize variables
            double averageOfRatios_GST  = 0.0;
            double averageOfRatios_WCST = 0.0;

            double ratioOfAverages_GST_NUMERATOR   = 0.0;
            double ratioOfAverages_GST_DENOMINATOR  = 0.0;
            double ratioOfAverages_WCST_NUMERATOR   = 0.0;
            double ratioOfAverages_WCST_DENOMINATOR  = 0.0;

            // Loop through loci
            for (uint32_t SNPindex = 0 ; SNPindex < polymorphicLoci.size() ; SNPindex++)
            {
                // Compute meanAlleleFreq
                double meanAlleleFreq = 0.0;
                for (int patch_index : patchesToConsiderIndices)
                {
                    meanAlleleFreq += alleleFreqs[patch_index][SNPindex];
                }
                meanAlleleFreq /= patchesToConsiderIndices.size();

                // Compute SS (Sum Of Square, the numerator of a variance calculation)
                double SS = 0.0;
                for (int patch_index : patchesToConsiderIndices)
                {
                    SS += pow(alleleFreqs[patch_index][SNPindex] - meanAlleleFreq,2);
                }
                

                // Handy variables
                int r = 2;    // Number of sampled patches
                double pbar = meanAlleleFreq;
                double GSTnum = SS / r;
                double GSTden = pbar - pow(pbar,2);
                double s_square = SS / (r-1);
                double WCSTnum = s_square;
                double WCSTden = GSTden + s_square / r;
            
                
            
                // GST
                ratioOfAverages_GST_NUMERATOR += GSTnum;
                ratioOfAverages_GST_DENOMINATOR += GSTden;
                averageOfRatios_GST += GSTnum / GSTden;

                // WCST
                ratioOfAverages_WCST_NUMERATOR += WCSTnum;
                ratioOfAverages_WCST_DENOMINATOR += WCSTden;
                averageOfRatios_WCST += WCSTnum / WCSTden;


                /*          
                std::cout << "d = " << d << "\n";            
                std::cout << "r = " << r << "\n";            
                std::cout << "s_square = " << s_square << "\n";            
                std::cout << "pbar = " << pbar << "\n";

                std::cout << "(SS / r - 1) = " << (SS / r - 1) << "\n";
                std::cout << "pbar - pow(pbar,2) = " << pbar - pow(pbar,2)<< "\n";
                std::cout << "s_square / r = " << s_square / r << "\n";            

                std::cout << "meanGST = " << meanGST << "\n";
                std::cout << "meanWCST = " << meanWCST << "\n";
                */
            }

            averageOfRatios_GST  /= polymorphicLoci.size();
            averageOfRatios_WCST /= polymorphicLoci.size();
            
            double ratioOfAverages_GST = ratioOfAverages_GST_NUMERATOR / ratioOfAverages_GST_DENOMINATOR;
            double ratioOfAverages_WCST = ratioOfAverages_WCST_NUMERATOR / ratioOfAverages_WCST_DENOMINATOR;

            s += s_tab + std::to_string(averageOfRatios_GST) + s_tab + std::to_string(averageOfRatios_WCST) + s_tab + std::to_string(ratioOfAverages_GST) + s_tab + std::to_string(ratioOfAverages_WCST);

        } while (std::prev_permutation(patchesToConsider.begin(), patchesToConsider.end()));
    }

    s += std::string("\n");
    file.write(s);
    file.close();

}


template<typename ntrlit_t, typename selit_t>
void OutputWriter::WriteOutputs_sampleSequences_process(OutputFile& file, std::vector<uint32_t>& T4data, ntrlit_t it_ntrlBegin, selit_t it_selBegin, ntrlit_t it_ntrlEnd, selit_t it_selEnd, SampleSequenceData& SSD, Haplotype& haplo)
{
    //std::cout << "Enters in WriteOutputs_sampleSequences_process\n";
    const std::vector<uint32_t> forRef;
    //const std::vector<uint32_t>& T4data = SSP->Gmap.T4_nbLoci ? SSP->T4Tree.getMutationsOfID(haplo.T4ID) : forRef;
    std::vector<uint32_t>::const_iterator T4it = SSP->Gmap.T4_nbLoci ? T4data.cbegin() : forRef.cbegin();
    std::vector<uint32_t>::const_iterator T4itEnd = SSP->Gmap.T4_nbLoci ? T4data.cend() : forRef.cbegin();

    ntrlit_t ntrlit = it_ntrlBegin;
    selit_t selit = it_selBegin;

    compressedString s("_");
    SSD.start();
    s.reserve(SSD.nbLoci() + 10);
    while (SSD.isMoreLoci())
    {
        auto locusD = SSP->Gmap.getLocusTypeAndItsIndex(SSD.nextLocus());
        //std::cout << "locusD.type = " << (size_t) locusD.type << "\n";
        //std::cout << "locusD.locusInType = " << locusD.locusInType << "\n";

        switch (locusD.type)
        {
            case 1:
                s.add(std::to_string(haplo.getT1_Allele(locusD.locusInType)));
                break;

            case 2:
                s.add(std::to_string(haplo.getT2_Allele(locusD.locusInType)));
                break;

            case 3:           
                s.add(std::to_string(haplo.getT3_Allele(locusD.locusInType)));
                break;

            case 4:
                for ( ; T4it < T4itEnd && *T4it < locusD.locusInType; ++T4it);
                s.add(T4it != T4itEnd ? std::to_string( *T4it == locusD.locusInType ) : "0");
                break;

            case 50:           
                for ( ; ntrlit != it_ntrlEnd && *ntrlit < locusD.locusInType; ++ntrlit);
                s.add(ntrlit != it_ntrlEnd ? std::to_string(*ntrlit == locusD.locusInType) : "0");
                break;

            case 51:           
                for ( ; selit != it_selEnd && *selit < locusD.locusInType; ++selit);
                s.add(selit != it_selEnd ? std::to_string(*selit == locusD.locusInType) : "0");
                break;

            default:
                std::cout << "Internal error in 'OutputWriter::WriteOutputs_sampleSequences(Pop& pop, OutputFile& file)'. Received locus of unknown type '" << (int) locusD.type << "'\n";
                abort();
        }
    }
    s.add("\n");
    file.write(s.toString());
}


void OutputWriter::WriteOutputs_sampleSequences(std::vector<std::vector<std::vector<uint32_t>>>& T4data, Pop& pop, OutputFile& file, size_t mutPlacingIndex)
{
    if (SSP->whenDidExtinctionOccur != -1)
    {
        file.openT4(mutPlacingIndex);
        std::string s("Species is extinct\n");
        file.write(s);
        file.close();
        return;
    }


    SSP->sampleSequenceDataContainer.start();

    /*if (SSP->Gmap.T4_nbLoci)
    {
        SSP->T4Tree.simplify(pop);
        SSP->T4Tree.placeMutations(pop);
    } */

    file.openT4(mutPlacingIndex);

    while (SSP->sampleSequenceDataContainer.areThereMoreSequences())
    {
        SampleSequenceData& SSD = SSP->sampleSequenceDataContainer.getNextSequence();

        std::string s = "P" + std::to_string(SSD.patch) +
            "I" + std::to_string(SSD.ind) +
            "H" + std::to_string(SSD.haplo) + ":";
        file.write(s);

        if (SSD.patch < GP->PatchNumber && SSD.ind < SSP->patchSize[SSD.patch])
        {
            Haplotype& haplo = pop.getPatch(SSD.patch).getInd(SSD.ind).getHaplo(SSD.haplo);

            if (SSP->Gmap.isT56ntrlCompress)
            {
                if (SSP->Gmap.isT56selCompress)
                {
                    WriteOutputs_sampleSequences_process(file,
                        T4data[SSD.patch][SSD.ind * 2 * SSD.haplo], 
                        haplo.T6ntrl_AllelesBegin(), 
                        haplo.T6sel_AllelesBegin(), 
                        haplo.T6ntrl_AllelesEnd(), 
                        haplo.T6sel_AllelesEnd(),
                        SSD,
                        haplo
                    );
                } else
                {
                    WriteOutputs_sampleSequences_process(file, 
                        T4data[SSD.patch][SSD.ind * 2 * SSD.haplo],
                        haplo.T6ntrl_AllelesBegin(), 
                        haplo.T5sel_AllelesBegin(), 
                        haplo.T6ntrl_AllelesEnd(), 
                        haplo.T5sel_AllelesEnd(),
                        SSD,
                        haplo
                    );
                }
            } else
            {
                if (SSP->Gmap.isT56selCompress)
                {
                    WriteOutputs_sampleSequences_process(file,
                        T4data[SSD.patch][SSD.ind * 2 * SSD.haplo], 
                        haplo.T5ntrl_AllelesBegin(), 
                        haplo.T6sel_AllelesBegin(), 
                        haplo.T5ntrl_AllelesEnd(), 
                        haplo.T6sel_AllelesEnd(),
                        SSD,
                        haplo
                    );
                } else
                {
                    WriteOutputs_sampleSequences_process(file, 
                        T4data[SSD.patch][SSD.ind * 2 * SSD.haplo],
                        haplo.T5ntrl_AllelesBegin(), 
                        haplo.T5sel_AllelesBegin(), 
                        haplo.T5ntrl_AllelesEnd(), 
                        haplo.T5sel_AllelesEnd(),
                        SSD,
                        haplo
                    );
                }
            }
        }
    }


    file.close();
}


void OutputWriter::WriteOutputs_T4_paintedHaploSegmentsDiversity_file_header(OutputFile& file)
{
    file.open();
    std::string s;
    if (SSP->SegDiversityFile_includeMainColor)
    {
        s = "paintedGeneration\tobservedGeneration\tleft\tright\tpatchSampled\tnbHaplotypesSampled\tnbColors\tcolorHighestFrequency\tdiversity\n";
    } else
    {
        s = "paintedGeneration\tobservedGeneration\tleft\tright\tpatchSampled\tnbHaplotypesSampled\tnbColors\tdiversity\n";
    }
        
    file.write(s);
    file.close();
}

void OutputWriter::WriteOutputs_T4_paintedHaploSegmentsDiversity_file(Pop& pop, OutputFile& file)
{
    if (SSP->whenDidExtinctionOccur != -1)
    {
        file.open();
        std::string s("Species is extinct\n");
        file.write(s);
        file.close();
        return;
    }

    SSP->T4Tree.assertIsFullySimplified();
    //SSP->T4Tree.simplify(pop);

    //std::cout << "file.T4_paintedHaplo_information.size() = " << file.T4_paintedHaplo_information.size() << "\n";
    for (auto& info : file.T4_paintedHaplo_information )
    {
        assert(info.observedGeneration > info.paintedGeneration);
        
        if (info.observedGeneration == GP->CurrentGeneration || (KillOnDemand::justAboutToKill && info.observedGeneration >= GP->CurrentGeneration ))
        {

            // Create T4IDs for current focal haplotypes
            std::vector<uint32_t> T4IDs;
            std::vector<uint32_t> T4IDs_patches;
            assert(info.patch_indices.size() == info.nbHaplotypesPerPatch.size());
            for (size_t i = 0 ; i < info.patch_indices.size() ; ++i)
            {
                const auto& patch_index  = info.patch_indices[i];
                const auto& nbHaplotypes = info.nbHaplotypesPerPatch[i];

                if (i != 0) assert(patch_index >= info.patch_indices[i-1]);

                if (GP->PatchNumber > patch_index) // This 'if' is for KillOnDemand::justAboutToKill
                {
                    auto& patch = pop.getPatch(patch_index);

                    for (size_t sampledHaploIndex = 0 ; sampledHaploIndex < nbHaplotypes ; ++sampledHaploIndex)
                    {
                        size_t ind_index = sampledHaploIndex / 2;
                        size_t haplo_index = sampledHaploIndex % 2;

                        if (ind_index < SSP->patchSize[patch_index])
                        {
                            T4IDs.push_back(patch.getInd(ind_index).getHaplo(haplo_index).T4ID);
                            T4IDs_patches.push_back(patch_index);
                        }
                    }
                }
            }


            // Compute the string (one string will be several lines of input)
            if (T4IDs.size())
            {
                SSP->T4Tree.writePaintedHaplotypesDiversity(T4IDs, T4IDs_patches, info.paintedGeneration, GP->CurrentGeneration, file);
            }
        }
    }
}


void OutputWriter::WriteOutputs_T4_paintedHaplo_file(Pop& pop, OutputFile& file)
{
    if (SSP->whenDidExtinctionOccur != -1)
    {
        file.open();
        std::string s("Species is extinct\n");
        file.write(s);
        file.close();
        return;
    }
    
    SSP->T4Tree.assertIsFullySimplified();
    bool gotAtLeastOneOutput = false; // for assertion only
    assert(file.T4_paintedHaplo_information.size());

    //std::cout << "file.T4_paintedHaplo_information.size() = " << file.T4_paintedHaplo_information.size() << "\n";
    for (auto& info : file.T4_paintedHaplo_information )
    {
        assert(info.observedGeneration > info.paintedGeneration);
        //std::cout << "info.observedGeneration = " << info.observedGeneration << "\n";
        //std::cout << "info.paintedGeneration = " << info.paintedGeneration << "\n";
        if (info.observedGeneration == GP->CurrentGeneration)
        {
            // Assertions
            gotAtLeastOneOutput = true; // for assertion only


            // Create T4IDs for current focal haplotypes
            std::vector<uint32_t> T4IDs;
            std::vector<uint32_t> T4IDs_patches;
            assert(info.patch_indices.size() == info.nbHaplotypesPerPatch.size());
            for (size_t i = 0 ; i < info.patch_indices.size() ; ++i)
            {
                const auto& patch_index  = info.patch_indices[i];
                const auto& nbHaplotypes = info.nbHaplotypesPerPatch[i];

                auto& patch = pop.getPatch(patch_index);

                for (size_t sampledHaploIndex = 0 ; sampledHaploIndex < nbHaplotypes ; ++sampledHaploIndex)
                {
                    size_t ind_index = sampledHaploIndex / 2;
                    size_t haplo_index = sampledHaploIndex % 2;

                    if (ind_index < SSP->patchSize[patch_index])
                    {
                        T4IDs.push_back(patch.getInd(ind_index).getHaplo(haplo_index).T4ID);
                        T4IDs_patches.push_back(patch_index);
                    }
                }
            }


            // Compute the string (one string will be several lines of input)
            SSP->T4Tree.writePaintedHaplotypes(T4IDs, T4IDs_patches, info.paintedGeneration, info.observedGeneration, file);
        }
    }
    assert(gotAtLeastOneOutput);    
}


void OutputWriter::WriteOutputs_fitnessStats_header(OutputFile& file)
{
    file.open();

    std::string s;
    
    s += std::string("Generation");
    for ( int patch_index = 0 ; patch_index < GP->maxEverPatchNumber ; ++patch_index )
    {
        s += std::string("\tP") + std::to_string(patch_index) + std::string("_meanFit");
        s += std::string("\tP") + std::to_string(patch_index) + std::string("_varFit");
    }
    s += std::string("\n");
    
    file.write(s);
    file.close();
}


void OutputWriter::WriteOutputs_fitnessStats(Pop& pop, OutputFile& file)
{
    file.open();
    file.write(std::to_string(GP->CurrentGeneration));

    std::string s;
    const std::string s_tab("\t");

    for ( int patch_index = 0 ; patch_index < GP->maxEverPatchNumber ; ++patch_index )
    {
        if (shouldNABePrinted(patch_index))
        {
            s += s_tab + "NA" + s_tab + "NA";
        } else
        {
            std::vector<double> fits(SSP->patchSize[patch_index]);
            double meanFit = 0.0;
            for (int ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ind_index++)
            {
                fits[ind_index] = pop.getPatch(patch_index).getInd(ind_index, patch_index).CalculateFitness(patch_index);
                meanFit += fits[ind_index];
            }
            meanFit /= SSP->patchSize[patch_index];
            double varFit = 0.0;
            for (int ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ind_index++)
            {
                varFit += pow(meanFit - fits[ind_index],2);
            }
            varFit /= SSP->patchSize[patch_index];

            s += s_tab + std::to_string(meanFit) + s_tab + std::to_string(varFit);
        }
    }
    s += std::string("\n");
    
    file.write(s);
    file.close();
}

void OutputWriter::WriteOutputs_patchSize_header(OutputFile& file)
{
    file.open();

    std::string s;
    
    s += std::string("Generation");
    for ( int patch_index = 0 ; patch_index < GP->maxEverPatchNumber ; ++patch_index )
    {
        s += std::string("\tP") + std::to_string(patch_index);
    }
    s += std::string("\n");
    
    file.write(s);
    file.close();
}

void OutputWriter::WriteOutputs_patchSize(OutputFile& file)
{
    file.open();
    file.write(std::to_string(GP->CurrentGeneration));

    std::string s;
    const std::string s_tab("\t");

    for ( int patch_index = 0 ; patch_index < GP->maxEverPatchNumber ; ++patch_index )
    {
        if (shouldNABePrinted(patch_index))
        {
            assert(patch_index >= SSP->patchSize.size());
            s += s_tab + "NA";   
        } else
        {
            assert(patch_index < SSP->patchSize.size());
            s += s_tab + std::to_string(SSP->patchSize[patch_index]);   
        }
    }
    s += std::string("\n");
    
    file.write(s);
    file.close();
}


void OutputWriter::WriteOutputs_fitnessSubsetLoci_header(OutputFile& file)
{
    file.open();
    assert(SSP->subsetT1LociForfitnessSubsetLoci_file.size() == SSP->subsetT2LociForfitnessSubsetLoci_file.size());
    assert(SSP->subsetT1LociForfitnessSubsetLoci_file.size() == SSP->subsetT3LociForfitnessSubsetLoci_file.size());
    assert(SSP->subsetT1LociForfitnessSubsetLoci_file.size() == SSP->subsetT1epistasisLociForfitnessSubsetLoci_file.size());

    std::string s;
    
    s += std::string("Generation");
    for ( int patch_index = 0 ; patch_index < GP->maxEverPatchNumber ; ++patch_index )
    {
        for (int ind_index = 0 ; ind_index < SSP->maxEverpatchCapacity[patch_index] ; ind_index++)
        {
            for (int lociSetIndex = 0 ; lociSetIndex < SSP->subsetT1LociForfitnessSubsetLoci_file.size() ; lociSetIndex++)
            {
                s += std::string("\tP") + std::to_string(patch_index) + std::string("_I") + std::to_string(ind_index) + "_LociSet" + std::to_string(lociSetIndex) + "_totalFit";
                s += std::string("\tP") + std::to_string(patch_index) + std::string("_I") + std::to_string(ind_index) + "_LociSet" + std::to_string(lociSetIndex) + "_T1Fit";
                s += std::string("\tP") + std::to_string(patch_index) + std::string("_I") + std::to_string(ind_index) + "_LociSet" + std::to_string(lociSetIndex) + "_T1EpistasisFit";
                s += std::string("\tP") + std::to_string(patch_index) + std::string("_I") + std::to_string(ind_index) + "_LociSet" + std::to_string(lociSetIndex) + "_T2Fit";
                s += std::string("\tP") + std::to_string(patch_index) + std::string("_I") + std::to_string(ind_index) + "_LociSet" + std::to_string(lociSetIndex) + "_T3Fit";
                s += std::string("\tP") + std::to_string(patch_index) + std::string("_I") + std::to_string(ind_index) + "_LociSet" + std::to_string(lociSetIndex) + "_T5Fit";
            }
        }
    }
    s += std::string("\n");
    
    file.write(s);
    file.close();
}

void OutputWriter::WriteOutputs_fitnessSubsetLoci(Pop& pop, OutputFile& file)
{
    file.open();
    file.write(std::to_string(GP->CurrentGeneration));

    std::string s;
    const std::string s_tab("\t");

    for ( int patch_index = 0 ; patch_index < GP->maxEverPatchNumber ; ++patch_index )
    {
        for (int ind_index = 0 ; ind_index < SSP->maxEverpatchCapacity[patch_index] ; ind_index++)
        {
            for (int lociSetIndex = 0 ; lociSetIndex < SSP->subsetT1LociForfitnessSubsetLoci_file.size() ; lociSetIndex++)
            {
                if (shouldNABePrinted(patch_index,ind_index))
                {
                    s += s_tab + "NA" + s_tab + "NA" + s_tab + "NA" + s_tab + "NA" + s_tab + "NA" + s_tab + "NA";
                } else
                {
                    auto fitnessComponents = pop.getPatch(patch_index).getInd(ind_index, patch_index).CalculateFitnessComponentsOnSubsetOfLoci(SSP->Habitats[patch_index], lociSetIndex);
                    assert(fitnessComponents.size() == 6);

                    std::string T1Fit =             std::to_string(fitnessComponents[0]);
                    std::string T2Fit =             std::to_string(fitnessComponents[1]);
                    std::string T1EpistasisFit =    std::to_string(fitnessComponents[2]);
                    std::string T3Fit =             std::to_string(fitnessComponents[3]);
                    std::string T56Fit =             std::to_string(fitnessComponents[4]);

                    std::string totalFit;
                    if (SSP->additiveEffectAmongLoci)
                    {
                        double totFit = 1 - fitnessComponents[0] - fitnessComponents[1] - fitnessComponents[2] - fitnessComponents[3] - fitnessComponents[4];
                        if (totFit < 0.0) totFit = 0.0;                        
                        totalFit = std::to_string(totFit);

                    } else
                    {
                        totalFit = std::to_string(fitnessComponents[0] * fitnessComponents[1] * fitnessComponents[2] * fitnessComponents[3] * fitnessComponents[4]);
                    }
                    

                    s += s_tab + totalFit + s_tab + T1Fit + s_tab + T1EpistasisFit + s_tab + T2Fit + s_tab + T3Fit + s_tab + T56Fit;
                }
            }
        }
    }
    s += std::string("\n");
    
    file.write(s);
    file.close();
}

void OutputWriter::WriteOutputs_fitness_header(OutputFile& file)
{
    file.open();

    std::string s;
    
    s += std::string("Generation");
    for ( int patch_index = 0 ; patch_index < GP->maxEverPatchNumber ; ++patch_index )
    {
        for (int ind_index = 0 ; ind_index < SSP->maxEverpatchCapacity[patch_index] ; ind_index++)
        {
            s += std::string("\tP") + std::to_string(patch_index) + std::string("_I") + std::to_string(ind_index) + "_totalFit";
            s += std::string("\tP") + std::to_string(patch_index) + std::string("_I") + std::to_string(ind_index) + "_T1Fit";
            s += std::string("\tP") + std::to_string(patch_index) + std::string("_I") + std::to_string(ind_index) + "_T1EpistasisFit";
            s += std::string("\tP") + std::to_string(patch_index) + std::string("_I") + std::to_string(ind_index) + "_T2Fit";
            s += std::string("\tP") + std::to_string(patch_index) + std::string("_I") + std::to_string(ind_index) + "_T3Fit";
            s += std::string("\tP") + std::to_string(patch_index) + std::string("_I") + std::to_string(ind_index) + "_T5Fit";
        }
    }
    s += std::string("\n");
    
    file.write(s);
    file.close();
}

void OutputWriter::WriteOutputs_fitness(Pop& pop, OutputFile& file)
{
    file.open();
    file.write(std::to_string(GP->CurrentGeneration));

    std::string s;
    const std::string s_tab("\t");

    for ( int patch_index = 0 ; patch_index < GP->maxEverPatchNumber ; ++patch_index )
    {
        for (int ind_index = 0 ; ind_index < SSP->maxEverpatchCapacity[patch_index] ; ind_index++)
        {
            if (shouldNABePrinted(patch_index,ind_index))
            {
                s += s_tab + "NA" + s_tab + "NA" + s_tab + "NA" + s_tab + "NA" + s_tab + "NA" + s_tab + "NA";
            } else
            {
                auto fitnessComponents = pop.getPatch(patch_index).getInd(ind_index, patch_index).CalculateFitnessComponents(SSP->Habitats[patch_index]);
                assert(fitnessComponents.size() == 6);

                std::string T1Fit =             std::to_string(fitnessComponents[0]);
                std::string T2Fit =             std::to_string(fitnessComponents[1]);
                std::string T1EpistasisFit =    std::to_string(fitnessComponents[2]);
                std::string T3Fit =             std::to_string(fitnessComponents[3]);
                std::string T5Fit =             std::to_string(fitnessComponents[4]);
                std::string totalFit;
                if (SSP->additiveEffectAmongLoci)
                {
                    double totFit = 1 - fitnessComponents[0] - fitnessComponents[1] - fitnessComponents[2] - fitnessComponents[3] - fitnessComponents[4];
                    if (totFit < 0.0) totFit = 0.0;                        
                    totalFit = std::to_string(totFit);
                } else
                {
                    totalFit = std::to_string(fitnessComponents[0] * fitnessComponents[1] * fitnessComponents[2] * fitnessComponents[3] * fitnessComponents[4]);
                }

                s += s_tab + totalFit + s_tab + T1Fit + s_tab + T1EpistasisFit + s_tab + T2Fit + s_tab + T3Fit + s_tab + T5Fit;
            }
        }
    }
    s += std::string("\n");
    
    file.write(s);
    file.close();
}

void OutputWriter::WriteOutputs_T1_AlleleFreq_header(OutputFile& file)
{
    file.open();

    std::string s;
    
    s += std::string("Generation");
    for ( int patch_index = 0 ; patch_index < GP->maxEverPatchNumber ; ++patch_index )
    {
        for (auto& T1_locus : file.getSubsetToConsider(SSP->speciesIndex))
        {
            s += std::string("\tP") + std::to_string(patch_index) + std::string("_L") + std::to_string(T1_locus.locus);
        }
    }
    s += std::string("\n");
    
    file.write(s);
    file.close();
}

void OutputWriter::WriteOutputs_T56_AlleleFreq_header(OutputFile& file)
{
    file.open();
    assert(SSP->Gmap.T5_nbLoci == SSP->Gmap.T5sel_nbLoci + SSP->Gmap.T5ntrl_nbLoci);
    assert(SSP->Gmap.T6_nbLoci == SSP->Gmap.T6sel_nbLoci + SSP->Gmap.T6ntrl_nbLoci);

    std::string s;
    
    s += std::string("Generation");
    for ( int patch_index = 0 ; patch_index < GP->maxEverPatchNumber ; ++patch_index )
    {
        for (uint32_t T5_locus = 0 ; T5_locus < SSP->Gmap.T56_nbLoci ; T5_locus++)
        {
            s += std::string("\tP") + std::to_string(patch_index) + std::string("_L") + std::to_string(T5_locus);
        }
    }
    s += std::string("\n");
    
    file.write(s);
    file.close();
}


void OutputWriter::WriteOutputs_T1_AlleleFreq(Pop& pop, OutputFile& file)
{
    file.open();

    file.write(std::to_string(GP->CurrentGeneration));

    const std::string s_tab("\t");
    
    for ( int patch_index = 0 ; patch_index < GP->maxEverPatchNumber ; ++patch_index )
    {
        for (auto& T1_locus : file.getSubsetToConsider(SSP->speciesIndex))
        {
            std::string s;
            if (shouldNABePrinted(patch_index))
            {
                s += s_tab + "NA";
            } else
            {
                s.reserve(SSP->patchSize[patch_index] * 2 + 40);
                double nbOnes = 0;
                for (int ind_index=0;ind_index<SSP->patchSize[patch_index];++ind_index)
                {
                    for (int haplo_index=0;haplo_index<SSP->ploidy;haplo_index++)
                    { 
                        nbOnes += pop.getPatch(patch_index).getInd(ind_index).getHaplo(haplo_index).getT1_Allele(T1_locus.byte_index,T1_locus.bit_index);
                    }
                }
                s += s_tab + std::to_string(nbOnes / (SSP->patchSize[patch_index] * SSP->ploidy));
            }
                
            file.write(s);
        }
    }
    file.write(std::string("\n"));
    file.close();
}


void OutputWriter::WriteOutputs_T56_AlleleFreq(Pop& pop, OutputFile& file)
{
    file.open();

    file.write(std::to_string(GP->CurrentGeneration));

    std::vector<std::vector<unsigned>> obsFreqs = pop.computePatchSpecificT56Frequencies();

    
    for ( int patch_index = 0 ; patch_index < GP->maxEverPatchNumber ; ++patch_index )
    {
        for (int locus = 0 ; locus < SSP->Gmap.T56_nbLoci ; locus++)
        {
            std::string s("\t");
            if (shouldNABePrinted(patch_index))
            {
                s += "NA";
            } else
            {           
                s += std::to_string((double) obsFreqs[patch_index][locus] / (2.0 * (double)SSP->TotalpatchSize));
            }
                
            file.write(s);
        }
    }
    file.write(std::string("\n"));
    file.close();
}

void OutputWriter::WriteOutputs_T1_MeanLD_header(OutputFile& file)
{
    if (SSP->Gmap.T1_nbLoci > 100)
    {
        std::cout <<  "Warning, you are asking for meanLD outputs with " << SSP->Gmap.T1_nbLoci << " loci. This represent a lot of computation! It might take a lot of time to compute (the computational time will also depends greatly on the total population size)." << std::endl;
    }

    int nbChroms = SSP->ChromosomeBoundaries.size() + 1;
    
    std::string s;
    s += std::string("Generation");
    for ( int patch_index = 0 ; patch_index < GP->maxEverPatchNumber ; ++patch_index )
    {
        for (int chrom_index = 0 ; chrom_index < nbChroms ; chrom_index++)
        {
            s += std::string("\tD_P")  + std::to_string(patch_index) + std::string("_WithinChrom") + std::to_string(chrom_index);
            s += std::string("\tDp_P") + std::to_string(patch_index) + std::string("_WithinChrom") + std::to_string(chrom_index);
            s += std::string("\tr2_P") + std::to_string(patch_index) + std::string("_WithinChrom") + std::to_string(chrom_index);
        }

        for (int chromA_index = 0 ; chromA_index < (nbChroms - 1) ; chromA_index++)
        {
            for (int chromB_index = chromA_index + 1 ; chromB_index < nbChroms ; chromB_index++)
            {
                assert(nbChroms > 1);
                s += std::string("\tD_P")  + std::to_string(patch_index) + std::string("_AmongChroms") + std::to_string(chromA_index) + std::to_string(chromB_index);
                s += std::string("\tDp_P") + std::to_string(patch_index) + std::string("_AmongChroms") + std::to_string(chromA_index) + std::to_string(chromB_index);
                s += std::string("\tr2_P") + std::to_string(patch_index) + std::string("_AmongChroms") + std::to_string(chromA_index) + std::to_string(chromB_index);
            }
        }
    }

    s += std::string("\n");
    file.open();
    file.write(s);
    file.close();
}

void OutputWriter::WriteOutputs_T1_MeanLD(Pop& pop, OutputFile& file)
{
    std::string s;
    const std::string s_tab = "\t";
    std::string s_NA = "NA";
    
    s += std::to_string(GP->CurrentGeneration);
    s.reserve(8 * GP->PatchNumber);
    
    std::vector<int> T1_AlleleChromBoundaries;
    T1_AlleleChromBoundaries.push_back(0);
    for (int b : SSP->ChromosomeBoundaries)
    {
        T1_AlleleChromBoundaries.push_back( SSP->Gmap.FromLocusToNextT1Locus(b) );
    }
    assert(T1_AlleleChromBoundaries.back() < SSP->Gmap.T1_nbLoci);
    T1_AlleleChromBoundaries.push_back(SSP->Gmap.T1_nbLoci);
    assert((T1_AlleleChromBoundaries.size()-2) == SSP->ChromosomeBoundaries.size()); // nbChroms = SSP->ChromosomeBoundaries.size()

    for ( int patch_index = 0 ; patch_index < GP->maxEverPatchNumber ; ++patch_index )
    {
        for ( int IndependentHaplotypeIndex = 0 ; IndependentHaplotypeIndex < T1_AlleleChromBoundaries.size()-1 ; IndependentHaplotypeIndex++)
        {
            if (shouldNABePrinted(patch_index))
            {
                s += s_tab + "NA";
                s += s_tab + "NA";
                s += s_tab + "NA";
            } else
            {
                int startSite   = T1_AlleleChromBoundaries[IndependentHaplotypeIndex];
                int endSite     = T1_AlleleChromBoundaries[IndependentHaplotypeIndex + 1];
                assert(startSite < endSite);
                // explore range [startSite, endSite)

                double D_mean  = 0.0;
                double Dp_mean = 0.0;
                double r2_mean = 0.0;
                double averageOverHowMany = 0.0; // Used double for future division but it may not be very efficient
                
                for ( int site_a = startSite ; site_a < endSite - 1 ; site_a++ )
                {
                    if (file.isLocusInSubset(site_a, SSP->speciesIndex))
                    {
                        for ( int site_b = site_a + 1 ; site_b < endSite ; site_b++ )
                        {
                            if (file.isLocusInSubset(site_b, SSP->speciesIndex))
                            {
                                int byte_a = ceil(site_a / EIGHT);
                                int bit_a  = site_a % EIGHT;
                                int byte_b = ceil(site_b / EIGHT);
                                int bit_b  = site_b % EIGHT;

                                int NbOneOne = 0;
                                int NbaOne   = 0;
                                int NbbOne   = 0;
                                for ( int ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ++ind_index )
                                {
                                    for ( int haplo_index = 0 ; haplo_index < SSP->ploidy ; haplo_index++ )
                                    {
                                        Haplotype& chrom = pop.getPatch(patch_index).getInd(ind_index).getHaplo(haplo_index);
                                        NbaOne   += chrom.getT1_Allele(byte_a, bit_a);
                                        NbbOne   += chrom.getT1_Allele(byte_b, bit_b);
                                        NbOneOne += chrom.getT1_Allele(byte_a, bit_a) && chrom.getT1_Allele(byte_b, bit_b);
                                    }
                                }

                                double totalNbAlleles = (double) SSP->patchSize[patch_index] * SSP->ploidy;
                                // px is the frequency of the allele x (called a and b here but they could be called 0 and 1 to match the rest of the naming system). qx = 1 - px
                                double pa = NbaOne / totalNbAlleles;
                                double pb = NbbOne / totalNbAlleles;
                                //assert(pa>=0.0 && pa <= 1.0);
                                //assert(pb>=0.0 && pb <= 1.0);
                                double qa = 1 - pa;
                                double qb = 1 - pb;
                                if (pa != 0.0 && qa != 0.0 && pb != 0.0 && qb != 0.0) // This is not optimal. I should only loop over polymorphic site in the region but...hey I'm lazy
                                {
                                    // Compute stats for this pair of loci
                                    averageOverHowMany++;
                                    double D = NbOneOne/totalNbAlleles - pa*pb;
                                    //assert(D >= -1.0 && D <= 1.0);
                                    double r2 = pow(D,2) / ( pa * qa * pb * qb );

                                    // Because of round-off error I must do
                                    if (r2 > 1.0)
                                    {
                                        //assert(r2 < 1.00001);
                                        {
                                            r2 = 1.0;
                                        }
                                    }
                                    if (r2 < -1.0)
                                    {
                                        //assert(r2 > -1.00001);
                                        {
                                            r2 = -1.0;
                                        }
                                    }

                                    //assert(r2 >= -1.0 && r2 <= 1.0);
                                    double Dp;
                                    if (D > 0.0)
                                    {
                                        double Dmax = std::min(pa*qb, qa*pb);
                                        Dp = D / Dmax;
                                    } else if (D < 0.0)
                                    {
                                        double Dmin = std::max(-pa*pb, -qa*qb);
                                        Dp = D / Dmin;
                                    } else
                                    {
                                        //assert(D==0.0);
                                        Dp = 0.0;
                                    }
                                    if (Dp < -1.0)
                                    {
                                        //assert(Dp > -1.1);
                                        Dp=-1.0;
                                    }
                                    if (Dp > 1.0)
                                    {
                                        //assert(Dp < 1.1);   
                                        Dp=1.0;
                                    }

                                    // Add statistics to the mean objects
                                    D_mean += D;
                                    Dp_mean += Dp;
                                    r2_mean += r2;
                                }
                            }
                        }
                    }
                }
                if (averageOverHowMany > 0)
                {
                    D_mean  /= averageOverHowMany;
                    Dp_mean /= averageOverHowMany;
                    r2_mean /= averageOverHowMany;

                    //assert(D_mean >= -1.0 && D_mean <= 1.0);
                    //assert(Dp_mean >= -1.0 && Dp_mean <= 1.0);
                    //assert(r2_mean >= -1.0 && r2_mean <= 1.0);
                } else
                {
                    D_mean  = std::numeric_limits<double>::quiet_NaN();
                    Dp_mean = std::numeric_limits<double>::quiet_NaN();
                    r2_mean = std::numeric_limits<double>::quiet_NaN();

                }            

                s += s_tab + std::to_string(D_mean);
                s += s_tab + std::to_string(Dp_mean);
                s += s_tab + std::to_string(r2_mean);
            }
        }
        

        // Between Haplotypes meanLD        
        for ( int IndependentHaplotypeIndex_A = 0 ; IndependentHaplotypeIndex_A < T1_AlleleChromBoundaries.size()-2 ; ++IndependentHaplotypeIndex_A)
        {

            int startSite_A   = T1_AlleleChromBoundaries[IndependentHaplotypeIndex_A];
            int endSite_A     = T1_AlleleChromBoundaries[IndependentHaplotypeIndex_A + 1];                    
            assert(startSite_A < endSite_A);
            // explore range [startSite_A, endSite_A)

            for ( int IndependentHaplotypeIndex_B = IndependentHaplotypeIndex_A + 1 ; IndependentHaplotypeIndex_B < T1_AlleleChromBoundaries.size()-1 ; ++IndependentHaplotypeIndex_B)
            {
                if (shouldNABePrinted(patch_index))
                {
                    s += s_tab + "NA";
                    s += s_tab + "NA";
                    s += s_tab + "NA";
                } else
                {
                    int startSite_B   = T1_AlleleChromBoundaries[IndependentHaplotypeIndex_B];
                    int endSite_B     = T1_AlleleChromBoundaries[IndependentHaplotypeIndex_B + 1];                    
                    assert(startSite_B < endSite_B);
                    // explore range [startSite_B, endSite_B)

                    double D_mean  = 0.0;
                    double Dp_mean = 0.0;
                    double r2_mean = 0.0;
                    double averageOverHowMany = 0.0; // Used double for future division but it may not be very efficient

                    for (int site_a = startSite_A ; site_a < endSite_A ; ++site_a)
                    {
                        if (file.isLocusInSubset(site_a, SSP->speciesIndex))
                        {
                            for (int site_b = startSite_B ; site_b < endSite_B ; ++site_b)
                            {
                                if (file.isLocusInSubset(site_b, SSP->speciesIndex))
                                {
                                    int byte_a = ceil(site_a / EIGHT);
                                    int bit_a  = site_a % EIGHT;
                                    int byte_b = ceil(site_b / EIGHT);
                                    int bit_b  = site_b % EIGHT;

                                    int NbOneOne = 0;
                                    int NbaOne   = 0;
                                    int NbbOne   = 0;
                                    for ( int ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ++ind_index )
                                    {
                                        for ( int haplo_index = 0 ; haplo_index < SSP->ploidy ; haplo_index++ )
                                        {
                                            Haplotype& chrom = pop.getPatch(patch_index).getInd(ind_index).getHaplo(haplo_index);
                                            NbaOne   += chrom.getT1_Allele(byte_a, bit_a);
                                            NbbOne   += chrom.getT1_Allele(byte_b, bit_b);
                                            NbOneOne += chrom.getT1_Allele(byte_a, bit_a) && chrom.getT1_Allele(byte_b, bit_b);
                                        }
                                    }

                                    // Compute stats for this pair of loci
                                    double totalNbAlleles = (double) SSP->patchSize[patch_index] * SSP->ploidy;
                                    double pa = NbaOne / totalNbAlleles;
                                    double pb = NbbOne / totalNbAlleles;
                                    //assert(pa>=0.0 && pa <= 1.0);
                                    //assert(pb>=0.0 && pb <= 1.0);
                                    double qa = 1 - pa;
                                    double qb = 1 - pb;

                                    if (pa != 0.0 && qa != 0.0 && pb != 0.0 && qb != 0.0)
                                    {
                                        averageOverHowMany++;
                                        double totalNbAlleles = (double) SSP->patchSize[patch_index] * SSP->ploidy;
                                        // px is the frequency of the allele x (called a and b here but they could be called 0 and 1 to match the rest of the naming system). qx = 1 - px
                                        
                                        double D = NbOneOne/totalNbAlleles - pa*pb;
                                        //assert(D<=1 && D>=-1);
                                        double r2 = pow(D,2) / ( pa * qa * pb * qb );

                                        // Because of round-off error I must do
                                        if (r2 > 1.0)
                                        {
                                            //assert(r2 < 1.00001);
                                            {
                                                r2 = 1.0;
                                            }
                                        }
                                        if (r2 < -1.0)
                                        {
                                            //assert(r2 > -1.00001);
                                            {
                                                r2 = -1.0;
                                            }
                                        }

                                        double Dp;
                                        if (D > 0.0)
                                        {
                                            double Dmax = std::min(-pa*pb,-qa*qb);
                                            Dp = D / Dmax;
                                        } else if (D < 0.0)
                                        {
                                            double Dmin = std::max(pa*qb,qa*pb);
                                            Dp = D / Dmin;
                                        } else
                                        {
                                            //assert(D==0.0);
                                            Dp = 0.0;
                                        }
                                        //assert(Dp >= -1.0 && Dp <= 1.0);


                                        // Add statistics to the mean objects
                                        D_mean += D;
                                        Dp_mean += Dp;
                                        r2_mean += r2;
                                    }
                                }
                            }
                        }
                    }
                    if (averageOverHowMany > 0)
                    {
                        D_mean  /= averageOverHowMany;
                        Dp_mean /= averageOverHowMany;
                        r2_mean /= averageOverHowMany;

                        //assert(D_mean >= -1.0 && D_mean <= 1.0);
                        //assert(Dp_mean >= -1.0 && Dp_mean <= 1.0);
                        //assert(r2_mean >= -1.0 && r2_mean <= 1.0);
                    } else
                    {
                        D_mean  = std::numeric_limits<double>::quiet_NaN();
                        Dp_mean = std::numeric_limits<double>::quiet_NaN();
                        r2_mean = std::numeric_limits<double>::quiet_NaN();

                    }            

                    s += s_tab + std::to_string(D_mean);
                    s += s_tab + std::to_string(Dp_mean);
                    s += s_tab + std::to_string(r2_mean);   
                }
            }
        }
    }
    s +=std::string("\n");

    file.open();
    file.write(s);
    file.close();
}

void OutputWriter::WriteOutputs_T1_LongestRun_header(OutputFile& file)
{
    file.open();
    
    std::string s;
    s += std::string("Generation");
    for ( int patch_index = 0 ; patch_index < GP->maxEverPatchNumber ; ++patch_index )
    {
        for (int ind_index = 0 ; ind_index < SSP->maxEverpatchCapacity[patch_index] ; ind_index++)
        {
            for (int haplo_index = 0 ; haplo_index < SSP->ploidy ; haplo_index++ )
            {
                s + std::string("\tP") + std::to_string(patch_index) + std::string("_I") + std::to_string(ind_index) + std::string("_H") + std::to_string(haplo_index) + std::string("_A0");
                
                s + std::string("\tP") + std::to_string(patch_index) + std::string("_I") + std::to_string(ind_index) + std::string("_H") + std::to_string(haplo_index) + std::string("_A1");
            }
        }
    }
    file.write(s);
    file.close();
}

void OutputWriter::WriteOutputs_T1_LongestRun(Pop& pop, OutputFile& file)
{
    file.open();

    const std::string s_tab = "\t";
    
    file.write(std::to_string(GP->CurrentGeneration));
    
    int RunZeros;
    int RunOnes;
    int LongestRunZeros;
    int LongestRunOnes;
    int PreviousT1_Allele;
    double CurrentT1_Allele;
    
    for ( int patch_index = 0 ; patch_index < GP->maxEverPatchNumber ; ++patch_index )
    {
        std::string s;
        if (SSP->patchSize.size() > patch_index)
        {
            s.reserve(10 * SSP->patchSize[patch_index]);
        }
            
        for (int ind_index = 0 ; ind_index < SSP->maxEverpatchCapacity[patch_index] ; ind_index++)
        {
            for (int haplo_index=0;haplo_index<SSP->ploidy;haplo_index++)
            {
                if (shouldNABePrinted(patch_index,ind_index))
                {
                    s += s_tab + "NA" + s_tab + "NA";
                } else
                {
                    LongestRunZeros = 0;
                    LongestRunOnes  = 0;
                    RunZeros        = 0;
                    RunOnes         = 0;
                    PreviousT1_Allele  =-1;
                    for (auto& T1_locus : file.getSubsetToConsider(SSP->speciesIndex))
                    {
                        CurrentT1_Allele = pop.getPatch(patch_index).getInd(ind_index).getHaplo(haplo_index).getT1_Allele(T1_locus.byte_index,T1_locus.bit_index);
                        if (CurrentT1_Allele)
                        {
                            RunZeros=0;
                            RunOnes++;
                        } else
                        {
                            RunOnes=0;
                            RunZeros++;
                        }
                        LongestRunOnes  = std::max(LongestRunOnes,RunOnes);
                        LongestRunZeros = std::max(LongestRunZeros,RunZeros);
                    }
                    assert(LongestRunZeros>=0);
                    assert(LongestRunOnes>=0);
                    assert((LongestRunZeros+LongestRunOnes)>0);
                    s += s_tab + std::to_string(LongestRunZeros) + s_tab + std::to_string(LongestRunOnes);
                }
            }
        }
        file.write(s);
    }

    file.write(std::string("\n"));
    file.close();
}


void OutputWriter::WriteOutputs_T2_LargeOutput_header(OutputFile& file)
{
    file.open();

    const std::string s_tabP = "\tP";
    std::string s__I = "_I";
    std::string s__H = "_H";
    std::string s__L = "_L";
    
    file.write(std::string("Generation"));
    for ( int patch_index = 0 ; patch_index < GP->maxEverPatchNumber ; ++patch_index )
    {
        std::string s;
        s.reserve(16 * SSP->patchCapacity[patch_index] * SSP->Gmap.T2_nbLoci);
        for (int ind_index = 0 ; ind_index < SSP->maxEverpatchCapacity[patch_index] ; ind_index++)
        {
            for (int haplo_index = 0 ; haplo_index < SSP->ploidy ; ++haplo_index)
            {
                for (int T2_char_index = 0 ; T2_char_index < SSP->Gmap.T2_nbLoci ; ++T2_char_index)
                {
                    s += s_tabP + std::to_string(patch_index) + s__I + std::to_string(ind_index) + s__H + std::to_string(haplo_index) + s__L + std::to_string(T2_char_index);
                }
            }
        }
        file.write(s);
    }
    file.write(std::string("\n"));
    file.close();
}

void OutputWriter::WriteOutputs_T2_LargeOutput(Pop& pop, OutputFile& file)
{
    file.open();
    
    const std::string s_tab = "\t";
    std::string s;
    s.reserve(SSP->TotalpatchCapacity * SSP->Gmap.T1_nbLoci * 2.5);
    
    s += std::to_string( GP->CurrentGeneration);
    for ( int patch_index = 0 ; patch_index < GP->maxEverPatchNumber ; ++patch_index )
    {
        std::string s;
        for (int ind_index = 0 ; ind_index < SSP->maxEverpatchCapacity[patch_index] ; ind_index++)
        {
            for (int haplo_index = 0 ; haplo_index < SSP->ploidy ; haplo_index++)
            {
                for (int T2_char_index = 0 ; T2_char_index < SSP->Gmap.T2_nbLoci ; T2_char_index++)
                {
                    if (shouldNABePrinted(patch_index,ind_index))
                    {
                        s += s_tab + "NA";
                    } else
                    {
                        s += s_tab + std::to_string(pop.getPatch(patch_index).getInd(ind_index).getHaplo(haplo_index).getT2_Allele(T2_char_index));
                    }
                }
            }
        }
    }
    s +=std::string("\n");
    
    file.write(s);
    file.close();
}

void OutputWriter::WriteOutputs_T56_LargeOutput_header(OutputFile& file)
{
    file.open();
    
    std::string s("Generation");
    for ( int patch_index = 0 ; patch_index < GP->maxEverPatchNumber ; ++patch_index )
    {
        for (int ind_index = 0 ; ind_index < SSP->maxEverpatchCapacity[patch_index] ; ind_index++)
        {
            for (int haplo_index = 0 ; haplo_index < SSP->ploidy ; haplo_index++)
            {
                for (int locus = 0 ; locus < SSP->Gmap.T56_nbLoci ; locus++)
                {
                    s += std::string("\tP") + std::to_string(patch_index) + std::string("_I") + std::to_string(ind_index) + std::string("_H") + std::to_string(haplo_index) + std::string("_L") + std::to_string(locus);
                }
            }
        }
    }
    s += std::string("\n");
    file.write(s);
    file.close();
}

template<typename ntrlIterator, typename selIterator>
void OutputWriter::WriteOutputs_T56_LargeOutput_writeData(ntrlIterator ntrlIt, selIterator selIt, ntrlIterator ntrlItEnd, selIterator selItEnd, OutputFile& file, bool printNA)
{
    std::string s;
    const std::string s_tab = "\t";
    for (int T5locus = 0 ; T5locus < SSP->Gmap.T56_nbLoci ; T5locus++) // not the most clever solution
    {
        if (printNA)
        {
            s += s_tab + "NA";
        } else
        {
            auto T5genderLocus = SSP->Gmap.FromT56LocusToT56genderLocus(T5locus);
            if (T5genderLocus.isNtrl)
            {
                auto T5ntrlLocus = T5genderLocus.locusInGender;
                assert(SSP->Gmap.isT56neutral(T5ntrlLocus));

                if (ntrlIt != ntrlItEnd && *ntrlIt == T5ntrlLocus)
                {
                    ++ntrlIt;
                    s += s_tab + "1";
                } else
                {
                    s += s_tab + "0";
                }
            } else
            {
                auto T5selLocus = T5genderLocus.locusInGender;
                assert(!SSP->Gmap.isT56neutral(T5selLocus));

                if (selIt != selItEnd && *selIt == T5selLocus)
                {
                    ++selIt;
                    s += s_tab + "1";
                } else
                {
                    s += s_tab + "0";
                }
            }
        }                    
    }
    assert(ntrlIt == ntrlItEnd);
    assert(selIt == selItEnd);
    file.write(s);
}


void OutputWriter::WriteOutputs_T56_LargeOutput(Pop& pop, OutputFile& file)
{
    std::string s;
    s.reserve(SSP->TotalpatchCapacity * SSP->Gmap.T1_nbLoci * 3);
    
    file.open();
    file.write(std::to_string( GP->CurrentGeneration));
    for ( int patch_index = 0 ; patch_index < GP->maxEverPatchNumber ; ++patch_index )
    {
        for (int ind_index = 0 ; ind_index < SSP->maxEverpatchCapacity[patch_index] ; ind_index++)
        {
            for (int haplo_index=0;haplo_index<SSP->ploidy;haplo_index++)
            {
                Haplotype& haplo = pop.getPatch(patch_index).getInd(ind_index).getHaplo(haplo_index);

                //std::cout << "SSP->Gmap.isT56ntrlCompress = " << SSP->Gmap.isT56ntrlCompress << "\n";
                //std::cout << "SSP->Gmap.isT56selCompress = " << SSP->Gmap.isT56selCompress << "\n";


                if (SSP->Gmap.isT56ntrlCompress)
                {
                    // T6ntrl
                    if (SSP->Gmap.isT56selCompress)
                    {
                        // T6sel
                        ZipIterator<CompressedSortedDeque, CompressedSortedDeque::iterator> ntrlIt;
                        ZipIterator<CompressedSortedDeque, CompressedSortedDeque::iterator> ntrlItEnd;
                        CompressedSortedDeque::iterator selIt;
                        CompressedSortedDeque::iterator selItEnd;

                        if (SSP->Gmap.T6ntrl_nbLoci)
                        {
                            ntrlIt = haplo.T6ntrl_AllelesBegin();;
                            ntrlItEnd = haplo.T6ntrl_AllelesEnd();;
                        }

                        if (SSP->Gmap.T6sel_nbLoci)
                        {
                            selIt = haplo.T6sel_AllelesBegin();;
                            selItEnd = haplo.T6sel_AllelesEnd();;
                        }

                        WriteOutputs_T56_LargeOutput_writeData(ntrlIt, selIt, ntrlItEnd, selItEnd, file, shouldNABePrinted(patch_index,ind_index));
                    } else
                    {
                        // T5sel
                        ZipIterator<CompressedSortedDeque, CompressedSortedDeque::iterator> ntrlIt;
                        ZipIterator<CompressedSortedDeque, CompressedSortedDeque::iterator> ntrlItEnd;
                        std::vector<uint32_t>::iterator selIt;
                        std::vector<uint32_t>::iterator selItEnd;

                        if (SSP->Gmap.T6ntrl_nbLoci)
                        {
                            ntrlIt = haplo.T6ntrl_AllelesBegin();;
                            ntrlItEnd = haplo.T6ntrl_AllelesEnd();;
                        }

                        if (SSP->Gmap.T5sel_nbLoci)
                        {
                            selIt = haplo.T5sel_AllelesBegin();;
                            selItEnd = haplo.T5sel_AllelesEnd();;
                        }

                        WriteOutputs_T56_LargeOutput_writeData(ntrlIt, selIt, ntrlItEnd, selItEnd, file, shouldNABePrinted(patch_index,ind_index));
                    }
                } else
                {
                    // T5ntrl
                    if (SSP->Gmap.isT56selCompress)
                    {
                        //T6sel
                        ZipIterator<std::vector<uint32_t>, std::vector<uint32_t>::iterator> ntrlIt;
                        ZipIterator<std::vector<uint32_t>, std::vector<uint32_t>::iterator> ntrlItEnd;
                        CompressedSortedDeque::iterator selIt;
                        CompressedSortedDeque::iterator selItEnd;

                        if (SSP->Gmap.T5ntrl_nbLoci)
                        {
                            ntrlIt = haplo.T5ntrl_AllelesBegin();;
                            ntrlItEnd = haplo.T5ntrl_AllelesEnd();;
                        }

                        if (SSP->Gmap.T6sel_nbLoci)
                        {
                            selIt = haplo.T6sel_AllelesBegin();;
                            selItEnd = haplo.T6sel_AllelesEnd();;
                        }
                        WriteOutputs_T56_LargeOutput_writeData(ntrlIt, selIt, ntrlItEnd, selItEnd, file, shouldNABePrinted(patch_index,ind_index));
                    } else
                    {
                        //T5sel
                        ZipIterator<std::vector<uint32_t>, std::vector<uint32_t>::iterator> ntrlIt;
                        ZipIterator<std::vector<uint32_t>, std::vector<uint32_t>::iterator> ntrlItEnd;
                        std::vector<uint32_t>::iterator selIt;
                        std::vector<uint32_t>::iterator selItEnd;

                        if (SSP->Gmap.T5ntrl_nbLoci)
                        {
                            ntrlIt = haplo.T5ntrl_AllelesBegin();;
                            ntrlItEnd = haplo.T5ntrl_AllelesEnd();;
                        }

                        if (SSP->Gmap.T5sel_nbLoci)
                        {
                            selIt = haplo.T5sel_AllelesBegin();;
                            selItEnd = haplo.T5sel_AllelesEnd();;
                        }

                        
                        WriteOutputs_T56_LargeOutput_writeData(ntrlIt, selIt, ntrlItEnd, selItEnd, file, shouldNABePrinted(patch_index,ind_index));
                    }
                }
            }
        }
    }
    file.write(std::string("\n"));
    file.close();
}



void OutputWriter::WriteOutputs_T1_LargeOutput_header(OutputFile& file)
{
    file.open();
    
    std::string s("Generation");
    for ( int patch_index = 0 ; patch_index < GP->maxEverPatchNumber ; ++patch_index )
    {
        for (int ind_index = 0 ; ind_index < SSP->maxEverpatchCapacity[patch_index] ; ind_index++)
        {
            for (int haplo_index = 0 ; haplo_index < SSP->ploidy ; haplo_index++)
            {
                for (auto& T1_locus : file.getSubsetToConsider(SSP->speciesIndex))
                {
                    s += std::string("\tP") + std::to_string(patch_index) + std::string("_I") + std::to_string(ind_index) + std::string("_H") + std::to_string(haplo_index) + std::string("_L") + std::to_string(T1_locus.locus);
                }
            }
        }
    }
    s += std::string("\n");
    file.write(s);
    file.close();
}


void OutputWriter::WriteOutputs_T1_LargeOutput(Pop& pop, OutputFile& file)
{
    std::string s;
    const std::string s_tab = "\t";
    s.reserve(SSP->TotalpatchCapacity * SSP->Gmap.T1_nbLoci * 3);
    
    s += std::to_string( GP->CurrentGeneration);
    for ( int patch_index = 0 ; patch_index < GP->maxEverPatchNumber ; ++patch_index )
    {
        for (int ind_index = 0 ; ind_index < SSP->maxEverpatchCapacity[patch_index] ; ind_index++)
        {
            for (int haplo_index=0;haplo_index<SSP->ploidy;haplo_index++)
            {
                for (auto& T1_locus : file.getSubsetToConsider(SSP->speciesIndex))
                {
                    if (shouldNABePrinted(patch_index,ind_index))
                    {
                        s += s_tab + "NA";
                    } else
                    {
                        s += s_tab + std::to_string(pop.getPatch(patch_index).getInd(ind_index).getHaplo(haplo_index).getT1_Allele(T1_locus.byte_index,T1_locus.bit_index));
                    }
                }
            }
        }
    }
    s +=std::string("\n");
    file.open();
    file.write(s);
    file.close();
}


void OutputWriter::WriteOutputs_T4_LargeOutput_header(OutputFile& file, size_t mutPlacingIndex)
{
    file.openT4(mutPlacingIndex);
    
    std::string s("Generation");
    for ( int patch_index = 0 ; patch_index < GP->maxEverPatchNumber ; ++patch_index )
    {
        for (int ind_index = 0 ; ind_index < SSP->maxEverpatchCapacity[patch_index] ; ind_index++)
        {
            for (int haplo_index = 0 ; haplo_index < SSP->ploidy ; haplo_index++)
            {
                for (int T4_locus = 0 ; T4_locus < SSP->Gmap.T4_nbLoci ; T4_locus++)
                {
                    s += std::string("\tP") + std::to_string(patch_index) + std::string("_I") + std::to_string(ind_index) + std::string("_H") + std::to_string(haplo_index) + std::string("_L") + std::to_string(T4_locus);
                }
            }
        }
    }
    s += std::string("\n");
    file.write(s);
    file.close();
}

void OutputWriter::WriteOutputs_T4_LargeOutput(OutputFile& file, std::vector<std::vector<std::vector<uint32_t>>>& data, size_t mutPlacingIndex)
{
    std::string s;
    const std::string s_tab = "\t";
    s.reserve(SSP->TotalpatchCapacity * SSP->Gmap.T4_nbLoci * 3);

    assert(data.size() == GP->PatchNumber);
    
    s += std::to_string( GP->CurrentGeneration);
    for ( int patch_index = 0 ; patch_index < GP->maxEverPatchNumber ; ++patch_index )
    {
        assert(data[patch_index].size() == SSP->patchSize[patch_index]);
        for (int ind_index = 0 ; ind_index < SSP->maxEverpatchCapacity[patch_index] ; ind_index++)
        {
            for (int haplo_index=0;haplo_index<SSP->ploidy;haplo_index++)
            {   
                if (shouldNABePrinted(patch_index,ind_index))
                {
                    for (int T4_locus = 0 ; T4_locus < SSP->Gmap.T4_nbLoci ; T4_locus++)   s += s_tab + "NA";                  
                } else
                {
                    auto it = data[patch_index][ind_index*2+haplo_index].begin();
                    for (int T4_locus = 0 ; T4_locus < SSP->Gmap.T4_nbLoci ; T4_locus++)
                    {
                        if (*it == T4_locus)
                        {
                            s += s_tab + "1";
                            ++it;
                        } else
                        {
                            s += s_tab + "0";
                        }
                    }
                    assert(it == data[patch_index][ind_index*2+haplo_index].end());
                }
            }
        }
    }
    s +=std::string("\n");
    file.openT4(mutPlacingIndex);
    file.write(s);
    file.close();
}




void OutputWriter::WriteOutputs_T1_HybridIndex_header(OutputFile& file)
{   
    std::string s("Generation");
    for ( int patch_index = 0 ; patch_index < GP->maxEverPatchNumber ; ++patch_index )
    {
        for (int ind_index = 0 ; ind_index < SSP->maxEverpatchCapacity[patch_index] ; ind_index++)
        {
            s += std::string("\tP") + std::to_string( patch_index) + std::string("_I" ) + std::to_string(ind_index);
        }
    }

    file.open();
    s += std::string("\n");
    file.write(s);
    file.close();
}


void OutputWriter::WriteOutputs_T1_AverageHybridIndex_header(OutputFile& file)
{
    std::string s("Generation");
    for ( int patch_index = 0 ; patch_index < GP->maxEverPatchNumber ; ++patch_index )
    {
        s += std::string("\tP") + std::to_string( patch_index);
    }

    file.open();
    s += std::string("\n");
    file.write(s);
    file.close();   
}

void OutputWriter::WriteOutputs_T1_HybridIndex(Pop& pop, OutputFile& file)
{
    std::string s;
    const std::string s_tab = "\t";
    
    s += std::to_string( GP->CurrentGeneration);
    s.reserve(SSP->TotalpatchCapacity * 3);
    
    
    for ( int patch_index = 0 ; patch_index < GP->maxEverPatchNumber ; ++patch_index )
    {
        for (int ind_index = 0 ; ind_index < SSP->maxEverpatchCapacity[patch_index] ; ind_index++)
        {
            if (shouldNABePrinted(patch_index,ind_index))
            {
                s += s_tab + "NA";
            } else
            {
                double nbOnes = 0.0;
                for (int haplo_index=0;haplo_index<SSP->ploidy;haplo_index++)
                {
                    for (auto& T1_locus : file.getSubsetToConsider(SSP->speciesIndex))
                    {
                        nbOnes+=pop.getPatch(patch_index).getInd(ind_index).getHaplo(haplo_index).getT1_Allele(T1_locus.byte_index,T1_locus.bit_index);
                    }   
                }
                double r =  nbOnes / (double)(file.getSubsetToConsider(SSP->speciesIndex).size() * SSP->ploidy);
                assert(r>=0 && r<=1);
                s += s_tab + std::to_string(r);
            }   
        }
    }
    s +=std::string("\n");
    file.open();
    file.write(s);
    file.close();
}


void OutputWriter::WriteOutputs_T1_AverageHybridIndex(Pop& pop, OutputFile& file)
{
    std::string s;
    const std::string s_tab = "\t";
    
    s += std::to_string( GP->CurrentGeneration);
    s.reserve(SSP->TotalpatchCapacity * 3);
    
    
    for ( int patch_index = 0 ; patch_index < GP->maxEverPatchNumber ; ++patch_index )
    {
        double meanr = 0.0;
        uint32_t divideBy = 0;
        for (int ind_index = 0 ; ind_index < SSP->maxEverpatchCapacity[patch_index] ; ind_index++)
        {
            if (shouldNABePrinted(patch_index,ind_index))
            {
                // nothing to do
            } else
            {
                double nbOnes = 0.0;
                for (int haplo_index=0;haplo_index<SSP->ploidy;haplo_index++)
                {
                    for (auto& T1_locus : file.getSubsetToConsider(SSP->speciesIndex))
                    {
                        nbOnes+=pop.getPatch(patch_index).getInd(ind_index).getHaplo(haplo_index).getT1_Allele(T1_locus.byte_index,T1_locus.bit_index);
                    }   
                }
                double r =  nbOnes / (double)(file.getSubsetToConsider(SSP->speciesIndex).size() * SSP->ploidy);
                assert(r>=0 && r<=1);
                meanr += r;
                ++divideBy;
            }   
        }
        meanr /= divideBy;
        assert(meanr>=0 && meanr<=1);
        s += s_tab + std::to_string(meanr);
    }
    s += std::string("\n");
    file.open();
    file.write(s);
    file.close();
}

void OutputWriter::WriteOutputs_T1_ExpectiMinRec_header(OutputFile& file)
{
    std::string s("Generation");
    for ( int patch_index = 0 ; patch_index < GP->maxEverPatchNumber ; ++patch_index )
    {
        s += std::string("\tP") + std::to_string(patch_index);
    }

    s += std::string("\n");
    file.open();
    file.write(s);
    file.close();
}

void OutputWriter::WriteOutputs_T1_ExpectiMinRec(Pop& pop, OutputFile& file)
{
    std::string s;
    const std::string s_tab = "\t";
    s.reserve(GP->PatchNumber * 5);
    
    s += std::to_string(GP->CurrentGeneration);
    
    for ( int patch_index = 0 ; patch_index < GP->maxEverPatchNumber ; ++patch_index )
    {
        double MinimNbRecs = 0.0;
        for (int ind_index=0;ind_index<SSP->patchSize[patch_index];++ind_index)
        {
            for (int haplo_index=0;haplo_index<SSP->ploidy;haplo_index++)
            {
                bool PreviousT1_Allele = pop.getPatch(patch_index).getInd(ind_index).getHaplo(haplo_index).getT1_Allele(0,0);
                for (auto& T1_locus : file.getSubsetToConsider(SSP->speciesIndex))
                {
                    bool CurrentT1_Allele = pop.getPatch(patch_index).getInd(ind_index).getHaplo(haplo_index).getT1_Allele(T1_locus.byte_index,T1_locus.bit_index);
                    if (CurrentT1_Allele != PreviousT1_Allele)
                    {
                        MinimNbRecs++;
                        PreviousT1_Allele = CurrentT1_Allele;
                    }
                }
                assert(MinimNbRecs>=0);
            }
        }
        if (shouldNABePrinted(patch_index))
        {
            s += s_tab +  "NA";
        } else
        {
            s += s_tab +  std::to_string(MinimNbRecs / (double)(SSP->patchSize[patch_index] * SSP->ploidy));
        }
    }
    s +=std::string("\n");
    file.open();
    file.write(s);
    file.close();
}


void OutputWriter::WriteOutputs_T1_vcf(Pop& pop, OutputFile& file)
{
    //std::cout << "In OutputWriter::WriteOutputs_T1_vcf at generation " <<GP->CurrentGeneration<<"\n";

    // Find Polymorphic Sites
    auto polymorphicLoci = file.removeSitesWeDontWant(pop.listT1PolymorphicLoci(), SSP->speciesIndex);


    const std::string s_tab = "\t";
    const std::string s_P   = "P";
    const std::string s_dotI = ".I";

    file.open();

    // Write header
    {
        std::string s;

        s += std::string("##fileformat=VCFv4.2\n##Generation=") + std::to_string( GP->CurrentGeneration) + std::string("\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
        s.reserve(s.size() + SSP->TotalpatchCapacity * 14);
        
        for (int patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
        {
            for (int ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ++ind_index)
            {
                s += s_tab + s_P + std::to_string(patch_index) + s_dotI + std::to_string(ind_index);
            }
        }
        s += std::string("\n");
        //std::cout << "writing: " << s << "\n";
        file.write(s);

    }
           
    if (polymorphicLoci.size() > SSP->Gmap.T1_nbLoci)
    {
        std::cout << "Internal bug: OutputWriter line 1929 polymorphicLoci.size() = "<<polymorphicLoci.size()<<"  SSP->Gmap.T1_nbLoci = "<<SSP->Gmap.T1_nbLoci<<"\n";
        for (int i = 0 ; i < polymorphicLoci.size(); i++)
        {
            std::cout << i << ": " << polymorphicLoci[i].locus << " ";
        }
        std::cout << "\n";
        assert(polymorphicLoci.size() <= SSP->Gmap.T1_nbLoci);
        abort();
    }

    
    // Write Data
    std::string s_vert = "|";
    for (auto& TrackedMutation : polymorphicLoci)
    {

        int byte_index  = TrackedMutation.byte_index;
        int bit_index   = TrackedMutation.bit_index;
        int locus       = TrackedMutation.locus;
        if (locus != byte_index * EIGHT + bit_index)
        {
            /*
            std::cout << "locus = " << locus << "\n";
            std::cout << "byte_index = " << byte_index << "\n";
            std::cout << "bit_index = " << bit_index << "\n";
            */
            assert(locus == byte_index * EIGHT + bit_index);
        }
            

        //std::cout << "Print locus " << locus << "(byte_index" << byte_index << ", bit_index " << bit_index << ").\n";
        
        std::vector<int>::iterator it = std::lower_bound(SSP->ChromosomeBoundaries.begin(), SSP->ChromosomeBoundaries.end(), locus);
        
        int WhichIndependentSegment = it - SSP->ChromosomeBoundaries.begin();
        
        std::string s;
        s += std::to_string(WhichIndependentSegment + 1) + std::string("\t") + std::to_string(locus) + std::string("\t.\tA\tT\t.\t.\t.\tGT"); // arbitrarily decide that 0 is A and 1 is T. I chose 'A' for ancestral as 0 is often the non-mutated allele although depends upon starting condtions asked by the user.
        s.reserve(SSP->TotalpatchCapacity * 4 + 10);
        
        for (int patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
        {
            for (int ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ++ind_index)
            {
                s += s_tab
                +
                std::to_string(pop.getPatch(patch_index).getInd(ind_index).getHaplo(0).getT1_Allele(byte_index,bit_index))
                +
                s_vert
                +
                std::to_string(pop.getPatch(patch_index).getInd(ind_index).getHaplo(1).getT1_Allele(byte_index,bit_index));
            }
        }
        s += "\n";
        file.write(s);
    }
    
    file.close();
}

void OutputWriter::WriteOutputs_Tx_SNPfreq_file_header(OutputFile& file, size_t mutPlacingIndex)
{
    file.openT4(mutPlacingIndex);
    std::string s("generation\tpatch\tlocus\tfrequency\n");
    file.write(s);
    file.close();
}

/*
void OutputWriter::mergeT4countsWithOtherCounts(T4TreeRec::PatchCountData& T4freqs, std::vector<uint32_t>& counts)
{
    if (T4freqs.isVec)
    {
        assert(T4freqs.vecFreqs.size() == SSP->Gmap.T4_nbLoci);
        for (size_t T4locus = 0 ; T4locus < SSP->Gmap.T4_nbLoci ; ++T4locus)
        {
            if (T4freqs.vecFreqs[T4locus] != 0)
                counts[SSP->Gmap.FromT4LocusToLocus(T4locus)] = T4freqs.vecFreqs[T4locus];
        }
    } else
    {
        for (auto it = T4freqs.mapFreqs.begin() ; it != T4freqs.mapFreqs.end() ; ++it)
        {
            assert(it->second > 0);
            auto T4locus = it->first;
            counts[SSP->Gmap.FromT4LocusToLocus(T4locus)] = it->second;
        }
    }
}

void OutputWriter::mergeT4countsWithOtherCounts(T4TreeRec::PatchCountData& T4freqs, std::map<uint32_t,uint32_t>& counts)
{
    if (T4freqs.isVec)
    {
        assert(T4freqs.vecFreqs.size() == SSP->Gmap.T4_nbLoci);
        for (size_t T4locus = 0 ; T4locus < SSP->Gmap.T4_nbLoci ; ++T4locus)
        {
            if (T4freqs.vecFreqs[T4locus] != 0)
                counts[SSP->Gmap.FromT4LocusToLocus(T4locus)] = T4freqs.vecFreqs[T4locus];
        }
    } else
    {
        counts.insert(T4freqs.mapFreqs.begin(), T4freqs.mapFreqs.end());
    }
}*/

template<typename CountContainerType> void OutputWriter::WriteOutputs_SNPfreq_computeFreq(CountContainerType& counts, Pop& pop, size_t patch_index)
{
    //double maxCount = 2 * SSP->patchSize[patch_index];
    for (uint32_t ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ++ind_index)
    {
        for (uint32_t haplo_index = 0 ; haplo_index < 2 ; ++haplo_index)
        {
            auto& haplo = pop.getPatch(patch_index).getInd(ind_index).getHaplo(haplo_index);

            // T1
                for (size_t T1_locus = 0 ; T1_locus < SSP->Gmap.T1_nbLoci ; ++T1_locus)
                    if (haplo.getT1_Allele(T1_locus))
                        ++counts[SSP->Gmap.FromT1LocusToLocus(T1_locus)];

            // T2
                for (size_t T2_locus = 0 ; T2_locus < SSP->Gmap.T2_nbLoci ; ++T2_locus)
                    if (haplo.getT2_Allele(T2_locus)) // If any value great than 0
                        ++counts[SSP->Gmap.FromT2LocusToLocus(T2_locus)];

            // T4 Nothing to do                

            // T56ntrl
            {
                if (SSP->Gmap.T56ntrl_nbLoci)
                {
                    if (SSP->Gmap.isT56ntrlCompress)
                    {
                        auto itEnd = haplo.T6ntrl_AllelesEnd();
                        for (auto it = haplo.T6ntrl_AllelesBegin() ; it != itEnd ; ++it)
                            ++counts[SSP->Gmap.FromT56ntrlLocusToLocus(*it)];
                    } else
                    {
                        auto itEnd = haplo.T5ntrl_AllelesEnd();
                        for (auto it = haplo.T5ntrl_AllelesBegin() ; it != itEnd ; ++it)
                            ++counts[SSP->Gmap.FromT56ntrlLocusToLocus(*it)];
                    }
                }
            }

            // T56sel
            {
                if (SSP->Gmap.T56sel_nbLoci)
                {
                    if (SSP->Gmap.isT56ntrlCompress)
                    {
                        auto itEnd = haplo.T6sel_AllelesEnd();
                        for (auto it = haplo.T6sel_AllelesBegin() ; it != itEnd ; ++it)
                            ++counts[SSP->Gmap.FromT56selLocusToLocus(*it)];
                    } else
                    {
                        auto itEnd = haplo.T5sel_AllelesEnd();
                        for (auto it = haplo.T5sel_AllelesBegin() ; it != itEnd ; ++it)
                            ++counts[SSP->Gmap.FromT56selLocusToLocus(*it)];
                    }
                }
            }
        }
    }
}

void OutputWriter::WriteOutputs_Tx_SNPfreq_file(OutputFile& file, std::vector<T4TreeRec::PatchCountData>& T4freqs, size_t mutPlacingIndex, Pop& pop)
{
    bool printInfo = false;

    if (printInfo) std::cout << "Enters in 'WriteOutputs_Tx_SNPfreq_file'\n";
    if (SSP->whenDidExtinctionOccur != -1)
    {
        assert(SSP->fecundityForFitnessOfOne != -1);
        file.openT4(mutPlacingIndex);
        std::string s("Species is extinct\n");
        file.write(s);
        file.close();
        return;
    }
    //SSP->T4Tree.simplify(pop);
    //SSP->T4Tree.placeMutations(pop);
    assert(T4freqs.size() == GP->PatchNumber);
    file.openT4(mutPlacingIndex);

    for (uint32_t patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
    {
        if (SSP->patchSize[patch_index] == 0) continue;
        double maxCount = 2 * SSP->patchSize[patch_index];



        if (SSP->TxSNPfreq_assumeMuchDiversity)
        {
            // Tx is map
            std::map<uint32_t, uint32_t> counts;
            WriteOutputs_SNPfreq_computeFreq(counts, pop, patch_index);


            if (T4freqs[patch_index].isVec)
            {
                // T4 is vector
                assert(T4freqs[patch_index].vecFreqs.size() == SSP->Gmap.T4_nbLoci);

                size_t T4locus = 0;
                
                for (auto it = counts.begin(); it != counts.end(); ++it)
                {
                    auto locusOfT4 = SSP->Gmap.FromT4LocusToLocus(T4locus);
                    auto locusOfTx = it->first;
                    assert(locusOfTx != locusOfT4);

                    for (;T4locus < SSP->Gmap.T4_nbLoci && locusOfTx < locusOfT4; ++T4locus)
                    {
                        if (T4freqs[patch_index].vecFreqs[T4locus] != 0)
                        {
                            auto count = T4freqs[patch_index].vecFreqs[T4locus];
                            double freq = (double)count / maxCount;
                            auto locus = SSP->Gmap.FromT4LocusToLocus(T4locus);
                            std::string s = std::to_string(GP->CurrentGeneration) + "\t" + std::to_string(patch_index) + "\t" + std::to_string(locus) + "\t" + std::to_string(freq) + "\n";
                            file.write(s);
                        }
                    }

                    assert(it->second != 0);
                    auto locus = it->first;
                    double freq = (double)it->second / maxCount ;
                    std::string s = std::to_string(GP->CurrentGeneration) + "\t" + std::to_string(patch_index) + "\t" + std::to_string(locus) + "\t" + std::to_string(freq) + "\n";
                    file.write(s);
                }

                for (;T4locus < SSP->Gmap.T4_nbLoci; ++T4locus)
                {
                    if (T4freqs[patch_index].vecFreqs[T4locus] != 0)
                    {
                        auto count = T4freqs[patch_index].vecFreqs[T4locus];
                        double freq = (double)count / maxCount;
                        auto locus = SSP->Gmap.FromT4LocusToLocus(T4locus);
                        std::string s = std::to_string(GP->CurrentGeneration) + "\t" + std::to_string(patch_index) + "\t" + std::to_string(locus) + "\t" + std::to_string(freq) + "\n";
                        file.write(s);
                    }
                }

            } else
            {
                // T4 is map
                // Merge the two maps
                counts.insert(T4freqs[patch_index].mapFreqs.begin(), T4freqs[patch_index].mapFreqs.end());

                // write counts
                for (auto it = counts.begin() ; it != counts.end() ; ++it)
                {
                    assert(it->second != 0);
                    auto locus = it->first;
                    double freq = (double)it->second / maxCount ;
                    std::string s = std::to_string(GP->CurrentGeneration) + "\t" + std::to_string(patch_index) + "\t" + std::to_string(locus) + "\t" + std::to_string(freq) + "\n";
                    file.write(s);
                }
            }


        } else
        {
            // Tx is vector
            std::vector<uint32_t> counts(SSP->Gmap.TotalNbLoci, 0);
            WriteOutputs_SNPfreq_computeFreq(counts, pop, patch_index);

            
            // Just merge T4 (whatever format) into Tx
            if (T4freqs[patch_index].isVec)
            {
                // T4 is vector
                assert(T4freqs[patch_index].vecFreqs.size() == SSP->Gmap.T4_nbLoci);
                for (size_t T4locus = 0 ; T4locus < SSP->Gmap.T4_nbLoci ; ++T4locus)
                {
                    if (T4freqs[patch_index].vecFreqs[T4locus] != 0)
                    {
                        auto count = T4freqs[patch_index].vecFreqs[T4locus];
                        auto locus = SSP->Gmap.FromT4LocusToLocus(T4locus);
                        assert(counts[locus] == 0);
                        counts[locus] = count;
                    }
                }
            } else
            {
                // T4 is map
                for (auto it = T4freqs[patch_index].mapFreqs.begin() ; it != T4freqs[patch_index].mapFreqs.end() ; ++it)
                {
                    assert(it->second > 0);
                    auto count = it->second;
                    auto locus = SSP->Gmap.FromT4LocusToLocus(it->first);
                    assert(counts[locus] == 0);
                    counts[locus] = count;
                }
            }

            // write
            for (size_t locus = 0 ; locus < SSP->Gmap.TotalNbLoci ; ++locus)
            {
                if (counts[locus] != 0)
                {
                    double freq = (double)counts[locus] / maxCount ;
                    std::string s = std::to_string(GP->CurrentGeneration) + "\t" + std::to_string(patch_index) + "\t" + std::to_string(locus) + "\t" + std::to_string(freq) + "\n";
                    file.write(s);
                }
            }
        }            
    }
    
    file.close();
}


void OutputWriter::WriteOutputs_T4_SNPfreq_file_header(OutputFile& file, size_t mutPlacingIndex)
{
    file.openT4(mutPlacingIndex);
    std::string s("generation\tpatch\tlocus\tfrequency\n");
    file.write(s);
    file.close();
}


void OutputWriter::WriteOutputs_T4_SNPfreq_file(OutputFile& file, std::vector<T4TreeRec::PatchCountData>& T4freqs, size_t mutPlacingIndex)
{
    if (SSP->whenDidExtinctionOccur != -1)
    {
        assert(SSP->fecundityForFitnessOfOne != -1);
        file.openT4(mutPlacingIndex);
        std::string s("Species is extinct\n");
        file.write(s);
        file.close();
        return;
    }
    //SSP->T4Tree.simplify(pop);
    //SSP->T4Tree.placeMutations(pop);

    file.openT4(mutPlacingIndex);

    assert(T4freqs.size() == GP->PatchNumber);
    for (uint32_t patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
    {
        // Compute counts

        /*
        //std::map<uint32_t, uint32_t> counts;
        std::vector<uint32_t> counts(SSP->Gmap.T4_nbLoci,0);

        for (uint32_t ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ++ind_index)
        {
            for (uint32_t haplo_index = 0 ; haplo_index < 2 ; ++haplo_index)
            {
                for (auto& mut : data[patch_index][ind_index * 2 * haplo_index]) // SSP->T4Tree.getMutationsOfID( pop.getPatch(patch_index).getInd(ind_index).getHaplo(haplo_index).T4ID ) )
                {
                    ++counts[mut];
                }
            }
        }
        */

        
        // Write
        double maxCount = 2 * SSP->patchSize[patch_index];

        if (T4freqs[patch_index].isVec)
        {
            assert(T4freqs[patch_index].vecFreqs.size() == SSP->Gmap.T4_nbLoci);
            for (size_t locus = 0 ; locus < SSP->Gmap.T4_nbLoci ; ++locus)
            {
                if (T4freqs[patch_index].vecFreqs[locus] != 0)
                {
                    double freq = (double) T4freqs[patch_index].vecFreqs[locus]  / maxCount;
                    std::string s = std::to_string(GP->CurrentGeneration) + "\t" + std::to_string(patch_index) + "\t" + std::to_string(locus) + "\t" + std::to_string(freq) + "\n";
                    file.write(s);
                }
            }
        } else
        {
            for (auto it = T4freqs[patch_index].mapFreqs.begin() ; it != T4freqs[patch_index].mapFreqs.end() ; ++it)
            {
                assert(it->second > 0);
                double freq = (double) it->second  / maxCount;
                auto locus = it->first;
                std::string s = std::to_string(GP->CurrentGeneration) + "\t" + std::to_string(patch_index) + "\t" + std::to_string(locus) + "\t" + std::to_string(freq) + "\n";
                file.write(s);
            }
        }
    }
    
    file.close();
}



void OutputWriter::WriteOutputs_T4_vcf(OutputFile& file, std::vector<std::vector<std::vector<uint32_t>>>& data, size_t mutPlacingIndex)
{
    //std::cout << "In OutputWriter::WriteOutputs_T1_vcf at generation " <<GP->CurrentGeneration<<"\n";

    // Get current states
    if (SSP->whenDidExtinctionOccur != -1)
    {
        file.openT4(mutPlacingIndex);
        std::string s("Species is extinct\n");
        file.write(s);
        file.close();
        return;
    }
    
   /* std::cout << "data:\n";
    for (uint32_t patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
    {
        std::cout << "P" << patch_index << "\n";
        for (uint32_t ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ++ind_index)
        {
            std::cout << "\tI" << ind_index << "\n\t\tH1: ";
            for (auto& mut : data[patch_index][ind_index*2])
                std::cout << mut << " ";
            std::cout << "\n\t\tH2: ";
            for (auto& mut : data[patch_index][ind_index*2+1])
                std::cout << mut << " ";
            std::cout << "\n";
        }
    }*/
    
    // Find polymorphic loci
    std::vector<double> alleleFreqs(SSP->Gmap.T4_nbLoci, 0.0);
    assert(data.size() == GP->PatchNumber);
    for (uint32_t patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
    {
        assert(data[patch_index].size() == 2*SSP->patchSize[patch_index]);
        for (uint32_t ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ++ind_index)
        {
            for (uint32_t haplo_index = 0 ; haplo_index < 2 ; ++haplo_index)
            {
                for (auto& mut : data[patch_index][ind_index*2+haplo_index])
                {
                    ++alleleFreqs[mut];
                }
            }
        }
    }
    std::vector<uint32_t> polymorphicLoci;
    for (uint32_t locus = 0 ; locus < SSP->Gmap.T4_nbLoci ; ++locus)
    {
        alleleFreqs[locus] /= SSP->TotalpatchSize * 2;
        //std::cout << "alleleFreqs["<<locus<<"] = " << alleleFreqs[locus] << "\n";
        assert(alleleFreqs[locus] >= 0.0 && alleleFreqs[locus] <= 1.0);
        if (alleleFreqs[locus] != 0.0 ) // && alleleFreqs[locus] != 1.0)
        {
            //std::cout << "alleleFreqs["<<locus<<"] = " << alleleFreqs[locus] << "\n";
            polymorphicLoci.push_back(locus);
        }
    }


    const std::string s_tab = "\t";
    const std::string s_P   = "P";
    const std::string s_dotI = ".I";

    file.openT4(mutPlacingIndex);

    std::vector<std::vector<std::vector<uint32_t>::const_iterator>> iterators(GP->PatchNumber); // iterators[patch_index][ind_haplo_index]
    // Write header and set iterators
    {
        std::string s;

        s += std::string("##fileformat=VCFv4.2\n##Generation=") + std::to_string( GP->CurrentGeneration) + std::string("\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
        s.reserve(s.size() + SSP->TotalpatchCapacity * 4);
        
        for (int patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
        {
            iterators[patch_index].reserve(SSP->patchSize[patch_index] * 2);
            for (int ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ++ind_index)
            {
                s += s_tab + s_P + std::to_string(patch_index) + s_dotI + std::to_string(ind_index);
                iterators[patch_index].push_back(data[patch_index][ind_index*2].cbegin());
                iterators[patch_index].push_back(data[patch_index][ind_index*2+1].cbegin());
            }
        }
        s += std::string("\n");
        //std::cout << "writing: " << s << "\n";
        file.write(s);

    }
      
    if (polymorphicLoci.size() > SSP->Gmap.T4_nbLoci)
    {
        std::cout << "Internal bug: OutputWriter polymorphicLoci.size() = "<<polymorphicLoci.size()<<"  SSP->Gmap.T4_nbLoci = "<<SSP->Gmap.T4_nbLoci<<"\n";
    }
 
    // Write Data
    std::string s_vert = "|";
    for (auto& locus : polymorphicLoci)
    {
        std::string s;
        s.reserve(SSP->TotalpatchCapacity * 6);
        
        // independent segments (chromosome breaks)
        {
            auto it = std::lower_bound(SSP->ChromosomeBoundaries.begin(), SSP->ChromosomeBoundaries.end(), locus);
            int WhichIndependentSegment = it - SSP->ChromosomeBoundaries.begin();
            s += std::to_string(WhichIndependentSegment + 1) + std::string("\t") + std::to_string(locus + 1) + std::string("\t.\tA\tT\t.\t.\t.\tGT"); // arbitrarily decide that 0 is A and 1 is T. I chose 'A' for ancestral as 0 is often the non-mutated allele although depends upon starting condtions asked by the user.
        }
        
        for (int patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
        {
            for (int ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ++ind_index)
            {
                bool a0 = false;
                bool a1 = false;

                auto& it0 = iterators[patch_index][ind_index*2];
                auto& it1 = iterators[patch_index][ind_index*2+1];


                while (it0 != data[patch_index][ind_index*2].end() && *it0 < locus)
                {
                    ++it0;
                }
                if (it0 != data[patch_index][ind_index*2].end() && *it0 == locus)
                {
                    a0 = true;
                }

                while (it1 != data[patch_index][ind_index*2+1].end() && *it1 < locus)
                {
                    ++it1;
                }
                if (it1 != data[patch_index][ind_index*2+1].end() && *it1 == locus)
                {
                    a1 = true;
                }
                

                s += s_tab
                +
                std::to_string(a0)
                +
                s_vert
                +
                std::to_string(a1);
            }
        }
        s += "\n";
        file.write(s);
    }
    
    file.close();
}

/*template<typename ntrlIteratorType, typename selIteratorType>
void OutputWriter::write_T56vcf_forGroupOfLoci(std::vector<uint32_t>& groupOfLoci, std::vector<double>& relFreqsForGroupOfLoci, Pop& pop, ntrlIteratorType& ntrlIteratorTypeInfo, selIteratorType& selIteratorTypeInfo, OutputFile& file)
{
    const std::string s_tab("\t");
    const std::string s_vert("|");
    std::vector<std::string> ss(groupOfLoci.size());

    // Compute in which chromosome each locus is found and Initialize strings
    for (uint32_t locusIndex = 0 ; locusIndex < groupOfLoci.size() ; ++locusIndex)
    {
        uint32_t locus = groupOfLoci[locusIndex];
        std::vector<int>::iterator it = std::lower_bound(SSP->ChromosomeBoundaries.begin(), SSP->ChromosomeBoundaries.end(), locus);
        int WhichIndependentSegment = it - SSP->ChromosomeBoundaries.begin();

        ss[locusIndex] = std::to_string(WhichIndependentSegment + 1) + std::string("\t") + std::to_string(locus + 1) + std::string("\t.\tA\tT\t.\t.\t.\tGT");
        ss[locusIndex].reserve(ss[locusIndex].size() + SSP->TotalpatchCapacity * 4 + 20);
    }


    // Loop
    for (int patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
    {
        for (int ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ++ind_index)
        {
            Haplotype& haplo0 =  pop.getPatch(patch_index).getInd(ind_index).getHaplo(0);
            Haplotype& haplo1 =  pop.getPatch(patch_index).getInd(ind_index).getHaplo(1);

            for (uint32_t locusIndex = 0 ; locusIndex < groupOfLoci.size() ; ++locusIndex)
            {
                uint32_t locus = groupOfLoci[locusIndex];
                auto T5locusGender = SSP->FromT56LocusToT56genderLocus[locus];

                bool a0 = false;
                bool a1 = false;

                if (T5locusGender.first) // ntrl
                {
                    // Take iterators
                    auto T5ntrlIt0 = haplo0.T56ntrl_AllelesIterator(ntrlIteratorTypeInfo, locus);
                    auto T5ntrlIt1 = haplo1.T56ntrl_AllelesIterator(ntrlIteratorTypeInfo, locus);
                    auto T5ntrlIt0End = haplo0.T56ntrl_AllelesEnd(ntrlIteratorTypeInfo);
                    auto T5ntrlIt1End = haplo1.T56ntrl_AllelesEnd(ntrlIteratorTypeInfo);
                    

                    // Get allelic states
                    //std::cout << "line 2648\n";
                    while (T5ntrlIt0 != T5ntrlIt0End && *T5ntrlIt0 < T5locusGender.second)
                    {
                        ++T5ntrlIt0;
                    }
                    if (T5ntrlIt0 != T5ntrlIt0End && *T5ntrlIt0 == T5locusGender.second)
                    {
                        a0 = true;
                    }
                    //std::cout << "line 2694\n";
                    while (T5ntrlIt1 != T5ntrlIt1End && *T5ntrlIt1 < T5locusGender.second)
                    {
                        ++T5ntrlIt1;
                    }
                    if (T5ntrlIt1 != T5ntrlIt1End && *T5ntrlIt1 == T5locusGender.second)
                    {
                        a1 = true;
                    }
                    //std::cout << "line 2704\n";

                } else // sel
                {
                    auto T5selIt0 = haplo0.T56sel_AllelesIterator(selIteratorTypeInfo, locus);
                    auto T5selIt1 = haplo1.T56sel_AllelesIterator(selIteratorTypeInfo, locus);
                    auto T5selIt0End = haplo0.T56sel_AllelesEnd(selIteratorTypeInfo);
                    auto T5selIt1End = haplo1.T56sel_AllelesEnd(selIteratorTypeInfo);

                    // Get allelic states
                    //std::cout << "line 2725\n";
                    while (T5selIt0 != T5selIt0End && *T5selIt0 < T5locusGender.second)
                    {
                        ++T5selIt0;
                    }
                    if (T5selIt0 != T5selIt0End && *T5selIt0 == T5locusGender.second)
                    {
                        a0 = true;
                    }
                    //std::cout << "line 2735\n";

                    while (T5selIt1 != T5selIt1End && *T5selIt1 < T5locusGender.second)
                    {
                        ++T5selIt1;
                    }
                    if (T5selIt1 != T5selIt1End && *T5selIt1 == T5locusGender.second)
                    {
                        a1 = true;
                    }
                    //std::cout << "line 2746\n";
                }



                // Expand string
                ss[locusIndex] += s_tab
                +
                std::to_string(a0)
                +
                s_vert
                +
                std::to_string(a1);    
            }      
        }
    }

    // Write on file
    for (uint32_t locusIndex = 0 ; locusIndex < groupOfLoci.size() ; ++locusIndex)
    {
        ss[locusIndex] += "\n";
        file.write(ss[locusIndex]);
    }
}*/


/*
void OutputWriter::WriteOutputs_T56_vcf_writeData(Pop& pop, std::vector<unsigned>& obsFreqs, OutputFile& file)
{       
    // Estimate Memory Used
    double memoryUsed = 0.0; // per ind, in "quadriBytes"
    for (uint32_t locus = 0 ; locus < SSP->Gmap.T56_nbLoci ; ++locus)
    {
        double relFreq = (double) obsFreqs[locus] / (2.0 * SSP->TotalpatchSize);
        memoryUsed += relFreq; // Assuming not compressed
    }
    const double fractionOfRAMincrease = 0.001;
    double memoryForBlockOfLoci = 0.0; // per ind, in "quadriBytes"



    // Make group of loci and write them
    std::vector<uint32_t> groupOfLoci;
    std::vector<double> relFreqsForGroupOfLoci;
    //uint32_t nbBlocks= 0;
    for (int locus = 0 ; locus < SSP->Gmap.T56_nbLoci ; ++locus)
    {
        //std::cout << "locus = "<< locus << "\n";
        double relFreq = (double) obsFreqs[locus] / (2.0 * SSP->TotalpatchSize);
        assert(relFreq >= 0.0 && relFreq <= 1.0);
        //std::cout << "obsFreqs["<<locus<<"] = " << obsFreqs[locus] << "\n";


        if (relFreq > 0.0 && relFreq < 1.0)
        {
            groupOfLoci.push_back(locus);
            relFreqsForGroupOfLoci.push_back(relFreq);
            memoryForBlockOfLoci += relFreq; // Neglecting the memory for the start of the line
        }

        if (locus == SSP->Gmap.T56_nbLoci - 1 || memoryForBlockOfLoci > fractionOfRAMincrease * memoryUsed)
        {
            
            //++nbBlocks;
            memoryForBlockOfLoci = 0;

            int patch_index = 0;
            for ( ; patch_index < GP->PatchNumber ; ++patch_index)
            {
                if (SSP->patchSize[patch_index])
                {
                    break;
                }
            }

            if (patch_index != GP->PatchNumber) // if the population is not completely empty
            {
                auto& haplo = pop.getPatch(patch_index).getInd(0).getHaplo(0); 
                if (SSP->Gmap.isT56ntrlCompress)
                {
                    // T6ntrl
                    if (SSP->Gmap.isT56selCompress)
                    {
                        // T6sel
                        auto ntrlIt     = haplo.T6ntrl_AllelesBegin();
                        auto selIt      = haplo.T6sel_AllelesBegin();
                        //auto ntrlItEnd  = haplo.T6ntrl_AllelesEnd();
                        //auto selItEnd   = haplo.T6sel_AllelesEnd();
                        write_T56vcf_forGroupOfLoci(groupOfLoci, relFreqsForGroupOfLoci, pop, ntrlIt, selIt, file);
                    } else
                    {
                        // T5sel
                        auto ntrlIt     = haplo.T6ntrl_AllelesBegin();
                        auto selIt      = haplo.T5sel_AllelesBegin();
                        //auto ntrlItEnd  = haplo.T6ntrl_AllelesEnd();
                        //auto selItEnd   = haplo.T5sel_AllelesEnd();
                        write_T56vcf_forGroupOfLoci(groupOfLoci, relFreqsForGroupOfLoci, pop, ntrlIt, selIt, file);
                    }
                } else
                {
                    // T5ntrl
                    if (SSP->Gmap.isT56selCompress)
                    {
                        //T6sel
                        auto ntrlIt     = haplo.T5ntrl_AllelesBegin();
                        auto selIt      = haplo.T6sel_AllelesBegin();
                        //auto ntrlItEnd  = haplo.T5ntrl_AllelesEnd();
                        //auto selItEnd   = haplo.T6sel_AllelesEnd();
                        write_T56vcf_forGroupOfLoci(groupOfLoci, relFreqsForGroupOfLoci, pop, ntrlIt, selIt, file);
                    } else
                    {
                        //T5sel
                        auto ntrlIt             = haplo.T5ntrl_AllelesBegin();
                        auto selIt              = haplo.T5sel_AllelesBegin();
                        //ZipIterator<std::vector<uint32_t>, std::vector<uint32_t>::iterator> ntrlItEnd          = haplo.T5ntrl_AllelesEnd();
                        //ZipIterator<std::vector<uint32_t>, std::vector<uint32_t>::iterator> selItEnd           = haplo.T5sel_AllelesEnd();
                        write_T56vcf_forGroupOfLoci(groupOfLoci, relFreqsForGroupOfLoci, pop, ntrlIt, selIt, file);
                    }
                }


                groupOfLoci.clear();
                relFreqsForGroupOfLoci.clear();
            }
        }
    }
    //std::cout << "nbBlocks = " << nbBlocks<<"\n";
}*/


/*void OutputWriter::WriteOutputs_T56_vcf(Pop& pop, OutputFile& file)
{
    //std::cout << "In OutputWriter::WriteOutputs_T56_vcf at generation " <<GP->CurrentGeneration<<"\n"
    const std::string s_tab = "\t";
    const std::string s_P   = "P";
    const std::string s_dotI = ".I";


    file.open();

    // Write header
    {
        std::string s;

        s += std::string("##fileformat=VCFv4.2\n##Generation=") + std::to_string( GP->CurrentGeneration) + std::string("\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
        s.reserve(s.size() + SSP->TotalpatchCapacity * 14);
        
        for (int patch_index = 0 ; patch_index < GP->PatchNumber ; ++patch_index)
        {
            for (int ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ++ind_index)
            {
                s += s_tab + s_P + std::to_string(patch_index) + s_dotI + std::to_string(ind_index);
            }
        }
        s += std::string("\n");
        //std::cout << "writing: " << s << "\n";
        file.write(s);
    }

    if (SSP->TotalpatchSize)
    {
        // Compute frequencies
        std::vector<unsigned> obsFreqs = pop.computeT56Frequencies();

          
        // Write Data
        WriteOutputs_T56_vcf_writeData(pop, obsFreqs, file);
    }
    file.close();
}*/

void OutputWriter::WriteOutputs_T56_vcf(Pop& pop, OutputFile& file)
{
    //std::cout << "In OutputWriter::WriteOutputs_T56_vcf at generation " <<GP->CurrentGeneration<<"\n"
    const std::string s_tab = "\t";
    const std::string s_P   = "P";
    const std::string s_dotI = ".I";
    const std::string s_vert = "|";


    auto alleleFrequencies = pop.computePatchSpecificT56Frequencies();

    // One file per patch
    for (int patch_index = 0 ; patch_index < GP->PatchNumber ; patch_index++)
    {
        file.openPatchSpecific(patch_index);

        ////////////////////
        /// Write header ///
        ////////////////////
        {
            std::string s;

            s += std::string("##fileformat=VCFv4.2\n##Generation=") + std::to_string( GP->CurrentGeneration) + std::string("\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
            s.reserve(s.size() + SSP->patchSize[patch_index] * 15);
            
            for (int ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ++ind_index)
            {
                s += s_tab + s_P + std::to_string(patch_index) + s_dotI + std::to_string(ind_index);
            }
            
            s += std::string("\n");
            //std::cout << "writing: " << s << "\n";
            file.write(s);
        }


        std::vector<uint32_t> polymorphicLoci;
        polymorphicLoci.reserve(10);
        assert(alleleFrequencies[patch_index].size() == SSP->Gmap.T56_nbLoci);
        for (uint32_t locus = 0 ; locus < SSP->Gmap.T56_nbLoci ; ++locus)
        {
            double x = alleleFrequencies[patch_index][locus];
            if (x != 0.0 )// && x != 1.0)
            {
                polymorphicLoci.push_back(locus);                   
            }
        }



        ////////////////////////////////////////////
        /// Stupid repetitive code to write data ///
        ////////////////////////////////////////////
        
        if (SSP->patchSize[patch_index] && polymorphicLoci.size())
        {
            // initalize s
            // s will be cleared at the end of each line and reused for the next line because I think the compiler is too stupid to do it on his own.
            std::string s;
            s.reserve(
                4 * SSP->TotalpatchSize // main data
                + 80 // intro line that follows + security
            );



            Patch& patch = pop.getPatch(patch_index);

            if (SSP->Gmap.isT56ntrlCompress)
            {
                // T6ntrl
                if (SSP->Gmap.isT56selCompress)
                {
                    // T6sel
                    for (auto locus : polymorphicLoci)
                    {
                        // Beginning of line
                        s += std::to_string(1) + std::string("\t") + std::to_string(locus + 1) + std::string("\t.\tA\tT\t.\t.\t.\tGT");

                        // write data
                        auto locusInGender = SSP->Gmap.FromT56LocusToT56genderLocus(locus).locusInGender;
                        for (uint32_t ind_index = 0 ; ind_index < SSP->patchSize[patch_index]; ++ind_index)
                        {
                            auto ntrlIt0     = patch.getInd(ind_index).getHaplo(0).T6ntrl_AllelesBegin();
                            auto selIt0      = patch.getInd(ind_index).getHaplo(0).T6sel_AllelesBegin();
                            auto ntrlIt1     = patch.getInd(ind_index).getHaplo(1).T6ntrl_AllelesBegin();
                            auto selIt1      = patch.getInd(ind_index).getHaplo(1).T6sel_AllelesBegin();

                            auto ntrlIt0End     = patch.getInd(ind_index).getHaplo(0).T6ntrl_AllelesEnd();
                            auto selIt0End      = patch.getInd(ind_index).getHaplo(0).T6sel_AllelesEnd();
                            auto ntrlIt1End     = patch.getInd(ind_index).getHaplo(1).T6ntrl_AllelesEnd();
                            auto selIt1End      = patch.getInd(ind_index).getHaplo(1).T6sel_AllelesEnd();

                            if (SSP->Gmap.FromT56LocusToT56genderLocus(locus).isNtrl)
                            {
                                // ntrl
                                while (ntrlIt0 != ntrlIt0End && *ntrlIt0 < locusInGender) ++ntrlIt0;
                                while (ntrlIt1 != ntrlIt1End && *ntrlIt1 < locusInGender) ++ntrlIt1;

                                s += s_tab + std::to_string(ntrlIt0 != ntrlIt0End && *ntrlIt0 == locusInGender) + s_vert + std::to_string(ntrlIt1 != ntrlIt1End && *ntrlIt1 == locusInGender);
                            } else
                            {
                                // sel
                                while (selIt0 != selIt0End && *selIt0 < locusInGender) ++selIt0;
                                while (selIt1 != selIt1End && *selIt1 < locusInGender) ++selIt1;

                                s += s_tab + std::to_string(selIt0 != selIt0End && *selIt0 == locusInGender) + s_vert + std::to_string(selIt1 != selIt1End && *selIt1 == locusInGender);
                            }
                        }

                         // write one line at a time
                        s += std::string("\n");
                        file.write(s);
                        s.resize(0);
                    }
                    
                } else
                {
                    // T5sel
                    for (auto locus : polymorphicLoci)
                    {
                        // Beginning of line
                        s += std::to_string(1) + std::string("\t") + std::to_string(locus + 1) + std::string("\t.\tA\tT\t.\t.\t.\tGT");

                        // write data
                        auto locusInGender = SSP->Gmap.FromT56LocusToT56genderLocus(locus).locusInGender;
                        for (uint32_t ind_index = 0 ; ind_index < SSP->patchSize[patch_index]; ++ind_index)
                        {
                            auto ntrlIt0     = patch.getInd(ind_index).getHaplo(0).T6ntrl_AllelesBegin();
                            auto selIt0      = patch.getInd(ind_index).getHaplo(0).T5sel_AllelesBegin();
                            auto ntrlIt1     = patch.getInd(ind_index).getHaplo(1).T6ntrl_AllelesBegin();
                            auto selIt1      = patch.getInd(ind_index).getHaplo(1).T5sel_AllelesBegin();

                            auto ntrlIt0End     = patch.getInd(ind_index).getHaplo(0).T6ntrl_AllelesEnd();
                            auto selIt0End      = patch.getInd(ind_index).getHaplo(0).T5sel_AllelesEnd();
                            auto ntrlIt1End     = patch.getInd(ind_index).getHaplo(1).T6ntrl_AllelesEnd();
                            auto selIt1End      = patch.getInd(ind_index).getHaplo(1).T5sel_AllelesEnd();

                            if (SSP->Gmap.FromT56LocusToT56genderLocus(locus).isNtrl)
                            {
                                // ntrl
                                while (ntrlIt0 != ntrlIt0End && *ntrlIt0 < locusInGender) ++ntrlIt0;
                                while (ntrlIt1 != ntrlIt1End && *ntrlIt1 < locusInGender) ++ntrlIt1;

                                s += s_tab + std::to_string(ntrlIt0 != ntrlIt0End && *ntrlIt0 == locusInGender) + s_vert + std::to_string(ntrlIt1 != ntrlIt1End && *ntrlIt1 == locusInGender);
                            } else
                            {
                                // sel
                                while (selIt0 != selIt0End && *selIt0 < locusInGender) ++selIt0;
                                while (selIt1 != selIt1End && *selIt1 < locusInGender) ++selIt1;

                                s += s_tab + std::to_string(selIt0 != selIt0End && *selIt0 == locusInGender) + s_vert + std::to_string(selIt1 != selIt1End && *selIt1 == locusInGender);
                            }
                        }
                        // write one line at a time
                        s += std::string("\n");
                        file.write(s);
                        s.resize(0);
                    }
                }
            } else
            {
                // T5ntrl
                if (SSP->Gmap.isT56selCompress)
                {
                    //T6sel
                    for (auto locus : polymorphicLoci)
                    {
                        // Beginning of line
                        s += std::to_string(1) + std::string("\t") + std::to_string(locus + 1) + std::string("\t.\tA\tT\t.\t.\t.\tGT");

                        // write data
                        auto locusInGender = SSP->Gmap.FromT56LocusToT56genderLocus(locus).locusInGender;
                        for (uint32_t ind_index = 0 ; ind_index < SSP->patchSize[patch_index]; ++ind_index)
                        {
                            auto ntrlIt0     = patch.getInd(ind_index).getHaplo(0).T5ntrl_AllelesBegin();
                            auto selIt0      = patch.getInd(ind_index).getHaplo(0).T6sel_AllelesBegin();
                            auto ntrlIt1     = patch.getInd(ind_index).getHaplo(1).T5ntrl_AllelesBegin();
                            auto selIt1      = patch.getInd(ind_index).getHaplo(1).T6sel_AllelesBegin();

                            auto ntrlIt0End     = patch.getInd(ind_index).getHaplo(0).T5ntrl_AllelesEnd();
                            auto selIt0End      = patch.getInd(ind_index).getHaplo(0).T6sel_AllelesEnd();
                            auto ntrlIt1End     = patch.getInd(ind_index).getHaplo(1).T5ntrl_AllelesEnd();
                            auto selIt1End      = patch.getInd(ind_index).getHaplo(1).T6sel_AllelesEnd();

                            if (SSP->Gmap.FromT56LocusToT56genderLocus(locus).isNtrl)
                            {
                                // ntrl
                                while (ntrlIt0 != ntrlIt0End && *ntrlIt0 < locusInGender) ++ntrlIt0;
                                while (ntrlIt1 != ntrlIt1End && *ntrlIt1 < locusInGender) ++ntrlIt1;

                                s += s_tab + std::to_string(ntrlIt0 != ntrlIt0End && *ntrlIt0 == locusInGender) + s_vert + std::to_string(ntrlIt1 != ntrlIt1End && *ntrlIt1 == locusInGender);
                            } else
                            {
                                // sel
                                while (selIt0 != selIt0End && *selIt0 < locusInGender) ++selIt0;
                                while (selIt1 != selIt1End && *selIt1 < locusInGender) ++selIt1;

                                s += s_tab + std::to_string(selIt0 != selIt0End && *selIt0 == locusInGender) + s_vert + std::to_string(selIt1 != selIt1End && *selIt1 == locusInGender);
                            }
                        }

                        // write one line at a time
                        s += std::string("\n");
                        file.write(s);
                        s.resize(0);
                    }
                 
                } else
                {
                    //T5sel
                    for (auto locus : polymorphicLoci)
                    {
                        // Beginning of line
                        s += std::to_string(1) + std::string("\t") + std::to_string(locus + 1) + std::string("\t.\tA\tT\t.\t.\t.\tGT");

                        // write data
                        auto locusInGender = SSP->Gmap.FromT56LocusToT56genderLocus(locus).locusInGender;
                        for (uint32_t ind_index = 0 ; ind_index < SSP->patchSize[patch_index]; ++ind_index)
                        {
                            auto ntrlIt0     = patch.getInd(ind_index).getHaplo(0).T5ntrl_AllelesBegin();
                            auto selIt0      = patch.getInd(ind_index).getHaplo(0).T5sel_AllelesBegin();
                            auto ntrlIt1     = patch.getInd(ind_index).getHaplo(1).T5ntrl_AllelesBegin();
                            auto selIt1      = patch.getInd(ind_index).getHaplo(1).T5sel_AllelesBegin();

                            auto ntrlIt0End     = patch.getInd(ind_index).getHaplo(0).T5ntrl_AllelesEnd();
                            auto selIt0End      = patch.getInd(ind_index).getHaplo(0).T5sel_AllelesEnd();
                            auto ntrlIt1End     = patch.getInd(ind_index).getHaplo(1).T5ntrl_AllelesEnd();
                            auto selIt1End      = patch.getInd(ind_index).getHaplo(1).T5sel_AllelesEnd();

                            if (SSP->Gmap.FromT56LocusToT56genderLocus(locus).isNtrl)
                            {
                                // ntrl
                                while (ntrlIt0 != ntrlIt0End && *ntrlIt0 < locusInGender) ++ntrlIt0;
                                while (ntrlIt1 != ntrlIt1End && *ntrlIt1 < locusInGender) ++ntrlIt1;

                                s += s_tab + std::to_string(ntrlIt0 != ntrlIt0End && *ntrlIt0 == locusInGender) + s_vert + std::to_string(ntrlIt1 != ntrlIt1End && *ntrlIt1 == locusInGender);
                            } else
                            {
                                // sel
                                while (selIt0 != selIt0End && *selIt0 < locusInGender) ++selIt0;
                                while (selIt1 != selIt1End && *selIt1 < locusInGender) ++selIt1;

                                s += s_tab + std::to_string(selIt0 != selIt0End && *selIt0 == locusInGender) + s_vert + std::to_string(selIt1 != selIt1End && *selIt1 == locusInGender);
                            }
                        }

                        // write one line at a time
                        s += std::string("\n");
                        file.write(s);
                        s.resize(0);
                    }
                }
            }

            s.resize(0);
        }


        file.close();
    }
}


void OutputWriter::WriteOutputs_T3_LargeOutput_header(OutputFile& file)
{
    
    std::string s("Generation");
    for ( int patch_index = 0 ; patch_index < GP->maxEverPatchNumber ; ++patch_index )
    {
        for (int ind_index = 0 ; ind_index < SSP->maxEverpatchCapacity[patch_index] ; ind_index++)
        {
            for (int haplo_index = 0 ; haplo_index < SSP->ploidy ; haplo_index++)
            {
                for (int T3_locus = 0 ; T3_locus < SSP->Gmap.T3_nbLoci; T3_locus++)
                {
                    s += std::string("\tP") + std::to_string(patch_index) + std::string("_I" ) + std::to_string(ind_index) + std::string("_H" ) + std::to_string(haplo_index) + std::string("_L" ) + std::to_string(T3_locus);
                }
            }
        }
    }
    s += std::string("\n");
    
    file.open();
    file.write(s);
    file.close();
}


void OutputWriter::WriteOutputs_T3_LargeOutput(Pop& pop, OutputFile& file)
{
    const std::string s_tab("\t");
    std::string s;
    s.reserve(SSP->TotalpatchCapacity * SSP->Gmap.T3_nbLoci * 8 );

    s += std::to_string(GP->CurrentGeneration);
    for ( int patch_index = 0 ; patch_index < GP->maxEverPatchNumber ; ++patch_index )
    {
        for (int ind_index = 0 ; ind_index < SSP->maxEverpatchCapacity[patch_index] ; ind_index++)
        {
            for (int haplo_index = 0 ; haplo_index < SSP->ploidy ; haplo_index++)
            {
                for (int T3_locus = 0 ; T3_locus < SSP->Gmap.T3_nbLoci; T3_locus++)
                {
                    if (shouldNABePrinted(patch_index, ind_index))
                    {
                        s += s_tab + "NA";
                    } else
                    {
                        s += s_tab + std::to_string(pop.getPatch(patch_index).getInd(ind_index).getHaplo(haplo_index).getT3_Allele(T3_locus));
                    }
                }
            }
        }
    }

    s += std::string("\n");
    file.open();
    file.write(s);
    file.close();
}


void OutputWriter::WriteOutputs_T3_MeanVar_header(OutputFile& file)
{
    std::string s("Generation");
    for ( int patch_index = 0 ; patch_index < GP->maxEverPatchNumber ; ++patch_index )
    {
        s += std::string("\tP") + std::to_string(patch_index) + std::string("_meanT3fitness");
        s += std::string("\tP") + std::to_string(patch_index) + std::string("_varT3fitness");
        for (int dim = 0 ; dim < SSP->T3_PhenoNbDimensions ; dim++)
        {
            s += std::string("\tP") + std::to_string(patch_index) + std::string("_dim") + std::to_string(dim) + std::string("_meanPheno");
            s += std::string("\tP") + std::to_string(patch_index) + std::string("_dim") + std::to_string(dim) + std::string("_varPheno");
        }
    }
    s += std::string("\n");

    file.open();
    file.write(s);
    file.close();
}


void OutputWriter::WriteOutputs_T3_MeanVar(Pop& pop, OutputFile& file)
{
    const std::string s_tab("\t");
    std::string s;
    s.reserve(GP->PatchNumber * 16 );

    std::vector<double> T3fits;

    s += std::to_string(GP->CurrentGeneration) + s_tab;
    for ( int patch_index = 0 ; patch_index < GP->maxEverPatchNumber ; ++patch_index )
    {
        if (shouldNABePrinted(patch_index))
        {
            s += s_tab + "NA" + s_tab + "NA";
            for (int dim = 0 ; dim < SSP->T3_PhenoNbDimensions; dim++)
            {
                s += s_tab + "NA" + s_tab + "NA";
            }
        } else
        {
            const int& Habitat = SSP->Habitats[patch_index];
            T3fits.resize(SSP->patchSize[patch_index]);
            Patch& CurrentPatch = pop.getPatch(patch_index);

            // Set what we want to calculate to 0.0 for phneotypes. Note that phenotypes are have one value per dimension
            std::vector<double> MeanP;
            std::vector<double> VarP;
            std::vector<std::vector<double>> Phenotypes;
            for (int dim = 0 ; dim < SSP->T3_PhenoNbDimensions; dim++)
            {
                MeanP.push_back(0.0);
                VarP.push_back(0.0);

                std::vector<double> v(SSP->patchSize[patch_index]);
                Phenotypes.push_back(v);
            }

            // Set what we want to calculate to 0.0 for fitnesses. Note that fitnesses are have one value overall (one value per patch, not one value per patch per dimension)
            double MeanFit = 0.0; // Only fitness regarding phenotype (T3 loci)
            double VarFit = 0.0;  // Only fitness regarding phenotype (T3 loci)
            std::vector<double> Fitnesses(SSP->patchSize[patch_index]);

            // calculate all phenotypes and fitnesses
            for (int ind_index=0;ind_index<SSP->patchSize[patch_index];++ind_index)
            {
                Individual::resetT3phenotype();
                CurrentPatch.getInd(ind_index).CalculateT3Phenotype(Habitat); // Save Phenotype in the static 'T3_IndPhenotype'
                double W;
                if (SSP->T3_isSelection)
                {
                    W = Individual::CalculateT3Fitness(Habitat);  
                } else
                {
                    W = 1.0;
                }
                Fitnesses[ind_index] = W;
                MeanFit += W;

                for (int dim = 0 ; dim < SSP->T3_PhenoNbDimensions; dim++)
                {
                    Phenotypes[dim][ind_index] = Individual::T3_IndPhenotype[dim];
                    MeanP[dim] += Individual::T3_IndPhenotype[dim];
                }
            }

            // Get means and variances for fitness
            MeanFit /= SSP->patchSize[patch_index];
            for (int ind_index=0;ind_index<SSP->patchSize[patch_index];++ind_index)
            {
                VarFit += pow(Fitnesses[ind_index] - MeanFit,2);
            }
            VarFit /= SSP->patchSize[patch_index];

            s += s_tab + std::to_string(MeanFit) + s_tab + std::to_string(VarFit);

            // Get means and variances for phenotypes (one value per dimension)
            for (int dim = 0 ; dim < SSP->T3_PhenoNbDimensions; dim++)
            {
                MeanP[dim] /= SSP->patchSize[patch_index];
                for (int ind_index=0;ind_index<SSP->patchSize[patch_index];++ind_index)
                {
                    VarP[dim] += pow(Phenotypes[dim][ind_index] - MeanP[dim],2);
                }
                VarP[dim] /= SSP->patchSize[patch_index];

                s += s_tab + std::to_string(MeanP[dim]) + s_tab + std::to_string(VarP[dim]);
            }
        }
    }

    s += std::string("\n");

    file.open();
    file.write(s);
    file.close();
}


void OutputWriter::WriteOutputs_extinction(OutputFile& file)
{
    assert(GP->nbGenerations == GP->CurrentGeneration);

    SpeciesSpecificParameters* oldSSP = SSP;

    std::string s;
    s += "Species\tGeneration\n";
    for (int speciesIndex = 0 ; speciesIndex < GP->nbSpecies ; speciesIndex++)
    {
        SSP = allParameters.getSSPaddress(speciesIndex);
        int g = SSP->whenDidExtinctionOccur;
        std::string& speciesName = SSP->speciesName;
        if (g != -1)
        {
            assert(g >= 0 && g <= GP->nbGenerations);
            s += speciesName + "\t" + std::to_string(g) + "\n";
        } else
        {
            s += speciesName + std::string("\tNA\n");
        }
    }

    SSP = oldSSP;

    file.open();
    file.write(s);
    file.close();
}


void OutputWriter::WriteOutputs_T1or4or56SFS_header(OutputFile& file, size_t mutPlacingIndex)
{
    std::string s;
    
    s += "Generation";

    for (uint32_t patch_index = 0 ; patch_index < GP->PatchNumber ; patch_index++)
    {
        for (double binSize : SSP->outputSFSbinSizes)
        {
            uint32_t nbBins = (uint32_t) std::ceil(1 / binSize);
            double from = 0.0;
            for (uint32_t bin_index = 0 ; bin_index < nbBins ; bin_index++)
            {
                double to = from + binSize;
                if (to > 1.0)
                    to = 1.0;

                s += "\tP" + std::to_string(patch_index) + "_" + std::to_string(from) + "-" + std::to_string(to);

                from = to;
            }
        }
    }

    for (double binSize : SSP->outputSFSbinSizes)
    {
        uint32_t nbBins = (uint32_t) std::ceil(1 / binSize);
        double from = 0.0;
        for (uint32_t bin_index = 0 ; bin_index < nbBins ; bin_index++)
        {
            double to = from + binSize;
            if (to > 1.0)
                to = 1.0;

            s += "\tallPatches" + std::to_string(from) + "-" + std::to_string(to);

            from = to;
        }
    }

    s += "\tnbSNPs\n";

    if (mutPlacingIndex == 0)
    {
        file.openT4(mutPlacingIndex);
        file.write(s);
        file.close();
    } else
    {
        file.open();
        file.write(s);
        file.close();
    }
}


void OutputWriter::WriteOutputs_T4SFS_header(OutputFile& file, size_t mutPlacingIndex)
{
    WriteOutputs_T1or4or56SFS_header(file, mutPlacingIndex);
}

/*void OutputWriter::WriteOutputs_T4CoalescenceFst_header(OutputFile& file)
{
    std::string s("Generation\tTt\tTs\tTb\tGst\tWCst\n");
    file.open();
    file.write(s);
    file.close();
}

void OutputWriter::WriteOutputs_T4CoalescenceFst(OutputFile& file)
{
    std::string s;
    std::string s_tab("\t");
    s += std::to_string(GP->CurrentGeneration) + s_tab;

    std::vector<double> v = SSP->T4Tree.computeCoalescenceFstStatistics();
    assert(v.size() == 5);
    s += std::to_string(v[0]) + s_tab + std::to_string(v[1]) + s_tab + std::to_string(v[2]) + s_tab + std::to_string(v[3]) + s_tab + std::to_string(v[4]) + std::string("\n");

    file.open();
    file.write(s);
    file.close();   
}*/

void OutputWriter::WriteOutputs_T1SFS_header(OutputFile& file)
{
    WriteOutputs_T1or4or56SFS_header(file);
}

void OutputWriter::WriteOutputs_T56SFS_header(OutputFile& file)
{
    WriteOutputs_T1or4or56SFS_header(file);
}



void OutputWriter::WriteOutputs_T4SFS(OutputFile& file, std::vector<std::vector<std::vector<uint32_t>>>& data, size_t mutPlacingIndex)
{
    if (SSP->whenDidExtinctionOccur != -1)
    {
        file.openT4(mutPlacingIndex);
        std::string s("Species is extinct\n");
        file.write(s);
        file.close();
        return;
    }

    //SSP->T4Tree.simplify(pop);
    //auto data = SSP->T4Tree.placeMutations(pop);
    std::vector<std::vector<double>> obsFreqs(GP->PatchNumber);
    for (uint32_t patch_index = 0 ; patch_index < GP->PatchNumber ; patch_index++)
    {
        obsFreqs[patch_index].resize(SSP->patchSize[patch_index], 0.0);
        assert(data[patch_index].size() == 2*SSP->patchSize[patch_index]);
        for (uint32_t ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ind_index++)
        {
            for (uint32_t haplo_index = 0 ; haplo_index < 2 ; haplo_index++)
            {
                for (auto& mut : data[patch_index][ind_index*2+haplo_index])
                    ++obsFreqs[patch_index][mut];
            }
        }
    }

    for (uint32_t patch_index = 0 ; patch_index < GP->PatchNumber ; patch_index++)
    {
        for (uint32_t locus = 0 ; locus < SSP->Gmap.T4_nbLoci ; ++locus)
        {
            obsFreqs[patch_index][locus] /= 2*SSP->patchSize[patch_index];
            assert(obsFreqs[patch_index][locus] >= 0.0 && obsFreqs[patch_index][locus] <= 1.0);
        }
    }

    WriteOutputs_T1or4or5SFS(obsFreqs, file, mutPlacingIndex);
}

void OutputWriter::WriteOutputs_T1SFS(Pop& pop, OutputFile& file)
{
    // Build obsFreqs
    std::vector<std::vector<double>> obsFreqs(GP->PatchNumber); // obsFreqs[patch_index][locus]
    for (uint32_t patch_index = 0 ; patch_index < GP->PatchNumber ; patch_index++)
    {
        // count
        obsFreqs[patch_index].resize(SSP->Gmap.T1_nbLoci, 0.0);
        for (uint32_t ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ind_index++)
        {
            for (uint32_t haplo_index = 0 ; haplo_index < 2 ; haplo_index++)
            {
                for (uint32_t locus = 0 ; locus < SSP->Gmap.T1_nbLoci ; locus++)
                {
                    if (pop.getPatch(patch_index).getInd(ind_index).getHaplo(haplo_index).getT1_Allele(locus / 8, locus % 8))
                        obsFreqs[patch_index][locus]++;
                
                }
            }
        }

        // divide
        for (uint32_t locus = 0 ; locus < SSP->Gmap.T1_nbLoci ; locus++)
        {
            obsFreqs[patch_index][locus] /= 2 * SSP->patchSize[patch_index];
            assert(obsFreqs[patch_index][locus] >= 0.0 && obsFreqs[patch_index][locus] <= 1.0);
        }
    }
    
    
    WriteOutputs_T1or4or5SFS(obsFreqs, file);
}


void OutputWriter::WriteOutputs_T1_haplotypeFreqs(Pop& pop, OutputFile& file)
{
    std::vector<std::string> uniqueHaplotypes;
    std::vector<uint32_t> haplotypeFrequencies;

    auto& loci = file.subset[SSP->speciesIndex];

    for (uint32_t patch_index = 0 ; patch_index < GP->PatchNumber; ++patch_index )
    {
        auto& patch = pop.getPatch(patch_index);
        for (uint32_t ind_index = 0 ; ind_index < SSP->patchSize[patch_index]; ++ind_index )
        {
            auto& ind = patch.getInd(ind_index);
            for (uint32_t haplo_index = 0 ; haplo_index < 2; ++haplo_index )
            {
                auto& haplo = ind.getHaplo(haplo_index);

                std::string haplo_s;
                haplo_s.reserve(loci.size());
                for (auto& locus : loci)
                {
                    haplo_s += std::to_string(haplo.getT1_Allele(locus.byte_index, locus.bit_index));
                }

                auto it = std::lower_bound(uniqueHaplotypes.begin(), uniqueHaplotypes.end(), haplo_s);
                uint32_t index = it - uniqueHaplotypes.begin();
                if (it != uniqueHaplotypes.end() && *it == haplo_s)
                {
                    ++haplotypeFrequencies[index];
                } else
                {
                    uniqueHaplotypes.insert(it, haplo_s);
                    haplotypeFrequencies.insert(haplotypeFrequencies.begin() + index, 1);
                }
            }
        }
    }
    assert(uniqueHaplotypes.size() > 0);
    assert(haplotypeFrequencies.size() == uniqueHaplotypes.size());
    assert(std::accumulate(haplotypeFrequencies.begin(), haplotypeFrequencies.end(), 0) == 2*SSP->TotalpatchSize);

    std::string s;
    s.reserve(uniqueHaplotypes.size() * uniqueHaplotypes[0].size() * 3);
    s += "Generation\t";
    std::string s_tab("\t");
    for (uint32_t index = 0 ; index < uniqueHaplotypes.size() ; ++index)
    {
        auto& haplo_s = uniqueHaplotypes[index];
        s += haplo_s + s_tab;
    }
    s += "\n" + std::to_string(GP->CurrentGeneration) + s_tab;

    for (uint32_t index = 0 ; index < uniqueHaplotypes.size() ; ++index)
    {
        double freq = (double) haplotypeFrequencies[index] / (2*SSP->TotalpatchSize);
        s += std::to_string(freq) + s_tab;
    }
    s += "\n";

    file.open();
    file.write(s);
    file.close();
}



void OutputWriter::WriteOutputs_T56SFS(Pop& pop, OutputFile& file)
{
    std::vector<std::vector<unsigned>> obsFreqs = pop.computePatchSpecificT56Frequencies();
    std::vector<std::vector<double>> obsRelFreqs(obsFreqs.size());
    
    // divide
    for (int patch_index = 0 ; patch_index < GP->PatchNumber ; patch_index++)
    {
        obsRelFreqs[patch_index].resize(SSP->Gmap.T56_nbLoci,0.0);
        for (uint32_t locus = 0 ; locus < SSP->Gmap.T56_nbLoci ; locus++)
        {
            obsRelFreqs[patch_index][locus] = (double) obsFreqs[patch_index][locus] / (2.0 * SSP->patchSize[patch_index]);
            assert(obsRelFreqs[patch_index][locus] >= 0.0 && obsRelFreqs[patch_index][locus] <= 1.0);
        }
    }

    WriteOutputs_T1or4or5SFS(obsRelFreqs, file);
}


void OutputWriter::WriteOutputs_T1or4or5SFS(std::vector<std::vector<double>>& obsFreqs, OutputFile& file, size_t mutPlacingIndex)
{
    // Format is obsFreqs[patch_index][locus]

    std::string s;
    s += std::to_string(GP->CurrentGeneration);

    uint32_t nbLoci = obsFreqs[0].size();
    std::vector<double> obsRelFreqsPop(nbLoci, 0.0); // obsRelFreqsPop[locus] // for the entire population



    // Build s
    for (uint32_t patch_index = 0 ; patch_index < GP->PatchNumber ; patch_index++)
    {
        assert(obsFreqs[patch_index].size() == nbLoci);

        // order loci by frequency. I could merge sorting with the next loop for performance
        std::sort(obsFreqs[patch_index].begin(), obsFreqs[patch_index].end());

        // set obsRelFreqsPop
        for (uint32_t locus = 0 ; locus < nbLoci ; locus++)
            obsRelFreqsPop[locus] += obsFreqs[patch_index][locus]; // division by patchNumber later

        // build s
        for (double binSize : SSP->outputSFSbinSizes)
        {
            uint32_t TotalNbLociFound = 0; // For security

            uint32_t nbBins = (uint32_t) std::ceil(1 / binSize);
            double from = 0.0;
            for (uint32_t bin_index = 0 ; bin_index < nbBins ; bin_index++)
            {
                double to = from + binSize;
                if (to == 1.0)
                    to = 1.01; // To include the frequencies of exactly 1.0

                // as usual, 'from' is included and 'to' is excluded. But there is exception for the frequency of exactly 1.0, where 'to' is included
                
                auto itFrom = std::lower_bound(obsFreqs[patch_index].begin(), obsFreqs[patch_index].end(), from);
                auto itTo   = std::lower_bound(itFrom, obsFreqs[patch_index].end(), to);

                uint32_t nbLociFound = itTo - itFrom;
                TotalNbLociFound += nbLociFound;
                double freqInBin = (double)nbLociFound / (double)nbLoci;
                    
                s += "\t" + std::to_string(freqInBin);

                from = to;
            }
            assert(TotalNbLociFound == nbLoci);
        }
    }

    // set obsRelFreqsPop
    uint32_t nbSNPs = 0;
    for (uint32_t locus = 0 ; locus < nbLoci ; locus++)
    {
        obsRelFreqsPop[locus] /= GP->PatchNumber;
        
        if (obsRelFreqsPop[locus] == 0.0 || obsRelFreqsPop[locus] == 1.0 )
        {
            nbSNPs++;
        } else
        {
            assert(obsRelFreqsPop[locus] > 0.0 && obsRelFreqsPop[locus] < 1.0);
        }
    }

    // For the whole population
    for (double binSize : SSP->outputSFSbinSizes)
    {
        uint32_t TotalNbLociFound = 0; // For security

        uint32_t nbBins = (uint32_t) std::ceil(1 / binSize);
        double from = 0.0;
        for (uint32_t bin_index = 0 ; bin_index < nbBins ; bin_index++)
        {
            double to = from + binSize;
            if (to == 1.0)
                to = 1.01; // To include the frequencies of exactly 1.0

            // as usual, 'from' is included and 'to' is excluded. But there is exception for the frequency of exactly 1.0, where 'to' is included
            
            auto itFrom = std::lower_bound(obsRelFreqsPop.begin(), obsRelFreqsPop.end(), from);
            auto itTo   = std::lower_bound(itFrom, obsRelFreqsPop.end(), to);

            uint32_t nbLociFound = itTo - itFrom;
            TotalNbLociFound += nbLociFound;
            double freqInBin = (double)nbLociFound / (double)nbLoci;
                
            s += "\t" + std::to_string(freqInBin);

            from = to;
        }
        assert(TotalNbLociFound == nbLoci);
    }

    s += "\t" + std::to_string(nbSNPs) + "\n";

    if (mutPlacingIndex == 0)
    {
        file.openT4(mutPlacingIndex);
        file.write(s);
        file.close();
    } else
    {
        file.open();
        file.write(s);
        file.close();
    }
}




void OutputWriter::imitateSequencingError(Pop& pop)
{
    assert(SSP != NULL);
    assert(SSP->patchSize.size() == GP->PatchNumber);
    for (int patch_index = 0 ; patch_index < GP->PatchNumber ; patch_index++)
    {
        Patch& patch = pop.getPatch(patch_index);
        for (int ind_index = 0 ; ind_index < SSP->patchSize[patch_index] ; ind_index++)
        {
            for (int haplo_index = 0  ; haplo_index < 2 ; haplo_index++)
            {
                imitateSequencingError(patch.getInd(ind_index).getHaplo(haplo_index));
            }
        }
    }
}

void OutputWriter::imitateSequencingError(Haplotype& TransmittedChrom)
{
    // This piece of code should be very related to Mutate_T1 from LifeCycle
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'OutputWriter::imitateSequencingError'\n";
#endif
    int Habitat = 0;

    // T1
    if (SSP->Gmap.T1_nbLoci)
    {
        TransmittedChrom.setAllW_T1(-1);
        std::poisson_distribution<uint32_t> d_nbMuts(GP->sequencingErrorRate * SSP->Gmap.T1_nbLoci);
        auto nbMuts = d_nbMuts(GP->rngw.getRNG());

        for (uint32_t i = 0 ; i < nbMuts ; i++)
        {
            auto MutPosition = GP->rngw.uniform_int_distribution(SSP->Gmap.T1_nbLoci);
            TransmittedChrom.mutateT1_Allele(MutPosition, Habitat);
        }    
    }

    // T2
    if (SSP->Gmap.T2_nbLoci)
    {
        TransmittedChrom.setAllW_T2(-1);
        std::poisson_distribution<uint32_t> d_nbMuts(GP->sequencingErrorRate * SSP->Gmap.T2_nbLoci);
        auto nbMuts = d_nbMuts(GP->rngw.getRNG());

        for (uint32_t i = 0 ; i < nbMuts ; i++)
        {
            auto MutPosition = GP->rngw.uniform_int_distribution(SSP->Gmap.T2_nbLoci);
            TransmittedChrom.AddMutT2_Allele(MutPosition, Habitat);
        }    
    }

    // T3
    if (SSP->Gmap.T3_nbLoci)
    {
        std::poisson_distribution<uint32_t> d_nbMuts(GP->sequencingErrorRate * SSP->Gmap.T3_nbLoci);
        auto nbMuts = d_nbMuts(GP->rngw.getRNG());

        for (uint32_t i = 0 ; i < nbMuts ; i++)
        {
            auto MutPosition = GP->rngw.uniform_int_distribution(SSP->Gmap.T3_nbLoci);
            TransmittedChrom.mutateT3_Allele(MutPosition);
        }    
    }

    // T4 Sorry can't do it

    // T56
    if (SSP->Gmap.T56_nbLoci)
    {
        std::poisson_distribution<uint32_t> d_nbMuts(GP->sequencingErrorRate * SSP->Gmap.T56_nbLoci);
        auto nbMuts = d_nbMuts(GP->rngw.getRNG());

        for (uint32_t i = 0 ; i < nbMuts ; i++)
        {
            auto MutPosition = GP->rngw.uniform_int_distribution(SSP->Gmap.T56_nbLoci);

            auto locusGender = SSP->Gmap.FromT56LocusToT56genderLocus(MutPosition);
            if (locusGender.isNtrl)
            {
                TransmittedChrom.mutateT56ntrl_Allele(locusGender.locusInGender);
            } else
            {
                TransmittedChrom.mutateT56sel_Allele(locusGender.locusInGender, Habitat);
            }
        }    
    }
}

void OutputWriter::WriteOutputs(Pop& realPop)
{
    /*
    std::cout << "GP->CurrentGeneration = " << GP->CurrentGeneration << "\n";
    std::cout << "this->isTime() = " << this->isTime() << "\n";
    */    

    if (this->isTime())
    { 
        // output the pop
        //std::cout << "About to enter WriteOutputs_forDefinedPop\n";
        WriteOutputs_forDefinedPop(realPop);

        // output the population with sequencing errors
        if (GP->sequencingErrorRate > 0.0)
        {
            // Maybe not the best solution (will much increase the RAM usage) but by far the easiest to implement!

            // Copy the entire pop.
            Pop seqErrorPop = realPop;

            // Make the sequencing error
            imitateSequencingError(seqErrorPop);

            // Produce outputs
            OutputFile::sequencingErrorStringToAddToFilnames = "_sequencingError";
            WriteOutputs_forDefinedPop(seqErrorPop);
            OutputFile::sequencingErrorStringToAddToFilnames = "";
        }

    }




    /*
     ##################################
     #### genealogy -> coalescence ####
     ##################################
     */
    if (SSP->genealogy.isCoalesce())
    {
        #ifdef DEBUG
        std::cout << "Write coalescence\n";
        #endif
        if (GP->CurrentGeneration == GP->startAtGeneration)
        {
            assert(this->isFile(genealogy));
            auto& file = this->get_OutputFiles(genealogy)[0];
            SSP->genealogy.removeFilesAtStart(file);
        }
        if ( SSP->genealogy.isTimeToCoalesce() )
        {
            assert(this->isFile(genealogy));
            auto& file = this->get_OutputFiles(genealogy)[0];
            SSP->genealogy.coalesce(file);
        }
    }

    if (SSP->genealogy.isTimeToMerge())
    {
        #ifdef DEBUG
        std::cout << "Merge coalescence\n";
        #endif
        assert(this->isFile(genealogy));
        auto& file = this->get_OutputFiles(genealogy)[0];
        SSP->genealogy.mergeFiles(file);
    }


}

void OutputWriter::WriteOutputs_forDefinedPop_general(Pop& pop)
{
    /*
     ##################
     ### extinction ###
     ##################
     */
    //std::cout << "line 2674\n";
    
    if (SSP->speciesIndex == 0 && this->isFile(extinction) && GP->CurrentGeneration == GP->nbGenerations)
    {
        OutputFile& file = this->get_OutputFiles(extinction)[0];
        this->WriteOutputs_extinction(file);
    }
    
    /*
     ###############
     ### fitness ###
     ###############
     */
    //std::cout << "line 2687\n";
    if (this->isFile(fitness))
    {
        OutputFile& file = this->get_OutputFiles(fitness)[0];
        if ( GP->CurrentGeneration == GP->startAtGeneration)
        {
            this->WriteOutputs_fitness_header(file);
        }
        if (file.isTime())
        {
            #ifdef DEBUG
            std::cout << "Write fitness\n";
            #endif
            this->WriteOutputs_fitness(pop, file);
        }
    }

    /*
     #########################
     ### fitnessSubsetLoci ###
     #########################
     */
    //std::cout << "line 2709\n";
    if (this->isFile(fitnessSubsetLoci))
    {
        OutputFile& file = this->get_OutputFiles(fitnessSubsetLoci)[0];
        if ( GP->CurrentGeneration == GP->startAtGeneration)
        {
            this->WriteOutputs_fitnessSubsetLoci_header(file);
        }
        if (file.isTime())
        {
            #ifdef DEBUG
            std::cout << "Write fitnessSubsetLoci\n";
            #endif
            this->WriteOutputs_fitnessSubsetLoci(pop, file);
        }
    }

    /*
     ####################
     ### fitnessStats ###
     ####################
     */
    //std::cout << "line 2731\n";
    if (this->isFile(fitnessStats))
    {
        OutputFile& file = this->get_OutputFiles(fitnessStats)[0];
        if ( GP->CurrentGeneration == GP->startAtGeneration)
        {
            this->WriteOutputs_fitnessStats_header(file);
        }
        if (file.isTime())
        {
            #ifdef DEBUG
            std::cout << "Write fitnessStats\n";
            #endif
            this->WriteOutputs_fitnessStats(pop, file);
        }
    }

    

    /*
     ####################
     ##### patchSize ####
     ####################
     */
    //std::cout << "line 3088\n";
    if (this->isFile(patchSize))
    {
        auto& files = this->get_OutputFiles(patchSize);
        for (auto& file : files)
        {
            if ( GP->CurrentGeneration == GP->startAtGeneration)
            {
                this->WriteOutputs_patchSize_header(file);
            }
            if (file.isTime())
            {
                #ifdef DEBUG
                std::cout << "Write patchSize\n";
                #endif
                this->WriteOutputs_patchSize(file); // It does not need pop. It uses SSP
            }
        }
    }

    /*
     ##########################
     #### extraGeneticInfo ####
     ##########################
     */
    //std::cout << "line 3142\n";
    if (SSP->speciesIndex == 0 && this->isFile(extraGeneticInfo))
    {
        if (SSP->Gmap.T1_nbLoci || SSP->Gmap.T2_nbLoci)
        {
            if (SSP->Gmap.T3_nbLoci) std::cout << "WARNING: extraGeneticInfo will be computed ignoring loci of type 3!\n";
            std::vector<OutputFile>& files = this->get_OutputFiles(extraGeneticInfo);
            for (auto& file : files)
            {
                if (file.isTime())
                {
                    #ifdef DEBUG
                    std::cout << "Write extraGeneticInfo\n";
                    #endif
                    assert(GP->CurrentGeneration == GP->startAtGeneration);
                    this->WriteOutputs_extraGeneticInfo(file);
                }
            }
        }
    }
}

void OutputWriter::WriteOutputs_forDefinedPop_T1(Pop& pop)
{
    /*
     ###########################
     ### T1_Allele Frequency ###
     ###########################
     */
    //std::cout << "line 2756\n";
    if (this->isFile(T1_AlleleFreqFile))
    {
        if (SSP->Gmap.T1_nbLoci)
        {
            std::vector<OutputFile>& files = this->get_OutputFiles(T1_AlleleFreqFile);
            for (auto& file : files)
            {
                if ( GP->CurrentGeneration == GP->startAtGeneration)
                {
                    this->WriteOutputs_T1_AlleleFreq_header(file);
                }
                if (file.isTime())
                {
                    #ifdef DEBUG
                    std::cout << "Write T1_AlleleFreqFile\n";
                    #endif
                    this->WriteOutputs_T1_AlleleFreq(pop, file);
                }
            }
        }

    }


    /*
     ####################
     #### T1_MeanLD  ####
     ####################
     */
    //std::cout << "line 2786\n";
    if (this->isFile(MeanLDFile))
    {
        if (SSP->Gmap.T1_nbLoci)
        {
            std::vector<OutputFile>& files = this->get_OutputFiles(MeanLDFile);
            for (auto& file : files)
            {
                if ( GP->CurrentGeneration == GP->startAtGeneration)
                {
                    this->WriteOutputs_T1_MeanLD_header(file);
                }
                if (file.isTime())
                {
                    #ifdef DEBUG
                    std::cout << "Write MeanLDFile\n";
                    #endif
                    this->WriteOutputs_T1_MeanLD(pop, file);
                }
            }
        }
    }


    


    /*
     ########################
     #### T1_LongestRun  ####
     ########################
     */
    //std::cout << "line 2814\n";
    if (this->isFile(LongestRunFile))
    {
        if (SSP->Gmap.T1_nbLoci)
        {
            std::vector<OutputFile>& files = this->get_OutputFiles(LongestRunFile);
            for (auto& file : files)
            {
                if ( GP->CurrentGeneration == GP->startAtGeneration)
                {
                    this->WriteOutputs_T1_LongestRun_header(file);
                }
                if (file.isTime())
                {
                    #ifdef DEBUG
                    std::cout << "Write LongestRunFile\n";
                    #endif
                    this->WriteOutputs_T1_LongestRun(pop, file);
                }
            }
        }
    }


    /*
     ##########################
     #### T1_LargeOutput  #####
     ##########################
     */
    //std::cout << "line 2767\n";
    if (this->isFile(T1_LargeOutputFile))
    {
        if (SSP->Gmap.T1_nbLoci)
        {
            std::vector<OutputFile>& files = this->get_OutputFiles(T1_LargeOutputFile);
            for (auto& file : files)
            {
                if ( GP->CurrentGeneration == GP->startAtGeneration)
                {
                    this->WriteOutputs_T1_LargeOutput_header(file);
                }
                if (file.isTime())
                {
                    #ifdef DEBUG
                    std::cout << "Write T1_LargeOutputFile\n";
                    #endif
                    this->WriteOutputs_T1_LargeOutput(pop, file);
                }
            }
        }
    }


    /*
     ##########################
     #### T1 Hybrid Index  ####
     ##########################
     */
    //std::cout << "line 2925\n";
    if (this->isFile(HybridIndexFile))
    {
        if (SSP->Gmap.T1_nbLoci)
        {
            std::vector<OutputFile>& files = this->get_OutputFiles(HybridIndexFile);
            for (auto& file : files)
            {
                if ( GP->CurrentGeneration == GP->startAtGeneration)
                {
                    this->WriteOutputs_T1_HybridIndex_header(file);
                }
                if (file.isTime())
                {
                    #ifdef DEBUG
                    std::cout << "Write HybridIndexFile\n";
                    #endif
                    this->WriteOutputs_T1_HybridIndex(pop, file);
                }
            }
        }
    }


    /*
     ##################################
     #### T1 Average Hybrid Index  ####
     ##################################
     */
    //std::cout << "line 2925\n";
    if (this->isFile(T1_AverageHybridIndexFile))
    {
        if (SSP->Gmap.T1_nbLoci)
        {
            std::vector<OutputFile>& files = this->get_OutputFiles(T1_AverageHybridIndexFile);
            for (auto& file : files)
            {
                if ( GP->CurrentGeneration == GP->startAtGeneration)
                {
                    this->WriteOutputs_T1_AverageHybridIndex_header(file);
                }
                if (file.isTime())
                {
                    #ifdef DEBUG
                    std::cout << "Write T1_AverageHybridIndexFile\n";
                    #endif
                    this->WriteOutputs_T1_AverageHybridIndex(pop, file);
                }
            }
        }
    }

    /*
     ##########################
     #### T1 ExpectiMinRec ####
     ##########################
     */
    //std::cout << "line 2952\n";
    if (this->isFile(ExpectiMinRecFile))
    {
        if (SSP->Gmap.T1_nbLoci)
        {
            std::vector<OutputFile>& files = this->get_OutputFiles(ExpectiMinRecFile);
            for (auto& file : files)
            {
                if ( GP->CurrentGeneration == GP->startAtGeneration)
                {
                    this->WriteOutputs_T1_ExpectiMinRec_header(file);
                }
                if (file.isTime())
                {
                    #ifdef DEBUG
                    std::cout << "Write ExpectiMinRecFile\n";
                    #endif
                    this->WriteOutputs_T1_ExpectiMinRec(pop, file);
                }
            }
        }
    }

    /*
     ################
     #### T1 vcf ####
     ################
     */
    //std::cout << "line 2981\n";
    if (this->isFile(T1_vcfFile))
    {
        if (SSP->Gmap.T1_nbLoci)
        {
            std::vector<OutputFile>& files = this->get_OutputFiles(T1_vcfFile);
            for (auto& file : files)
            {
                if (file.isTime())
                {
                    #ifdef DEBUG
                    std::cout << "Write T1_vcfFile \n";
                    #endif
                    this->WriteOutputs_T1_vcf(pop, file);
                }
            }
        }
    }



     /*
     ##########################
     ### T1_haplotypeFreqs  ###
     ##########################
     */
    
    if (this->isFile(T1_haplotypeFreqs_file))
    {
        if (SSP->Gmap.T1_nbLoci)
        {
            std::vector<OutputFile>& files = this->get_OutputFiles(T1_haplotypeFreqs_file);
            for (auto& file : files)
            {
                if (file.isTime())
                {
                    #ifdef DEBUG
                    std::cout << "Write T1_haplotypeFreqs_file\n";
                    #endif
                    this->WriteOutputs_T1_haplotypeFreqs(pop, file);
                }
            }
        }
    }


    /*
     ################
     #### T1 FST ####
     ################
     */
    //std::cout << "line 3114\n";
    if (this->isFile(T1_FST))
    {
        if (SSP->Gmap.T1_nbLoci)
        {
            std::vector<OutputFile>& files = this->get_OutputFiles(T1_FST);
            for (auto& file : files)
            {
                if ( GP->CurrentGeneration == GP->startAtGeneration)
                {
                    this->WriteOutputs_T1_FST_header(file);
                }
                if (file.isTime())
                {
                    #ifdef DEBUG
                    std::cout << "Write T1_FST\n";
                    #endif
                    this->WriteOutputs_T1_FST(pop, file);
                }
            }
        }
    }

    /*
     ###############
     #### T1SFS ####
     ###############
     */
    if (this->isFile(T1_SFS_file))
    {
        if (SSP->Gmap.T1_nbLoci)
        {
            std::vector<OutputFile>& files = this->get_OutputFiles(T1_SFS_file);
            for (auto& file : files)
            {
                if ( GP->CurrentGeneration == GP->startAtGeneration)
                {
                    this->WriteOutputs_T1SFS_header(file);
                }
                if (file.isTime())
                {
                    #ifdef DEBUG
                    std::cout << "Write T1SFS\n";
                    #endif
                    this->WriteOutputs_T1SFS(pop, file);
                }
            }
        }
    }
}

void OutputWriter::WriteOutputs_forDefinedPop_T2(Pop& pop)
{
    /*
     #########################
     #### T2_LargeOutput  ####
     #########################
     */
    //std::cout << "line 2842\n";
    if (this->isFile(T2_LargeOutputFile))
    {
        if (SSP->Gmap.T2_nbLoci)
        {
            OutputFile& file = this->get_OutputFiles(T2_LargeOutputFile)[0];
            if ( GP->CurrentGeneration == GP->startAtGeneration)
            {
                this->WriteOutputs_T2_LargeOutput_header(file);
            }
            if (file.isTime())
            {
                #ifdef DEBUG
                std::cout << "Write T2_LargeOutputFile\n";
                #endif
                this->WriteOutputs_T2_LargeOutput(pop, file);
            }
        }
    }
}

void OutputWriter::WriteOutputs_forDefinedPop_T3(Pop& pop)
{
     /*
     ##########################
     #### T3 Large Outputs ####
     ##########################
     */
    //std::cout << "line 3032\n";
    if (this->isFile(T3_LargeOutputFile))
    {
        if (SSP->Gmap.T3_nbLoci)
        {
            std::vector<OutputFile>& files = this->get_OutputFiles(T3_LargeOutputFile);
            for (auto& file : files)
            {
                if ( GP->CurrentGeneration == GP->startAtGeneration)
                {
                    this->WriteOutputs_T3_LargeOutput_header(file);
                }
                if (file.isTime())
                {
                    #ifdef DEBUG
                    std::cout << "Write T3_LargeOutputFile\n";
                    #endif
                    this->WriteOutputs_T3_LargeOutput(pop, file);
                }
            }
        }
    }

    /*
     ####################
     #### T3 VarMean ####
     ####################
     */
    //std::cout << "line 3036\n";
    if (this->isFile(T3_MeanVarFile))
    {
        if (SSP->Gmap.T3_nbLoci)
        {
            std::vector<OutputFile>& files = this->get_OutputFiles(T3_MeanVarFile);
            for (auto& file : files)
            {
                if ( GP->CurrentGeneration == GP->startAtGeneration)
                {
                    this->WriteOutputs_T3_MeanVar_header(file);
                }
                if (file.isTime())
                {
                    #ifdef DEBUG
                    std::cout << "Write T3_MeanVarFile\n";
                    #endif
                    this->WriteOutputs_T3_MeanVar(pop, file);
                }
            }
        }
    }
}

void OutputWriter::WriteOutputs_forDefinedPop_T4AndSampledSeq(Pop& pop)
{
    //std::cout << "Enters 'OutputWriter::WriteOutputs_forDefinedPop_T4AndSampledSeq'\n";
    /*
        ---------------------
        - ################# -
        - #### Headers #### -
        - ################# -
        ---------------------
    */

    if (SSP->whenDidExtinctionOccur == -1)
    {
        if ( GP->CurrentGeneration == GP->startAtGeneration)
            {
                /*
                 ########################################
                 #### T4 painted haplotype diversity ####
                 ########################################
                */
                
                if (this->isFile(T4_paintedHaploSegmentsDiversity_file))
                {
                    if (SSP->Gmap.T4_nbLoci)
                    {
                        std::vector<OutputFile>& files = this->get_OutputFiles(T4_paintedHaploSegmentsDiversity_file);
                        for (auto& file : files)
                        {
                            this->WriteOutputs_T4_paintedHaploSegmentsDiversity_file_header(file);
                        }
                    }
                }

                for (size_t mutPlacingIndex = 0 ; mutPlacingIndex < SSP->T4_nbMutationPlacingsPerOutput ; ++mutPlacingIndex)
                {
                    /*
                     ##########################
                     #### T4_LargeOutput  #####
                     ##########################
                     */
                    //std::cout << "line 2896\n";
                    if (this->isFile(T4_LargeOutputFile))
                    {
                        if (SSP->Gmap.T4_nbLoci)
                        {
                            std::vector<OutputFile>& files = this->get_OutputFiles(T4_LargeOutputFile);
                            for (auto& file : files)
                            {
                                this->WriteOutputs_T4_LargeOutput_header(file, mutPlacingIndex);
                            }
                        }
                    }

                    /*
                     ##############################
                     #### T4 frequency of SNPs ####
                     ##############################
                     */
                    
                    if (this->isFile(T4_SNPfreq_file))
                    {
                        if (SSP->Gmap.T4_nbLoci)
                        {
                            std::vector<OutputFile>& files = this->get_OutputFiles(T4_SNPfreq_file);
                            for (auto& file : files)
                            {
                                this->WriteOutputs_T4_SNPfreq_file_header(file, mutPlacingIndex);
                            }
                        }
                    }



                    /*
                     ###############
                     #### T4SFS ####
                     ###############
                     */
                    if (this->isFile(T4_SFS_file))
                    {
                        if (SSP->Gmap.T4_nbLoci)
                        {
                            std::vector<OutputFile>& files = this->get_OutputFiles(T4_SFS_file);
                            for (auto& file : files)
                            {
                                this->WriteOutputs_T4SFS_header(file, mutPlacingIndex);
                            }
                        }
                    }

                    /*
                     ##############################
                     #### Tx frequency of SNPs ####
                     ##############################
                     */
                    
                    if (this->isFile(Tx_SNPfreq_file))
                    {
                        std::vector<OutputFile>& files = this->get_OutputFiles(Tx_SNPfreq_file);
                        for (auto& file : files)
                        {
                            //std::cout << "About to enter 'WriteOutputs_Tx_SNPfreq_file_header'\n";
                            this->WriteOutputs_Tx_SNPfreq_file_header(file, mutPlacingIndex);
                        }
                    }
                }
            }



        /*
            ---------------------------------------------------------------------------------------
            - ################################################################################### -
            - #### Figure if simplification or simplification and mutation placing necessary #### -
            - ################################################################################### -
            ---------------------------------------------------------------------------------------
        */


        bool isAnyOutputWithoutMutPlacing = false;
        bool isAnyOutputWithMutPlacing = false;
        bool isAnyOutputWithMutPlacing_needFreqs = false;
        bool isAnyOutputWithMutPlacing_needAllData = false;

        if (this->isFile(T4_LargeOutputFile))
            for (auto& file : this->get_OutputFiles(T4_LargeOutputFile))
                if (file.isTime())
                {
                    isAnyOutputWithMutPlacing = true;
                    isAnyOutputWithMutPlacing_needAllData = true;
                }

        if (this->isFile(T4_vcfFile))
            for (auto& file : this->get_OutputFiles(T4_vcfFile))
                if (file.isTime())
                {
                    isAnyOutputWithMutPlacing = true;
                    isAnyOutputWithMutPlacing_needAllData = true;
                }

        if (this->isFile(T4_SNPfreq_file))
            for (auto& file : this->get_OutputFiles(T4_SNPfreq_file))
                if (file.isTime())
                {
                    isAnyOutputWithMutPlacing = true;
                    isAnyOutputWithMutPlacing_needFreqs = true;
                }

        if (this->isFile(T4_SFS_file))
            for (auto& file : this->get_OutputFiles(T4_SFS_file))
                if (file.isTime())
                {
                    isAnyOutputWithMutPlacing = true;
                    isAnyOutputWithMutPlacing_needAllData = true;
                }

        if (this->isFile(Tx_SNPfreq_file))
            for (auto& file : this->get_OutputFiles(Tx_SNPfreq_file))
                if (file.isTime())
                {
                    isAnyOutputWithMutPlacing = true;
                    isAnyOutputWithMutPlacing_needFreqs = true;
                }

        if (this->isFile(T4_paintedHaplo_file))
            for (auto& file : this->get_OutputFiles(T4_paintedHaplo_file))
                if (file.isTime())
                {
                    isAnyOutputWithoutMutPlacing = true;
                }


        if (this->isFile(T4_paintedHaploSegmentsDiversity_file))
            for (auto& file : this->get_OutputFiles(T4_paintedHaploSegmentsDiversity_file))
                if (file.isTime())
                {
                    isAnyOutputWithoutMutPlacing = true;
                }


        // If not output, I can exit

        //std::cout << "line 5356: GP->CurrentGeneration = " << GP->CurrentGeneration << "\n";
        //std::cout << "isAnyOutputWithoutMutPlacing = " << isAnyOutputWithoutMutPlacing << "\n";
        //std::cout << "isAnyOutputWithMutPlacing = " << isAnyOutputWitMutPlacing << "\n";

        if (!(isAnyOutputWithoutMutPlacing || isAnyOutputWithMutPlacing))
            return;
        

        //std::cout << "About to simplify (SSP->Gmap.T4_nbLoci = "<< SSP->Gmap.T4_nbLoci <<")\n";
        if (SSP->Gmap.T4_nbLoci) SSP->T4Tree.simplify(pop);

        /*
            ------------------------------------------------------
            - ################################################## -
            - #### Outputs that don't need mutation placing #### -
            - ################################################## -
            ------------------------------------------------------
        */
        if (isAnyOutputWithoutMutPlacing)
        {
            /*
             ##############################
             #### T4 painted haplotype ####
             ##############################
             */
            
            if (this->isFile(T4_paintedHaplo_file))
            {
                if (SSP->Gmap.T4_nbLoci)
                {
                    std::vector<OutputFile>& files = this->get_OutputFiles(T4_paintedHaplo_file);
                    for (auto& file : files)
                    {
                        if (file.isTime())
                        {
                            #ifdef DEBUG
                            std::cout << "Write T4_paintedHaplo_file \n";
                            #endif
                            this->WriteOutputs_T4_paintedHaplo_file(pop, file);
                        }
                    }
                }
            }


             /*
             ########################################
             #### T4 painted haplotype diversity ####
             ########################################
             */
            
            if (this->isFile(T4_paintedHaploSegmentsDiversity_file))
            {
                if (SSP->Gmap.T4_nbLoci)
                {
                    std::vector<OutputFile>& files = this->get_OutputFiles(T4_paintedHaploSegmentsDiversity_file);
                    for (auto& file : files)
                    {
                        if (file.isTime())
                        {
                            #ifdef DEBUG
                            std::cout << "Write T4_paintedHaploSegmentsDiversity_file \n";
                            #endif
                            this->WriteOutputs_T4_paintedHaploSegmentsDiversity_file(pop, file);
                        }
                    }
                }
            }

            /*
             ##########################
             #### T4CoalescenceFst ####
             ##########################
             */
            /*if (this->isFile(T4CoalescenceFst))
            {
                if (SSP->Gmap.T4_nbLoci)
                {
                    std::vector<OutputFile>& files = this->get_OutputFiles(T4CoalescenceFst);
                    for (auto& file : files)
                    {
                        if ( GP->CurrentGeneration == GP->startAtGeneration)
                        {
                            this->WriteOutputs_T4CoalescenceFst_header(file);
                        }
                        if (file.isTime())
                        {
                            #ifdef DEBUG
                            std::cout << "Write T4CoalescenceFst\n";
                            #endif
                            this->WriteOutputs_T4CoalescenceFst(file); // No need to give Pop as T4Tree is in SSP
                        }
                    }
                }
            }*/
        }




        /*
            ------------------------------------------------
            - ############################################ -
            - #### Outputs that need mutation placing #### -
            - ############################################ -
            ------------------------------------------------
        */

        if (isAnyOutputWithMutPlacing)
        {
            auto rng_save = GP->rngw.getRNGState();
            assert(SSP->T4_nbMutationPlacingsPerOutput > 0);
            for (size_t mutPlacingIndex = 0 ; mutPlacingIndex < SSP->T4_nbMutationPlacingsPerOutput ; ++mutPlacingIndex)
            {
                /*
                {
                    std::cout << "sleep for 10 seconds before placing mutations" << std::endl;
                    //std::chrono::seconds dura(10);
                    //std::this_thread::sleep_for( dura );
                    std::cout << "end of sleeping" << std::endl;
                }
                */
                

                std::pair< std::vector<std::vector<std::vector<uint32_t>>>, std::vector<T4TreeRec::PatchCountData> > allData;

                if (SSP->Gmap.T4_nbLoci)
                {
                    unsigned char placeMutationOutputsType; 
                    assert(isAnyOutputWithMutPlacing_needFreqs || isAnyOutputWithMutPlacing_needAllData);

                    if (isAnyOutputWithMutPlacing_needFreqs && isAnyOutputWithMutPlacing_needAllData) placeMutationOutputsType = 'b';
                    if (isAnyOutputWithMutPlacing_needFreqs && !isAnyOutputWithMutPlacing_needAllData) placeMutationOutputsType = 'f';
                    if (!isAnyOutputWithMutPlacing_needFreqs && isAnyOutputWithMutPlacing_needAllData) placeMutationOutputsType = 'r';
                    
                    if (GP->CurrentGeneration == GP->nbGenerations && mutPlacingIndex == SSP->T4_nbMutationPlacingsPerOutput - 1)
                    {
                        //std::cout << "placeMutations: Completely Delete tree\n";
                        allData = SSP->T4Tree.placeMutations(pop, false, true, placeMutationOutputsType);
                    } else
                    {
                        if (SSP->T4_nbMutationPlacingsPerOutput == 1 && SSP->T4_respectPreviousMutationsOutputs)
                        {
                            // delete tree
                            //std::cout << "placeMutations: Partially delete tree\n";
                            allData = SSP->T4Tree.placeMutations(pop, true, false, placeMutationOutputsType);
                        } else
                        {
                            // keep the tree
                            //std::cout << "placeMutations: Keep tree\n";
                            assert(!SSP->T4_respectPreviousMutationsOutputs);
                            allData = SSP->T4Tree.placeMutations(pop, false, false, placeMutationOutputsType);
                        }
                    }
                    //std::cout << "Mutations placed\n";
                }

                auto& data = allData.first;
                auto& T4freqs = allData.second;

                /*
                {
                    std::cout << "sleep for 10 seconds after placing mutations" << std::endl;
                    //std::chrono::seconds dura(10);
                    //std::this_thread::sleep_for( dura );
                    std::cout << "end of sleeping" << std::endl;
                }
                */
                

                /*
                 #######################
                 ### Sample sequence ###
                 #######################
                 */
                
                if (this->isFile(sampleSeq_file))
                {
                    OutputFile& file = this->get_OutputFiles(sampleSeq_file)[0];
                    if (file.isTime())
                    {
                        #ifdef DEBUG
                        std::cout << "Write sampleSeq\n";
                        #endif
                        this->WriteOutputs_sampleSequences(data, pop, file, mutPlacingIndex);
                    }
                }


                /*
                 ##########################
                 #### T4_LargeOutput  #####
                 ##########################
                 */
                //std::cout << "line 2896\n";
                if (this->isFile(T4_LargeOutputFile))
                {
                    if (SSP->Gmap.T4_nbLoci)
                    {
                        std::vector<OutputFile>& files = this->get_OutputFiles(T4_LargeOutputFile);
                        for (auto& file : files)
                        {
                            if (file.isTime())
                            {
                                #ifdef DEBUG
                                std::cout << "Write T4_LargeOutputFile\n";
                                #endif
                                this->WriteOutputs_T4_LargeOutput(file, data, mutPlacingIndex);
                            }
                        }
                    }
                }

                /*
                 ################
                 #### T4 vcf ####
                 ################
                 */
                //std::cout << "line 3006\n";
                if (this->isFile(T4_vcfFile))
                {
                    if (SSP->Gmap.T4_nbLoci)
                    {
                        std::vector<OutputFile>& files = this->get_OutputFiles(T4_vcfFile);
                        for (auto& file : files)
                        {
                            if (file.isTime())
                            {
                                #ifdef DEBUG
                                std::cout << "Write T4_vcfFile \n";
                                #endif
                                this->WriteOutputs_T4_vcf(file, data, mutPlacingIndex);
                            }
                        }
                    }
                }


                  /*
                 ##############################
                 #### T4 frequency of SNPs ####
                 ##############################
                 */
                
                if (this->isFile(T4_SNPfreq_file))
                {
                    if (SSP->Gmap.T4_nbLoci)
                    {
                        std::vector<OutputFile>& files = this->get_OutputFiles(T4_SNPfreq_file);
                        for (auto& file : files)
                        {
                            if (file.isTime())
                            {
                                #ifdef DEBUG
                                std::cout << "Write T4_SNPfreq_file \n";
                                #endif
                                this->WriteOutputs_T4_SNPfreq_file(file, T4freqs, mutPlacingIndex);
                            }
                        }
                    }
                }



                /*
                 ###############
                 #### T4SFS ####
                 ###############
                 */
                if (this->isFile(T4_SFS_file))
                {
                    if (SSP->Gmap.T4_nbLoci)
                    {
                        std::vector<OutputFile>& files = this->get_OutputFiles(T4_SFS_file);
                        for (auto& file : files)
                        {
                            if (file.isTime())
                            {
                                #ifdef DEBUG
                                std::cout << "Write T4SFS\n";
                                #endif
                                this->WriteOutputs_T4SFS(file, data, mutPlacingIndex);
                            }
                        }
                    }
                }


                  /*
                 ##############################
                 #### Tx frequency of SNPs ####
                 ##############################
                 */
                
                if (this->isFile(Tx_SNPfreq_file))
                {
                    std::vector<OutputFile>& files = this->get_OutputFiles(Tx_SNPfreq_file);
                    for (auto& file : files)
                    {
                        if (file.isTime())
                        {
                            #ifdef DEBUG
                            std::cout << "Write Tx_SNPfreq_file \n";
                            #endif
                            this->WriteOutputs_Tx_SNPfreq_file(file, T4freqs, mutPlacingIndex, pop);
                        }
                    }
                }

            
                /* WATCH OUT:
                    WriteOutputs_Tx_SNPfreq_file EMPTIES T4freqs FOR PERFORMANCE REASONS. DO NOT PLACE ANY OUTPUT AFTER WriteOutputs_Tx_SNPfreq_file()

                */
            }
            GP->rngw.setRNG(rng_save);
        }
    }
}

void OutputWriter::WriteOutputs_forDefinedPop_T56(Pop& pop)
{
    /*
     ###########################
     ### T5_Allele Frequency ###
     ###########################
     */
    //std::cout << "line 2756\n";
    if (this->isFile(T56_AlleleFreqFile))
    {
        if (SSP->Gmap.T56_nbLoci)
        {
            std::vector<OutputFile>& files = this->get_OutputFiles(T56_AlleleFreqFile);
            for (auto& file : files)
            {
                if ( GP->CurrentGeneration == GP->startAtGeneration)
                {
                    this->WriteOutputs_T56_AlleleFreq_header(file);
                }
                if (file.isTime())
                {
                    #ifdef DEBUG
                    std::cout << "Write T5_AlleleFreqFile\n";
                    #endif
                    this->WriteOutputs_T56_AlleleFreq(pop, file);
                }
            }
        }
    }

    /*
     ##########################
     #### T5_LargeOutput  #####
     ##########################
     */
    if (this->isFile(T56_LargeOutputFile))
    {
        if (SSP->Gmap.T56_nbLoci)
        {
            std::vector<OutputFile>& files = this->get_OutputFiles(T56_LargeOutputFile);
            for (auto& file : files)
            {
                if ( GP->CurrentGeneration == GP->startAtGeneration)
                {
                    this->WriteOutputs_T56_LargeOutput_header(file);
                }
                if (file.isTime())
                {
                    #ifdef DEBUG
                    std::cout << "Write T5_LargeOutputFile\n";
                    #endif
                    this->WriteOutputs_T56_LargeOutput(pop, file);
                }
            }
        }
    }


    /*
     ################
     #### T5 vcf ####
     ################
     */
    //std::cout << "line 2981\n";
    if (this->isFile(T56_vcfFile))
    {
        if (SSP->Gmap.T56_nbLoci)
        {
            std::vector<OutputFile>& files = this->get_OutputFiles(T56_vcfFile);
            for (auto& file : files)
            {
                if (file.isTime())
                {
                    #ifdef DEBUG
                    std::cout << "Write T5_vcfFile \n";
                    #endif
                    this->WriteOutputs_T56_vcf(pop, file);
                }
            }
        }
    }

    /*
     ###############
     #### T5SFS ####
     ###############
     */
    if (this->isFile(T56_SFS_file))
    {
        if (SSP->Gmap.T56_nbLoci)
        {
            std::vector<OutputFile>& files = this->get_OutputFiles(T56_SFS_file);
            for (auto& file : files)
            {
                if ( GP->CurrentGeneration == GP->startAtGeneration)
                {
                    this->WriteOutputs_T56SFS_header(file);
                }
                if (file.isTime())
                {
                    #ifdef DEBUG
                    std::cout << "Write T5SFS\n";
                    #endif
                    this->WriteOutputs_T56SFS(pop,file);
                }
            }
        }
    }
}


void OutputWriter::WriteOutputs_forDefinedPop(Pop& pop)
{
#ifdef CALLENTRANCEFUNCTIONS
std::cout << "Enters in 'WriteOutputs_forDefinedPop'\n";
#endif 

    std::cout.precision(8);


    WriteOutputs_forDefinedPop_general(pop);
    WriteOutputs_forDefinedPop_T56(pop);
    if (GP->CurrentGeneration == GP->nbGenerations)
    {
        pop.freeT56Memory();
    }
    WriteOutputs_forDefinedPop_T1(pop);
    WriteOutputs_forDefinedPop_T2(pop);
    WriteOutputs_forDefinedPop_T3(pop);
    WriteOutputs_forDefinedPop_T4AndSampledSeq(pop);

}


bool OutputWriter::shouldNABePrinted(int patch_index)
{
    return patch_index >= GP->PatchNumber;
}

bool OutputWriter::shouldNABePrinted(int patch_index, int ind_index)
{
    if (shouldNABePrinted(patch_index))
    {
        return true;
    } else
    {
        assert(SSP->patchSize.size() > patch_index);
        if (ind_index >= SSP->patchSize[patch_index])
        {
            return true;
        }
    }
    return false;
}
