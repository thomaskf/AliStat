/*
 *
 * UserOptions.cpp
 * AliStat
 * 
 * CSIRO Open Source Software License Agreement (variation of the BSD / MIT License)
 * Copyright (c) 2014, Commonwealth Scientific and Industrial Research Organisation (CSIRO) ABN 41 687 119 230.
 * 
 * All rights reserved. CSIRO is willing to grant you a license to this AliStat on the following terms, except 
 * where otherwise indicated for third party material.
 * 
 * Redistribution and use of this software in source and binary forms, with or without modification, are
 * permitted provided that the following conditions are met:
 * 
 * •	Redistributions of source code must retain the above copyright notice, this list of conditions and the
 *      following disclaimer.
 * •	Redistributions in binary form must reproduce the above copyright notice, this list of conditions and
 *      the following disclaimer in the documentation and/or other materials provided with the distribution.
 * •	Neither the name of CSIRO nor the names of its contributors may be used to endorse or promote products
 *      derived from this software without specific prior written permission of CSIRO.
 * 
 * EXCEPT AS EXPRESSLY STATED IN THIS AGREEMENT AND TO THE FULL EXTENT PERMITTED BY APPLICABLE LAW, THE SOFTWARE
 * IS PROVIDED "AS-IS". CSIRO MAKES NO REPRESENTATIONS, WARRANTIES OR CONDITIONS OF ANY KIND, EXPRESS OR IMPLIED,
 * INCLUDING BUT NOT LIMITED TO ANY REPRESENTATIONS, WARRANTIES OR CONDITIONS REGARDING THE CONTENTS OR ACCURACY
 * OF THE SOFTWARE, OR OF TITLE, MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, NON-INFRINGEMENT, THE ABSENCE
 * OF LATENT OR OTHER DEFECTS, OR THE PRESENCE OR ABSENCE OF ERRORS, WHETHER OR NOT DISCOVERABLE.
 * 
 * TO THE FULL EXTENT PERMITTED BY APPLICABLE LAW, IN NO EVENT SHALL CSIRO BE LIABLE ON ANY LEGAL THEORY (INCLUDING,
 * WITHOUT LIMITATION, IN AN ACTION FOR BREACH OF CONTRACT, NEGLIGENCE OR OTHERWISE) FOR ANY CLAIM, LOSS, DAMAGES
 * OR OTHER LIABILITY HOWSOEVER INCURRED.  WITHOUT LIMITING THE SCOPE OF THE PREVIOUS SENTENCE THE EXCLUSION OF
 * LIABILITY SHALL INCLUDE: LOSS OF PRODUCTION OR OPERATION TIME, LOSS, DAMAGE OR CORRUPTION OF DATA OR RECORDS;
 * OR LOSS OF ANTICIPATED SAVINGS, OPPORTUNITY, REVENUE, PROFIT OR GOODWILL, OR OTHER ECONOMIC LOSS; OR ANY SPECIAL,
 * INCIDENTAL, INDIRECT, CONSEQUENTIAL, PUNITIVE OR EXEMPLARY DAMAGES, ARISING OUT OF OR IN CONNECTION WITH THIS
 * AGREEMENT, ACCESS OF THE SOFTWARE OR ANY OTHER DEALINGS WITH THE SOFTWARE, EVEN IF CSIRO HAS BEEN ADVISED OF
 * THE POSSIBILITY OF SUCH CLAIM, LOSS, DAMAGES OR OTHER LIABILITY.
 * 
 * APPLICABLE LEGISLATION SUCH AS THE AUSTRALIAN CONSUMER LAW MAY APPLY REPRESENTATIONS, WARRANTIES, OR CONDITIONS,
 * OR IMPOSES OBLIGATIONS OR LIABILITY ON CSIRO THAT CANNOT BE EXCLUDED, RESTRICTED OR MODIFIED TO THE FULL EXTENT
 * SET OUT IN THE EXPRESS TERMS OF THIS CLAUSE ABOVE "CONSUMER GUARANTEES".  TO THE EXTENT THAT SUCH CONSUMER
 * GUARANTEES CONTINUE TO APPLY, THEN TO THE FULL EXTENT PERMITTED BY THE APPLICABLE LEGISLATION, THE LIABILITY
 * OF CSIRO UNDER THE RELEVANT CONSUMER GUARANTEE IS LIMITED (WHERE PERMITTED AT CSIRO’S OPTION) TO ONE OF
 * FOLLOWING REMEDIES OR SUBSTANTIALLY EQUIVALENT REMEDIES:
 * 
 * (a) THE REPLACEMENT OF THE SOFTWARE, THE SUPPLY OF EQUIVALENT SOFTWARE, OR SUPPLYING RELEVANT SERVICES AGAIN;
 * (b) THE REPAIR OF THE SOFTWARE;
 * (c) THE PAYMENT OF THE COST OF REPLACING THE SOFTWARE, OF ACQUIRING EQUIVALENT SOFTWARE, HAVING THE RELEVANT
 *     SERVICES SUPPLIED AGAIN, OR HAVING THE SOFTWARE REPAIRED.
 * 
 * IN THIS CLAUSE, CSIRO INCLUDES ANY THIRD PARTY AUTHOR OR OWNER OF ANY PART OF THE SOFTWARE OR MATERIAL
 * DISTRIBUTED WITH IT.  CSIRO MAY ENFORCE ANY RIGHTS ON BEHALF OF THE RELEVANT THIRD PARTY.
 *
 */


#include "UserOptions.h"

// read the arguments and collect the user options
int readArguments(int argc, char** argv, UserOptions* user_options, string* errMsg) {
    
    // return:
    //    0 if the user options are valid
    //    1 if there exists error message
    //    2 if a help menu is called
    
    
    int i,j;
    bool duplicateOption = false;
    bool emptyOption = false;
	bool used_option_m = false;
    vector<string> token;
    
    *errMsg = "";
    
    if (argc < 3) {
        // invoke the help menu
        return 2;
    }
    
    // set the user options as initial values
    user_options->reset();
    user_options->alignFile = argv[1];
    user_options->dataType = atoi(argv[2]);
    
    // check the value of dataType
    if (user_options->dataType < 1 || user_options->dataType > 3) {
        *errMsg = "Error! <data type> has to be 1, 2 or 3.";
        return 1;
    }
    
    
    // get the other options if there is any
    if (argc > 3) {
        GetOptions options;
        options.read(argc, argv);
        
        for (i=0; i<options.size(); i++) {
            
            char flag = options.flags[i];
            string value = options.values[i];
            
            switch (flag) {
                    
                case 'c':
                    if (user_options->codeType!=-1)
                        duplicateOption = true;
                    else if (value.length()==0)
                        emptyOption = true;
                    else if (user_options->dataType!=1)
                        *errMsg = "Error! The option '-c' is designed for data type = 1 only";
                    else {
                        user_options->codeType = atoi(value.c_str());
                        if (user_options->codeType < 0 || user_options->codeType > 13) {
                            *errMsg = "Error! <coding_type> has to be between 1 and 13.";
                        }
                    }
                    break;
                    
                case 'o':
                    if (user_options->prefixOut!="")
                        duplicateOption = true;
                    else if (value.length()==0)
                        emptyOption = true;
                    else
                        user_options->prefixOut = value;
                    break;
                    
                case 'n':
                    if (user_options->seqNameFile!="")
                        duplicateOption = true;
                    else if (value.length()==0)
                        emptyOption = true;
                    else
                        user_options->seqNameFile = value;
                    break;
                    
                case 'p':
                    if (user_options->partitionFile!="")
                        duplicateOption = true;
                    else if (value.length()==0)
                        emptyOption = true;
                    else if (user_options->slideWindowSize!=-1 ||
                             user_options->slideWindowStep!=-1)
                        *errMsg = "Error! The options '-p' and '-s' cannot be used at the same time";
                    else
                        user_options->partitionFile = value;
                    break;
                    
                case 's':
                    if (user_options->slideWindowSize!=-1 ||
                        user_options->slideWindowStep!=-1)
                        duplicateOption = true;
                    else if (value.length()==0)
                        emptyOption = true;
                    else if (user_options->partitionFile!="")
                        *errMsg = "Error! The options '-p' and '-s' cannot be used at the same time";
                    else {
                        tokenizer(value, ",", &token);
                        if (token.size()!=2) {
                            *errMsg = "Error! Regarding the option '-s', there should be 2 integers";
                        } else {
                            user_options->slideWindowSize = atoi(token[0].c_str());
                            user_options->slideWindowStep = atoi(token[1].c_str());
                            if (user_options->slideWindowSize <= 0)
                                *errMsg = "Error! Sliding window size should be > 0";
                            else if (user_options->slideWindowStep <= 0)
                                *errMsg = "Error! Sliding window step size should be > 0";
                        }
                    }
                    break;
                    
                case 't':
                    int tableSelected[TOTAL_TABLES_NUM];
                    if (user_options->outTables.size() > 0)
                        duplicateOption = true;
                    else if (value.length()==0) {
                        for (j=0; j<TOTAL_TABLES_NUM; j++)
                            tableSelected[j] = 1;
                    } else {
                        for (j=0; j<TOTAL_TABLES_NUM; j++)
                            tableSelected[j] = 0;
                        tokenizer(value, ",", &token);
                        for (j=0; j<token.size(); j++) {
                            int tableID = atoi(token[j].c_str());
                            if (tableID >= 1 && tableID <= TOTAL_TABLES_NUM) {
                                if (!tableSelected[tableID-1]) {
                                    tableSelected[tableID-1] = 1;
                                } else {
                                    *errMsg = "Error! Regarding the option '-t', there are duplicate table ID";
                                    break;
                                }
                            } else {
                                *errMsg = "Error! Regarding the option '-t', some table IDs are not valid";
                                break;
                            }
                        }
                    }
                    if (*errMsg == "") {
                        for (j=0; j<TOTAL_TABLES_NUM; j++) {
                            if (tableSelected[j]) {
                                user_options->outTables.push_back(j+1);
                            }
                        }
                    }
                    break;
                    
                case 'r':
                    if (user_options->reorder!=-1)
                        duplicateOption = true;
                    else if (value.length()==0)
                        emptyOption = true;
                    else if (value == "row")
                        user_options->reorder = 1;
                    else if (value == "col")
                        user_options->reorder = 2;
                    else if (value == "both")
                        user_options->reorder = 3;
                    else
                        *errMsg = "Error! Regarding the option '-r', the value should be either 'row', 'col' or 'both'.";
                    break;
                
                case 'i':
                    if (user_options->makeHeatMap!=0)
                        duplicateOption = true;
                    else if (value.length()==0)
                        user_options->makeHeatMap = 3;
					else {
						user_options->makeHeatMap = atoi(value.c_str());
                        if (user_options->makeHeatMap < 1 || user_options->makeHeatMap > 3) {
                            *errMsg = "Error! Regarding the option '-i', the value should be 1 - 3.";
                        }
					}
                    break;

                case 'm':
                    if (used_option_m) {
                        duplicateOption = true;
                    } else {
					    used_option_m = true;
						user_options->outputAlign = 1;
					    if (value=="") {
                            user_options->maskThres = -1.0;
                        } else {
                            user_options->maskThres = atof(value.c_str());
                            if (user_options->maskThres < 0.0 || user_options->maskThres > 1.0) {
                                *errMsg = "Error! Regarding the option '-m', the value should be 0.0 - 1.0.";
                            }
                        }
					}
                    break;
                    
                case 'd':
                    if (user_options->computePDist > 0)
                        duplicateOption = true;
                    else
                        user_options->computePDist = 1;
                    break;

				case 'h':
                    return 2;
                    break;
                    
                default:
                    *errMsg = "Unknown option '-" + string(1,flag);
                    break;
                    
            } // case
            
            if (duplicateOption) {
                *errMsg = "Error! Duplicate option -" + string(1,flag);
            } else if (emptyOption) {
                *errMsg = "Error! Empty value for the option -" + string(1,flag);
            }
            
            if (*errMsg != "") {
                return 1;
            }
            
        }
        
        if (user_options->dataType==3 && user_options->partitionFile=="") {
            // if dataType is 3, then partition file has to be supplied
            *errMsg = "Error! Partifion file has to be supplied when <data type> = 3";
            return 1;
        }
    
    }
    
    // set the default value if necessary
    if (user_options->codeType==-1) {
        user_options->codeType = 0;
    }
    if (user_options->prefixOut == "") {
        user_options->setDefaultPrefixOut();
    }
    if (user_options->reorder==-1) {
        user_options->reorder = 0;
    }
    
    return 0;
}

// reset
void UserOptions::reset() {
    
    alignFile = "";
    dataType = 0;
    prefixOut = "";
    codeType = -1;
    seqNameFile = "";
    partitionFile = "";
    slideWindowSize = -1;
    slideWindowStep = -1;
    outTables.clear();
    reorder = -1;
    namesConsidered.clear();
    maskThres = -1.0;
    computePDist = 0;
	makeHeatMap = 0;
	outputAlign = 0;
}

// by default, prefixOut = <alignment file> w/o .ext
void UserOptions::setDefaultPrefixOut() {
    
    prefixOut = "";
    
    /*
    // become obsoleted on 13/05/2014
    size_t pos = alignFile.find_last_of(".");
    if (pos != string::npos && pos > 0) {
        // found
        prefixOut = alignFile.substr(0,(int)pos);
    } else {
        prefixOut = alignFile;
    }
     */
    
}

// by default, outTables = { 1, 2, 3, 4, 5, 6 }
void UserOptions::setDefaultOutTables() {
    
    int i;
    int totTables = 6;
    for (i=1; i<=totTables; i++)
        outTables.push_back(i);
    
}

// load the file including the sequence names considered
void UserOptions::loadSeqNameFile() {
    
    ifstream fin;
    fin.open(seqNameFile.c_str());
    string aline;
    while (getline(fin, aline)) {
        trimEnd(aline);
        if (aline.length() > 0) {
            namesConsidered.insert(aline);
        }
    }
    fin.close();
    
}

// load the partition file
void UserOptions::loadPartitionFile(int seqLen) {
    
    set<string> partitionNames;
    ifstream fin;
    fin.open(partitionFile.c_str());
    string aline;
    size_t pos;
    string partDataType;
    int partDataTypeID;
    int isCodon;
    string name;
    string posStr;
    int frPos, toPos;
    string errMsg = "";
    vector<string> token0, token1, token2;
    PartitionPos* newPart;
    int i,j;
    
    while (getline(fin,aline)) {
        trimEnd(aline);
        if (aline.length() > 0) {
            pos = aline.find_first_of('=');
            if (pos == string::npos) {
                // no "="
                errMsg = "no '='";
                break;
            } else {
                
                // get partition name
                // and the data type of the partition (if there is)
                
                name = "";
                partDataTypeID = 0;
                isCodon = 0;
                posStr = "";
                
                if (pos > 0) {
                    tokenizer(aline.substr(0,(int)pos), ",", &token0);
                    if (token0.size() > 1) {
                        partDataType = token0[0];
                        trimAll(partDataType);
                        if (partDataType == "NU")
                            partDataTypeID = 1;
                        else if (partDataType == "AA")
                            partDataTypeID = 2;
                        else if (partDataType.length() > 0) {
                            errMsg = "data type of the partition is unknown";
                            break;
                        }
                        name = token0[1];
                    } else {
                        name = token0[0];
                    }
                    trimAll(name);
                }
                
                // check whether the data type of the partition is correct
                if (this->dataType==1 && partDataTypeID==2) {
                    errMsg = "<data type> = 1, but amino acid partition exists";
                    break;
                } else if (this->dataType==2 && partDataTypeID==1) {
                    errMsg = "<data type> = 2, but nucleotide partition exists";
                    break;
                } else if (this->dataType==3 && partDataTypeID==0) {
                    errMsg = "<data type> = 3, but data type of the partition is unknown";
                    break;
                }
                
                // check whether the partition name is empty
                if (name.length() == 0) {
                    errMsg = "empty partition name";
                    break;
                }
                
                // check whether the same partition name appears before
                if (partitionNames.find(name)==partitionNames.end()) {
                    partitionNames.insert(name);
                } else {
                    errMsg = "duplicate partition name";
                    break;
                }
                
                // check whether it is a codon region
                if (pos+1 < aline.length()) {
                    posStr = aline.substr(pos+1);
                    removeAllSpaces(posStr);
                    
                    if (posStr.length() > 2 && posStr.substr(posStr.length()-2)=="\\3") {
                        isCodon = 1;
                        posStr.resize(posStr.length()-2);
                    }
                    
                    if (isCodon && (this->dataType==2 || partDataTypeID==2)) {
                        errMsg = "codon partition cannot be applied for amino acid data type";
                        break;
                    }
                }
                
                // creat a new partition
                newPart = new PartitionPos(name, seqLen, partDataTypeID, isCodon);

                // get the set of positions
                if (posStr.length() > 0) {
                    tokenizer(posStr, ",", &token1);
                    if (token1.size() > 0) {
                        for (i=0; i<token1.size(); i++) {
                            tokenizer(token1[i], "-", &token2);
                            if (token2.size() == 1) {
                                frPos = toPos = atoi(token2[0].c_str());
                                newPart->insertPair(frPos, toPos);
                            } else if (token2.size() > 1) {
                                frPos = atoi(token2[0].c_str());
                                toPos = atoi(token2[1].c_str());
                                if (isCodon) {
                                    for (j=frPos; j<=toPos; j+=3)
                                        newPart->insertPair(j,j);
                                } else {
                                    newPart->insertPair(frPos, toPos);
                                }
                            }
                        }
                    }
                }
                
                // check whether the new partition is empty
                if (newPart->size() == 0) {
                    errMsg = "no (valid) position";
                    break;
                }
                
                //save the new partition
                partitionPosList.push_back(newPart);
            }
        }
    }

    fin.close();
    
    if (errMsg == "" && partitionPosList.size() == 0) {
        cerr << "Error! No partition information inside the file " << partitionFile << endl;
        exit(1);
    }

    if (errMsg!="") {
        cerr << "Error! Inside the partition file " << partitionFile << "," << endl;
        cerr << errMsg << " in the line:" << endl;
        cerr << "     " << aline << endl;
        cerr << endl;
        if (dataType==3) {
            cerr << "When <data type> = 3, the format of the partition file should be:" << endl;
            cerr << "-----------------------------------------------------------------" << endl;
            cerr << "<NU/AA>, <partition name> = <start pos> - <end pos>, ..." << endl;
            cerr << "For example:" << endl;
            cerr << "   NU,part1=1-50,60-100" << endl;
            cerr << "   NU,part2=101-200" << endl;
            cerr << "   NU,part3=201-220,240-260,270-300" << endl << endl;
            
            cerr << "If considering the specific codon positions, the format should be:" << endl;
            cerr << "<NU/AA>, <partition name> = <start pos> - <end pos>, ...\\3" << endl;
            cerr << "For example:" << endl;
            cerr << "   NU,codon1=1-7\\3" << endl;
            cerr << "   NU,codon2=2-8\\3" << endl;
            cerr << "   NU,codon3=3-9\\3" << endl << endl;
            
            cerr << "(NA - nucleotide; AA - amino acid)" << endl;
            cerr << "(enumeration starts with 1, and all partition names are unique)" << endl;
        } else {
            cerr << "The format of the partition file should be:" << endl;
            cerr << "-------------------------------------------" << endl;
            cerr << "<partition name 1> = <start pos 1> - <end pos 1>, <start 2> - <end 2>, ..." << endl;
            cerr << "<partition name 2> = ..." << endl;
            cerr << "For example:" << endl;
            cerr << "   part1=1-50,60-100" << endl;
            cerr << "   part2=101-200" << endl;
            cerr << "   part3=201-220,240-260,270-300" << endl << endl;

            if (dataType==1) {
                cerr << "If considering the specific codon positions, the format should be:" << endl;
                cerr << "<partition name 1> = <start pos 1> - <end pos 1>, ...\\3" << endl;
                cerr << "<partition name 2> = ...\\3" << endl;
                cerr << "For example:" << endl;
                cerr << "   codon1=1-7\\3" << endl;
                cerr << "   codon2=2-8\\3" << endl;
                cerr << "   codon3=3-9\\3" << endl << endl;
            }

            cerr << "(enumeration starts with 1, and all partition names are unique)" << endl;
        }


        exit(1);
    }
    
    // consolidate all partitions
    for (i=0; i<partitionPosList.size(); i++) {
        partitionPosList[i]->consolidate();
    }
    
    // print the partition information
    cout << "================================================================================" << endl;
    cout << "Partition information:" << endl;
    cout << "--------------------------------------------------------------------------------" << endl;
    
    cout << "Name\tPositions";
    if (dataType==3)
        cout << "\tData type";
    cout << endl;
    for (i=0; i<partitionPosList.size(); i++) {
        partitionPosList[i]->print(dataType==3);
    }
    cout << "================================================================================" << endl;
    
}

// create partitions from sliding windows
void UserOptions::createPartFrSlideWin(int seqLen) {
    
    int i;
    PartitionPos* newPart;
    string name;
    int frPos, toPos;
    
    for (i=0; i<=seqLen-slideWindowSize; i+=slideWindowStep) {
        frPos = i+1;
        toPos = frPos+slideWindowSize-1;
        name = intToStr(frPos) + "-" + intToStr(toPos);
        // creat a new partition
        newPart = new PartitionPos(name, seqLen, 0, 0);
        newPart->insertPair(frPos, toPos);
        //save the new partition
        partitionPosList.push_back(newPart);
    }
    
    // check whether the last partition covers the last part of the alignment
    if (toPos < seqLen) {
        toPos = seqLen;
        frPos = toPos - slideWindowSize + 1;
        name = intToStr(frPos) + "-" + intToStr(toPos);
        // creat a new partition
        newPart = new PartitionPos(name, seqLen, 0, 0);
        newPart->insertPair(frPos, toPos);
        //save the new partition
        partitionPosList.push_back(newPart);
    }

    // consolidate all partitions
    for (i=0; i<partitionPosList.size(); i++) {
        partitionPosList[i]->consolidate();
    }
    
    // print the partition information
    cout << "================================================================================" << endl;
    cout << "Partition information:" << endl;
    cout << "--------------------------------------------------------------------------------" << endl;
    cout << "Index\tPositions" << endl;
    for (i=0; i<partitionPosList.size(); i++) {
        cout << i << "\t";
        partitionPosList[i]->printPosOnly();
    }
    cout << "================================================================================" << endl;
}


// constructor
PartitionPos::PartitionPos(string name, int seqLen, int partDataTypeID, int isCodon) {
    this->name = name;
    this->seqLen = seqLen;
    partLen = 0;
    dataType = partDataTypeID;
    this->isCodon = isCodon;
}

int PartitionPos::size() {
    return (int) pair_list.size();
}

void PartitionPos::insertPair(int frPos, int toPos) {
    // the enumeration of frPos and toPos start from 1
    if (frPos > 0 && toPos > 0 && toPos >= frPos && toPos <= seqLen) {
        pair_list.push_back(pair<int,int>(frPos-1, toPos-1));
    }
}

void PartitionPos::consolidate() {
    int i,j;
    int frPos, toPos;
    int* selectPos = new int[seqLen];
    
    // indicate whether the previous selected position has been inputted into the new list
    bool posConsidered = true;
    
    memset(selectPos, 0, seqLen*sizeof(int));
    for (i=0; i<pair_list.size(); i++) {
        frPos = pair_list[i].first;
        toPos = pair_list[i].second;
        for (j=frPos; j<=toPos; j++) {
            selectPos[j] = 1;
        }
    }
    
    partLen = 0;
    for (i=0; i<seqLen; i++)
        partLen += selectPos[i];
    
    vector<pair<int,int> > new_list;
    for (i=0; i<seqLen; i++) {
        if (!selectPos[i] && !posConsidered) {
            // input the position pair into the new list
            new_list.push_back(pair<int,int>(frPos,i-1));
            posConsidered = true;
        } else if (selectPos[i] && posConsidered) {
            frPos = i;
            posConsidered = false;
        }
    }
    if (!posConsidered) {
        new_list.push_back(pair<int,int>(frPos,seqLen-1));
    }
    pair_list.assign(new_list.begin(), new_list.end());
    delete[] selectPos;
}

void PartitionPos::print(bool showDataType) {
    int i;
    cout << name << "\t";
    for (i=0; i<pair_list.size(); i++) {
        if (i > 0)
            cout << ", ";
        if (pair_list[i].first==pair_list[i].second)
            cout << pair_list[i].first+1;
        else
            cout << pair_list[i].first+1 << "-" << pair_list[i].second+1;
    }
    if (showDataType) {
        cout << "\t";
        switch(dataType) {
            case 1:
                cout << "nucleotide";
                break;
            case 2:
                cout << "amino acid";
                break;
        }
    }
    if (isCodon) {
        cout << "\t[codon site]";
    }
    cout << endl;
}

void PartitionPos::printPosOnly() {
    int i;
    for (i=0; i<pair_list.size(); i++) {
        if (i > 0)
            cout << ", ";
        if (pair_list[i].first==pair_list[i].second)
            cout << pair_list[i].first+1;
        else
            cout << pair_list[i].first+1 << "-" << pair_list[i].second+1;
    }
    cout << endl;
}
