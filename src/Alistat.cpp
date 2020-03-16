/*
 *
 * Alistat.cpp
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


#include "Alistat.h"


int main(int argc, char** argv) {
    
    vector<string> seqNames;
    vector<string> seqs;
    int seqLen,seqNum;
    char* validCharArr = NULL;
    char* validSNChars = NULL;   // valid characters for single nucleotide (i.e. A, C, G, T)
    char* validDNChars = NULL;   // valid characters for di-nucleotide (i.e. 1 - 16)
    char* validCDChars = NULL;   // valid characters for codon (i.e. 1 - 64)
    char* valid10GTChars = NULL; // valid characters for 10GT (i.e. ACGTUWSMKRY)
    char* valid14GTChars = NULL; // valid characters for 14GT (i.e. ACGTUWSMKRYBDHV)
    char* validAAChars = NULL;   // valid characters for amino acid (i.e. ARNDCEQGHILKMFPSTWYV)
    UserOptions user_options;
    string errMsg;
    set<string>* seqNameSelected = NULL;
    Partition *partition;
    int i;
	char* recodeMatrix = NULL;
	char* recodeSNMatrix = NULL;
    char* recodeDNMatrix = NULL;
    char* recodeCDMatrix = NULL;
    char* recode10GTMatrix = NULL;
    char* recode14GTMatrix = NULL;
	char* recodeAAMatrix = NULL;
	string partitionName;
	int showHeader;
	int showStatus;

    // read the arguments
    int status = readArguments(argc, argv, &user_options, &errMsg);
    
    if (status==1) {
        // error
        cerr << errMsg << endl;
        exit(1);
    }
    
    if (status==2) {
        // invoke the help page
        outputUsage(argv[0]);
        exit(1);
    }
    
    // To create an array to indicate which characters are valid
    if (user_options.dataType == 7) {
        // Mixture of nucleotides and amino acids (NA)
        validSNChars = validCharArray(1,0); // default DNA characters
        validDNChars = validCharArray(2,0);
        validCDChars = validCharArray(3,0);
        valid10GTChars = validCharArray(4,0);
        valid14GTChars = validCharArray(5,0);
        validAAChars = validCharArray(6,0); // default Amino acid characters
    } else {
        validCharArr = validCharArray(user_options.dataType, user_options.codeType);
    }
	
	// To create a recoding matrix if necessary
	if (user_options.computePDist == 1) {
		if (user_options.dataType == 7) {
            // Mixture of nucleotides and amino acids (NA)
			recodeSNMatrix = makeRecodeMatrix(1,0); // default DNA recoding matrix
            recodeDNMatrix = makeRecodeMatrix(2,0);
            recodeCDMatrix = makeRecodeMatrix(3,0);
            recode10GTMatrix = makeRecodeMatrix(4,0);
            recode14GTMatrix = makeRecodeMatrix(5,0);
			recodeAAMatrix = makeRecodeMatrix(6,0); // default Amino acid recoding matrix
		} else {
			recodeMatrix = makeRecodeMatrix(user_options.dataType, user_options.codeType);
		}
	}
    
    // Output the start message
    outputStartMessage(user_options, validCharArr, validSNChars, validDNChars, validCDChars, valid10GTChars, valid14GTChars, validAAChars);
    
    // get the selected sequence names
    if (user_options.seqNameFile != "") {
        user_options.loadSeqNameFile();
        seqNameSelected = &(user_options.namesConsidered);
    }
    
    // To read the alignments from a file with fasta format
    // input:
    // (1) fileName - the file name
    // (2) seqNameSelected - the array of sequence names which should be considered only
    // output:
    // (1) seqNames - a set of sequence names, and
    // (2) seqs - a set of sequences
    // (3) SeqNum - number of sequences
    // (4) seqLen - the length of all the sequences (i.e. all the sequences should have the same length)
    readFasta(user_options.alignFile, seqNameSelected, &seqNames, &seqs, seqNum, seqLen);
    
    partition = new Partition(&seqNames, &seqs, user_options.dataType, seqNum, seqLen, validCharArr,
                              &user_options);
    
    if (user_options.slideWindowSize != -1 || user_options.partitionFile != "") {
        // sliding window is used OR A partition file is provided
        
        string term;
        if (user_options.slideWindowSize != -1) {
            // create partitions from sliding windows
            user_options.createPartFrSlideWin(seqLen);
            term = "partition with positions:";
        } else {
            // A partition file is provided
            user_options.loadPartitionFile(seqLen);
            term = "partition";
        }
        
        for (i=0; i<user_options.partitionPosList.size(); i++) {
            partitionName = user_options.partitionPosList[i]->name;
            
            partition->setPartitionPos(user_options.partitionPosList[i]);
            partition->buildUpdateSeqs();
            partition->setupArray();

            if (user_options.dataType == 7) {
                switch(user_options.partitionPosList[i]->dataType) {
                    case 1:
                        partition->setValidCharArray(validSNChars);
                        break;
                    case 2:
                        partition->setValidCharArray(validDNChars);
                        break;
                    case 3:
                        partition->setValidCharArray(validCDChars);
                        break;
                    case 4:
                        partition->setValidCharArray(valid10GTChars);
                        break;
                    case 5:
                        partition->setValidCharArray(valid14GTChars);
                        break;
                    case 6:
                        partition->setValidCharArray(validAAChars);
                        break;
                }
            }
            
			if (user_options.computePDist == 1) {
				if (user_options.dataType == 7) {
                    switch(user_options.partitionPosList[i]->dataType) {
                        case 1:
                            partition->setRecodeMatrix(recodeSNMatrix);
                            break;
                        case 2:
                            partition->setRecodeMatrix(recodeDNMatrix);
                            break;
                        case 3:
                            partition->setRecodeMatrix(recodeCDMatrix);
                            break;
                        case 4:
                            partition->setRecodeMatrix(recode10GTMatrix);
                            break;
                        case 5:
                            partition->setRecodeMatrix(recode14GTMatrix);
                            break;
                        case 6:
                            partition->setRecodeMatrix(recodeAAMatrix);
                            break;
                    }
				} else {
					partition->setRecodeMatrix(recodeMatrix);
				}
			}
            
            partition->dataType = user_options.partitionPosList[i]->dataType;
            partition->setDataTypeStr();

            if (user_options.briefOutput) {
				showStatus = 0;
				partition->computeAllFigures(showStatus);
				if (i==0)
					showHeader = 1;
				else
					showHeader = 0;
				partition->outSummaryToScreen(user_options.alignFile, partitionName, showHeader);
			} else {
				showStatus = 1;
				cout << "Performing analysis on " << term << " " << partitionName << endl;
				partition->computeAllFigures(showStatus);
				partition->reorderAlignment();
				partition->outputAllFigures(user_options.alignFile, fileName(user_options.prefixOut, partitionName));
				partition->outputAlignment(fileName(user_options.prefixOut, partitionName));
				partition->outputHeatMap(fileName(user_options.prefixOut, partitionName));
				// output the R script
				outputRScript(fileName(user_options.prefixOut, partitionName), user_options, partitionName);
				// Output the finish message
				outputFinishMessage(fileName(user_options.prefixOut, partitionName), seqNum, user_options.partitionPosList[i]->partLen, &user_options);
            }
        }
    } else {
        // No partition file is provided

        partition->buildUpdateSeqs();
        partition->setupArray();
        partition->dataType = user_options.dataType;
        partition->setDataTypeStr();

        if (user_options.briefOutput) {
			showStatus = 0;
        	partition->computeAllFigures(showStatus);
        	showHeader = 1;
        	partitionName = "";
			partition->outSummaryToScreen(user_options.alignFile, partitionName, showHeader);
        } else {
			showStatus = 1;
			if (user_options.computePDist == 1) {
				partition->setRecodeMatrix(recodeMatrix);
			}
			partition->computeAllFigures(showStatus);
			partition->reorderAlignment();
			partition->outputAllFigures(user_options.alignFile, user_options.prefixOut);
			partition->outputAlignment(user_options.prefixOut);
			partition->outputHeatMap(user_options.prefixOut);
			// output the R script
			outputRScript(user_options.prefixOut, user_options, "");
			// Output the finish message
			outputFinishMessage(user_options.prefixOut, seqNum, seqLen, &user_options);
        }
    }
    
    // output the message about how to use the R scripts generated if necessary
    outputRMessage(user_options);

    // release memory
    delete[] validCharArr;

    return 0;
}



// To read the alignments from a file with fasta format
// input:
// (1) fileName - the file name
// (2) seqNameSelected - the array of sequence names which should be considered only
// output:
// (1) seqNames - a set of sequence names, and
// (2) seqs - a set of sequences
// (3) SeqNum - number of sequences
// (4) seqLen - the length of all the sequences (i.e. all the sequences should have the same length)
void readFasta(string fileName, set<string>* seqNameSelected, vector<string>* seqNames,
               vector<string>* seqs, int& seqNum, int& seqLen) {

    ifstream fin;
    string aline;
    string seqName;
    size_t spacePos;
    int i;
    bool seqToIgnore = true;
    
    fin.open(fileName.c_str());
    while (getline(fin,aline)) {
        trimEnd(aline); // remove the space of invisible chars at the end of the sequence
        if (aline.length() == 0)
            continue;
        if (aline[0]=='>') {
            // get the sequence name
            spacePos = aline.find(' ');
            if (spacePos==string::npos) {
                // no space
                if (aline.length() > 1)
                    seqName = aline.substr(1);
                else
                    seqName = "";
            } else {
                // space inside the name
                if (spacePos > 1)
                    seqName = aline.substr(1,spacePos-1);
                else
                    seqName = "";
            }
            if (seqNameSelected!=NULL && seqNameSelected->find(seqName)==seqNameSelected->end()) {
                seqToIgnore = true;
            } else {
                seqToIgnore = false;
                seqNames->push_back(seqName);
                // create a new empty sequence
                seqs->push_back("");
            }
        } else if (!seqToIgnore) {
            for (i=0; i<aline.length(); i++)
                seqs->at(seqs->size()-1).append(1, toupper(aline[i]));
        }
    }
    fin.close();
    
    seqNum = (int) seqNames->size();
    
    if (seqNum == 0) {
        cerr << "Error! No sequence has been loaded" << endl;
        exit(1);
    }
    
    seqLen = (int) seqs->at(0).length();
    for (i=1; i<seqNum; i++) {
        if (seqs->at(i).length() != seqLen) {
            cerr << "Error! The length of the " << i+1 << "-th sequence does not match with that of the prevous ones" << endl;
            exit(1);
        }
    }

}


// To create an array to indicate which characters are valid
// input: (1) dataType - 1-5 for DNA; 6 for amino acid; (2) codeType
// output: char array X, where X[c] is 1 if character c is valid
char* validCharArray(int dataType, int codeType) {
    
    char* validCharArr = new char[256];
    memset(validCharArr, 0, 256 * sizeof(char));
    
    string validChars;
    int i;
    
    if (dataType == 1) {
        // SN
        validChars = "ACGTU";
    } else if (dataType == 2) {
        // DN
        validChars = string(16,0);
        for (i=0; i<16; i++) {
            validChars[i] = i+1;
        }
    } else if (dataType == 3) {
        // Codons (CD)
        validChars = string(64,0);
        for (i=0; i<64; i++) {
            validChars[i] = i+1;
        }
    } else if (dataType == 4) {
        // 10GT
        validChars = "ACGTUWSMKRY";
    } else if (dataType == 5) {
        // 14GT
        validChars = "ACGTUWSMKRYBDHV";
    } else if (dataType == 6) {
        // AA
        validChars = "ARNDCEQGHILKMFPSTWYV";
    } else {
        // Mixture of nucleotides and amino acids (NA)
        validChars = "";
    }
    
    string otherValidChars = "";
    if (dataType == 1 && codeType!=0) {
        switch(codeType) {
            case 1:
                otherValidChars = "R";
                break;
            case 2:
                otherValidChars = "Y";
                break;
            case 3:
                otherValidChars = "S";
                break;
            case 4:
                otherValidChars = "W";
                break;
            case 5:
                otherValidChars = "K";
                break;
            case 6:
                otherValidChars = "M";
                break;
            case 7:
                otherValidChars = "KM";
                break;
            case 8:
                otherValidChars = "SW";
                break;
            case 9:
                otherValidChars = "RY";
                break;
            case 10: // B (i.e. CGT)
				// also consider S (i.e. CG), Y (i.e. CT), K (i.e. GT)
                otherValidChars = "BSYK";
                break;
            case 11: // D (i.e. AGT)
				// also consider R (i.e. AG), W (i.e. AT), K (i.e. GT)
                otherValidChars = "DRWK";
                break;
            case 12: // H (i.e. ACT)
				// also consider M (i.e. AC), W (i.e. AT), Y (i.e. CT)
                otherValidChars = "HMWY";
                break;
            case 13: // V (i.e. ACG)
				// also consider M (i.e. AC), R (i.e. AG), S (i.e. CG)
                otherValidChars = "VMRS";
                break;
        }
    }
    
    for (i=0; i<validChars.length(); i++)
        validCharArr[(int)validChars[i]] = 1;
    for (i=0; i<otherValidChars.length(); i++)
        validCharArr[(int)otherValidChars[i]]=1;
    
    return validCharArr;
}

// To create the recoding matrix
// which indicates the mapping between the original character and the recoded character
char* makeRecodeMatrix(int dataType, int codeType) {
    char* recodeMatrix = new char[256];
	memset(recodeMatrix, 0, 256 * sizeof(char));
	
	if (dataType == 7) {
        // Mixture of nucleotides and amino acids (NA)
		return recodeMatrix;
	}
	
	int i;
	char* frChars = NULL;
	char* toChars = NULL;
	int numChars = 0;
	if (dataType == 1) {
		// DNA
        switch(codeType) {
			case 1: // A,G -> R
				frChars= (char*) "ACGTRU";
				toChars= (char*) "RCRTRU";
				numChars=6;
				break;
			case 2: // C,T -> Y
				frChars= (char*) "ACGTYU";
				toChars= (char*) "AYGYYY";
				numChars=6;
				break;
			case 3: // C,G -> S
				frChars= (char*) "ACGTSU";
				toChars= (char*) "ASSTSU";
				numChars=6;
				break;
			case 4: // A,T -> W
				frChars= (char*) "ACGTWU";
				toChars= (char*) "WCGWWW";
				numChars=6;
				break;
			case 5: // G,T -> K
				frChars= (char*) "ACGTKU";
				toChars= (char*) "ACKKKK";
				numChars=6;
				break;
			case 6: // A,C -> M
				frChars= (char*) "ACGTMU";
				toChars= (char*) "MMGTMU";
				numChars=6;
				break;
			case 7: // A,C -> M; G,T -> K
				frChars= (char*) "ACGTMKU";
				toChars= (char*) "MMKKMKK";
				numChars=7;
				break;
			case 8: // A,T -> W; C,G -> S
				frChars= (char*) "ACGTWSU";
				toChars= (char*) "WSSWWSW";
				numChars=7;
				break;
			case 9: // A,G -> R; C,T -> Y
				frChars= (char*) "ACGTRYU";
				toChars= (char*) "RYRYRYY";
				numChars=7;
				break;
			case 10: // C,G,T -> B
				// also consider S (i.e. CG), Y (i.e. CT), K (i.e. GT)
				frChars= (char*) "ACGTSYKBU";
				toChars= (char*) "ABBBBBBBB";
				numChars=9;
				break;
			case 11: // A,G,T -> D
				// also consider R (i.e. AG), W (i.e. AT), K (i.e. GT)
				frChars= (char*) "ACGTRWKDU";
				toChars= (char*) "DCDDDDDDD";
				numChars=9;
				break;
			case 12: // A,C,T -> H
				// also consider M (i.e. AC), W (i.e. AT), Y (i.e. CT)
				frChars= (char*) "ACGTMWYHU";
				toChars= (char*) "HHGHHHHHH";
				numChars=9;
				break;
			case 13: // A,C,G -> V
				// also consider M (i.e. AC), R (i.e. AG), S (i.e. CG)
				frChars= (char*) "ACGTMRSVU";
				toChars= (char*) "VVVTVVVVU";
				numChars=9;
				break;
			default: // 0 or other cases
				frChars= (char*) "ACGTU";
				toChars= (char*) "ACGTU";
				numChars=5;
				break;
		}
    /*
    } else if (dataType > 1 && dataType <= 5) {
        frChars= (char*) "ACGTU";
        toChars= (char*) "ACGTU";
        numChars=5;
    */
	} else if (dataType == 6) {
		// Amino acid
        frChars =  (char*) "ARNDCEQGHILKMFPSTWYV";
        toChars =  (char*) "ARNDCEQGHILKMFPSTWYV";
		numChars = 20;
	}
	
    if (frChars != NULL) {
        for (i=0; i<numChars; i++) {
            recodeMatrix[(int)frChars[i]]=toChars[i];
        }
    } else {
        for (i=0; i<256; i++) {
            recodeMatrix[i] = i;
        }
    }
	return recodeMatrix;
}


// To print out the sequence names and the binary sequences
void printNameAndBinarySeq(int seqNum, vector<string>* seqNames, vector<BoolArray*> * binSeqs) {
    int i;
    for (i=0; i<seqNum; i++) {
        cout << seqNames->at(i) << endl;
        binSeqs->at(i)->print();
    }
}

