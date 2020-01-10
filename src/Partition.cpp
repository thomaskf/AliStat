/*
 *
 * Partition.cpp
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


#include "Partition.h"

// constructor
Partition::Partition(vector<string>* seqNames, vector<string>* origSeqs, int dataType, int seqNum, int seqLen,
                     char* validCharArray, UserOptions* userOptions ) {
    
    int i;
    this->seqNames = seqNames;
    this->origSeqs = origSeqs;
    this->dataType = dataType;
    this->seqNum = seqNum;
    this->seqLen = seqLen;
    this->partLen = seqLen;
    this->partName = "";
    this->validCharArray = validCharArray;
    this->userOptions = userOptions;
    this->partPos = NULL;
    Cr = new int[seqNum];
    Cc = new int[seqLen];
    Cc_part = new int[seqLen];
    Cij = new int[seqNum*seqNum];
    Cr_max = Cr_min = -1;
    Cc_max = Cc_min = -1;
    Cij_max = Cij_min = -1;
    memset(Cc, 0, seqLen*sizeof(int));
    row_index = new int[seqNum];
    for (i=0; i<seqNum; i++)
        row_index[i] = i;
    col_index = new int[seqLen];
    for (i=0; i<seqLen; i++)
        col_index[i] = i;
	recodeMatrix = NULL;
    if (userOptions->computePDist == 1) {
        Pij = new int[seqNum*seqNum];
    } else {
        Pij = NULL;
    }
}

// update the array of "validCharArray"
void Partition::setValidCharArray(char* validCharArray) {
    this->validCharArray = validCharArray;
}

// update the array of "recodeMatrix"
void Partition::setRecodeMatrix(char* recodeMatrix) {
	this->recodeMatrix = recodeMatrix;
}
	
// update the array "selectPos"
void Partition::setPartitionPos(PartitionPos* partPos) {

    int i,j,k;
    int frPos,toPos;
    this->partPos = partPos;
    this->partLen = partPos->partLen;
    this->partName = partPos->name;
    
    i=0;
    for (j=0; j<partPos->pair_list.size(); j++) {
        frPos = partPos->pair_list[j].first;
        toPos = partPos->pair_list[j].second;
        for (k=frPos; k<=toPos; k++)
            col_index[i++] = k;
    }
    
}


// compute all figures
void Partition::computeAllFigures(int showStatus) {
    computeCaCr();
    computeCc();
	if (userOptions->computePDist == 1) {
		computePijCij();
	} else {
		computeCij(showStatus);
	}
}


// output all figures
void Partition::outputAllFigures(string seqFile, string prefixOut) {

    // Output the summary to the file <prefixOut>.summary.txt
    int i;
    outputSummary(seqFile, prefixOut, (dataType!=3?dataType:partPos->dataType), (partPos!=NULL?partPos->isCodon:0),validCharArray, seqNum, partLen,
                  Ca, Cr_max, Cr_min, Cc_max, Cc_min, Cij_max, Cij_min,
                  Pij_max, Pij_min, Pij_avg, Pij_exist,
                  partName, userOptions);

    
    for (i=0; i<userOptions->outTables.size(); i++) {
        switch(userOptions->outTables[i]) {
            case 1:
                // Output the C scores for individual sequences (Cr) (i.e. table 1)
                // to the file <prefixOut>.table1.csv
                outputTable1(prefixOut, Cr, row_index, seqNames, seqNum, partLen);
                break;
            case 2:
                // Output the Distribution of C scores for individual sites (Cc) (i.e. table 2)
                // to the file <prefixOut>.table2.csv
                // format 0 : "Index\tSite ID\tCc"; format 1 : "Site ID\tCc"
                outputTable2(prefixOut, Cc, col_index, seqNum, partLen, partPos==NULL);
                break;
            case 3:
                // Output the Summary of distribution of C scores for individual sites (Cc) (i.e. table 3)
                // to the file <prefixOut>.table3.csv
                outputTable3(prefixOut, Cc_part, Cc_min, Cc_max, seqNum, partLen);
                break;
            case 4:
                if (seqNum > 1) {
                    // Output the Distribution of C scores for individual pairs of sequences (Cij) (i.e. table 4)
                    // to the file <prefixOut>.table4.csv
                    outputTable4(prefixOut, seqNames, Cij, row_index, seqNum, partLen);
                }
                break;
            case 5:
                if (seqNum > 1) {
                    // Output the Incompleteness scores (Iij = 1.0 - Cij) (i.e. table 5)
                    // to the file <prefixOut>.table5.csv
                    outputTable5(prefixOut, seqNames, Cij, row_index, seqNum, partLen);
                }
                break;
            case 6:
                if (seqNum > 1) {
                    // Output the Completeness scores (Cij) and Incompleteness scores (Iij = 1.0 - Cij) (i.e. table 6)
                    // to the file <prefixOut>.table6.csv
                    outputTable6(prefixOut, seqNames, Cij, row_index, seqNum, partLen);
                }
                break;
        }
    }
	
	// output the p distance table if necessary
	if (userOptions->computePDist == 1) {
		// Output the p-distances (Pij)
		// to the file <prefixOut>.p-dist.csv
		outputTablePDist(prefixOut, seqNames, Pij, Cij, row_index, seqNum, partLen);
	}
}


// output summary of figures to screen
void Partition::outSummaryToScreen(string seqFile, string partName, int showHeader) {
	outSumToScreen(seqFile, seqNum, partLen, Ca, Cr_max, Cr_min, Cc_max, Cc_min, Cij_max, Cij_min, partName, showHeader);
}


// output heat map files
void Partition::outputHeatMap(string prefixOut) {
	if (userOptions->makeHeatMap==1 || userOptions->makeHeatMap==3) {
		// Output the triangular heatmap
		outputTriHeatmap(prefixOut, seqNames, Cij, row_index, seqNum, partLen);
	}
	if (userOptions->makeHeatMap==2 || userOptions->makeHeatMap==3) {
		// Output the full heatmap
		outputFullHeatmap(prefixOut, seqNames, Cij, row_index, seqNum, partLen);
	}
}

// mask and output the alignment
void Partition::outputAlignment(string prefixOut) {

#define maxColPerLine 60
#define maskedSign    'I'
#define discardedSign 'E'
#define maskName      "Mask"

    if (!userOptions->outputAlign)
        return;
    
    int i,j;
    int index,pos;
    string seqMask;
    string seqDisc;
    string seqSign;
    
    string outMaskFile = fileName(prefixOut, "Mask.fst");
    string outDiscFile = fileName(prefixOut, "Disc.fst");
    string outStatFile = fileName(prefixOut, "Stat.fst");
    
    ofstream foutMask;
    ofstream foutDisc;
    ofstream foutStat;
    
    foutStat.open(outStatFile.c_str());
    
    // build up the sequence to show whether the column is included in 'Mask.fa' file
    seqSign = "";
    seqSign.append(seqLen, discardedSign);

	// three cases:
	// case 1: userOptions->maskThres = -1, then output all columns
	// case 2: userOptions->maskThres = 0, then output the columns with Cc values > 0
	// case 3: userOptions->maskThres > 0, then output the columns with Cc values >= userOptions->maskThres
	if (userOptions->maskThres==-1.0) {
		for (j=0; j<partLen; j++) {
			pos = col_index[j];
			seqSign[pos] = maskedSign;
		}
	} else if (userOptions->maskThres==0.0) {
		for (j=0; j<partLen; j++) {
			pos = col_index[j];
			if (Cc[pos]>0)
				seqSign[pos] = maskedSign;
		}
	} else {
		for (j=0; j<partLen; j++) {
			pos = col_index[j];
			if ((double)Cc[pos]/seqNum>=userOptions->maskThres) {
				seqSign[pos] = maskedSign;
			}
		}
	}
    foutStat << ">" << maskName << endl;
    printSeq(seqSign, &foutStat, maxColPerLine);
    for (i=0; i<seqNum; i++) {
        foutStat << ">" + seqNames->at(i) << endl;
        printSeq(origSeqs->at(i), &foutStat, maxColPerLine);
    }
    
    foutStat.close();
    
    
    foutMask.open(outMaskFile.c_str());
    foutDisc.open(outDiscFile.c_str());

    for (i=0; i<seqNum; i++) {
        index = row_index[i];
        seqMask = "";
        seqDisc = "";
        for (j=0; j<partLen; j++) {
            pos = col_index[j];
			if (seqSign[pos] == maskedSign)
                seqMask.append(1, (origSeqs->at(index))[pos]);
            else
                seqDisc.append(1, (origSeqs->at(index))[pos]);
        }
        if (seqMask.length() > 0) {
            foutMask << ">" + seqNames->at(index) << endl;
            printSeq(seqMask, &foutMask, maxColPerLine);
        }
        if (seqDisc.length() > 0) {
            foutDisc << ">" + seqNames->at(index) << endl;
            printSeq(seqDisc, &foutDisc, maxColPerLine);
        }
    }
    
    foutMask.close();
    foutDisc.close();
    
}


// ===================================================
// Supplementary functions
// ===================================================

// compute all values of Ca, Cr, Cr_max, Cr_min
void Partition::computeCaCr() {
    int i,j;

    for (i=0; i<seqNum; i++) {
        Cr[i] = 0;
        for (j=0; j<partLen; j++) {
            Cr[i] += validCharArray[origSeqs->at(i)[col_index[j]]];
        }
    }
    
    // get Cr_max, Cr_min
    if (seqNum > 0) {
        Cr_max = Cr_min = Cr[0];
        for (i=1; i<seqNum; i++) {
            if (Cr[i] < Cr_min)
                Cr_min = Cr[i];
            else if (Cr[i] > Cr_max)
                Cr_max = Cr[i];
        }
    }
    
    // get Ca
    Ca = 0;
    for (i=0; i<seqNum; i++) {
        Ca += Cr[i];
    }
}

// compute all values of Cc, Cc_max, Cc_min
void Partition::computeCc() {
    int i,j;
    int pos;
    int first=true;
    
    for (i=0; i<partLen; i++) {
        pos = col_index[i];
        if (Cc[pos]==0) {
            for (j=0; j<seqNum; j++) {
                Cc[pos] += validCharArray[(origSeqs->at(j))[pos]];
            }
        }
        // get Cc_max, Cc_min
        if (first) {
            Cc_max = Cc_min = Cc[pos];
            first = false;
        } else if (Cc[pos] < Cc_min)
            Cc_min = Cc[pos];
        else if (Cc[pos] > Cc_max)
            Cc_max = Cc[pos];
    }
    
    // compute Cc_part
    for (i=0; i<partLen; i++) {
        Cc_part[i] = Cc[col_index[i]];
    }
}


// compute all values of Cij, Cij_max, Cij_min
void Partition::computeCij(int showStatus) {
    
    if (seqNum < 2)
        return;
    
    rebuildBinaries();

    int i,j;
    
    // for i<j
    if (showStatus)
    	cout << "| 0%                                                                      100% |" << endl << flush;
    int totBarLen = 80;
    int curBarLen = 0;
    int newBarLen;
    for (i=0; i<seqNum-1; i++) {
        
        // for showing status
	    if (showStatus) {
			newBarLen = ((i+1)*totBarLen)/(seqNum-1);
			for (j=curBarLen; j<newBarLen; j++)
				cout << "*" << flush;
			curBarLen = newBarLen;
        }
        
        for (j=i+1; j<seqNum; j++) {
            Cij[i*seqNum+j] = boolArray[i].countInter(&(boolArray[j]));
        }
    }
    if (showStatus)
    	cout << endl;
    
    // for i>j
    for (i=1; i<seqNum; i++)
        for (j=0; j<i; j++)
            Cij[i*seqNum+j] = Cij[j*seqNum+i];
    
    // for i==j
    for (i=0; i<seqNum; i++)
        Cij[i*seqNum+i] = partLen;
    
    // compute Cij_max and Cij_min
    Cij_max = Cij_min = Cij[1];
    for (i=0; i<seqNum-1; i++) {
        for (j=i+1; j<seqNum; j++) {
            if (Cij_min > Cij[i*seqNum+j])
                Cij_min = Cij[i*seqNum+j];
            else if (Cij_max < Cij[i*seqNum+j])
                Cij_max = Cij[i*seqNum+j];
        }
    }
}

// compute the p-distances and Cij between sequences
void Partition::computePijCij() {
	
    if (seqNum < 2) {
		Pij_exist = false;
        return;
	}
    
    rebuildRecodeSeqs();

    int i,j,k;
    
    // for i<j
    cout << "| 0%                                                                      100% |" << endl << flush;
    int totBarLen = 80;
    int curBarLen = 0;
    int newBarLen;
    for (i=0; i<seqNum-1; i++) {
        
        // for showing status
        newBarLen = ((i+1)*totBarLen)/(seqNum-1);
        for (j=curBarLen; j<newBarLen; j++)
            cout << "*" << flush;
        curBarLen = newBarLen;
        
        for (j=i+1; j<seqNum; j++) {
			Pij[i*seqNum+j] = 0;
			Cij[i*seqNum+j] = 0;
			for (k=0; k<partLen; k++) {
				if (recodeSeqs[i].at(k)!=0 && recodeSeqs[j].at(k)!=0) {
					Cij[i*seqNum+j]++;
					if (recodeSeqs[i].at(k) != recodeSeqs[j].at(k))
						Pij[i*seqNum+j]++;
				}
			}
        }
    }
    cout << endl;
    
    // for i>j
    for (i=1; i<seqNum; i++) {
        for (j=0; j<i; j++) {
            Cij[i*seqNum+j] = Cij[j*seqNum+i];
			Pij[i*seqNum+j] = Pij[j*seqNum+i];
		}
	}
    
    // for i==j
    for (i=0; i<seqNum; i++) {
        Cij[i*seqNum+i] = partLen;
        Pij[i*seqNum+i] = 0;
	}
    
    // compute Cij_max and Cij_min
    Cij_max = Cij_min = Cij[1];
    for (i=0; i<seqNum-1; i++) {
        for (j=i+1; j<seqNum; j++) {
            if (Cij_min > Cij[i*seqNum+j])
                Cij_min = Cij[i*seqNum+j];
            else if (Cij_max < Cij[i*seqNum+j])
                Cij_max = Cij[i*seqNum+j];
        }
    }
	
    // compute Pij_max, Pij_min and Pij_avg
	double Pij_sum = 0.0;
	int Pij_num = 0;
	bool isFirst = true;
	double Pij_curr;
    for (i=0; i<seqNum-1; i++) {
        for (j=i+1; j<seqNum; j++) {
			if (Cij[i*seqNum+j] > 0) {
				Pij_curr = (double)Pij[i*seqNum+j] / Cij[i*seqNum+j];
				if (isFirst) {
					isFirst = false;
					Pij_max = Pij_curr;
					Pij_min = Pij_curr;
				} else if (Pij_max < Pij_curr) {
					Pij_max = Pij_curr;
				} else if (Pij_min > Pij_curr) {
					Pij_min = Pij_curr;
				}
				Pij_sum += Pij_curr;
				Pij_num++;
			}
        }
    }
	if (Pij_num > 0) {
		Pij_exist = true;
		Pij_avg = Pij_sum / Pij_num;
	} else {
		Pij_exist = false;
	}
}

// rebuild the binary sequences
void Partition::rebuildBinaries() {
    int i,j;
    int pos;
    
    boolArray.resize(seqNum);
    for (i=0; i<seqNum; i++) {
        boolArray[i].clear();
        for (j=0; j<partLen; j++) {
            pos = col_index[j];
            boolArray[i].push_back(validCharArray[origSeqs->at(i)[pos]]);
        }
    }
}

// rebuild the recode sequences
void Partition::rebuildRecodeSeqs() {
	int i,j;
	int pos;
	
	recodeSeqs.resize(seqNum);
    for (i=0; i<seqNum; i++) {
        recodeSeqs[i].clear();
        for (j=0; j<partLen; j++) {
            pos = col_index[j];
            recodeSeqs[i].append(1, recodeMatrix[origSeqs->at(i)[pos]]);
        }
    }
}


// reorder the alignment
void Partition::reorderAlignment() {
    
#define descendingOrder 1
    
    int i;

    if (userOptions->reorder & 1) {
        // 1 or 3
        // reorder the rows
        sortUnsignInt(Cr, seqNum, descendingOrder, row_index);
    }
    
    if (userOptions->reorder & 2) {
        // 2 or 3
        // reorder the columns
        
        // sort the Cc values
        int* sortedIndex;
        if (partPos == NULL) {
            sortedIndex = col_index;
        } else {
            sortedIndex = new int[partLen];
        }
        sortUnsignInt(Cc_part, partLen, descendingOrder, sortedIndex);
        
        // get the corresponding column indices
        if (partPos != NULL) {
            for (i=0; i<partLen; i++) {
                sortedIndex[i] = col_index[sortedIndex[i]];
            }
            for (i=0; i<partLen; i++) {
                col_index[i] = sortedIndex[i];
            }
        }
        
        // release memory
        if (partPos != NULL) {
            delete[] sortedIndex;
        }
    }
}

