/*
 *
 * UserOptions.h
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


#ifndef __alistat__UserOptions__
#define __alistat__UserOptions__

#include <iostream>
#include <fstream>
#include <string>
#include <set>
#include <vector>
#include "tool_box.h"
#include "GetOptions.h"

#define TOTAL_TABLES_NUM 6

class PartitionPos {
public:
    string name;
    int seqLen;
    int partLen;
    vector<pair<int,int> > pair_list;
    int dataType; // 0: as specified in the argument; 1: nucleotide; 2: amino acid
    int isCodon; // 0: not a codon region; 1: is a codon region
    
    PartitionPos(string name, int seqLen, int partDataTypeID, int isCodon);
    int size();
    // the enumeration of frPos and toPos start from 1
    void insertPair(int frPos, int toPos);
    void consolidate();
    void print(bool showDataType);
    void printPosOnly();
};


class UserOptions {
public:
    
    // ======================================================
    // Variables
    // ======================================================
    
    string alignFile; // multiple alignment file with FASTA format
    
    int dataType;
    /*
     1  - Single nucleotides (SN);
     2  - Di-nucleotides (DN);
     3  - Codons (CD);
     4  - 10-state genotype data (10GT);
     5  - 14-state genotype data (14GT);
     6  - Amino acids (AA);
     7  - Mixture of nucleotides and amino acids (NA)
    */
    
    string dataTypeStr;
    
    int codeType;
    // 0  - A, C, G, T [default]
    // 1  - C, T, R (i.e. C, T, AG)
    // 2  - A, G, Y (i.e. A, G, CT)
    // 3  - A, T, S (i.e. A, T, CG)
    // 4  - C, G, W (i.e. C, G, AT)
    // 5  - A, C, K (i.e. A, C, GT)
    // 6  - G, T, M (i.e. G, T, AC)
    // 7  - K, M    (i.e. GT,   AC)
    // 8  - S, W    (i.e. GC,   AT)
    // 9  - R, Y    (i.e. AG,   CT)
    // 10 - A, B    (i.e. A,   CGT)
    // 11 - C, D    (i.e. C,   AGT)
    // 12 - G, H    (i.e. G,   ACT)
    // 13 - T, V    (i.e. T,   ACG)
    
    string prefixOut; // prefix for output files (default: <alignment file> w/o .ext)
    
    string seqNameFile; // FILE contains the names the sequences should only be considered
    
    string partitionFile; // FILE contains the partition information
    
    // sliding window analysis
    int slideWindowSize;
    int slideWindowStep;
    
    vector<int> outTables; // tables which should only be outputted
    
    int reorder; // 0 : no change; 1: reorder the rows; 2: reorder the columns; 3: reorder both rows and columns
    
    set<string> namesConsidered; // the sequence names to be considered
    
    vector<PartitionPos*> partitionPosList; // the list of partition positions
    
    double maskThres; // the threshold for Cc to mask the alignment
                      // -1.0 : no masking is necessary
	
	int outputAlign; // whether outputting the alignment file
	                 // 1: yes; 0: no (default: no)
					  
	int computePDist; // 1 : compute the p-dist; 0 : no (default : no)
	
	int makeHeatMap; // 0 : No; 1 : Triangular heat map;
				     // 2 : Rectangular heat map; 3: Both
				     
	bool briefOutput; // true if option -b is used
    
    int color_scheme; // 1, 2 or 3; default is 1
    
    // ======================================================
    // Functions
    // ======================================================
    
    // load the file including the sequence names considered
    void loadSeqNameFile();
    
    // load the partition file
    void loadPartitionFile(int seqLen);
    
    // create partitions from sliding windows
    void createPartFrSlideWin(int seqLen);
    
    // reset
    void reset();
    
    // (obsoleted) by default, prefixOut = <alignment file> w/o .ext
    // (changed on 12/05/2014) by default, prefixOut = empty string
    void setDefaultPrefixOut();
    
    // by default, outTables = { 1, 2, 3, 4, 5 }
    void setDefaultOutTables();

    // get the data type description
    void getDataTypeStr();
};

#endif /* defined(__alistat__UserOptions__) */
