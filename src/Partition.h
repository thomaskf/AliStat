/*
 *
 * Partition.h
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


#ifndef __alistat__Partition__
#define __alistat__Partition__

#include <iostream>
#include "UserOptions.h"
#include "BoolArray.h"
#include "Output.h"
#include "tool_box.h"

class Partition {

public:
    
    // constructor
    Partition(vector<string>* seqNames, vector<string>* origSeqs, int dataType, int seqNum, int seqLen, char* validCharArray, UserOptions* userOption);
    
    // update the array of "validCharArray"
    void setValidCharArray(char* validCharArray);
    
    // update the array of "recodeMatrix"
    void setRecodeMatrix(char* recodeMatrix);

    // update the array "selectPos"
    void setPartitionPos(PartitionPos* partPos);
    
    // compute all figures
    void computeAllFigures();
    
    // reorder the alignment
    void reorderAlignment();
    
    // output all figures
    void outputAllFigures(string seqFile, string prefixOut);
    
    // mask and output the alignment
    void outputAlignment(string prefixOut);

	// output heat map files
	void outputHeatMap(string prefixOut);
	
private:
    
    // original sequences information
    vector<string>* seqNames;
    vector<string>* origSeqs;
    int dataType;
    int seqNum;
    int seqLen;
    int partLen;
    string partName;
    
    // valid character set
    char* validCharArray;
	
	// recoding matrix
	char* recodeMatrix;
    
    // user options
    UserOptions* userOptions;

    // partition pos information
    PartitionPos* partPos;

    // binary sequences;
    vector<BoolArray> boolArray;
	
	// recode sequences;
	vector<string> recodeSeqs;
    
    // index for row and column
    int* row_index; // sorted row index according to Cr if necessary
    int* col_index; // its size = partLen
                    // selected for those belong to the partition
                    // sorted if necessary

    // statistics
    int Ca;
    int* Cr; int Cr_max; int Cr_min;
    int* Cc; int Cc_max; int Cc_min;
    int* Cc_part; // a set of Cc values according to the positions inside the partition
    int* Cij; int Cij_max; int Cij_min;
	// p distances between sequences
	int* Pij; double Pij_max; double Pij_min; double Pij_avg; bool Pij_exist;
    
    // compute all values of Ca, Cr, Cr_max, Cr_min
    void computeCaCr();
    
    // compute all values of Cc, Cc_max, Cc_min
    void computeCc();
    
    // compute all values of Cij, Cij_max, Cij_min
    void computeCij();
    
    // rebuild the binary sequences
    void rebuildBinaries();
	
	// rebuild the recode sequences
	void rebuildRecodeSeqs();
	
	// compute the p-distances and Cij between sequences
	void computePijCij();
    
};


#endif /* defined(__alistat__Partition__) */
