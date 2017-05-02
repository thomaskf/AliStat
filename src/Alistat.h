/*
 *
 * Alistat.h
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


#ifndef Alistat_alistat_h
#define Alistat_alistat_h

#include <vector>
#include <set>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include "BoolArray.h"
#include "UserOptions.h"
#include "tool_box.h"
#include "Output.h"
#include "Partition.h"

using namespace std;


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
               vector<string>* seqs, int& seqNum, int& seqLen);



// read the arguments and collect the user options
// return:
//    0 if the user options are valid
//    1 if there exists error message
//    2 if a help menu is called
int readArguments(int argc, char** argv, UserOptions* user_options, string* errMsg);



// To create an array to indicate which characters are valid
// input: (1) dataType - 1 for DNA; 2 for amino acid; (2) codeType
// output: char array X, where X[c] is 1 if character c is valid
char* validCharArray(int dataType, int codeType);


// To create the recoding matrix
// which indicates the mapping between the original character and the recoded character
char* makeRecodeMatrix(int dataType, int codeType);


// To print out the sequence names and the binary sequences
void printNameAndBinarySeq(int seqNum, vector<string>* seqNames, vector<BoolArray*> * binSeqs);


#endif
