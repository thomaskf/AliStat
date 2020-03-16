/*
 *
 * tool_box.h
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


#ifndef __Alistat__tool_box__
#define __Alistat__tool_box__

#include <iostream>
#include <string>
#include <cstring>
#include <vector>
#include <map>
#include <stdlib.h>
#include <time.h>
#include <math.h>

using namespace std;


// split the seq into many parts by "separators"
// note: the vector<string> *result cannot be NULL
static void tokenizer(string seq, string separators, vector<string>* result) {

    result->clear();
    int startpos = (int) seq.find_first_not_of(separators);
    while (startpos != string::npos) {
        int endpos = (int) seq.find_first_of(separators, startpos);
        if (endpos != string::npos) {
            result->push_back(seq.substr(startpos, endpos-startpos));
            startpos = (int) seq.find_first_not_of(separators, endpos);
        } else {
            result->push_back(seq.substr(startpos));
            break;
        }
    }
}

// return true if c is space or non-printable character
static bool isSpaceOrNonprintableChar(char c) {
    return (c <= 32 || c >= 127);
}


// remove the non-printable characters or space at the end of the string
static void trimEnd(string& str) {

    int i;
    for (i=(int)str.length()-1; i>=0; i--) {
        if (!isSpaceOrNonprintableChar(str[i]))
            break;
    }
    if (i==str.length()-1) {
        // do nothing
    } else if (i>=0) {
        // resize the string
        str.resize(i+1);
    } else {
        // empty the string
        str = "";
    }
}

// remove the non-printable characters or space in the front of the string
static void trimFront(string& str) {
    
    int i;
    for (i=0; i<(int)str.length(); i++) {
        if (!isSpaceOrNonprintableChar(str[i]))
            break;
    }
    if (i==0) {
        // do nothing
    } else if (i<(int)str.length()) {
        // resize the string
        str=str.substr(i);
    } else {
        // empty the string
        str = "";
    }
}

// remove the non-printable charactes or space in the front or at the end of the string
static void trimAll(string& str) {
    trimFront(str);
    trimEnd(str);
}

// remove all the non-printable characters or space in the string
static void removeAllSpaces(string& str) {

    int i,j;
    j=0;
    for (i=0; i<(int)str.length(); i++) {
        if (!isSpaceOrNonprintableChar(str[i])) {
            // not a space / nonprintable character
            if (i > j)
                str[j] = str[i];
            j++;
        }
    }
    if (j==0)
        str = "";
    else if (j < str.length())
        str.resize(j);
}

// get the distribution
// input: (1) an integer array; (2) array size; (3) min value; and (4) max value
// output: the distribution of the values inside the integer array
static int* getDist(int* intArrays, int size, int min, int max) {
    
    int i;
    if (min > max)
        return NULL;
    
    int* distributions = new int[max-min+1];
    memset(distributions, 0, (max-min+1)*sizeof(int));
    for (i=0; i<size; i++) {
        distributions[intArrays[i]-min]++;
    }
    return distributions;
}

// sort the integers (by radix sort)
// output: the sorted index ID
// order: 0 - ascending; 1 - descending
static void sortUnsignInt(int* intArrays, int size, int order, int* sortedIndexID) {
    int digits = 8; // either 8 or 16
    int iterationNum = 32/digits;
    int bucketNum = (int) pow(2,8);
    int mask = 0x000000FF;
    int* buffer[2];
    int i,k;
    
    //initialization
    buffer[0] = new int[size];
    buffer[1] = new int[size];
    for (int i=0; i<size; i++)
        (buffer[0])[i] = i;
    
    int* numEachBucket = new int[bucketNum];
    int* startPosEachBucket = new int[bucketNum];
    int* numEachBucketSoFar = new int[bucketNum];
    
    int currBufferID = 0;
    int* currBuffer = buffer[0];
    int* preBuffer = buffer[1];
    
    for (i=0; i<iterationNum; i++) {
        
        // update the pointers of the buffers
        preBuffer = currBuffer;
        currBufferID = (currBufferID+1)%2;
        currBuffer = buffer[currBufferID];
        
        // compute the starting positions for each bucket
        for (k=0; k<bucketNum; k++)
            numEachBucket[k] = 0;
        for (k=0; k<size; k++)
            numEachBucket[(intArrays[k]>>(i*digits))&mask]++;
        if (order) {
            // descending order
            startPosEachBucket[bucketNum-1]=0;
            for (k=bucketNum-2; k>=0; k--) {
                startPosEachBucket[k]=startPosEachBucket[k+1]+numEachBucket[k+1];
            }
        } else {
            // ascending order
            startPosEachBucket[0]=0;
            for (k=1; k<bucketNum; k++) {
                startPosEachBucket[k]=startPosEachBucket[k-1]+numEachBucket[k-1];
            }
        }
        
        // initialize the array of storing the number of characters coming across
        for (k=0; k<bucketNum; k++) {
            numEachBucketSoFar[k]=0;
        }
        
        // perform classification according to the corresponding value
        for (k=0; k<size; k++) {
            int index = preBuffer[k];
            int value = (intArrays[index]>>(i*digits))&mask;
            int pos = startPosEachBucket[value]+numEachBucketSoFar[value];
            currBuffer[pos] = index;
            numEachBucketSoFar[value]++;
        }
        
    }
    
    for (i=0; i<size; i++)
        sortedIndexID[i] = currBuffer[i];
    
    delete[] preBuffer;
    delete[] currBuffer;
    delete[] numEachBucket;
    delete[] startPosEachBucket;
    delete[] numEachBucketSoFar;
}

static void printSeq(string& str, ofstream *fout, int maxCols) {
    // print out the sequence with at most <maxCols> columns per line
    int i = 0;
    while (i < str.length()) {
        if (i > 0 && i%maxCols==0)
            (*fout) << endl;
        (*fout) << str[i];
        i++;
    }
    (*fout) << endl;
}

static void printSeq(string& str, ofstream *fout, int maxCols, vector<int>& pos_select, int base_unit) {
    // print out the sequence with at most <maxCols> columns per line
    int i = 0;
    int j,k;
    for (j=0; j<pos_select.size(); j++) {
        for (k=0; k<base_unit; k++) {
            if (i > 0 && i%maxCols==0)
                (*fout) << endl;
            (*fout) << str[pos_select[j]+k];
            i++;
        }
    }
    (*fout) << endl;
}

static string intToStr(int i) {
    // integer -> string
    
    string intToChr = "0123456789";
    
    if (i<0)
        return "-" + intToStr(-i);
    else if (i<10)
        return string(1,intToChr[i]);
    else
        return intToStr(i/10)+intToChr[i%10];
}

static string packedStr(string str, int len) {
	// return the string with length eactly "len"
	if (str.length() == len)
		return str;
	else if (str.length() > len)
		return str.substr(0,len);
	else {
		string s = str;
		s.append(len+1, ' ');
		return s.substr(0,len);
	}
		
}

#endif /* defined(__Alistat__tool_box__) */
