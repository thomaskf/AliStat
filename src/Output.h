/*
 *
 * Output.h
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


#ifndef __alistat__Output__
#define __alistat__Output__

#define SEPARATOR ","
#define SEPARATOR_TAB_4_5 "\t" // the separator used for table 4 and table 5

#define HEATMAP_BLOCK_DIM 5 // the dimension of the squares for heatmap
#define HEATMAP_BLOCK_DIM_FEW_SEQ 10 // the dimension of the squares for heatmap if the number of sequences is "FEW"
#define HEATMAP_FEW_SEQ_NUM 10 // the number of sequences is regarded as "FEW" if the number of sequences <= 10
#define HEATMAP_FONT "Arial" // Before: "Bitstream Vera Sans"
#define HEATMAP_FONT_SIZE 4 // the font size of the text for heatmap
#define HEATMAP_FONT_SIZE_FEW_SEQ 8 // the font size of the text for heatmap if the number of sequences is "FEW"
#define HEATMAP_TEXT_ALIGN "end" // the alignment of the text for heatmap
#define HEATMAP_LEGEND_BLOCK_DIM 10 // the dimension of the squares for the legend of the heatmap
#define HEATMAP_LEGEND_FONT "Arial" // Before: "Bitstream Vera Sans"
#define HEATMAP_LEGEND_FONT_SIZE 8 // the font size of the text for the legend of the heatmap
#define HEATMAP_LEGEND_TEXT_ALIGN "front" // the alignment of the text for the legend of the heatmap

#define HEATMAP_GRAY_COLOR "#cccccc" // gray color code for heatmap
#define HEATMAP_NAME_ANGLE 315 // angle of the sequence name for heatmap
#define HEATMAP_LEGEND_DESC_ANGLE 0 // angle of the description for the legend of the heatmap

#define HEATMAP_TOP_BORDER 5
#define HEATMAP_BOTTOM_BORDER 5
#define HEATMAP_LEFT_BORDER 5
#define HEATMAP_RIGHT_BORDER 5
#define HEATMAP_GAP_BW_TEXT 3
#define HEATMAP_GAP_BW_LEGEND 20
#define HEATMAP_LEGEND_GAP_BW_DESC 5

#define R_MESSAGE "To generate the charts by using our R scripts: \n\t $ R CMD BATCH <R script file>"

#include <iostream>
#include <string>
#include <vector>
#include "UserOptions.h"
#include "Version.h"

using namespace std;
static int HEATMAP_COLOR_NUM = 11;
static string HEATMAP_COLOR_MAP[] = {"#ef3b2c","#08306B","#08519C","#2171B5","#4292C6","#6BAED6","#9ECAE1","#C6DBEF","#DEEBF7","#F7FBFF","#FFFFFF"};
// static string HEATMAP_COLOR_MAP[] = {"#000000","#08306B","#08519C","#2171B5","#4292C6","#6BAED6","#9ECAE1","#C6DBEF","#DEEBF7","#F7FBFF","#FFFFFF"};
// static string HEATMAP_COLOR_MAP[] = {"#000000","#800026","#BD0026","#E31A1C","#FC4E2A","#FD8D3C","#FEB24C","#FED976","#FFEDA0","#FFFFCC","#FFFFFF"};
static string HEATMAP_COLOR_DESC[] = {"= 0.0", "&lt; 0.1", "&lt; 0.2", "&lt; 0.3", "&lt; 0.4", "&lt; 0.5", "&lt; 0.6", "&lt; 0.7", "&lt; 0.8", "&lt; 0.9", "&#x2264; 1.0"};

// output the usage of this program
void outputUsage(char* progName);

// Output the summary to the screen
void outSumToScreen(string seqFile, int seqNum, int seqLen, double Ca, double Cr_max, double Cr_min,
                   double Cc_max, double Cc_min, double Cij_max, double Cij_min, string partName, int showHeader);
                   
// Output the summary to the file <seqFile>.summary.txt
void outputSummary(string seqFile, string prefixOut, int dataType, int isCodon, char* validCharArr,
                   int seqNum, int seqLen, double Ca, double Cr_max, double Cr_min, double Cc_max,
                   double Cc_min, double Cij_max, double Cij_min,
                   double Pij_max, double Pij_min, double Pij_avg, bool Pij_exist,
                   string partName, UserOptions* userOptions);

// Output the C scores for individual sequences (Cr) (i.e. table 1)
// to the file <prefixOut>.table1.csv
void outputTable1(string prefixOut, int* Cr, int* row_index, vector<string>* seqNames, int seqNum, int partLen);

// Output the Distribution of C scores for individual sites (Cc) (i.e. table 2)
// to the file <prefixOut>.table2.csv
// format 0 : "Index\tSite ID\tCc"; format 1 : "Site ID\tCc"
void outputTable2(string prefixOut, int* Cc, int* col_index, int seqNum, int partLen, int format);

// Output the Summary of distribution of C scores for individual sites (Cc) (i.e. table 3)
// to the file <prefixOut>.table3.csv
void outputTable3(string prefixOut, int* Cc_part, int Cc_min, int Cc_max, int seqNum, int partLen);

// Output the Distribution of C scores for individual pairs of sequences (Cij) (i.e. table 4)
// to the file <prefixOut>.table4.csv
void outputTable4(string prefixOut, vector<string>* seqNames, int* Cij, int* row_index, int seqNum, int seqLen);

// Output the Incompleteness scores (Iij = 1.0 - Cij) (i.e. table 5)
// to the file <prefixOut>.table5.csv
void outputTable5(string prefixOut, vector<string>* seqNames, int* Cij, int* row_index, int seqNum, int seqLen);

// Output the Completeness scores (Cij) and Incompleteness scores (Iij = 1.0 - Cij) (i.e. table 6)
// to the file <prefixOut>.table6.csv
void outputTable6(string prefixOut, vector<string>* seqNames, int* Cij, int* row_index, int seqNum, int seqLen);

// Output the p-distances (Pij)
// to the file <prefixOut>.p-dist.csv
void outputTablePDist(string prefixOut, vector<string>* seqNames, int* Pij, int* Cij, int* row_index, int seqNum, int seqLen);

// Output the R script
void outputRScript(string prefixOut, UserOptions& userOption, string partName);

// Output the start message
void outputStartMessage(UserOptions &user_options, char* validCharArr, char* validNUChars, char* validAAChars);

// Output the finish message
void outputFinishMessage(string prefixOut, int seqNum, int seqLen, UserOptions* userOptions);

// Output the triangular heatmap
void outputTriHeatmap(string prefixOut, vector<string>* seqNames, int* Cij, int* row_index, int seqNum, int seqLen);

// Output the full heatmap
void outputFullHeatmap(string prefixOut, vector<string>* seqNames, int* Cij, int* row_index, int seqNum, int seqLen);

// The corresponding file name
string fileName(string prefixOut, string file);

// output the message about how to use the R scripts generated
void outputRMessage(UserOptions &user_options);

#endif /* defined(__alistat__Output__) */
