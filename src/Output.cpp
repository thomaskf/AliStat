/*
 *
 * Output.cpp
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


#include "Output.h"

// the updated usage of this program
void outputUsage(char* progName) {
    cout << "================================================================================" << endl;
    cout << "                        Welcome to Alistat - Version " << VERSION << endl << endl;
    
    cout << "Syntax: " << progName << " <alignment file> <data type> [other options]" << endl;
    cout << "        " << progName << " -h" << endl << endl;

    cout << "  <alignment file> : Multiple alignment file in FASTA format" << endl;
    
    /*
    cout << "  <data type>      : 1  - Nucleotide;" << endl;
    cout << "                     2  - Amino acid;" << endl;
    cout << "                     3  - Mixture of nucleotide and amino acid" << endl;
    */
    
    cout << "  <data type>      : 1  - Single nucleotides (SN);" << endl;
    cout << "                     2  - Di-nucleotides (DN);" << endl;
    cout << "                     3  - Codons (CD);" << endl;
    cout << "                     4  - 10-state genotype data (10GT);" << endl;
    cout << "                     5  - 14-state genotype data (14GT);" << endl;
    cout << "                     6  - Amino acids (AA);" << endl;
    cout << "                     7  - Mixture of nucleotides and amino acids (NA)" << endl;
    cout << "                          (User has to specify the data type for each partition" << endl;
    cout << "                           inside the partition file)" << endl;

    cout << "other options:" << endl;
    cout << "  -b               : Report the brief summary of figures to the screen" << endl;
    cout << "                     Output format:" << endl;
    cout << "                     File,#seqs,#sites,Ca,Cr_max,Cr_min,Cc_max,Cc_min,Cij_max,Cij_min" << endl;
    cout << "                     (this option cannot work with the options: -o,-t,-r,-m,-i,-d)" << endl;
    cout << "  -c <coding_type> : 0  - A, C, G, T [default];" << endl;
    cout << "                     1  - C, T, R (i.e. C, T, AG);" << endl;
    cout << "                     2  - A, G, Y (i.e. A, G, CT);" << endl;
    cout << "                     3  - A, T, S (i.e. A, T, CG);" << endl;
    cout << "                     4  - C, G, W (i.e. C, G, AT);" << endl;
    cout << "                     5  - A, C, K (i.e. A, C, GT);" << endl;
    cout << "                     6  - G, T, M (i.e. G, T, AC);" << endl;
    cout << "                     7  - K, M    (i.e. GT,   AC);" << endl;
    cout << "                     8  - S, W    (i.e. GC,   AT);" << endl;
    cout << "                     9  - R, Y    (i.e. AG,   CT);" << endl;
    cout << "                     10 - A, B    (i.e. A,   CGT);" << endl;
    cout << "                     11 - C, D    (i.e. C,   AGT);" << endl;
    cout << "                     12 - G, H    (i.e. G,   ACT);" << endl;
    cout << "                     13 - T, V    (i.e. T,   ACG);" << endl;
    cout << "                     (this option is only valid for <data type> = 1)" << endl;
    cout << "  -o <FILE>        : Prefix for output files" << endl;
    cout << "                     (default: <alignment file> w/o .ext)" << endl;
    cout << "  -n <FILE>        : Only consider sequences with names listed in FILE" << endl;
    cout << "  -p <FILE>        : Specify the partitions" << endl;
    cout << "                     For <data type> = 1 - 6, partition file format:" << endl;
    cout << "                       \"<partition name 1>=<start pos>-<end pos>, ...\"" << endl;
    cout << "                       \"<partition name 2>= ...\"" << endl;
    cout << "                       Example:" << endl;
    cout << "                           part1=1-50,60-100" << endl;
    cout << "                           part2=101-200" << endl;
    cout << "                       (enumeration starts with 1)" << endl;
    cout << "                     For <data type> = 7, partition file format:" << endl;
    cout << "                       \"<SN/DN/CD/10GT/14GT/AA>, <partition name 1>=<start pos>-<end pos>, ...\"" << endl;
    cout << "                       \"<SN/DN/CD/10GT/14GT/AA>, <partition name 2>= ...\"" << endl;
    cout << "                       Example:" << endl;
    cout << "                           SN,part1=1-50,60-100" << endl;
    cout << "                           AA,part2=101-200" << endl;
    cout << "                           CD,part3=201-231" << endl;
//    cout << "                     For codon data, in order to consider specific codon" << endl;
//    cout << "                     position:" << endl;
//    cout << "                       \"<partition name>=<start pos>-<end pos>, ... \\3\"" << endl;
//    cout << "                        or" << endl;
//    cout << "                       \"<CD>, <partition name>=<start pos>-<end pos>, ... \\3\"" << endl;
//    cout << "                        Example:" << endl;
//    cout << "                           CD,codon1=1-7\\3" << endl;
//    cout << "                           CD,codon2=2-8\\3" << endl;
//    cout << "                           CD,codon3=3-9\\3" << endl;
    cout << "                     (this option cannot be used with '-s' at the same time)" << endl;
    cout << "  -s <n1,n2>       : Sliding window analysis: window size = n1; step size = n2" << endl;
    cout << "                     (this option cannot be used with '-p' at the same time)" << endl;
    cout << "  -t <n3,n4,...>   : Only output the tables n3, n4, ..." << endl;
    cout << "                     1 - C scores for individual sequences (Cr)" << endl;
    cout << "                     2 - C scores for individual sites (Cc)" << endl;
    cout << "                     3 - Distribution of C scores for individual sites (Cc)" << endl;
    cout << "                     4 - Matrix with C scores for pairs of sequences (Cij)" << endl;
    cout << "                     5 - Matrix with incompleteness scores for pairs of" << endl;
	cout << "                         sequences (Iij = 1 - Cij)" << endl;
    cout << "                     6 - Table with C score and incompleteness scores for pairs" << endl;
    cout << "                         of sequences (Cij & Iij)" << endl;
    cout << "                     (default: the program does not output any tables)" << endl;
    cout << "                     If \"-t\" option is used but no <n3,n4,...>, then the program" << endl;
    cout << "                     outputs all tables" << endl;
    cout << "  -r <row|col|both>: Reorder the rows/columns (or both) of the alignment" << endl;
    cout << "                     according to the Cr/Cc scores" << endl;
    cout << "                     All the tables are displayed according to the reordered" << endl;
    cout << "                     alignment" << endl;
    cout << "                     To output the reordered alignment, please also use the" << endl;
	cout << "                     option -m" << endl;
	/*
    cout << "  -m <criteria>    : Keep the rows and/or the columns of alignment according to" << endl;
    cout << "                     <criteria>, and output the resulting alignment in the file" << endl;
	cout << "                     with extension '.mask.fst', and output the discarded alignment" << endl;
    cout << "                     in the file with extension '.disc.fst'" << endl;
	cout << "                     Example of <criteria>:" << endl;
	cout << "                       \"Cr > 0.8\"            : keep the seqs with Cr > 0.8" << endl;
	cout << "                       \"Cc >= 0.4\"           : keep the sites with Cr >= 0.4" << endl;
	cout << "                       \"Cr >= 0.9, Cc > 0.5\" : keep the seqs with Cr >= 0.9 and" << endl;
	cout << "                                               the sites with Cc > 0.5" << endl;
    cout << "                     (if no <criteria>, the whole alignment is outputted)" << endl;
	cout << "                     An algorithm for selecting sequences and sites based on Ca" << endl;
	cout << "                     score has been implemented. To use the algorithm, the user" << endl;
	cout << "                     should specify the Ca value and the weight on keeping the" << endl;
	cout << "                     sequences or sites (W, where 0 <= W <= 1). W = 0.5 means" << endl;
	cout << "                     equal weight assigned to seqs and sites, while greater W" << endl;
	cout << "                     implies higher weight on keeping the sequences." << endl; 
	cout << "                     For example:" << endl;
	cout << "                       \"Ca >= 0.6, W = 0.8\"  : Get a sub-alignment (with higher" << endl;
	cout << "                                               weight to keep seqs) such that" << endl;
	cout << "                                               Ca >= 0.6" << endl;
	cout << "                     If no W is given, the program maximizes the size of the" << endl;
	cout << "                     sub-alignment in order to meet the Ca requirement" << endl;
	cout << "                     For example:" << endl;
	cout << "                       \"Ca >= 0.7\"           : Get a sub-alignment such that" << endl;
	cout << "                                               Ca >= 0.7 and its size is" << endl;
	cout << "                                               maximal." << endl;
	 */
    cout << "  -m <n5>          : Mask the alignment; 0 <= n5 <= 1" << endl;
    cout << "                     Output (1) the alignment with columns Cc >= n5 in the file" << endl;
    cout << "                     'Mask.fst', (2) the alignment with columns Cc < n5 in the" << endl;
    cout << "                     file 'Disc.fst', and (3) the alignment with an extra row in" << endl;
    cout << "                     the first line to indicate whether the column is masked in" << endl;
    cout << "                     the file 'Stat.fst'." << endl;
    cout << "                     (Special case: if no <n5>, whole alignment is outputted in" << endl;
    cout << "                      the file 'Mask.fst'; if <n5> is 0, the alignment with" << endl;
	cout << "                      columns Cc > 0 is outputted in the file 'Mask.fst')" << endl;
	cout << "  -i <1|2|3>       : Generate heat map image for Cij scores of sequence pairs" << endl;
	cout << "                     1 - Triangular heat map" << endl;
	cout << "                     2 - Rectangular heat map" << endl;
	cout << "                     3 - Both" << endl;
	cout << "                     (if no number, then both triangular & rectangular heat map" << endl;
	cout << "                      files are outputted)" << endl;
	cout << "  -d               : Report the p distances between sequences" << endl;
	cout << "                     in the file with extension '.p-dist.csv'" << endl;
    cout << "                     (default: disabled)" << endl;
	cout << "                     Note: computation of p distances may take long time for" << endl;
	cout << "                          large number of sequences" << endl;
    cout << "  -u               : Color scheme of the heatmaps" << endl;
    cout << "                     1 - Default color scheme, suitable for color-blind persons" << endl;
    cout << "                     2 - Another color scheme" << endl;
	cout << "  -h               : This help page" << endl;
    cout << "================================================================================" << endl;
}


string fileName(string prefixOut, string file) {
    if (prefixOut.length()>0)
        return prefixOut + "." + file;
    else
        return file;
}

// Output the summary to the screen
void outSumToScreen(string seqFile, int seqNum, int seqLen, double Ca, double Cr_max, double Cr_min,
                   double Cc_max, double Cc_min, double Cij_max, double Cij_min, string partName, int showHeader) {

	if (showHeader) {
		if (partName.length() == 0)
			cout << "File_name, #sequences, #sites, Ca, Cr_max, Cr_min, Cc_max, Cc_min, Cij_max, Cij_min" << endl;
		else
			cout << "File_name, Partition_name, #sequences, #sites, Ca, Cr_max, Cr_min, Cc_max, Cc_min, Cij_max, Cij_min" << endl;
	}
	                   
    cout << seqFile;
    
    if (partName.length() > 0)
    	cout << ", " << partName;

    cout << ", " << seqNum;
    cout << ", " << seqLen;
    cout << ", " << (double)Ca/(seqNum*seqLen);
    cout << ", " << (double)Cr_max/seqLen;
    cout << ", " << (double)Cr_min/seqLen;
    cout << ", " << (double)Cc_max/seqNum;
    cout << ", " << (double)Cc_min/seqNum;
    cout << ", " << (double)Cij_max/seqLen;
    cout << ", " << (double)Cij_min/seqLen;
    cout << endl;
                   
}

// Output the summary to the file <prefixOut>.summary.txt
void outputSummary(string seqFile, string prefixOut, int dataType, int isCodon, char* validCharArr,
                   int seqNum, int seqLen, double Ca, double Cr_max, double Cr_min, double Cc_max,
                   double Cc_min, double Cij_max, double Cij_min,
                   double Pij_max, double Pij_min, double Pij_avg, bool Pij_exist,
                   string partName, UserOptions* userOptions, string dataTypeStr) {
    
    int i;
    
    string outFileName = fileName(prefixOut, "Summary.txt");
    ofstream fout;
    string term = "alignment   ";

    fout.open(outFileName.c_str());
    fout << "Name of input file ............................................ " << seqFile << endl;
    if (userOptions->slideWindowSize != -1) {
        fout << "Sliding window positions ...................................... " << partName << endl;
        term = "slide window";
    } else if (partName != "") {
        fout << "Partition name ................................................ " << partName << endl;
        term = "partition   ";
    }

    fout << "Type of data .................................................. " << dataTypeStr << endl;
    if (dataType == 1 || dataType == 4 || dataType == 5 || dataType == 6) {
        fout << "Characters used ............................................... ";
        for (i=0; i<256; i++) {
            if (validCharArr[i]) {
                fout << (char) i;
            }
        }
        fout << endl;
    }
    
    fout << "Number of sequences in the "<< term << " ....................... " << seqNum << endl;
    if (dataType==3)
    	fout << "Number of codon sites in the " << term << " ..................... " << seqLen << endl;
    else
    	fout << "Number of sites in the " << term << " ........................... " << seqLen << endl;
    fout << "Pairs of sequences ............................................ " << seqNum * (seqNum-1) / 2 << endl;
    fout << "Completeness (C) score for the " << term << " (Ca) .............. " << (double)Ca/(seqNum*seqLen) << endl;
    fout << "Maximum C-score for individual sequences (Cr_max) ............. " << (double)Cr_max/seqLen << endl;
    fout << "Minimum C-score for individual sequences (Cr_min) ............. " << (double)Cr_min/seqLen << endl;
    fout << "Maximum C-score for individual sites (Cc_max) ................. " << (double)Cc_max/seqNum << endl;
    fout << "Minimum C-score for individual sites (Cc_min) ................. " << (double)Cc_min/seqNum << endl;
    
    if (Cij_min != -1 && Cij_max != -1) {
        
        fout << "Maximum C-score for pairs of sequences (Cij_max, i!=j) ........ " << (double)Cij_max/seqLen << endl;
        fout << "Minimum C-score for pairs of sequences (Cij_min, i!=j) ........ " << (double)Cij_min/seqLen << endl;
        
    }

	if (userOptions->computePDist==1) {
        
        if (Pij_exist) {
            fout << "Maximum P-dist between sequences (Pij_max) .................... " << Pij_max << endl;
            fout << "Minimum P-dist between sequences (Pij_min) .................... " << Pij_min << endl;
            fout << "Average P-dist between sequences (Pij_avg) .................... " << Pij_avg << endl;
        } else {
            fout << "Maximum P-dist between sequences (Pij_max) .................... NA" << endl;
            fout << "Minimum P-dist between sequences (Pij_min) .................... NA" << endl;
            fout << "Average P-dist between sequences (Pij_avg) .................... NA" << endl;
        }
        
    }
    
    fout << endl;
    
    for (i=0; i<userOptions->outTables.size(); i++) {
        switch(userOptions->outTables[i]) {
            case 1:
                fout << "C scores for individual sequences (Cr) stored in .............. " << fileName(prefixOut, "Table_1.csv") << endl;
                break;
            case 2:
                fout << "C scores for individual sites (Cc) stored in .................. " << fileName(prefixOut, "Table_2.csv") << endl;
                break;
            case 3:
                fout << "Distribution of C scores for individual sites (Cc) in ......... " << fileName(prefixOut, "Table_3.csv") << endl;
                break;
            case 4:
                if (Cij_min != -1 && Cij_max != -1) {
                    fout << "Matrix with C scores for individual pairs of seqs (Cij) in .... " << fileName(prefixOut, "Table_4.dis") << endl;
                }
                break;
            case 5:
                if (Cij_min != -1 && Cij_max != -1) {
                    fout << "Matrix with incompleteness scores (Iij = 1.0 - Cij) in ........ " << fileName(prefixOut, "Table_5.dis") << endl;
                }
                break;
            case 6:
                if (Cij_min != -1 && Cij_max != -1) {
                    fout << "Table with C scores (Cij) and incompleteness scores (Iij) in .. " << fileName(prefixOut, "Table_6.csv") << endl;
                }
                break;
        }
    }
	
	if (userOptions->makeHeatMap==1 || userOptions->makeHeatMap==3)
		fout << "Triangular heatmap for individual pairs of seqs (Cij) in ...... " << fileName(prefixOut, "Triangular.svg") << endl;
	
	if (userOptions->makeHeatMap==2 || userOptions->makeHeatMap==3)
        fout << "Rectangular heatmap for individual pairs of seqs (Cij) in ..... " << fileName(prefixOut, "Full.svg") << endl;
    
    if (userOptions->outputAlign) {
        fout << "The masked alignment (Cc >= " << userOptions->maskThres << ") in ........................... " << fileName(prefixOut, "Mask.fst") << endl;
        fout << "The alignment of columns with low Cc values (Cc < " << userOptions->maskThres << ") in ..... " << fileName(prefixOut, "Disc.fst") << endl;
        fout << "The file indicating which columns are masked in ............... " << fileName(prefixOut, "Stat.fst") << endl;
    }
	
	if (userOptions->computePDist==1) {
        fout << "P-distance between sequences in ............................... " << fileName(prefixOut, "P-dist.csv") << endl;
    }
	
    // showing the R-script files
    for (i=0; i<userOptions->outTables.size(); i++) {
        switch(userOptions->outTables[i]) {
            case 1:
                fout << "R script for generating histogram for Cr distribution ......... " << fileName(prefixOut, "Histogram_Cr.R") << endl;
                break;
            case 2:
                fout << "R script for generating histogram for Cc distribution ......... " << fileName(prefixOut, "Histogram_Cc.R") << endl;
                break;
            case 3:
                fout << "R script for generating cumulative curve for Cc distribution .. " << fileName(prefixOut, "Cumulative_Cc.R") << endl;
                break;
        }
    }

    fout.close();
}

// Output the C scores for individual sequences (Cr) (i.e. table 1)
// to the file <prefixOut>.table1.csv
void outputTable1(string prefixOut, int* Cr, int* row_index, vector<string>* seqNames, int seqNum, int partLen) {
    
    int i;
    string outFileName = fileName(prefixOut, "Table_1.csv");
    ofstream fout;
    fout.open(outFileName.c_str());
    /*
    fout << "TABLE 1: C scores for individual sequences (Cr)" << endl;
    fout << "(Proportion of sites with unambiguous characters in a given sequence)" << endl << endl;
     */ 
    fout << "Seq ID" << SEPARATOR << "Sequence" << SEPARATOR << "Valid sites" << SEPARATOR << "Cr" << endl;
    for (i=0; i<seqNum; i++) {
        fout << row_index[i]+1 << SEPARATOR << seqNames->at(row_index[i]) << SEPARATOR << Cr[row_index[i]] << SEPARATOR;
        fout << (double) Cr[row_index[i]] / partLen << endl;
    }
    fout.close();
    
}

// Output the Distribution of C scores for individual sites (Cc) (i.e. table 2)
// to the file <prefixOut>.table2.csv
// format 0 : "Index\tSite ID\tCc"; format 1 : "Site ID\tCc"
void outputTable2(string prefixOut, int* Cc, int* col_index, int seqNum, int partLen, int format) {
    
    int i;
    string outFileName = fileName(prefixOut, "Table_2.csv");
    ofstream fout;
    fout.open(outFileName.c_str());
    /*
    fout << "TABLE 2: C scores for individual sites (Cc)" << endl << endl;
     */ 
    if (format) {
        fout << "Site ID" << SEPARATOR << "Cc" << endl;
        for (i=0; i<partLen; i++) {
            fout << col_index[i]+1 << SEPARATOR << (double) Cc[col_index[i]] / seqNum << endl;
        }
    } else {
        fout << "Index" << SEPARATOR << "Site ID" << SEPARATOR << "Cc" << endl;
        for (i=0; i<partLen; i++) {
            fout << i << SEPARATOR << col_index[i]+1 << SEPARATOR << (double) Cc[col_index[i]] / seqNum << endl;
        }
    }
    fout.close();
    
}

// Output the Summary of distribution of C scores for individual sites (Cc) (i.e. table 3)
// to the file <prefixOut>.table3.csv
void outputTable3(string prefixOut, int* Cc_part, int Cc_min, int Cc_max, int seqNum, int partLen) {
    
    int i;

    // To get the distribution of Cc scores
    int* distribute;
    distribute = getDist(Cc_part, partLen, Cc_min, Cc_max);

    if (distribute != NULL) {
    
        string outFileName = fileName(prefixOut, "Table_3.csv");
        ofstream fout;
        fout.open(outFileName.c_str());
        /*
        fout << "TABLE 3: Distribution of C scores for individual sites (Cc)" << endl;
        fout << "(Proportion of sequences with unambiguous characters at a given site." << endl;
        fout << "Here X represents the number of ambiguous characters at a given site)" << endl << endl;
         */
        fout << "X" << SEPARATOR << "Cc" << SEPARATOR << "Sites" << SEPARATOR << "Proportion" << SEPARATOR << "Cumulative" << endl;
        double cum = 0.0;
        for (i=seqNum; i>=Cc_max+1; i--) {
            fout << seqNum-i << SEPARATOR << (double)i/seqNum  << SEPARATOR << 0 << SEPARATOR << 0 << SEPARATOR << 0 << endl;
        }
        for (i=Cc_max; i>=Cc_min; i--) {
            cum += (double) distribute[i-Cc_min]/partLen;
            fout << seqNum-i << SEPARATOR << (double)i/seqNum  << SEPARATOR << distribute[i-Cc_min] << SEPARATOR << (double) distribute[i-Cc_min]/partLen << SEPARATOR << cum << endl;
        }
        for (i=Cc_min-1; i>=0; i--) {
            fout << seqNum-i << SEPARATOR << (double)i/seqNum  << SEPARATOR << 0 << SEPARATOR << 0 << SEPARATOR << cum << endl;
        }
        fout.close();
        
    }
    
    delete[] distribute;
    
}

// Output the Distribution of C scores for individual pairs of sequences (Cij) (i.e. table 4)
// to the file <prefixOut>.table4.csv
void outputTable4(string prefixOut, vector<string>* seqNames, int* Cij, int* row_index, int seqNum, int seqLen) {
    
    int i,j;
    int i_index,j_index;
    string outFileName = fileName(prefixOut, "Table_4.dis");
    ofstream fout;
    fout.open(outFileName.c_str());
    
    /*
     fout << "TABLE 4: C scores for individual pairs of sequences (Cij)" << endl;
     fout << "(Proportion of homologous sites with unambiguous characters" << endl;
     fout << "for pairs of sequences)" << endl << endl;
     
     fout << "Seq. 1.\tSeq. 2.\tSites\tCij" << endl;
     for (i=0; i<seqNum-1; i++)
     for (j=i+1; j<seqNum; j++)
     fout << seqNames->at(i) << "\t" << seqNames->at(j) << "\t" << CijMatrix->at(i*seqNum+j) << "\t" << (double)CijMatrix->at(i*seqNum+j)/seqLen << endl;
     */
    
    fout << seqNum << endl;
    
    for (i=0; i<seqNum; i++) {
        i_index = row_index[i];
        fout << packedStr(seqNames->at(i_index), 10);
        for (j=0; j<seqNum; j++) {
            j_index = row_index[j];
            fout << SEPARATOR_TAB_4_5 << (double)Cij[i_index*seqNum+j_index]/seqLen;
        }
        fout << endl;
    }
    
    fout.close();
    
}

// Output the Incompleteness scores (Iij = 1.0 - Cij) (i.e. table 5)
// to the file <prefixOut>.table5.csv
void outputTable5(string prefixOut, vector<string>* seqNames, int* Cij, int* row_index, int seqNum, int seqLen) {
    
    int i,j;
    int i_index,j_index;
    string outFileName = fileName(prefixOut, "Table_5.dis");
    ofstream fout;
    fout.open(outFileName.c_str());
    // fout << "TABLE 5: Incompleteness scores (Iij = 1.0 - Cij)" << endl;
    
    fout << seqNum << endl;
    
    for (i=0; i<seqNum; i++) {
        i_index = row_index[i];
        fout << packedStr(seqNames->at(i_index), 10);
        for (j=0; j<seqNum; j++) {
            j_index = row_index[j];
            fout << SEPARATOR_TAB_4_5 << 1.0 - (double)Cij[i_index*seqNum+j_index]/seqLen;
        }
        fout << endl;
    }
    
    fout.close();
    
}

// Output the completeness scores (Cij) and incompleteness scores (Iij = 1.0 - Cij) (i.e. table 6)
// to the file <prefixOut>.table6.csv
void outputTable6(string prefixOut, vector<string>* seqNames, int* Cij, int* row_index, int seqNum, int seqLen) {
    
    int i,j;
    int i_index,j_index;
    string outFileName = fileName(prefixOut, "Table_6.csv");
    ofstream fout;
    fout.open(outFileName.c_str());
    double cij;
    // fout << "TABLE 6: completeness scores (Cij) and incompleteness scores (Iij = 1.0 - Cij)" << endl;

    fout << "taxon 1" << SEPARATOR << "taxon 2" << SEPARATOR << "Cij" << SEPARATOR <<"Iij" << endl;
    for (i=0; i<seqNum-1; i++) {
        i_index = row_index[i];
        for (j=i+1; j<seqNum; j++) {
            j_index = row_index[j];
            cij = (double)Cij[i_index*seqNum+j_index]/seqLen;
            fout << seqNames->at(i_index) << SEPARATOR << seqNames->at(j_index) << SEPARATOR << cij << SEPARATOR << 1.0-cij << endl;
        }
    }
    
    fout.close();
    
}

// Output the p-distances (Pij)
// to the file <prefixOut>.p-dist.csv
void outputTablePDist(string prefixOut, vector<string>* seqNames, int* Pij, int* Cij, int* row_index, int seqNum, int seqLen) {
    
    int i,j;
    int i_index,j_index;
    string outFileName = fileName(prefixOut, "P-dist.csv");
    ofstream fout;
    fout.open(outFileName.c_str());
    // fout << "P Distances between sequences" << endl;
    
    fout << seqNum << endl;
    
    for (i=0; i<seqNum; i++) {
        i_index = row_index[i];
        fout << seqNames->at(i_index);
        for (j=0; j<seqNum; j++) {
            j_index = row_index[j];
			if (Cij[i_index*seqNum+j_index] > 0)
				fout << SEPARATOR << (double)Pij[i_index*seqNum+j_index]/Cij[i_index*seqNum+j_index];
			else
				fout << SEPARATOR << "NA";
        }
        fout << endl;
    }
    
    fout.close();
    
}


// Output the R script
void outputRScript(string prefixOut, UserOptions& user_options, string partName) {
    
    int i;
    for (i=0; i<user_options.outTables.size(); i++) {
        
        if (user_options.outTables[i] == 1) {
            // generate R script for Table_1.csv
            string outFileName = fileName(prefixOut, "Histogram_Cr.R");
            string tableFileName = fileName(prefixOut, "Table_1.csv");
            string outPdfName = fileName(prefixOut, "Histogram_Cr.pdf");
            string chartTitle;
            if (partName.length() >0) {
                chartTitle = "Distribution of Cr\\nin partition '" + partName + "'";
            } else {
                chartTitle = "Distribution of Cr";
            }
            ofstream fout;
            fout.open(outFileName.c_str());
            fout << "d <- read.csv(\"" << tableFileName << "\", header=T)" << endl;
            fout << "head(d)" << endl;
            fout << "pdf(height=4, width=3.5, file=\"" << outPdfName << "\")" << endl;
            fout << "cbks <- seq(0,1,0.05)" << endl;
            fout << "hist(d$Cr, breaks = cbks, col=\"lightblue\", freq = TRUE, xlim=c(0,1), xlab=\"Cr\", main=\"" << chartTitle << "\")" << endl;
            fout << "dev.off()" << endl;
            fout.close();
            
        } else if (user_options.outTables[i] == 2) {
            // generate R script for Table_2.csv
            string outFileName = fileName(prefixOut, "Histogram_Cc.R");
            string tableFileName = fileName(prefixOut, "Table_2.csv");
            string outPdfName = fileName(prefixOut, "Histogram_Cc.pdf");
            string chartTitle;
            if (partName.length() >0) {
                chartTitle = "Distribution of Cc\\nin partition '" + partName + "'";
            } else {
                chartTitle = "Distribution of Cc";
            }
            ofstream fout;
            fout.open(outFileName.c_str());
            fout << "d <- read.csv(\"" << tableFileName << "\", header=T)" << endl;
            fout << "head(d)" << endl;
            fout << "pdf(height=4, width=3.5, file=\"" << outPdfName << "\")" << endl;
            fout << "cbks <- seq(0,1,0.05)" << endl;
            fout << "hist(d$Cc, breaks = cbks, col=\"lightblue\", freq = TRUE, xlim=c(0,1), xlab=\"Cc\", main=\"" << chartTitle << "\")" << endl;
            fout << "dev.off()" << endl;
            fout.close();
            
        } else if (user_options.outTables[i] == 3) {
            // generate R script for Table_3.csv
            string outFileName = fileName(prefixOut, "Cumulative_Cc.R");
            string tableFileName = fileName(prefixOut, "Table_3.csv");
            string outPdfName = fileName(prefixOut, "Cumulative_Cc.pdf");
            string chartTitle;
            if (partName.length() >0) {
                chartTitle = "Cumulative distribution of Cc\\nin partition '" + partName + "'";
            } else {
                chartTitle = "Cumulative distribution of Cc";
            }
            ofstream fout;
            fout.open(outFileName.c_str());
            fout << "d <- read.csv(\"" << tableFileName << "\", header=T)" << endl;
            fout << "head(d)" << endl;
            fout << "pdf(height=4, width=3.5, file=\"" << outPdfName << "\")" << endl;
            fout << "plot(d$Cc, d$Cumulative, ylab=\"Cumulative proportion\", xlab=\"Cc\", main=\"" << chartTitle << "\", xlim=c(0,1), ylim=c(0,1), pch=20, col=\"black\", cex=0.8)" << endl;
            fout << "abline(h=1.00, lty=2, lwd=0.2)" << endl;
            fout << "abline(h=0.90, lty=2, lwd=0.2)" << endl;
            fout << "abline(h=0.75, lty=2, lwd=0.2)" << endl;
            fout << "abline(h=0.50, lty=2, lwd=0.2)" << endl;
            fout << "abline(h=0.25, lty=2, lwd=0.2)" << endl;
            fout << "abline(h=0.10, lty=2, lwd=0.2)" << endl;
            fout << "abline(h=0.00, lty=2, lwd=0.2)" << endl;
            fout << "dev.off()" << endl;
            fout.close();
        }
    }
}

// Output the start message
void outputStartMessage(UserOptions &user_options, char* validCharArr, char* validSNChars, char* validDNChars, char* validCDChars, char* valid10GTChars, char* valid14GTChars, char* validAAChars) {
    
    int i;
    
    cout << "================================================================================" << endl;
    cout << "                        Welcome to Alistat - Version " << VERSION << endl << endl;
    
    cout << "Name of input file ............................................ " << user_options.alignFile << endl;
    
    cout << "Type of data .................................................. " << user_options.dataTypeStr << endl;
    
    if (user_options.dataType==7) {
        cout << "Characters used for single nucleotide data .................... ";
        for (i=0; i<256; i++) {
            if (validSNChars[i]) {
                cout << (char) i;
            }
        }
        cout << endl;
        cout << "Characters used for 10-state genotype data .................... ";
        for (i=0; i<256; i++) {
            if (valid10GTChars[i]) {
                cout << (char) i;
            }
        }
        cout << endl;
        cout << "Characters used for 14-state genotype data .................... ";
        for (i=0; i<256; i++) {
            if (valid14GTChars[i]) {
                cout << (char) i;
            }
        }
        cout << endl;
        cout << "Characters used for amino acid data ........................... ";
        for (i=0; i<256; i++) {
            if (validAAChars[i]) {
                cout << (char) i;
            }
        }
        cout << endl;
    } else if (user_options.dataType < 2 || user_options.dataType > 3) {
        cout << "Characters used ............................................... ";
        for (i=0; i<256; i++) {
            if (validCharArr[i]) {
                cout << (char) i;
            }
        }
        cout << endl;
    } else if (user_options.dataType == 3) {
        cout << "Characters used ............................................... AAA,...,TTT/UUU" << endl;
    } else {
        cout << "Characters used ............................................... ACGTU" << endl;
    }

    if (user_options.prefixOut.length()>0)
        cout << "Prefix for output files ....................................... " << user_options.prefixOut << endl;
    
    if (user_options.seqNameFile!="") {
        cout << "Only consider the sequences with names inside the file ........ " << user_options.seqNameFile << endl;
    }
    
    if (user_options.partitionFile!="") {
        cout << "Partition file ................................................ " << user_options.partitionFile << endl;
    }
    
    if (user_options.slideWindowSize>0 && user_options.slideWindowStep>0) {
        cout << "Sliding window size ........................................... " << user_options.slideWindowSize << endl;
        cout << "Sliding window step size ...................................... " << user_options.slideWindowStep << endl;
    }
    
    if (user_options.outTables.size()>0) {
        cout << "Output the tables ............................................. ";
        for (i=0; i<user_options.outTables.size(); i++) {
            if (i>0)
                cout << ", ";
            cout << user_options.outTables[i];
        }
        cout << endl;
    }
    
    if (user_options.reorder>0) {
        cout << "Rearrangement of the alignment ................................ ";
        switch(user_options.reorder) {
            case 1:
                cout << "Row";
                break;
            case 2:
                cout << "Column";
                break;
            case 3:
                cout << "Row & column";
                break;
        }
        cout << endl;
    }
    
    if (user_options.maskThres>=0.0) {
        cout << "Cc Threshold for masking of the alignment...................... " << user_options.maskThres << endl;
    }
    
	cout << "Computation of P-distance between sequences ................... ";
	switch (user_options.computePDist) {
	    case 0:
			cout << "Not enabled";
			break;
		case 1:
			cout << "Enabled";
			break;
	}
	cout << endl;
    
    cout << "================================================================================" << endl;
    
}


// Output the finish message
void outputFinishMessage(string prefixOut, int seqNum, int seqLen, UserOptions* userOptions) {
    
    int i;
    
    string term = "alignment ...";
    
    if (userOptions->slideWindowSize != -1) {
        term = "slide window ";
    } else if (userOptions->partitionFile != "") {
        term = "partition ...";
    }

    cout << "================================================================================" << endl;
    cout << "Number of sequences in the " << term << "....................... " << seqNum << endl;
    if (userOptions->dataType == 3)
    	cout << "Number of codon sites in the " << term << "..................... " << seqLen/3 << endl;
    else
    	cout << "Number of sites in the " << term << "........................... " << seqLen << endl;
    cout << "Summary of statistics stored in ............................... " << fileName(prefixOut, "Summary.txt") << endl;
    
    for (i=0; i<userOptions->outTables.size(); i++) {
        switch(userOptions->outTables[i]) {
            case 1:
                cout << "C scores for individual sequences (Cr) stored in .............. " << fileName(prefixOut, "Table_1.csv") << endl;
                break;
            case 2:
                cout << "C scores for individual sites (Cc) stored in .................. " << fileName(prefixOut, "Table_2.csv") << endl;
                break;
            case 3:
                cout << "Distribution of C scores for individual sites (Cc) in ......... " << fileName(prefixOut, "Table_3.csv") << endl;
                break;
            case 4:
                if (seqNum > 1) {
                    cout << "Matrix with C scores for individual pairs of seqs (Cij) in .... " << fileName(prefixOut, "Table_4.dis") << endl;
                }
                break;
            case 5:
                if (seqNum > 1) {
                    cout << "Matrix with incompleteness scores (Iij = 1.0 - Cij) in ........ " << fileName(prefixOut, "Table_5.dis") << endl;
                }
                break;
            case 6:
                if (seqNum > 1) {
                    cout << "Table with C scores (Cij) and incompleteness scores (Iij) in .. " << fileName(prefixOut, "Table_6.csv") << endl;
                }
                break;
        }
    }

    if (userOptions->makeHeatMap==1 || userOptions->makeHeatMap==3)
        cout << "Triangular heatmap for individual pairs of seqs (Cij) in ...... " << fileName(prefixOut, "Triangular.svg") << endl;

    if (userOptions->makeHeatMap==2 || userOptions->makeHeatMap==3)
        cout << "Rectangular heatmap for individual pairs of seqs (Cij) in ..... " << fileName(prefixOut, "Full.svg") << endl;

    if (userOptions->outputAlign) {
        cout << "The masked alignment in ....................................... " << fileName(prefixOut, "Mask.fst") << endl;
        cout << "The alignment of columns with low Cc values in ................ " << fileName(prefixOut, "Disc.fst") << endl;
        cout << "The file indicating which columns are masked in ............... " << fileName(prefixOut, "Stat.fst") << endl;
    }
    
    if (userOptions->computePDist==1) {
        cout << "The p-distances between sequences in .......................... " << fileName(prefixOut, "P-dist.csv") << endl;
    }
    
    // showing the R-script files
    for (i=0; i<userOptions->outTables.size(); i++) {
        switch(userOptions->outTables[i]) {
            case 1:
                cout << "R script for generating histogram for Cr distribution ......... " << fileName(prefixOut, "Histogram_Cr.R") << endl;
                break;
            case 2:
                cout << "R script for generating histogram for Cc distribution ......... " << fileName(prefixOut, "Histogram_Cc.R") << endl;
                break;
            case 3:
                cout << "R script for generating cumulative curve for Cc distribution .. " << fileName(prefixOut, "Cumulative_Cc.R") << endl;
                break;
        }
    }

    cout << "================================================================================" << endl;
    
}

// Output the SVG Header
void outputSVGHeader(ofstream& fout, int width, int height) {
    fout << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>" << endl;
    fout << "<!-- created by Alistat version " << VERSION << " -->" << endl;
    fout << "<svg" << endl;
    fout << "\txmlns:svg=\"http://www.w3.org/2000/svg\"" << endl;
    fout << "\txmlns=\"http://www.w3.org/2000/svg\"" << endl;
    fout << "\tversion=\"1.0\"" << endl;
    fout << "\twidth=\"" << width << "\"" << endl;
    fout << "\theight=\"" << height << "\"" << endl;
    fout << "\tid=\"svg2\">" << endl;
    fout << "\t<defs id=\"defs4\" />" << endl;
}


void outputSVGTail(ofstream& fout) {
    fout << "</svg>" << endl;
}
// Output the square
void outputSquare(ofstream& fout, double x, double y, string color, int size) {
    fout << "<rect" << endl;
    fout << "\twidth=\"" << size << "\"" << endl;
    fout << "\theight=\"" << size << "\"" << endl;
    fout << "\tx=\"" << x << "\"" << endl;
    fout << "\ty=\"" << y << "\"" << endl;
//    fout << "\tstyle=\"fill:" << color << ";fill-opacity:1;stroke:black;stroke-width:0.1;stroke-miterlimit:10;stroke-dasharray:none;stroke-opacity:1\" />" << endl;
    // LS Jermiin (17 Feb 2026) - replaced the line above with the line below
    fout << "\tstyle=\"fill:" << color << ";stroke-width:0.1;stroke:black\" />" << endl;
}

// output the sequence name
void outputSeqName(ofstream& fout, double x, double y, int angle, string txt, string font, int fontSize, string text_align) {
    fout << "<text" << endl;
    fout << "\tx=\"" << x << "\" y=\"" << y << "\"" << endl;
//    fout << "\tstyle=\"font-size:" << fontSize << "px;font-style:normal;font-weight:normal;text-align:" << text_align << ";text-anchor:" << text_align << ";fill:black;fill-opacity:1;stroke:none;stroke-width:1px;stroke-linecap:butt;stroke-linejoin:miter;stroke-opacity:1;font-family:" << font << "\"" << endl;
    // LS Jermiin (17 Feb 2026) - replaced the line above with the four lines below
    fout << "\tfont-family=\"monospace\"" << endl;
    if (angle > 0)
        fout << "\ttransform=\"rotate(" << angle << " " << x << " " << y << ")\"" << endl;
    fout << "\t>" << txt << "</text>" << endl;
}

// get the corresponding color code according to the value
string getColor(double value, int color_scheme) {
    // 0 : 0.0
    // 1 : (0.0, 0.1)
    // 2 : [0.1, 0.2)
    // 3 : [0.2, 0.3)
    // 4 : [0.3, 0.4)
    // 5 : [0.4, 0.5)
    // 6 : [0.5, 0.6)
    // 7 : [0.6, 0.7)
    // 8 : [0.7, 0.8)
    // 9 : [0.8, 0.9)
    // 10: [0.9, 1.0]
    
    switch (color_scheme) {
        case 2:
            if (value <= 0.0)
                return HEATMAP_COLOR_MAP2[0];
            else if (value >= 1.0)
                return HEATMAP_COLOR_MAP2[10];
            else
                return HEATMAP_COLOR_MAP2[(int)(value*10.0) + 1];
            break;
        default:
            if (value <= 0.0)
                return HEATMAP_COLOR_MAP1[0];
            else if (value >= 1.0)
                return HEATMAP_COLOR_MAP1[10];
            else
                return HEATMAP_COLOR_MAP1[(int)(value*10.0) + 1];
            break;
    }
}

// get the longest length among the first "num" strings inside the array
int longestLenFirst(vector<string>* seqArray, int* row_index, int num) {
    if (seqArray->size() == 0)
        return 0;
    int max_len = seqArray->at(row_index[0]).length();
    int i;
    for (i=1; i<num && i<seqArray->size(); i++)
        if (seqArray->at(row_index[i]).length() > max_len)
            max_len = seqArray->at(row_index[i]).length();
    return max_len * LENGTH_RATIO;
}

// get the longest length among the last "num" strings inside the array
int longestLenLast(vector<string>* seqArray, int* row_index, int num) {
    if (seqArray->size() == 0)
        return 0;
    int max_len = seqArray->at(row_index[seqArray->size()-1]).length();
    int i;
    for (i=seqArray->size()-2; i>=seqArray->size()-num && i>=0; i--)
        if (seqArray->at(row_index[i]).length() > max_len)
            max_len = seqArray->at(row_index[i]).length();
    return max_len * LENGTH_RATIO;
}

// get the longest length among all the strings inside the array
int longestLen(vector<string>* seqArray) {
    if (seqArray->size() == 0)
        return 0;
    int max_len = seqArray->at(0).length();
    int i;
    for (i=1; i<seqArray->size(); i++)
        if (seqArray->at(i).length() > max_len)
            max_len = seqArray->at(i).length();
    return max_len * LENGTH_RATIO;
}

// output the legend
void outputLegend(ofstream& fout, int x, int y, int color_scheme) {

    double x_text, y_text;
    x_text = x + HEATMAP_LEGEND_BLOCK_DIM + HEATMAP_LEGEND_GAP_BW_DESC;
    y_text = y + HEATMAP_LEGEND_BLOCK_DIM*0.8;
    int i;
    for (i=0; i<HEATMAP_COLOR_NUM; i++) {
        // print out the square
        switch(color_scheme) {
            case 2:
                outputSquare(fout, x, y, HEATMAP_COLOR_MAP2[i], HEATMAP_LEGEND_BLOCK_DIM);
                break;
            default:
                outputSquare(fout, x, y, HEATMAP_COLOR_MAP1[i], HEATMAP_LEGEND_BLOCK_DIM);
                break;
        }
        // print out the description for each color
        outputSeqName(fout, x_text, y_text, HEATMAP_LEGEND_DESC_ANGLE, HEATMAP_COLOR_DESC[i], HEATMAP_LEGEND_FONT, HEATMAP_LEGEND_FONT_SIZE, HEATMAP_LEGEND_TEXT_ALIGN);
        y += HEATMAP_LEGEND_BLOCK_DIM;
        y_text += HEATMAP_LEGEND_BLOCK_DIM;
    }
}

// Output the triangular heatmap
void outputTriHeatmap(string prefixOut, vector<string>* seqNames, int* Cij, int* row_index, int seqNum, int seqLen, int color_scheme) {
    
    int block_size = HEATMAP_BLOCK_DIM;
    int font_size = HEATMAP_FONT_SIZE;
    if (seqNum <= HEATMAP_FEW_SEQ_NUM) {
        block_size = HEATMAP_BLOCK_DIM_FEW_SEQ;
        font_size = HEATMAP_FONT_SIZE_FEW_SEQ;
    }

    int longestFirstTextLen = longestLenFirst(seqNames, row_index, 3);
    int longestLastTextLen = longestLenLast(seqNames, row_index, 3);
    int width = HEATMAP_LEFT_BORDER + longestFirstTextLen*font_size/2.0 + HEATMAP_GAP_BW_TEXT + 
                block_size*(seqNum-1) + HEATMAP_GAP_BW_LEGEND + HEATMAP_LEGEND_BLOCK_DIM +
                HEATMAP_LEGEND_GAP_BW_DESC + HEATMAP_LEGEND_FONT_SIZE*2.5 + HEATMAP_RIGHT_BORDER;
    int height = HEATMAP_TOP_BORDER + block_size*seqNum + (longestLastTextLen*font_size/2.0 - block_size) +
                 HEATMAP_BOTTOM_BORDER;
    int height_legend = HEATMAP_TOP_BORDER + HEATMAP_LEGEND_BLOCK_DIM * HEATMAP_COLOR_NUM + HEATMAP_BOTTOM_BORDER;
    if (height < height_legend)
        height = height_legend;
    double y = HEATMAP_TOP_BORDER;
    double x;
    double y_coord, x_coord, cij;
    int row, col, i, j;
    
    string outFileName = fileName(prefixOut, "Triangular.svg");
    ofstream fout;
    fout.open(outFileName.c_str());

    // output the header
    outputSVGHeader(fout, width, height);
    
    // output the triangular matrix
    for (i=0; i<seqNum; i++) {
        
        row = row_index[i];
        x = HEATMAP_LEFT_BORDER + longestFirstTextLen*font_size/2.0 + HEATMAP_GAP_BW_TEXT + i*block_size;
        
        // print out the sequence name
        if (seqNames->at(row).length() > 0) {
            x_coord = x - HEATMAP_GAP_BW_TEXT;
            y_coord = y + block_size*0.8;
            outputSeqName(fout, x_coord, y_coord, HEATMAP_NAME_ANGLE, seqNames->at(row), HEATMAP_FONT, font_size, HEATMAP_TEXT_ALIGN);
        }
        
        for (j=i+1; j<seqNum; j++) {
            col = row_index[j];
            // print out the squares representing the Cij values
            cij = (double)Cij[row*seqNum+col]/seqLen;
            outputSquare(fout, x, y, getColor(cij, color_scheme), block_size);
            x += block_size;
        }
        
        y += block_size;
    }
    
    // output the legend
    y = HEATMAP_TOP_BORDER;
    x = HEATMAP_LEFT_BORDER + longestFirstTextLen*font_size/2.0 + HEATMAP_GAP_BW_TEXT + 
                block_size*(seqNum-1) + HEATMAP_GAP_BW_LEGEND;
    outputLegend(fout, x, y, color_scheme);
    
    // output the tail
    outputSVGTail(fout);

    fout.close();
}

// Output the full heatmap
void outputFullHeatmap(string prefixOut, vector<string>* seqNames, int* Cij, int* row_index, int seqNum, int seqLen, int color_scheme) {
    
    int block_size = HEATMAP_BLOCK_DIM;
    int font_size = HEATMAP_FONT_SIZE;
    if (seqNum <= HEATMAP_FEW_SEQ_NUM) {
        block_size = HEATMAP_BLOCK_DIM_FEW_SEQ;
        font_size = HEATMAP_FONT_SIZE_FEW_SEQ;
    }
    
    int longestTextLen = longestLen(seqNames);
    int width = HEATMAP_LEFT_BORDER + longestTextLen*font_size/2.0 + HEATMAP_GAP_BW_TEXT + 
                block_size*seqNum + HEATMAP_GAP_BW_LEGEND + HEATMAP_LEGEND_BLOCK_DIM +
                HEATMAP_LEGEND_GAP_BW_DESC + HEATMAP_LEGEND_FONT_SIZE*2.5 + HEATMAP_RIGHT_BORDER;
    int height = HEATMAP_TOP_BORDER + block_size*seqNum + HEATMAP_BOTTOM_BORDER;
    int height_legend = HEATMAP_TOP_BORDER + HEATMAP_LEGEND_BLOCK_DIM * HEATMAP_COLOR_NUM + HEATMAP_BOTTOM_BORDER;
    if (height < height_legend)
        height = height_legend;
    double y = HEATMAP_TOP_BORDER;
    double x;
    double y_coord, x_coord, cij;
    int row, col, i, j;
    
    string outFileName = fileName(prefixOut, "Full.svg");
    ofstream fout;
    fout.open(outFileName.c_str());

    // output the header
    outputSVGHeader(fout, width, height);
    
    // output the full matrix
    for (i=0; i<seqNum; i++) {
        
        row = row_index[i];
        x = HEATMAP_LEFT_BORDER + longestTextLen*font_size/2.0 + HEATMAP_GAP_BW_TEXT;
        
        // print out the sequence name
        if (seqNames->at(row).length() > 0) {
            x_coord = x - HEATMAP_GAP_BW_TEXT;
            y_coord = y + block_size*0.8;
            outputSeqName(fout, x_coord, y_coord, 0, seqNames->at(row), HEATMAP_FONT, font_size, HEATMAP_TEXT_ALIGN);
        }

        for (j=0; j<i; j++) {
            col = row_index[j];
            // print out the squares representing the Cij values
            cij = (double)Cij[row*seqNum+col]/seqLen;
            outputSquare(fout, x, y, getColor(cij, color_scheme), block_size);
            x += block_size;
        }

        // print out the gray square
        outputSquare(fout, x, y, HEATMAP_GRAY_COLOR, block_size);
        x += block_size;
        
        for (j=i+1; j<seqNum; j++) {
            col = row_index[j];
            // print out the squares representing the Cij values
            cij = (double)Cij[row*seqNum+col]/seqLen;
            outputSquare(fout, x, y, getColor(cij, color_scheme), block_size);
            x += block_size;
        }
        
        y += block_size;
    }
    
    // output the legend
    y = HEATMAP_TOP_BORDER;
    x = HEATMAP_LEFT_BORDER + longestTextLen*font_size/2.0 + HEATMAP_GAP_BW_TEXT + 
                block_size*seqNum + HEATMAP_GAP_BW_LEGEND;
    outputLegend(fout, x, y, color_scheme);
    
    // output the tail
    outputSVGTail(fout);

    fout.close();
}

// output the message about how to use the R scripts generated
void outputRMessage(UserOptions &user_options) {
    
    bool anyRScriptGenerated = false;
    int i;
    for (i=0; i<user_options.outTables.size(); i++) {
        switch(user_options.outTables[i]) {
            case 1:
            case 2:
            case 3:
                anyRScriptGenerated = true;
                break;
        }
    }
    if (anyRScriptGenerated) {
        cout << R_MESSAGE << endl;
    }
    
}
