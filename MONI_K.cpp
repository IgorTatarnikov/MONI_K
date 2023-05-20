#include <string>
#include <cstdlib>
#include <iostream>
#include <sdsl/wt_blcd.hpp>
#include <sdsl/construct.hpp>
#include <tuple>
#include <algorithm>
#include <vector>

#include "MONI_K.h"

bool printMEMs = false;
bool printTables = false;
bool useTestFile = false;
bool useTestPatterns = false;

int main(int argc, char** argv) {
    std::cout << "Print tables? (0/1)" << std::endl;
    std::cin >> printTables;

    std::cout << "Print MEMs? (0/1)" << std::endl;
    std::cin >> printMEMs;

    std::cout << "Use test file? (0/1)" << std::endl;
    std::cin >> useTestFile;

    std::string fileName;
    if (useTestFile) {
        fileName = "../test";
    } else {
        std::cout << "Please enter the file names" << std::endl;
        std::cin >> fileName;
    }

    std::cout << "Use test patterns? (0/1)" << std::endl;
    std::cin >> useTestPatterns;
    std::vector<std::string> testPatterns;

    if (useTestPatterns) {
        testPatterns = readTestStrings("test_strings.txt");
    } else {
        std::cout << "Please enter the pattern" << std::endl;
        std::string pattern;
        std::cin >> pattern;
        testPatterns = {pattern};
    }

    int n, r;
    bool kConstruction = KCONSTRUCTION;

    std::ifstream summaryFile(fileName + "_Summary_Bin");

    summaryFile >> n;
    summaryFile >> r;

    summaryFile.close();

    std::cout << "Reading the original text" << std::endl;
    std::ifstream textFile(fileName + ".txt");
    std::stringstream buffer;
    buffer << textFile.rdbuf();
    std::string text = buffer.str();

    textFile.close();

    std::string table2NameBin = fileName + "_MONITable_Bin";
    std::string table3NameBin = fileName + "_PhiTable_Bin";
    std::string table4NameBin = fileName + "_InversePhiTable_Bin";

    MONI_Table tableMONI(table2NameBin, r, kConstruction);
    Phi_Table tablePhi(table3NameBin, r);
    Inverse_Phi_Table tableInversePhi(table4NameBin, r);

    int k = 3;

    std::string BWTHead;

    for (int i = 0; i < r; i++) {
        BWTHead += tableMONI.BWT_head[i];
    }

    //Construct a balanced wavelet tree for BWT head for fast predecessor queries
    sdsl::wt_blcd<> wtBlcd;
    construct_im(wtBlcd, BWTHead, 1);

    std::string pattern;

    while (!testPatterns.empty()) {
        pattern = testPatterns.back();
        testPatterns.pop_back();
        unsigned int p_n = pattern.length();

        //Find a starting row such that the last character of the pattern is the same as BWT head
        unsigned int startj;
        for (startj = 0; startj < r; startj++) {
            if (tableMONI.BWT_head[startj] == pattern[p_n - 1]) {
                break;
            }
        }

        auto moniKTable = (unsigned int**) calloc(p_n + 1, sizeof(int*));
        onlineMONIK(moniKTable, tableMONI, tablePhi, tableInversePhi, startj, r, wtBlcd, pattern, text, k);

        if (printMEMs) {
            std::cout << "The k-MEMs of " << pattern << " online:" << std::endl;
            std::cout << pattern.substr(0, moniKTable[p_n][4]) << std::endl;

            //Start from the last row (index 0 of the pattern) and print out the pattern of correct length whenever
            //the min(li , Li) at i+1 is greater than at i. This is a k-MEM of length min(li, Li)
            for (unsigned int i = p_n - 1; i > 0; i--) {
                if (moniKTable[i][4] >= moniKTable[i + 1][4]) {
                    std::cout << pattern.substr(p_n - i, moniKTable[i][4]) << std::endl;
                }
            }
        }

        //If the MONI table was constructed with values for a precalculated k, print those values as well.
        if (kConstruction) {
            auto preCalcMONIKTable = (unsigned int**) calloc(p_n + 1, sizeof(int*));
            preCalcMONIK(preCalcMONIKTable, tableMONI, startj, r, wtBlcd, pattern, text, k, n);

            if (printMEMs) {
                std::cout << "The k-MEMs of " << pattern << " pre-computed:" << std::endl;
                std::cout << pattern.substr(0, preCalcMONIKTable[p_n][4]) << std::endl;

                for (unsigned int i = p_n - 1; i > 0; i--) {
                    if (preCalcMONIKTable[i][4] >= preCalcMONIKTable[i + 1][4]) {
                        std::cout << pattern.substr(p_n - i, preCalcMONIKTable[i][4]) << std::endl;
                    }
                }
            }

            for (int i = 1; i <= p_n; i++) {
                if (preCalcMONIKTable[i][4] != moniKTable[i][4]) {
                    std::cout << "ERROR: Precomputed k-MEMs do not match online k-MEMs for " << pattern << std::endl;
                    break;
                }
            }

            for (int i = 0; i <= p_n; i++) {
                free(preCalcMONIKTable[i]);
            }

            free(preCalcMONIKTable);
        }

        for (int i = 0; i <= p_n; i++) {
            free(moniKTable[i]);
        }

        free(moniKTable);

        if (!useTestPatterns) {
            std::cout << "Please enter another pattern or 0 to terminate" << std::endl;
            std::cin >> pattern;

            if (pattern == "0") {
                break;
            }

            testPatterns = {pattern};
        }

    }

    return 0;
}


void onlineMONIK(unsigned int** moniKTable, MONI_Table &tableMONI, Phi_Table &tablePhi, Inverse_Phi_Table &tableInversePhi, unsigned int startj, int r,
            sdsl::wt_blcd<> &wtBlcd, const std::string &pattern, const std::string &text, int k) {

    unsigned int p_n = pattern.length();
    int numLCP = 2 * k - 2;

    unsigned int curri = p_n;
    unsigned int currj = startj;
    unsigned int currl = 0;
    unsigned int currq = tableMONI.head[currj];
    unsigned int currSA = tableMONI.SA_head[currj];
    auto LCPArrays = (unsigned int**) calloc(p_n + 1, sizeof(int*));

    unsigned int tempSA;
    unsigned int inversePhiRow;
    unsigned int phiRow;
    int windowSize = k-1;

    for (int i = 0; i <= p_n; i++) {
        moniKTable[i] = (unsigned int*) calloc(5, sizeof(int));

        tempSA = currSA;
        //Find the correct starting row for inverse phi given the current SA value
        inversePhiRow = pred2D(tempSA, tableInversePhi.SA_tail, 0, r, r);

        //Use inverse phi to take k-1 steps back through the SA
        for (int j = 0; j < (k-1); j++) {
            tempSA = inversePhi(tableInversePhi, inversePhiRow, tempSA, r);
        }

        //Find the correct starting row for phi given the current SA value (k-1 back)
        phiRow = pred2D(tempSA, tablePhi.SA_head, 0, r, r);

        //Use phi to take 2 * k - 2 steps forward through the SA and calculate the LCP values
        LCPArrays[i] = (unsigned int*) calloc(numLCP, sizeof(unsigned int));
        for (int j = 0; j < numLCP; j++) {
            LCPArrays[i][numLCP - (j + 1)] = LCPStep(tablePhi, phiRow, tempSA);
            tempSA = phi(tablePhi, phiRow, tempSA, r);
        }

        //Find the maximum minimum LCP value given a k-1 window size
        unsigned int maxMin = 0;
        unsigned int currMin;
        for (int j = 0; j <= numLCP - windowSize; j++) {
            currMin = *std::min_element(LCPArrays[i] + j, LCPArrays[i] + j + windowSize);
            maxMin = currMin > maxMin ? currMin : maxMin;
        }

        //Store the maxMin (L_tail) value, and either the l given by MONI or the maxMin (L_tail) value, whichever is smaller
        moniKTable[i][3] = maxMin;
        moniKTable[i][4] = currl < maxMin ? currl : maxMin;

        //If the current character is not the first character of the pattern, check if the BWT head is the same as the previous character
        if (curri > 0) {
            if (tableMONI.BWT_head[currj] != pattern[curri - 1]) {
                //If not update the row as per MONI
                updateRow(text, pattern[curri - 1], tableMONI, currj, currl, currq, currSA, r);
            }
        }

        moniKTable[i][0] = currq;
        moniKTable[i][1] = currSA;
        moniKTable[i][2] = currl;

        //Code to print out the output table as per table 5 in the original paper
        if (printTables) {
            std::cout << currq << "\t";
            std::cout << currSA << "\t";
            std::cout << currl << "\t";
            std::cout << "[";

            for (int j = 0; j < numLCP; j++) {
                std::cout << LCPArrays[i][j] << ",";
            }

            std::cout << "]" << "\t" << maxMin << "\t" << moniKTable[i][4] << std::endl;
        }

        curri--;
        currl++;

        //Take a LF step
        LFStep(tableMONI, currj, currq, currSA, wtBlcd, r);
    }
}


//MONI-k given that the table has pre-calculated values for the given k
void preCalcMONIK(unsigned int** moniKTable, MONI_Table& tableMONI, unsigned int startj, int r,
                            sdsl::wt_blcd<>& wtBlcd, const std::string& pattern, const std::string& text, int k, int n) {
    unsigned int p_n = pattern.length();
    int moniTableNumColumns = 5;

    auto BWTArrays = (char**) calloc(p_n + 1, sizeof(char*));

    //Set the starting values for MONI
    unsigned int curri = p_n;
    unsigned int currj = startj;
    unsigned int currl = 0;
    unsigned int currq = tableMONI.head[currj];
    unsigned int currSA = tableMONI.SA_head[currj];
    unsigned int currL;
    unsigned int currs;
    unsigned int offset;

    //Calculate the L and offset values at the head of the start run
    currL = tableMONI.L_head[currj];
    offset = tableMONI.offset_head[currj];
    currs = currq - offset;

    moniKTable[0] = (unsigned int*) calloc(8, sizeof(int));
    moniKTable[0][0] = currq;
    moniKTable[0][1] = currl;
    moniKTable[0][2] = currL;
    moniKTable[0][3] = currs;

    if (printTables) {
        std::cout << currq << "\t";
        std::cout << currl << "\t";
        std::cout << " " << "\t";
        std::cout << " " << "\t";
        std::cout << " " << std::endl;
    }

    currl++;
    curri--;

    //Take the first LF step
    LFStep(tableMONI, currj, currq, currSA, wtBlcd, r);
    currs = currq - offset;

    //Loop through the pattern calculating the length of the k-MEM at each position
    for (int i = 1; i < p_n; i++) {
        BWTArrays[i] = (char*) calloc(k, sizeof(char));
        moniKTable[i] = (unsigned int*) calloc(moniTableNumColumns, sizeof(int));

        moniKTable[i][0] = currq;
        moniKTable[i][1] = currl;
        moniKTable[i][2] = currL;
        moniKTable[i][3] = currs;
        moniKTable[i][4] = currl < currL ? currl : currL;


        //Values used to extract the BWT[currs...currs+k-1]
        unsigned int LFj = currj;
        unsigned int tempS = currs;
        unsigned int tempj = currj;

        //Move the LFj row count back if tempS is smaller than the head of the run at LFj
        //This makes sure we start at the correct row in the MONI table to extract the BWT head values
        while (tableMONI.head[LFj] > tempS) {
            LFj--;
            tempj--;
        }

        //Extract each character of the BWT[currs...currs+k-1], increment LFj when tempS crosses a run boundary
        for (int j = 0; j <= n && j < k; j++) {
            if (LFj + 1 < r && tempS >= tableMONI.head[LFj + 1]) {
                LFj++;
            }
            BWTArrays[i][j] = tableMONI.BWT_head[LFj];
            tempS++;
        }

        auto prevLetter = pattern[curri - 1];
        //Case 1 - BWT[currs] == ... === BWT[currs+k-1] ==  prevLetter
        //Take a simple LF step back with both currq and currs and continue
        if (tableMONI.BWT_head[tempj] == prevLetter && (tempj == (r - 1) || currs + 2 < tableMONI.head[tempj + 1])) {
            unsigned int pi = wtBlcd.rank(tempj, tableMONI.BWT_head[tempj]) +
                              std::get<1>(wtBlcd.lex_smaller_count(r, tableMONI.BWT_head[tempj]));
            currs = tableMONI.mu[pi] + currs - tableMONI.head[tempj];

            currl++;
            currL++;
            curri--;

            LFStep(tableMONI, currj, currq, currSA, wtBlcd, r);
        } else {
            //Case 2 - BWT[currq] != prevLetter
            //Reset currq and currl as per MONI
            if (tableMONI.BWT_head[currj] != prevLetter) {
                updateRow(text, prevLetter, tableMONI, currj, currl, currq, currSA, r);
                moniKTable[i][1] = currl;
                moniKTable[i][0] = currq;
            } else {
                //Case 3 - BWT[currq - 1] is not the end of a run
                //Choose b in [currs, currs+k-1] such that BWT[b] == prevLetter and BWT[b+1] != BWT[b]
                unsigned int b = currs;

                while (BWTArrays[i][b-currs] != prevLetter || BWTArrays[i][b - currs] == BWTArrays[i][b - currs + 1]) {
                    b++;
                }

                int lce;

                //If b is smaller than the head of the current run, use the tail of the correct run for LCE
                if (b < tableMONI.head[currj]) {
                    lce = LCE(text, currSA, tableMONI.SA_tail[tempj]);
                } else {
                    lce = LCE(text, currSA, tableMONI.SA_tail[currj]);
                }

                currl = lce < currl ? lce : currl;
                currq = b;
            }

            currl++;
            curri--;

            //Grab the correct L and offset values for the next LF step
            //If currq is the head of the run, use the head values
            //If currq is smaller than the head of the run, use the tail values for the previous run
            //Otherwise use the tail values for the current run
            if (currq == tableMONI.head[currj]) {
                currL = tableMONI.L_head[currj];
                offset = tableMONI.offset_head[currj];
            } else if (currq < tableMONI.head[currj]) {
                currL = tableMONI.L_tail[tempj];
                offset = tableMONI.offset_tail[tempj];
            } else {
                currL = tableMONI.L_tail[currj];
                offset = tableMONI.offset_tail[currj];
            }

            LFStep(tableMONI, currj, currq, currSA, wtBlcd, r);

            currs = currq - offset;
        }

        if (printTables) {
            std::cout << moniKTable[i][0] << "\t";
            std::cout << currSA << "\t";
            std::cout << moniKTable[i][1] << "\t";
            std::cout << moniKTable[i][2] << "\t";
            std::cout << moniKTable[i][4] << "\t";
            std::cout << moniKTable[i][3] << "\t";

            std::cout << "[" << BWTArrays[i] << "]" << std::endl;
        }
    }

    //Record the final values into the table
    BWTArrays[p_n] = (char*) calloc(k, sizeof(char));
    moniKTable[p_n] = (unsigned int*) calloc(moniTableNumColumns, sizeof(int));

    moniKTable[p_n][0] = currq;
    moniKTable[p_n][1] = currl;
    moniKTable[p_n][2] = currL;
    moniKTable[p_n][3] = currs;
    moniKTable[p_n][4] = currl < currL ? currl : currL;

    if (printTables) {
        std::cout << moniKTable[p_n][0] << "\t";
        std::cout << currSA << "\t";
        std::cout << moniKTable[p_n][1] << "\t";
        std::cout << moniKTable[p_n][2] << "\t";
        std::cout << moniKTable[p_n][4] << "\t";
        std::cout << moniKTable[p_n][3] << std::endl;
    }
}

//Basic longest common extension (LCE) query (has to have the entire text loaded into memory)
int LCE(const std::string& text, unsigned int a, unsigned int b) {
    int lce = 0;

    while (text[a] == text[b] && text[a] != STRING_SEPARATOR) {
        lce++;
        a++;
        b++;
    }

    return lce;
}

//LF step as per the MONI
void LFStep(MONI_Table tableMONI, unsigned int& currj, unsigned int& currq, unsigned int& currSA, sdsl::wt_blcd<>& wtBlcd, int r) {
    unsigned int pi = wtBlcd.rank(currj, tableMONI.BWT_head[currj]) + std::get<1>(wtBlcd.lex_smaller_count(r, tableMONI.BWT_head[currj]));

    currq = tableMONI.mu[pi] + currq - tableMONI.head[currj];
    currSA = currSA - 1;

    currj = exponentialSearch(tableMONI.head, tableMONI.finger[pi], currq, r);
}

//Exponential search starting at a given location, double until the target is found or exceeded, then binary search
//the range for the predecessor
unsigned int exponentialSearch(unsigned int* table, unsigned int start, unsigned int target, int r) {
    unsigned int delta = 1;
    unsigned int max = start + delta;

    while (max < r && table[max] <= target) {
        delta *= 2;
        max = start + delta;
    }

    return pred2D(target, table, start + (delta / 2), max, r);
}

//Predecessor query for a given column in a 2D array
unsigned int pred2D(unsigned int target, const unsigned int* array, unsigned int low, unsigned int high, int r) {
    unsigned int mid = (low + high + 1) / 2;
    while (low != high && mid < r) {
        mid = (low + high + 1) / 2;

        if (array[mid] <= target) {
            low = mid;
        } else {
            high = mid - 1;
        }
    }

    return low;
}

//Update values when the previous letter in the pattern does not match BWT head, as per MONI
void updateRow(const std::string& text, char prevLetter, MONI_Table tableMONI, unsigned int& currj, unsigned int& currl, unsigned int& currq, unsigned int& currSA, int r) {
    unsigned int upperBoundary = currj-1;
    unsigned int lowerBoundary = currj+1;

    while (upperBoundary < lowerBoundary && tableMONI.BWT_head[upperBoundary] != prevLetter) {
        upperBoundary--;
    }

    while (lowerBoundary < r && tableMONI.BWT_head[lowerBoundary] != prevLetter) {
        lowerBoundary++;
    }

    int lUpper = (upperBoundary < lowerBoundary) ? LCE(text, currSA, tableMONI.SA_tail[upperBoundary]) : -1;
    int lLower = (lowerBoundary < r) ? LCE(text, currSA, tableMONI.SA_head[lowerBoundary]) : -1;

    if (lUpper > lLower) {
        currj = upperBoundary;
        currl = currl < lUpper ? currl : lUpper;
        currq = tableMONI.tail[currj];
        currSA = tableMONI.SA_tail[currj];
    } else {
        currj = lowerBoundary;
        currl = currl < lLower ? currl : lLower;
        currq = tableMONI.head[currj];
        currSA = tableMONI.SA_head[currj];
    }
}

//Use Phi to take a step forward in the suffix array
unsigned int phi(Phi_Table tablePhi, unsigned int& currRow, unsigned int currSA, int r) {
    unsigned int newSA =  tablePhi.SA_tail[currRow] + currSA - tablePhi.SA_head[currRow];
    currRow = exponentialSearch(tablePhi.SA_head, tablePhi.finger[currRow], newSA, r);

    return newSA;
}

//Use Inverse PHi to take a step back in the suffix array
unsigned int inversePhi(Inverse_Phi_Table tableInversePhi, unsigned int& currRow, unsigned int currSA, int r) {
    unsigned int newSA = tableInversePhi.SA_head[currRow] + currSA - tableInversePhi.SA_tail[currRow];
    currRow = exponentialSearch(tableInversePhi.SA_tail, tableInversePhi.finger[currRow], newSA, r);

    return newSA;
}

//Use the LCP formula described in the paper to step forward in the LCP array
unsigned int LCPStep(Phi_Table tablePhi, unsigned int currRow, unsigned int currSA) {
    unsigned int lcp = (tablePhi.LCP_head[currRow] + tablePhi.SA_head[currRow]) < currSA ? 0 : tablePhi.LCP_head[currRow] + tablePhi.SA_head[currRow] - currSA;
    return lcp;
}

//Reads strings from a file and returns a vector of strings
std::vector<std::string> readTestStrings(const std::string& fileName) {
    std::ifstream file(fileName);
    std::vector<std::string> strings;
    std::string line;

    while (std::getline(file, line)) {
        strings.push_back(line);
    }

    file.close();

    return strings;
}
