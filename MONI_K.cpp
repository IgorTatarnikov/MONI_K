#include <string>
#include <cstdlib>
#include <iostream>
#include <sdsl/wt_blcd.hpp>
#include <sdsl/construct.hpp>
#include <tuple>
#include <algorithm>

#include "MONI_K.h"

int main(int argc, char** argv) {
    int table2NumColumns = TABLE2NUMCOLUMNS;
    int table3NumColumns = TABLE3NUMCOLUMNS;
    int table4NumColumns = TABLE4NUMCOLUMNS;

    std::cout << "Please enter the file name" << std::endl;

    std::string fileName;
    std::getline(std::cin, fileName);

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

    std::cout << "Please enter a pattern:" << std::endl;
    std::string pattern;
    std::getline(std::cin, pattern);

    int k = 3;
    int numLCP = 2 * k - 2;

    int p_n = pattern.length();

    //Find a starting row such that the last character of the pattern is the same as BWT head
    unsigned int startj = 5;
//    for (startj = 0; startj < r; startj++) {
//        if (tableMONI.BWT_head[startj] == pattern[p_n - 1]) {
//            break;
//        }
//    }

    //Initialize values for MONI
    unsigned int curri = p_n;
    unsigned int currj = startj;
    unsigned int currl = 0;
    unsigned int currq = tableMONI.head[currj];
    unsigned int currSA = tableMONI.SA_head[currj];

    auto moniKTable = (unsigned int**) calloc(p_n + 1, sizeof(int*));
    auto LCPArrays = (unsigned int**) calloc(p_n + 1, sizeof(int*));

    std::string BWTHead;

    for (int i = 0; i < r; i++) {
        BWTHead += tableMONI.BWT_head[i];
    }

    //Construct a balanced wavelet tree for BWT head for fast predecessor queries
    sdsl::wt_blcd<> wtBlcd;
    construct_im(wtBlcd, BWTHead, 1);

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
        std::cout << currq <<"\t";
        std::cout << currSA <<"\t";
        std::cout << currl << "\t";
        std::cout << "[";

        for (int j = 0; j < numLCP; j++) {
            std::cout << LCPArrays[i][j] << ",";
        }

        std::cout << "]" << "\t" << maxMin << "\t" << moniKTable[i][4] << std::endl;
//        std::cout << moniKTable[i][3] << "\t";
//        std::cout << moniKTable[i][4] << std::endl;

        curri--;
        currl++;

        //Take a LF step
        LFStep(tableMONI, currj, currq, currSA, wtBlcd, r);
    }

    std::cout << "The k-MEMs of " << pattern << " online:" << std::endl;
    std::cout << pattern.substr(0, moniKTable[p_n][4]) << std::endl;

    //Start from the last row (index 0 of the pattern) and print out the pattern of correct length whenever
    //the min(li , Li) at i+1 is greater than at i. This is a k-MEM of length min(li, Li)
    for (int i = p_n - 1; i > 0; i--) {
        if (moniKTable[i][4] >= moniKTable[i+1][4]) {
            std::cout << pattern.substr(p_n - i, moniKTable[i][4]) << std::endl;
        }
    }

    //If the MONI table was constructed with values for a precalculated k, pring those values as well.
    if (kConstruction) {
        unsigned int **preCalcMONIk = preCalcMONI(tableMONI, startj, r, wtBlcd, pattern, text, k, n);

        std::cout << "The k-MEMs of " << pattern << " pre-computed:" << std::endl;
        std::cout << pattern.substr(0, preCalcMONIk[p_n][4]) << std::endl;

        for (int i = p_n - 2; i > 0; i--) {
            if (preCalcMONIk[i][4] >= preCalcMONIk[i + 1][4]) {
                std::cout << pattern.substr(p_n - i + 1, preCalcMONIk[i][4]) << std::endl;
            }
        }
    }

    return 0;
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
    int pi = wtBlcd.rank(currj, tableMONI.BWT_head[currj]) + std::get<1>(wtBlcd.lex_smaller_count(r, tableMONI.BWT_head[currj]));

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
unsigned int pred2D(unsigned int target, unsigned int* array, unsigned int low, unsigned int high, int r) {
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

unsigned int phi(Phi_Table tablePhi, unsigned int& currRow, unsigned int currSA, int r) {
    unsigned int newSA =  tablePhi.SA_tail[currRow] + currSA - tablePhi.SA_head[currRow];
    currRow = exponentialSearch(tablePhi.SA_head, tablePhi.finger[currRow], newSA, r);

    return newSA;
}

unsigned int inversePhi(Inverse_Phi_Table tableInversePhi, unsigned int& currRow, unsigned int currSA, int r) {
    unsigned int newSA = tableInversePhi.SA_head[currRow] + currSA - tableInversePhi.SA_tail[currRow];
    currRow = exponentialSearch(tableInversePhi.SA_tail, tableInversePhi.finger[currRow], newSA, r);

    return newSA;
}
unsigned int LCPStep(Phi_Table tablePhi, unsigned int currRow, unsigned int currSA) {
    unsigned int lcp = (tablePhi.LCP_head[currRow] + tablePhi.SA_head[currRow]) < currSA ? 0 : tablePhi.LCP_head[currRow] + tablePhi.SA_head[currRow] - currSA;
    return lcp;
}

unsigned int** preCalcMONI(MONI_Table tableMONI, unsigned int startj, int r, sdsl::wt_blcd<>& wtBlcd, const std::string& pattern, const std::string& text, int k, int n) {
    unsigned int p_n = pattern.length();

    auto moniKTable = (unsigned int**) calloc(p_n + 1, sizeof(int*));
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
    unsigned int oldl;
    bool MONIReset = false;

    currL = tableMONI.L_tail[currj];
    offset = tableMONI.offset_tail[currj];
    currs = currq - offset;

    moniKTable[0] = (unsigned int*) calloc(8, sizeof(int));
    moniKTable[0][0] = currq;
    moniKTable[0][1] = currl;
    moniKTable[0][2] = currL;
    moniKTable[0][3] = currs;

    std::cout << currq <<"\t";
    std::cout << currl <<"\t";
    std::cout << currL << "\t";
    std::cout << moniKTable[0][4] << "\t";
    std::cout << currs << std::endl;

    currl++;
    curri--;
    LFStep(tableMONI, currj, currq, currSA, wtBlcd, r);

    for (int i = 1; i < p_n; i++) {
        BWTArrays[i] = (char*) calloc(k, sizeof(char));
        moniKTable[i] = (unsigned int*) calloc(8, sizeof(int));
        unsigned int LFj = currj;
        unsigned int LFs = currs;

        for (int j = 0; j <= n && j < k; j++) {
            if (LFj + 1 < r && LFs >= tableMONI.head[LFj+1]) {
                LFj++;
            }
            BWTArrays[i][j] = tableMONI.BWT_head[LFj];
            LFs++;
        }

        auto prevLetter = pattern[curri-1];
        if (tableMONI.BWT_head[currj] == prevLetter && (currj == (r-1) || currs + 2 < tableMONI.head[currj+1])) {
            moniKTable[i][0] = currq;
            moniKTable[i][1] = currl;
            moniKTable[i][2] = currL;
            moniKTable[i][3] = currs;
            moniKTable[i][4] = currl < currL ? currl : currL;

            std::cout << currq <<"\t";
            std::cout << currl <<"\t";
            std::cout << currL << "\t";
            std::cout << moniKTable[i][4] << "\t";
            std::cout << currs << std::endl;

            int pi = wtBlcd.rank(currj, tableMONI.BWT_head[currj]) + std::get<1>(wtBlcd.lex_smaller_count(r, tableMONI.BWT_head[currj]));
            currs = tableMONI.mu[pi] + currs - tableMONI.head[currj];

            currl++;
            currL++;
            curri--;
            LFStep(tableMONI, currj, currq, currSA, wtBlcd, r);
        } else {
            oldl = currl;
            if (tableMONI.BWT_head[currj] != prevLetter) {
                moniKTable[i][4] = currl < currL ? currl : currL;
                updateRow(text, prevLetter, tableMONI, currj, currl, currq, currSA, r);
                MONIReset = true;
            } else {
                unsigned int b = currs;

                while (BWTArrays[i][b - currs] == BWTArrays[i][b - currs + 1]) {
                    b++;
                }

                int lce = LCE(text, currSA, tableMONI.SA_tail[currj-1]);
                currl = lce < currl ? lce : currl;
                moniKTable[i][4] = currl < currL ? currl : currL;
                currq = b;
            }

            currL = tableMONI.L_tail[currj];

            if (MONIReset) {
                moniKTable[i][4] = oldl < currL ? oldl : currL;
                MONIReset = false;
            }

            moniKTable[i][0] = currq;
            moniKTable[i][1] = currl;
            moniKTable[i][2] = currL;
            moniKTable[i][3] = currs;

            std::cout << currq <<"\t";
            std::cout << currl <<"\t";
            std::cout << currL << "\t";
            std::cout << moniKTable[i][4] << "\t";
            std::cout << currs << std::endl;

            currl++;
            curri--;

            offset = tableMONI.offset_tail[currj];

            LFStep(tableMONI, currj, currq, currSA, wtBlcd, r);

            currs = currq - offset;
        }
    }

    moniKTable[p_n] = (unsigned int*) calloc(8, sizeof(int));
    moniKTable[p_n][0] = currq;
    moniKTable[p_n][1] = currl;
    moniKTable[p_n][2] = currL;
    moniKTable[p_n][3] = currs;
    moniKTable[p_n][4] = currl < currL ? currl : currL;

//    std::cout << currq <<"\t";
//    std::cout << currl <<"\t";
//    std::cout << currL << "\t";
//    std::cout << moniKTable[p_n][4] << "\t";
//    std::cout << currs << std::endl;

    return moniKTable;
}
