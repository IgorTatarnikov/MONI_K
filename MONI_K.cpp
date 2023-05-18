#include <string>
#include <fstream>
#include <cstdlib>
#include <iostream>
#include <sdsl/wt_huff.hpp>
#include <sdsl/wt_blcd.hpp>
#include <sdsl/wt_rlmn.hpp>
#include <sdsl/construct.hpp>
#include <tuple>
#include <algorithm>
#include <vector>
#include "MONI_K.h"

int main(int argc, char** argv) {
    int table2NumColumns = TABLE2NUMCOLUMNS;
    int table3NumColumns = TABLE3NUMCOLUMNS;
    int table4NumColumns = TABLE4NUMCOLUMNS;

    std::cout << "Please enter the file name" << std::endl;

    std::string fileName;
    std::getline(std::cin, fileName);

    int n, r;

    std::ifstream summaryFile(fileName + "_Summary");

    summaryFile >> n;
    summaryFile >> r;

    summaryFile.close();

    std::cout << "Reading the original text" << std::endl;
    std::ifstream textFile(fileName + ".txt");
    std::stringstream buffer;
    buffer << textFile.rdbuf();
    std::string text = buffer.str();

    textFile.close();

    std::string table2Name = fileName + "_Table2MONI";
    std::string table3Name = fileName + "_Table3MONI";
    std::string table4Name = fileName + "_Table4MONI";

    unsigned int** table2 = readTable(table2Name, r, table2NumColumns);
    unsigned int** table3 = readTable(table3Name, r, table3NumColumns);
    unsigned int** table4 = readTable(table4Name, r, table4NumColumns);

    std::cout << "Please enter a pattern:" << std::endl;
    std::string pattern;
    std::getline(std::cin, pattern);

    int k = 3;
    int numLCP = 2 * k - 2;

    int p_n = pattern.length();

    unsigned int curri = p_n;
    unsigned int currj = 5;
    unsigned int currl = 0;
    unsigned int currq = table2[currj][0];
    unsigned int currSA = table2[currj][1];


    auto moniKTable = (unsigned int**) calloc(p_n + 1, sizeof(int*));
    auto LCPArrays = (unsigned int**) calloc(p_n + 1, sizeof(int*));

    std::string BWTHead;

    for (int i = 0; i < r; i++) {
        BWTHead += (char) table2[i][4];
    }

    sdsl::wt_blcd<> wtBlcd;

    construct_im(wtBlcd, BWTHead, 1);

    unsigned int tempSA;
    unsigned int inversePhiRow;
    unsigned int phiRow;
    int windowSize = k-1;

    for (int i = 0; i <= p_n; i++) {
        moniKTable[i] = (unsigned int*) calloc(5, sizeof(int));

        tempSA = currSA;
        inversePhiRow = pred2D(tempSA, table4, 1, 0, r, r);

        for (int j = 0; j < (k-1); j++) {
            tempSA = inversePhi(table4, inversePhiRow, tempSA, r);
        }

        phiRow = pred2D(tempSA, table3, 0, 0, r, r);

        LCPArrays[i] = (unsigned int*) calloc(numLCP, sizeof(unsigned int));

        for (int j = 0; j < numLCP; j++) {
            LCPArrays[i][numLCP - (j + 1)] = LCPStep(table3, phiRow, tempSA);
            tempSA = phi(table3, phiRow, tempSA, r);
        }

        unsigned int maxMin = 0;
        unsigned int currMin;

        for (int j = 0; j <= numLCP - windowSize; j++) {
            currMin = *std::min_element(LCPArrays[i] + j, LCPArrays[i] + j + windowSize);
            maxMin = currMin > maxMin ? currMin : maxMin;
        }

        moniKTable[i][3] = maxMin;
        moniKTable[i][4] = currl < maxMin ? currl : maxMin;

        if (curri > 0) {
            auto prevLetter = (unsigned int) pattern[curri - 1];
            if (table2[currj][4] != prevLetter) {
                updateRow(text, prevLetter, table2, currj, currl, currq, currSA, r);
            }
        }

        moniKTable[i][0] = currq;
        moniKTable[i][1] = currSA;
        moniKTable[i][2] = currl;

        std::cout << currq <<"\t";
//        std::cout << currSA <<"\t";
        std::cout << currl << "\t";
//        std::cout << "[";
//
//        for (int j = 0; j < numLCP; j++) {
//            std::cout << LCPArrays[i][j] << ",";
//        }

//        std::cout << "]" << "\t" << maxMin << "\t" << moniKTable[i][4] << std::endl;
        std::cout << moniKTable[i][3] << "\t";
        std::cout << moniKTable[i][4] << std::endl;

        curri--;
        currl++;
        LFStep(table2, currj, currq, currSA, wtBlcd, r);
    }

    unsigned int** preCalcMONIk = preCalcMONI(table2, table3, table4, r, wtBlcd, pattern, text, k, n);

    std::cout << "The k-MEMs of " << pattern << " online:" << std::endl;
    std::cout << pattern.substr(0, moniKTable[p_n][4]) << std::endl;

    for (int i = p_n - 1; i > 0; i--) {
        if (moniKTable[i][4] >= moniKTable[i+1][4]) {
            std::cout << pattern.substr(p_n - i, moniKTable[i][4]) << std::endl;
        }
    }

    std::cout << "The k-MEMs of " << pattern << " pre-computed:" << std::endl;
    std::cout << pattern.substr(0, preCalcMONIk[p_n][4]) << std::endl;

    for (int i = p_n - 1; i > 0; i--) {
        if (preCalcMONIk[i][4] >= preCalcMONIk[i+1][4]) {
            std::cout << pattern.substr(p_n - i, preCalcMONIk[i][4]) << std::endl;
        }
    }


    return 0;
}


unsigned int** readTable(const std::string& fileName, int r, int numColumns) {
    std::ifstream tableFile(fileName);
    auto** table = (unsigned int**) calloc(r,sizeof(int*));
    std::cout << "Reading " << fileName << std::endl;

    for (int i = 0; i < r; i++) {
        table[i] = (unsigned int*) calloc(numColumns, sizeof(unsigned int));

        for (int j = 0; j < numColumns; j++) {
            tableFile >> table[i][j];
        }
    }

    tableFile.close();

    return table;
}

int LCE(const std::string& text, unsigned int a, unsigned int b) {
    int lce = 0;

    while (text[a] == text[b] && text[a] != STRING_SEPARATOR) {
        lce++;
        a++;
        b++;
    }

    return lce;
}

void LFStep(unsigned int** table2, unsigned int& currj, unsigned int& currq, unsigned int& currSA, sdsl::wt_blcd<>& wtBlcd, int r) {
    int pi = wtBlcd.rank(currj, (char) table2[currj][4]) + std::get<1>(wtBlcd.lex_smaller_count(r, (char) table2[currj][4]));

    currq = table2[pi][5] + currq - table2[currj][0];
    currSA = currSA - 1;

    currj = exponentialSearch(table2, 0, table2[pi][6], currq, r);
}

unsigned int exponentialSearch(unsigned int** table, unsigned int column, unsigned int start, unsigned int target, int r) {
    unsigned int delta = 1;
    unsigned int max = start + delta;

    while (max < r && table[max][column] <= target) {
        delta *= 2;
        max = start + delta;
    }

    return pred2D(target, table, column, start + (delta / 2), max, r);
}

unsigned int pred2D(unsigned int target, unsigned int **array, unsigned int column, unsigned int low, unsigned int high, int r) {
    unsigned int mid = (low + high + 1) / 2;
    while (low != high && mid < r) {
        mid = (low + high + 1) / 2;

        if (array[mid][column] <= target) {
            low = mid;
        } else {
            high = mid - 1;
        }
    }

    return low;
}

void updateRow(const std::string& text, unsigned int prevLetter, unsigned int** table2, unsigned int& currj, unsigned int& currl, unsigned int& currq, unsigned int& currSA, int r) {
    unsigned int upperBoundary = currj-1;
    unsigned int lowerBoundary = currj+1;

    while (upperBoundary < lowerBoundary && table2[upperBoundary][4] != prevLetter) {
        upperBoundary--;
    }

    while (lowerBoundary < r && table2[lowerBoundary][4] != prevLetter) {
        lowerBoundary++;
    }

    int lUpper = (upperBoundary < lowerBoundary) ? LCE(text, currSA, table2[upperBoundary][3]) : -1;
    int lLower = (lowerBoundary < r) ? LCE(text, currSA, table2[lowerBoundary][1]) : -1;

    if (lUpper > lLower) {
        currj = upperBoundary;
        currl = currl < lUpper ? currl : lUpper;
        currq = table2[currj][2];
        currSA = table2[currj][3];
    } else {
        currj = lowerBoundary;
        currl = currl < lLower ? currl : lLower;
        currq = table2[currj][0];
        currSA = table2[currj][1];
    }
}

unsigned int phi(unsigned int** table3, unsigned int& currRow, unsigned int currSA, int r) {
    unsigned int newSA =  table3[currRow][1] + currSA - table3[currRow][0];
    currRow = exponentialSearch(table3, 0, table3[currRow][3], newSA, r);

    return newSA;
}

unsigned int inversePhi(unsigned int** table4, unsigned int& currRow, unsigned int currSA, int r) {
    unsigned int newSA = table4[currRow][0] + currSA - table4[currRow][1];
    currRow = exponentialSearch(table4, 1, table4[currRow][2], newSA, r);

    return newSA;
}
unsigned int LCPStep(unsigned int** table3, unsigned int currRow, unsigned int currSA) {
    unsigned int lcp = (table3[currRow][2] + table3[currRow][0]) < currSA ? 0 : table3[currRow][2] + table3[currRow][0] - currSA;
    return lcp;
}

unsigned int** preCalcMONI(unsigned int** table2, unsigned int** table3, unsigned int** table4, int r, sdsl::wt_blcd<>& wtBlcd, const std::string& pattern, const std::string& text, int k, int n) {
    unsigned int p_n = pattern.length();

    unsigned int startj;
    auto lastLetter = (unsigned int) pattern[p_n - 1];

    for (startj = 0; startj < r; startj++) {
        if (table2[startj][4] == lastLetter) {
            break;
        }
    }

    auto moniKTable = (unsigned int**) calloc(p_n + 1, sizeof(int*));
    auto BWTArrays = (char**) calloc(p_n + 1, sizeof(char*));

    unsigned int curri = p_n;
    unsigned int currj = startj;
    unsigned int currl = 0;
    unsigned int currq = table2[currj][0];
    unsigned int currSA = table2[currj][1];
    unsigned int currL;
    unsigned int currs;
    unsigned int offset;
    unsigned int oldl;
    bool MONIReset = false;

    currl++;
    curri--;
    LFStep(table2, currj, currq, currSA, wtBlcd, r);

    currL = table2[currj][10];
    offset = table2[currj][9];
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

    for (int i = 1; i < p_n; i++) {
        BWTArrays[i] = (char*) calloc(k, sizeof(char));
        moniKTable[i] = (unsigned int*) calloc(8, sizeof(int));
        unsigned int LFj = currj;
        unsigned int LFs = currs;

        for (int j = 0; j <= n && j < k; j++) {
            if (LFj + 1 < r && LFs >= table2[LFj+1][0]) {
                LFj++;
            }
            BWTArrays[i][j] = (char) table2[LFj][4];
            LFs++;
        }

        auto prevLetter = (unsigned int) pattern[curri-1];
        if (table2[currj][4] == prevLetter && (currj == (r-1) || currs + 2 < table2[currj+1][0])) {
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

            int pi = wtBlcd.rank(currj, (char) table2[currj][4]) + std::get<1>(wtBlcd.lex_smaller_count(r, (char) table2[currj][4]));
            currs = table2[pi][5] + currs - table2[currj][0];

            currl++;
            currL++;
            curri--;
            LFStep(table2, currj, currq, currSA, wtBlcd, r);
        } else {
            oldl = currl;
            if (table2[currj][4] != prevLetter) {
                moniKTable[i][4] = currl < currL ? currl : currL;
                updateRow(text, prevLetter, table2, currj, currl, currq, currSA, r);
//                MONIReset = true;
            } else {
                unsigned int b = currs;

                while (BWTArrays[i][b - currs] == BWTArrays[i][b - currs + 1]) {
                    b++;
                }

                int lce = LCE(text, currSA, table2[currj-1][3]);
                currl = lce < currl ? lce : currl;
                moniKTable[i][4] = currl < currL ? currl : currL;
                currq = b;
            }

            currL = table2[currj][10];

//            if (MONIReset) {
//                moniKTable[i][4] = oldl < currL ? oldl : currL;
//                MONIReset = false;
//            }

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

            offset = table2[currj][9];

            LFStep(table2, currj, currq, currSA, wtBlcd, r);

            currs = currq - offset;
        }
    }

    moniKTable[p_n] = (unsigned int*) calloc(8, sizeof(int));
    moniKTable[p_n][0] = currq;
    moniKTable[p_n][1] = currl;
    moniKTable[p_n][2] = currL;
    moniKTable[p_n][3] = currs;
    moniKTable[p_n][4] = currl < currL ? currl : currL;

    std::cout << currq <<"\t";
    std::cout << currl <<"\t";
    std::cout << currL << "\t";
    std::cout << moniKTable[p_n][4] << "\t";
    std::cout << currs << std::endl;

    return moniKTable;
}
