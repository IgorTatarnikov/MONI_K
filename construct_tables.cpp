#include <string>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <algorithm>

#include "MONI_K.h"
#include "construct_tables.h"

int main(int argc, char** argv) {
    std::cout << "Please enter the file name" << std::endl;

    std::string file_name;
    std::getline(std::cin, file_name);

    bool kConstruction = KCONSTRUCTION;
    int k = 3;

    int table2NumColumns = TABLE2NUMCOLUMNS - 3;
    int table3NumColumns = TABLE3NUMCOLUMNS;
    int table4NumColumns = TABLE4NUMCOLUMNS;

    sdsl::csa_wt<> csa;
    sdsl::lcp_bitcompressed<> lcp;

    //Construct the data structures given the file name
    unsigned int n = constructDataStructures(&csa, &lcp, file_name);

    std::cout << "Finding the run heads" << std::endl;
    std::vector<unsigned int> run_heads;

    int r = 0;
    for (int i = 0; i < n; i++) {
        if (i == 0 || csa.bwt[i] != csa.bwt[(i - 1 + n) % n]) {
            run_heads.push_back(i);
            r++;
        }
    }

    //Filling the MONI table based on the BWT and LCP data structures
    auto** tableMONI = (unsigned int**) calloc(table2NumColumns, sizeof(unsigned int*));
    auto* BWTHeads = (char*) calloc(r, sizeof(char));
    auto** preCalcK = (unsigned char**) calloc(4, sizeof(unsigned char*));

    for (int i = 0; i < table2NumColumns; i++) {
        tableMONI[i] = (unsigned int*) calloc(r, sizeof(*tableMONI[i]));
    }

    for (int i = 0; i < 4; i++) {
        preCalcK[i] = (unsigned char*) calloc(r, sizeof(*preCalcK[i]));
    }

    //Breakpoint used to keep track of the progress of the table construction
    unsigned int breakpoint = r / 20;

    std::cout << "Filling MONI table with the head and SA head values" << std::endl;
    for (int i = 0; i < r; i++) {
        tableMONI[0][i] = run_heads[i];
        tableMONI[1][i] = csa[run_heads[i]];
        BWTHeads[i] = csa.bwt[run_heads[i]];

        if (breakpoint && (i % breakpoint == 0)) {
            std::cout << "Done with " << ((float) i/r) * 100 << "%" << std::endl;
        }
    }

    std::cout << "Filling MONI table with tail, SA tail, and mu values" << std::endl;
    for (int i = 0; i < r-1; i++) {
        tableMONI[2][i] = (tableMONI[0][(i + r + 1) % r] + n - 1) % n;
        tableMONI[3][i] = csa[tableMONI[2][i]];
        tableMONI[4][i] = csa.lf[tableMONI[0][i]];

        if (breakpoint && (i % breakpoint == 0)) {
            std::cout << "Done with " << ((float) i/r) * 100 << "%" << std::endl;
        }
    }

    tableMONI[2][r - 1] = n;
    tableMONI[3][r - 1] = csa[tableMONI[2][r - 1]];
    tableMONI[4][r - 1] = csa.lf[tableMONI[0][r - 1]];

    //Calculating the offset_tail and L_tail values for the tail of each run only if kConstruction is true
    if (kConstruction) {
        unsigned int LCPMin;
        unsigned char L_head;
        unsigned char offset_head;
        unsigned int start_head;
        unsigned int max_head;

        unsigned char L_tail;
        unsigned char offset_tail;
        unsigned int start_tail;
        unsigned int max_tail;
        int totalSpan = 2 * k - 2;
        int windowSize = k - 1;
        auto* tempLCPStoreTail = (unsigned int*) calloc(totalSpan, sizeof(unsigned int));
        auto* tempLCPStoreHead = (unsigned int*) calloc(totalSpan, sizeof(unsigned int));

        for (int i = 0; i < r; i++) {
            // n+1 accounts for the $ sign added by sdsl
            start_head = (csa.lf[tableMONI[0][i]] >= (k - 2)) ? csa.lf[tableMONI[0][i]] - k + 2 : 0;
            max_head = (csa.lf[tableMONI[0][i]] + k) < n + 1 ? csa.lf[tableMONI[0][i]] + k : n + 1;
            start_tail = (csa.lf[tableMONI[2][i]] >= (k - 2)) ? csa.lf[tableMONI[2][i]] - k + 2 : 0;
            max_tail = (csa.lf[tableMONI[2][i]] + k) < n + 1 ? csa.lf[tableMONI[2][i]] + k : n + 1;

            unsigned int total_head = max_head - start_head;
            unsigned int total_tail = max_tail - start_tail;

            L_head = 0;
            offset_head = 0;
            L_tail = 0;
            offset_tail = 0;

            for (unsigned int j = start_head; j < max_head; j++) {
                tempLCPStoreHead[j - start_head] = lcp[j];
            }

            for (unsigned int j = start_tail; j < max_tail; j++) {
                tempLCPStoreTail[j - start_tail] = lcp[j];
            }

            //Finding the maximum minimum LCP value in a k-1 window across the range of (q - k + 2) to (q + k - 1)
            for (int j = 0; j < total_head - 1; j++) {
                LCPMin = *std::min_element(tempLCPStoreHead + j, tempLCPStoreHead + windowSize + j);
                if (LCPMin > L_head) {
                    L_head = LCPMin;
                    offset_head = windowSize - j;
                }
            }

            for (int j = 0; j < total_tail - 1; j++) {
                LCPMin = *std::min_element(tempLCPStoreTail + j, tempLCPStoreTail + windowSize + j);
                if (LCPMin > L_tail) {
                    L_tail = LCPMin;
                    offset_tail = windowSize - j;
                }
            }

            //Record the offset_tail and L_tail
            preCalcK[0][i] = offset_tail;
            preCalcK[1][i] = L_tail;
            preCalcK[2][i] = offset_head;
            preCalcK[3][i] = L_head;

            if (breakpoint && (i % breakpoint == 0)) {
                std::cout << "Done with " << ((float) i/r) * 100 << "%" << std::endl;
            }

//            free(tempLCPStoreHead);
//            free(tempLCPStoreTail);
        }
    }

    std::cout << "Sorting the mu values" << std::endl;
    qsort(tableMONI[4], r, sizeof(int), rowComp);

    std::cout << "Finding the finger values for the MONI table by running a predecessor query" << std::endl;
    for (int i = 0; i < r; i++) {
        tableMONI[5][i] = pred(tableMONI[4][i], tableMONI[0], r, 1);

        if (breakpoint && (i % breakpoint == 0)) {
            std::cout << "Done with " << ((float) i/r) * 100 << "%" << std::endl;
        }
    }

    std::cout << "Filling phi and inverse phi tables" << std::endl;
    auto** tablePhi = (unsigned int**) calloc(r, sizeof(unsigned int*));
    auto** tableInversePhi = (unsigned int**) calloc(r, sizeof(unsigned int*));

    for (int i = 0; i < r; i++) {
        tablePhi[i] = (unsigned int*) calloc(table3NumColumns, sizeof(*tablePhi[i]));
        tableInversePhi[i] = (unsigned int*) calloc(table4NumColumns, sizeof(*tableInversePhi[i]));

        //Grabbing the SA head value from the MONI table
        tablePhi[i][0] = tableMONI[1][i];
        tableInversePhi[i][0] = tableMONI[1][i];

        //Grabbing the SA tail value from the MONI table (the SA value immediately preceding the head in the MONI table)
        tablePhi[i][1] = tableMONI[3][(i + r - 1) % r];
        tableInversePhi[i][1] = tablePhi[i][1];

        tablePhi[i][2] = lcp[tableMONI[0][i]];

        if (breakpoint && (i % breakpoint == 0)) {
            std::cout << "Done with " << ((float) i/r) * 100 << "%" << std::endl;
        }
    }

    //Sorting the phi table based on the SA head values
    std::cout << "Sorting phi table" << std::endl;
    qsort(tablePhi, r, sizeof(*tablePhi), rowCompArrayFirstElem);

    //Sorting the inverse phi table based on the SA tail values
    std::cout << "Sorting inverse phi table" << std::endl;
    qsort(tableInversePhi, r, sizeof(*tableInversePhi), rowCompArraySecondElem);

    std::cout << "Predecessor to fill finger column in phi and inverse phi tables" << std::endl;
    for (int i = 0; i < r; i++) {
        tablePhi[i][3] = pred2D(tablePhi[i][1], tablePhi, 0, r);
        tableInversePhi[i][2] = pred2D(tableInversePhi[i][0], tableInversePhi, 1, r);

        if (breakpoint && (i % breakpoint == 0)) {
            std::cout << "Done with " << ((float) i/r) * 100 << "%" << std::endl;
        }
    }

    //Summary file just contains the n and r values as ASCII characters
    std::cout << "Writing to file" << std::endl;
    std::ofstream summaryFile(file_name.substr(0,file_name.length()-4) + "_Summary_Bin");
    auto tableMONIFile = std::fstream(file_name.substr(0, file_name.length() - 4) + "_MONITable_Bin", std::ios::out | std::ios::binary);
    auto tablePhiFile = std::fstream(file_name.substr(0, file_name.length() - 4) + "_PhiTable_Bin", std::ios::out | std::ios::binary);
    auto tableInversePhiFile = std::fstream(file_name.substr(0, file_name.length() - 4) + "_InversePhiTable_Bin", std::ios::out | std::ios::binary);

    summaryFile << n << "\t" << r << std::endl;

    summaryFile.close();

    for (int i = 0; i < table2NumColumns; i++) {
        tableMONIFile.write((char*)tableMONI[i], r * sizeof(int));
    }

    tableMONIFile.write((char*) BWTHeads, r);

    for (int i = 0; i < 4; i++) {
        tableMONIFile.write((char*) preCalcK[i], r * sizeof(*preCalcK[i]));
    }

    tableMONIFile.close();


    for (int j = 0; j < table3NumColumns; j++) {
        for (int i = 0; i < r; i++) {
            tablePhiFile.write((char *) &tablePhi[i][j], sizeof(int));
        }
    }

    for (int j = 0; j < table4NumColumns; j++) {
        for (int i = 0; i < r; i++) {
            tableInversePhiFile.write((char*)&tableInversePhi[i][j], sizeof(int));
        }
    }

    tablePhiFile.close();
    tableInversePhiFile.close();

    return 0;
}

unsigned int constructDataStructures(sdsl::csa_wt<>* csa, sdsl::lcp_bitcompressed<>* lcp, const std::string& file_name) {
    std::ifstream input_file(file_name);
    std::stringstream buffer;
    std::cout << "Reading the input file" << std::endl;
    buffer << input_file.rdbuf();

    std::string text = buffer.str();
    unsigned int n = text.length();

    sdsl::cache_config cc(true);

    std::cout << "Building the BWT" << std::endl;
    construct_im(*csa, text, 1);

    std::cout << "Building the LCP Array" << std::endl;
    construct_im(*lcp, text, 1);

    input_file.close();

    return n;
}

int rowComp(const void *a, const void *b) {
    return(*((int *) a) - *((int *) b));
}

int rowCompArrayFirstElem(const void *a, const void * b) {
    return *((int**)a)[0] - *((int**)b)[0];
}

int rowCompArraySecondElem(const void *a, const void * b) {
    return (*((int**) a))[1] - (*((int**) b))[1];
}

//Predecessor query
unsigned int pred(unsigned int target, const unsigned int *array, int num, int size) {
    int low = 0;
    int high = num - 1;

    while (low != high) {
        int mid = (low + high + 1) / 2;

        if (array[mid * size] <= target) {
            low = mid;
        } else {
            high = mid - 1;
        }
    }

    return low;
}

//Predecessor query based on a specific column in a 2D array
unsigned int pred2D(unsigned int target, unsigned int **array, int column, int num) {
    int low = 0;
    int high = num - 1;

    while (low != high) {
        int mid = (low + high + 1) / 2;

        if (array[mid][column] <= target) {
            low = mid;
        } else {
            high = mid - 1;
        }
    }

    return low;
}
