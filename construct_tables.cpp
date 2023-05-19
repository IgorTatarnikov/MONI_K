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

    auto** table2 = (unsigned int**) calloc(table2NumColumns, sizeof(unsigned int*));
    auto** table3 = (unsigned int**) calloc(r, sizeof(unsigned int*));
    auto** table4 = (unsigned int**) calloc(r, sizeof(unsigned int*));
    auto* BWTHeads = (char*) calloc(r, sizeof(char));
    auto** preCalcK = (unsigned short**) calloc(2, sizeof(unsigned short*));

    for (int i = 0; i < table2NumColumns; i++) {
        table2[i] = (unsigned int*) calloc(r, sizeof(*table2[i]));
    }

    preCalcK[0] = (unsigned short*) calloc(r, sizeof(*preCalcK[0]));
    preCalcK[1] = (unsigned short*) calloc(r, sizeof(*preCalcK[1]));

    unsigned int curr_head;

    unsigned int breakpoint = r / 20;

    std::cout << "Filling table 2 part 1" << std::endl;
    for (int i = 0; i < r; i++) {
        curr_head = run_heads[i];
        table2[0][i] = curr_head;
        table2[1][i] = csa[curr_head];
        BWTHeads[i] = csa.bwt[curr_head];

        if (breakpoint && (i % breakpoint == 0)) {
            std::cout << "Done with " << ((float) i/r) * 100 << "%" << std::endl;
        }
    }

    std::cout << "Filling table 2 part 2" << std::endl;
    for (int i = 0; i < r-1; i++) {
        table2[2][i] = (table2[0][(i + r + 1) % r] + n - 1) % n;
        table2[3][i] = csa[table2[2][i]];
        table2[4][i] = csa.lf[table2[0][i]];

        if (breakpoint && (i % breakpoint == 0)) {
            std::cout << "Done with " << ((float) i/r) * 100 << "%" << std::endl;
        }
    }

    table2[2][r-1] = n;
    table2[3][r-1] = csa[table2[2][r-1]];
    table2[4][r-1] = csa.lf[table2[0][r-1]];


    unsigned int LCPMin;
    unsigned short L_tail;
    unsigned short offset_tail;
    unsigned int start_tail;
    unsigned int max_tail;
    int totalSpan = 2 * k - 2;
    int windowSize = k - 1;
    auto* tempLCPStoreTail = (unsigned int*) calloc(totalSpan, sizeof(unsigned int));

    if (kConstruction) {
        for (int i = 0; i < r; i++) {
            // n+1 accounts for the $ sign added by sdsl
            start_tail = (csa.lf[table2[2][i]] >= (k - 2)) ? csa.lf[table2[2][i]] - k + 2 : 0;
            max_tail = (csa.lf[table2[2][i]] + k) < n + 1 ? csa.lf[table2[2][i]] + k : n + 1;

            unsigned int total_tail = max_tail - start_tail;

            L_tail = 0;
            offset_tail = 0;

            for (unsigned int j = start_tail; j < max_tail; j++) {
                tempLCPStoreTail[j - start_tail] = lcp[j];
            }

            for (int j = 0; j < total_tail - 1; j++) {
                LCPMin = *std::min_element(tempLCPStoreTail + j, tempLCPStoreTail + windowSize + j);
                if (LCPMin > L_tail) {
                    L_tail = LCPMin;
                    offset_tail = windowSize - j;
                }
            }

            preCalcK[0][i] = offset_tail;
            preCalcK[1][i] = L_tail;

            if (breakpoint && (i % breakpoint == 0)) {
                std::cout << "Done with " << ((float) i/r) * 100 << "%" << std::endl;
            }
        }
    }

    std::cout << "Sorting table 2" << std::endl;
    qsort(table2[4], r, sizeof(int), rowComp);

    std::cout << "Predecessor table 2" << std::endl;
    for (int i = 0; i < r; i++) {
        table2[5][i] = pred(table2[4][i], table2[0], r, 1);

        if (breakpoint && (i % breakpoint == 0)) {
            std::cout << "Done with " << ((float) i/r) * 100 << "%" << std::endl;
        }
    }

    std::cout << "Filling table 3 and 4" << std::endl;
    for (int i = 0; i < r; i++) {
        table3[i] = (unsigned int*) calloc(table3NumColumns, sizeof(*table3[i]));
        table3[i][0] = table2[1][i];
        table3[i][1] = table2[3][(i+r-1) % r];
        table3[i][2] = lcp[table2[0][i]];

        table4[i] = (unsigned int*) calloc(table4NumColumns, sizeof(*table4[i]));
        table4[i][0] = table3[i][0];
        table4[i][1] = table3[i][1];

        if (breakpoint && (i % breakpoint == 0)) {
            std::cout << "Done with " << ((float) i/r) * 100 << "%" << std::endl;
        }
    }

    std::cout << "Sorting table 3" << std::endl;
    qsort(table3, r, sizeof(*table3), rowCompArrayFirstElem);

    std::cout << "Sorting table 4" << std::endl;
    qsort(table4, r, sizeof(*table4), rowCompArraySecondElem);

    std::cout << "Predecessor table 3 and 4" << std::endl;
    for (int i = 0; i < r; i++) {
        table3[i][3] = pred2D(table3[i][1], table3, 0, r);
        table4[i][2] = pred2D(table4[i][0], table4, 1, r);

        if (breakpoint && (i % breakpoint == 0)) {
            std::cout << "Done with " << ((float) i/r) * 100 << "%" << std::endl;
        }
    }


    std::cout << "Writing to file" << std::endl;
    std::ofstream summaryFile(file_name.substr(0,file_name.length()-4) + "_Summary_Bin");
    auto table2File = std::fstream(file_name.substr(0,file_name.length()-4) + "_Table2MONI_Bin", std::ios::out | std::ios::binary);
    auto table3File = std::fstream(file_name.substr(0,file_name.length()-4) + "_Table3MONI_Bin", std::ios::out | std::ios::binary);
    auto table4File = std::fstream(file_name.substr(0,file_name.length()-4) + "_Table4MONI_Bin", std::ios::out | std::ios::binary);

    summaryFile << n << "\t" << r << std::endl;

    summaryFile.close();

    for (int i = 0; i < table2NumColumns; i++) {
        table2File.write((char*)table2[i], r * sizeof(int));
    }

    table2File.write((char*) BWTHeads, r);
    table2File.write((char*) preCalcK[0], r * sizeof(*preCalcK[0]));
    table2File.write((char*) preCalcK[1], r * sizeof(*preCalcK[1]));

    table2File.close();


    for (int j = 0; j < table3NumColumns; j++) {
        for (int i = 0; i < r; i++) {
            table3File.write((char *) &table3[i][j], sizeof(int));
        }
    }

    for (int j = 0; j < table4NumColumns; j++) {
        for (int i = 0; i < r; i++) {
            table4File.write((char*)&table4[i][j], sizeof(int));
        }
    }

    table3File.close();
    table4File.close();

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