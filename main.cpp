#include <iostream>
#include <fstream>
#include <string>
#include <sdsl/wavelet_trees.hpp>
#include "LocationSummary.h"

using namespace std;

void LFStep(LocationSummary* currentLocation, int** table2);

int FingerSearch(int startIndex, int** table2);

int LCE(int index1, int index2);

int main() {
    int n = 14;
    int n_pattern = 13;
    int table2NumColumns = 6;
    int table3NumColumns = 4;
    int table4NumColumns = 3;
    int table5NumColumns = 4;

    string table2FileName = "../Table2";
    string table3FileName = "../Table3";
    string table4FileName = "../Table4";
    string table5FileName = "../Table5";

    ifstream table2Input (table2FileName);
    ifstream table3Input (table3FileName);
    ifstream table4Input (table4FileName);
    ifstream table5Input (table5FileName);

    if (!table2Input.is_open() || !table3Input.is_open() || !table4Input.is_open() || !table5Input.is_open()) {
        cout << "Could not open the toy table files" << endl;
        return -1;
    }

    int** table2 = (int**) calloc(n, sizeof(int*));
    int** table3 = (int**) calloc(n, sizeof(int*));
    int** table4 = (int**) calloc(n, sizeof(int*));
    int** table5 = (int**) calloc(n_pattern, sizeof(int*));
    char* table2BWTHead = (char*) malloc(sizeof(char) * n);
    char* table5Pi = (char*) malloc(sizeof(char) * n_pattern);
    char* table5BWT = (char*) malloc(sizeof(char) * n_pattern);

    for (int i = 0; i < n; i++) {
        table2[i] = (int*) calloc(table2NumColumns, sizeof(int));
        table3[i] = (int*) calloc(table3NumColumns, sizeof(int));
        table4[i] = (int*) calloc(table4NumColumns, sizeof(int));

        for (int j = 0; j < table2NumColumns; j++) {
            table2Input >> table2[i][j];
        }

        table2Input >> table2BWTHead[i];

        for (int j = 0; j < table3NumColumns; j++) {
            table3Input >> table3[i][j];
        }

        for (int j = 0; j < table4NumColumns; j++) {
            table4Input >> table4[i][j];
        }
    }

    for (int i = 0; i < n_pattern; i++) {
        table5[i] = (int*) calloc(table5NumColumns, sizeof(int));

        table5Input >> table5[i][0];
        table5Input >> table5Pi[i];

        for (int j = 1; j < table5NumColumns; j++) {
            table5Input >> table5[i][j];
        }

        table5Input >> table5BWT[i];
    }

    table2Input.close();
    table3Input.close();
    table4Input.close();
    table5Input.close();

    string pattern = "TAGATTACATTA";

    auto* currentLocation = (LocationSummary*) malloc(sizeof(LocationSummary));

    currentLocation -> i = 2;
    currentLocation -> l = 8;
    currentLocation -> q = 29;
    currentLocation -> SA = 0;
    currentLocation -> j = 7;

    if (table2BWTHead[currentLocation -> j] != pattern[currentLocation -> i-1]) {
        int upper_boundary = currentLocation -> j-1;
        int lower_boundary = currentLocation -> j+1;

        while (table2BWTHead[upper_boundary] != pattern[currentLocation -> i-1]) {
            upper_boundary--;
        }

        while (table2BWTHead[lower_boundary] != pattern[currentLocation -> i-1]) {
            lower_boundary++;
        }

        if (LCE(currentLocation -> SA, table2[upper_boundary][3]) > LCE(currentLocation -> SA, table2[lower_boundary][1])) {
//            currentLocation -> l = LCE(        currentLocation -> SA, table2[upper_boundary][3]);
            currentLocation -> l = 3;
            currentLocation -> q = table2[upper_boundary][2];
            currentLocation -> SA = table2[upper_boundary][3];
            currentLocation -> j = upper_boundary;
        } else {
            currentLocation -> l = LCE(currentLocation -> SA, table2[lower_boundary][1]);
            currentLocation -> q = table2[lower_boundary][0];
            currentLocation -> SA = table2[lower_boundary][1];
            currentLocation -> j = lower_boundary;
        }

        cout << "Testing " << upper_boundary << " " << lower_boundary << endl;
        cout << currentLocation -> l << "\t";
        cout << currentLocation -> q << "\t";
        cout << currentLocation -> SA << "\t";
        cout << currentLocation -> j << endl;
    }

    LFStep(currentLocation, table2);

    cout << currentLocation -> l << "\t";
    cout << currentLocation -> q << "\t";
    cout << currentLocation -> SA << "\t";
    cout << currentLocation -> j << endl;

    return 0;
}

int LCE(int index1, int index2) {
    return index2 - index1;
}

void LFStep(LocationSummary* currentLocation, int** table2) {
    int pi = 5;
    currentLocation -> l = currentLocation -> l + 1;
    currentLocation -> q = table2[pi][4] + currentLocation -> q - table2[currentLocation -> j][0];
    currentLocation -> SA = currentLocation -> SA - 1;
    currentLocation -> j = FingerSearch(pi, table2);
}

int FingerSearch(int startIndex, int** table2) {


    return 0;
}
