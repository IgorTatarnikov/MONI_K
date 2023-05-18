//
// Created by igort on 18/05/2023.
//

#include <iostream>
#include <fstream>

#include "Tables.h"

MONI_Table::MONI_Table(const std::string &fileName, int r, bool kConstruction) {
    head = (unsigned int*) calloc(r, sizeof(int));
    SA_head = (unsigned int*) calloc(r, sizeof(int));
    tail = (unsigned int*) calloc(r, sizeof(int));
    SA_tail = (unsigned int*) calloc(r, sizeof(int));
    BWT_head = (unsigned int*) calloc(r, sizeof(int));
    mu = (unsigned int*) calloc(r, sizeof(int));
    finger = (unsigned int*) calloc(r, sizeof(int));

    std::ifstream input(fileName, std::ios::in | std::ios::binary);

    input.read((char*) head, r * sizeof(int));
    input.read((char*) SA_head, r * sizeof(int));
    input.read((char*) tail, r * sizeof(int));
    input.read((char*) SA_tail, r * sizeof(int));
    input.read((char*) BWT_head, r * sizeof(int));
    input.read((char*) mu, r * sizeof(int));
    input.read((char*) finger, r * sizeof(int));

    if (kConstruction) {
        offset = (unsigned int*) calloc(r, sizeof(int));
        L = (unsigned int*) calloc(r, sizeof(int));
        input.read((char*) offset, r * sizeof(int));
        input.read((char*) L, r * sizeof(int));
    }

    input.close();
}

Phi_Table::Phi_Table(const std::string & fileName, int r) {
    SA_head = (unsigned int*) calloc(r, sizeof(int));
    SA_tail = (unsigned int*) calloc(r, sizeof(int));
    LCP_head = (unsigned int*) calloc(r, sizeof(int));
    finger = (unsigned int*) calloc(r, sizeof(int));

    std::ifstream input(fileName, std::ios::in | std::ios::binary);

    input.read((char*) SA_head, r * sizeof(int));
    input.read((char*) SA_tail, r * sizeof(int));
    input.read((char*) LCP_head, r * sizeof(int));
    input.read((char*) finger, r * sizeof(int));

    input.close();
}

Inverse_Phi_Table::Inverse_Phi_Table(const std::string & fileName, int r) {
    SA_head = (unsigned int*) calloc(r, sizeof(int));
    SA_tail = (unsigned int*) calloc(r, sizeof(int));
    finger = (unsigned int*) calloc(r, sizeof(int));

    std::ifstream input(fileName, std::ios::in | std::ios::binary);

    input.read((char*) SA_head, r * sizeof(int));
    input.read((char*) SA_tail, r * sizeof(int));
    input.read((char*) finger, r * sizeof(int));

    input.close();
}
