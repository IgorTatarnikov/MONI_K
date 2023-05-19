#include <iostream>
#include <fstream>

#include "Tables.h"

MONI_Table::MONI_Table(const std::string &fileName, int r, bool kConstruction) {
    head = (unsigned int*) calloc(r, sizeof(*head));
    SA_head = (unsigned int*) calloc(r, sizeof(*SA_head));
    tail = (unsigned int*) calloc(r, sizeof(*tail));
    SA_tail = (unsigned int*) calloc(r, sizeof(*SA_tail));
    mu = (unsigned int*) calloc(r, sizeof(*mu));
    finger = (unsigned int*) calloc(r, sizeof(*finger));

    BWT_head = (char*) calloc(r, 1);

    std::ifstream input(fileName, std::ios::in | std::ios::binary);

    input.read((char*) head, r * sizeof(int));
    input.read((char*) SA_head, r * sizeof(int));
    input.read((char*) tail, r * sizeof(int));
    input.read((char*) SA_tail, r * sizeof(int));
    input.read((char*) mu, r * sizeof(int));
    input.read((char*) finger, r * sizeof(int));

    input.read((char*) BWT_head, r);

    if (kConstruction) {
        offset = (unsigned short*) calloc(r, sizeof(*offset));
        L = (unsigned short*) calloc(r, sizeof(*L));
        input.read((char*) offset, r * sizeof(*offset));
        input.read((char*) L, r * sizeof(*L));
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
