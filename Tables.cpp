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
        offset_tail = (unsigned short*) calloc(r, sizeof(*offset_tail));
        L_tail = (unsigned short*) calloc(r, sizeof(*L_tail));
        offset_head = (unsigned short*) calloc(r, sizeof(*offset_head));
        L_head = (unsigned short*) calloc(r, sizeof(*L_head));
        input.read((char*) offset_tail, r * sizeof(*offset_tail));
        input.read((char*) L_tail, r * sizeof(*L_tail));
        input.read((char*) offset_head, r * sizeof(*offset_head));
        input.read((char*) L_head, r * sizeof(*L_head));
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
