//
// Created by igort on 18/05/2023.
//

#ifndef MONI_K_TABLES_H
#define MONI_K_TABLES_H

#include <string>


class MONI_Table {
public:
    unsigned int* head;
    unsigned int* SA_head;
    unsigned int* tail;
    unsigned int* SA_tail;
    unsigned int* mu;
    unsigned int* finger;
    char* BWT_head;
    unsigned char* offset_tail;
    unsigned int* L_tail;
    unsigned char* offset_head;
    unsigned int* L_head;

    MONI_Table(const std::string&, int, bool);
};

class Phi_Table {
public:
    unsigned int* SA_head;
    unsigned int* SA_tail;
    unsigned int* LCP_head;
    unsigned int* finger;

    Phi_Table(const std::string&, int);
};

class Inverse_Phi_Table {
public:
    unsigned int* SA_head;
    unsigned int* SA_tail;
    unsigned int* finger;

    Inverse_Phi_Table(const std::string&, int);
};


#endif //MONI_K_TABLES_H
