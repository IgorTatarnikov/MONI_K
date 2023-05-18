#ifndef MONI_K_MONI_K_H
#define MONI_K_MONI_K_H

#include <string>
#include <sdsl/wt_blcd.hpp>

#include "Tables.h"

#define TABLE2NUMCOLUMNS 9
#define TABLE3NUMCOLUMNS 4
#define TABLE4NUMCOLUMNS 3
#define STRING_SEPARATOR '#'

unsigned int** readTable(const std::string&, int, int);
int LCE(const std::string&, unsigned int, unsigned int);
void updateRow(const std::string&, unsigned int, MONI_Table, unsigned int&, unsigned int&, unsigned int&, unsigned int&, int);
void LFStep(MONI_Table, unsigned int&, unsigned int&, unsigned int&, sdsl::wt_blcd<>&, int);
unsigned int exponentialSearch(unsigned int*, unsigned int, unsigned int, int);
unsigned int pred2D(unsigned int, unsigned int*, unsigned int, unsigned int, int);
unsigned int phi(Phi_Table, unsigned int&, unsigned int, int r);
unsigned int inversePhi(Inverse_Phi_Table, unsigned int&, unsigned int, int r);
unsigned int LCPStep(Phi_Table, unsigned int, unsigned int);
unsigned int** preCalcMONI(MONI_Table, Phi_Table, Inverse_Phi_Table, int, sdsl::wt_blcd<>&, const std::string&, const std::string&, int, int);

#endif //MONI_K_MONI_K_H
