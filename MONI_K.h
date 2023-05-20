#ifndef MONI_K_MONI_K_H
#define MONI_K_MONI_K_H

#include <string>
#include <sdsl/wt_blcd.hpp>
#include <vector>

#include "Tables.h"

#define TABLE2NUMCOLUMNS 9
#define TABLE3NUMCOLUMNS 4
#define TABLE4NUMCOLUMNS 3
#define STRING_SEPARATOR '#'
#define KCONSTRUCTION true

std::vector<std::string> readTestStrings(const std::string& fileName);
int LCE(const std::string&, unsigned int, unsigned int);
void updateRow(const std::string&, char, MONI_Table, unsigned int&, unsigned int&, unsigned int&, unsigned int&, int);
void LFStep(MONI_Table, unsigned int&, unsigned int&, unsigned int&, sdsl::wt_blcd<>&, int);
unsigned int exponentialSearch(unsigned int*, unsigned int, unsigned int, int);
unsigned int pred2D(unsigned int, unsigned int*, unsigned int, unsigned int, int);
unsigned int phi(Phi_Table, unsigned int&, unsigned int, int r);
unsigned int inversePhi(Inverse_Phi_Table, unsigned int&, unsigned int, int r);
unsigned int LCPStep(Phi_Table, unsigned int, unsigned int);
unsigned int** preCalcMONI(MONI_Table, unsigned int, int, sdsl::wt_blcd<>&, const std::string&, const std::string&, int, int);

#endif //MONI_K_MONI_K_H
