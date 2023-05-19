#ifndef MONI_K_CONSTRUCT_TABLES_H
#define MONI_K_CONSTRUCT_TABLES_H

#include <sdsl/lcp_byte.hpp>
#include <sdsl/lcp_bitcompressed.hpp>
#include <sdsl/suffix_arrays.hpp>
#include <string>

int rowComp(const void*, const void*);
int rowCompArrayFirstElem(const void*, const void*);
int rowCompArraySecondElem(const void*, const void*);
unsigned int pred(unsigned int, const unsigned int*, int, int);
unsigned int pred2D(unsigned int, unsigned int**, int, int);
unsigned int constructDataStructures(sdsl::csa_wt<> *, sdsl::lcp_bitcompressed<> *, const std::string&);

#endif //MONI_K_CONSTRUCT_TABLES_H
