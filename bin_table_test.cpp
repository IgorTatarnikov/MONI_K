#include <fstream>

#include "Tables.h"

int main(int argc, char** argv) {
    int n, r;

    std::ifstream summaryFile("../test_Summary");

    summaryFile >> n;
    summaryFile >> r;

    summaryFile.close();

    std::string table2Name = "../test_Table2MONI_Bin";
    std::string table3Name = "../test_Table3MONI_Bin";
    std::string table4Name = "../test_Table4MONI_Bin";

    MONI_Table tableMONI(table2Name, r);
    Phi_Table tablePhi(table3Name, r);
    Inverse_Phi_Table tableInversePhi(table4Name, r);
}
