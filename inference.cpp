#include <iostream>
#include "functions.hpp"
#include <time.h>
#include <fstream>

using namespace std;


int main(int argc, char* argv[])
{
    
    time_t debut = time(NULL);
    
    // store tree
    RecTree tree = RecTree::readTree(argv[2]);
    cout << "tree OK" << endl;

    // store presence/absence matrix
    PAmatrix mat(argv[3]);
    cout << "matrix OK, size " << mat.getMatrix().size() << "x" << mat.getMatrix()[0].size() << endl;

    // build Inference object
    Inference inf(tree, mat);
    cout << "inference initialization OK" << endl;
    
    // inference
    if (argv[4] == "7") // 7 or 9 parameters?
    {
        if (argv[1] != "0") // fixed or random starting point?
        {
            inf.optim7(stoi(argv[1]));
        } else {
            array<double,7> start = {stod(argv[7]), stod(argv[8]), stod(argv[9]), stod(argv[10]),
            stod(argv[11]), stod(argv[12]), stod(argv[13])};
            inf.optim7(start);
        }
    }
    else if (argv[4] == "9") {
        if (argv[1] != "0") // fixed or random starting point?
        {
            inf.optim9(stoi(argv[1]));
        } else {
            array<double,9> start = {stod(argv[7]), stod(argv[8]), stod(argv[9]), stod(argv[10]),
            stod(argv[11]), stod(argv[12]), stod(argv[13]), stod(argv[14]), stod(argv[15])};
            inf.optim9(start);
        }
    } else {
        cout << "incorrect number of parameters (valid values: 7 or 9)" << endl;
    }
    inf.writeParam(argv[1], argv[5]);
    inf.assignCat(argv[6]);
    inf.printPangenomeCompo();
    
    time_t fin = time(NULL);
    cout << "DONE :)    Time duration: " << fin-debut << " seconds" << endl;
    return 0;
}
