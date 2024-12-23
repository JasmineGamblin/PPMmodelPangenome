/*
Main file to perform Maximum Likelihood inference under the PPM model.

This code is distributed under the GNU GPL license.

Author: Jasmine Gamblin
*/

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
    if (stoi(argv[1]) != 0) // fixed or random starting point?
    {
        inf.optim(stoi(argv[1]));
    } else {
        array<double,8> start = {stod(argv[6]), stod(argv[7]), stod(argv[8]), stod(argv[9]),
        stod(argv[10]), stod(argv[11]), stod(argv[12]), stod(argv[13])};
        inf.optim(start);
    }


    // store results
    inf.writeParam(argv[1], argv[4]); // store ML estimates of parameters
    inf.assignCat(argv[5]); // store inferred gene categories
    inf.printPangenomeCompo(); // print expected pangenome composition
    

    time_t fin = time(NULL);
    cout << "DONE :)    Time duration: " << fin-debut << " seconds" << endl;
    return 0;
}
