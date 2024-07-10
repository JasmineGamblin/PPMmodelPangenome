#include <iostream>
#include "functions.hpp"
#include <time.h>
#include <fstream>

using namespace std;


int main()
{

    // test
    string tree_file = "/home/jasminegamblin/Thèse/Ecoli_Refseq/subtree_600.nwk";
    RecTree tree = RecTree::readTree(tree_file);
    cout << "tree OK" << endl;

    // define alignment
    string mat_file = "/home/jasminegamblin/Thèse/Ecoli_Refseq/rep_cat2_cpp.txt";
    PAmatrix mat(mat_file);
    cout << "matrix OK, size " << mat.getMatrix().size() << "x" << mat.getMatrix()[0].size() << endl;

    // build Inference object
    Inference inf(tree, mat);
    cout << "inference init OK" << endl;
    inf.optim_2(1);
    // cout << inf.verif({10,0.2,10,2,5,10,20,0.01}) << endl;
    // cout << inf.verif({10,10,2,5,10,20,0.01,0.005,0.008,0.01}) << endl;
    // cout << inf.nb1({36839, 60.8, 0.002}) << endl;
    // cout << inf.nb2({487386, 179, 9584, 0.002}) << endl;
    cout << "OK!" <<endl;

    
/*     time_t debut_init = time(NULL);
    // read tree
    string tree_file = "../../../EcoliComplete/core_tree_ESCO443_rooted.nwk";
    RecTree tree = RecTree::readTree(tree_file);
    cout << "tree OK" << endl;

    // define alignment
    string mat_file = "../../../EcoliComplete/repartitions_ESCO440_cpp.txt"; // use compressed matrix
    PAmatrix mat(mat_file);
    cout << "matrix OK, size " << mat.getMatrix().size() << "x" << mat.getMatrix()[0].size() << endl;


    // build Inference object
    Inference inf(tree, mat);
    cout << "inference init OK" << endl;
    time_t debut_calc = time(NULL);
    inf.optim();
    time_t fin = time(NULL);
    cout << "time to init: " << debut_calc-debut_init << " seconds" << endl;
    cout << "time to calc: " << fin-debut_calc << " seconds" << endl; */
    

    return 0;
}