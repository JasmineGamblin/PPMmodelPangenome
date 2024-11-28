/*
Function headers to perform Maximum Likelihood inference under the PPM model.

This code is distributed under the GNU GPL license.

Author: Jasmine Gamblin
*/

#include <vector>
#include <map>
#include "objects.hpp"
#include "nelder_mead_dfoptim.h"

using namespace std;


class Inference
{
    public:

    //initialization of the inference object
    static void addBl(RecTree* rt, vector<double>& bl);
    static bool hasAllLeaves(RecTree* rt, vector<string>& rep);
    static int computeMRCA(RecTree* rt, PAmatrix& m, int rep);
    static void depths(RecTree* rt, vector<double>& d);
    Inference(RecTree& tree, PAmatrix mat);

    // transition probabilities for Private genes
    static double pLost(double t, double l);
    static double pRemain(double t, double l);

    // likelihood functions for Private genes
    double likSubtree(int nodeID, double l, double e);
    double likUpper1(int gainNodeID, double l, double e);
    double likRep1(int rep, double l, double e);

    // auxiliary functions for Mobile gene calculations
    static double p01(double t, double g, double l);
    static double p0leaf(int obs, double t, double g, double l, double e1, double e2);

    // likelihood function for Mobile genes
    double likRep2(int rep, double g, double l, double e1, double e2);

    // complete likelihood
    double getPobs1(double l, double e);
    double getPobs2(double g, double l, double e1, double e2);
    double computeI2(double N0, double l0, double i1, double l1, double g2, double l2, vector<double> s);
    double likRepTot(int rep, double N0, double l0, double i1, double l1, double i2, double g2, double l2,
    vector<double> s);
    double logLik(double N0, double l0, double i1, double l1, double g2, double l2, vector<double> s);
    static bool anyAbove(vector<double> start, vector<double> upper);
    void optim(int seed);
    void optim(array<double,8>& start);

    // post-optim
    void writeParam(string seed, string file);
    void assignCat(string mat_file);
    double nb0(double N0, double l0, double s);
    double nb1(double i1, double l1, double s);
    double nb2(double i2, double g2, double l2, double s_gain, double s_loss);
    void printPangenomeCompo();



    private:

    RecTree* tree; // species tree
    PAmatrix mat; // presence/absence matrix
    TreeLikFunctions tLik; // tree structure allowing to store intermediate results in likelihood calculations
    vector<double> branchLength; // branch lengths (indexed by ID of bottom node)
    vector<int> mrca; // mrca of gene patterns
    vector<double> height; // node heights
    vector<int> order; // node IDs ordered by increasing height
    vector<int> freq; // gene frequencies
    vector<vector<int>> subtrees;
    vector<string> leavesID; // leaf labels
    int nRep; // number of identical gene patterns (if the PAmatrix is compressed)
    int nObs; // number of genes
    int nNodes; // number of nodes in the species tree
    int nLeaves; // number of genomes
    double meanGenes; // mean number of genes per genome
    double H; // tree height
    double L; // total branch length of the tree
    double MLEN0; // maximum likelihood estimates of N0
    double MLEl0;
    double MLEi1;
    double MLEl1;
    double MLEi2;
    double MLEg2;
    double MLEl2;
    vector<double> MLEs;
    double MLElik; // maximum value of the log-likelihood
};
