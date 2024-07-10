#include <vector>
#include <map>
#include "objects.hpp"
#include "nelder_mead_dfoptim.h"

using namespace std;


class Inference
{
    public:

    static void addBl(RecTree* rt, vector<double>& bl);
    static bool hasAllLeaves(RecTree* rt, vector<string>& rep);
    static int computeMRCA(RecTree* rt, PAmatrix& m, int rep);
    static void depths(RecTree* rt, vector<double>& d);
    Inference(RecTree& tree, PAmatrix mat);

    // transition probabilities for cat. 1 genes
    static double pLost(double t, double l);
    static double pRemain(double t, double l);

    // likelihood functions for cat. 1
    double likSubtree(int nodeID, double l, double e);
    double likUpper1(int gainNodeID, double l, double e);
    double likRep1(int rep, double l, double e);

    // auxiliary functions for cat. 2 calculations
    static double p01(double t, double g, double l);
    static double p0leaf(int obs, double t, double g, double l, double e1, double e2);

    // likelihood function for cat. 2
    double likRep2(int rep, double g, double l, double e1, double e2);

    // complete lik
    double getPobs1(double l, double e);
    double getPobs2(double g, double l, double e1, double e2);
    double computeI2(double N0, double l0, double i1, double l1, double g2, double l2, vector<double> eps);
    double likRepTot(int rep, double N0, double l0, double i1, double l1, double i2, double g2, double l2,
    vector<double> eps);
    double verif(double l0, double l1, double g2, double l2, vector<double> eps);
    double logLik(double N0, double l0, double i1, double l1, double g2, double l2, vector<double> eps);
    static bool anyAbove(vector<double> start, vector<double> upper);
    void optim7(int seed);
    void optim7(array<double,7>& start);
    void optim9(int seed);
    void optim9(array<double,9>& start);

    // post-optim
    void writeParam(string seed, string file);
    void assignCat(string mat_file);
    double nb0(double N0, double l0, double eps);
    double nb1(double i1, double l1, double eps);
    double nb2(double i2, double g2, double l2, double eps);
    void printPangenomeCompo();



    private:

    RecTree* tree;
    PAmatrix mat;
    TreeLikFunctions tLik;
    vector<double> branchLength;
    vector<int> mrca;
    vector<double> height;
    vector<int> order;
    vector<int> freq;
    vector<vector<int>> subtrees;
    vector<string> leavesID;
    int nRep;
    int nObs;
    int nNodes;
    int nLeaves;
    double meanGenes;
    double H;
    double L;
    double MLEN0;
    double MLEl0;
    double MLEi1;
    double MLEl1;
    double MLEi2;
    double MLEg2;
    double MLEl2;
    vector<double> MLEeps;
    double MLElik;
};
