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
    double computeI2(const array<double,7>& param);
    double computeI2_0(const array<double,2>& param);
    double computeI2_1(const array<double,2>& param);
    double computeI2_2(const array<double,3>& param);
    double likRepTot(int rep, vector<double>& param, array<double,4>& error);
    double likRepTot_0(int rep, vector<double>& param, double error);
    double likRepTot_1(int rep, vector<double>& param, double error);
    double likRepTot_2(int rep, vector<double>& param, array<double,2>& error);
    double verif(vector<double> param);
    double logLik(const array<double,7>& param);
    double logLik_0(const array<double,2>& param);
    double logLik_1(const array<double,2>& param);
    double logLik_2(const array<double,3>& param);
    void optim(int seed);
    void optim_0(int seed);
    void optim_1(int seed);
    void optim_2(int seed);

    // post-optim
    void writeParam(string seed, string file);
    // void writeParam(string seed, string file, double N0, double l0, double i1, double l1, double g2, double l2,
    // double e);
    // void writeParam_2(string seed, string file, double g2, double l2, double e);
    void assignCat(string mat_file);
    void computeProbaCat(string file, double N0, double l0, double i1, double l1, double i2, double g2, double l2,
    double eps0, double eps1, double eps2);
    void computeProbaPattern(string file, double N0, double l0, double i1, double l1, double i2, double g2, double l2,
    double eps0, double eps1, double eps2);

    // others
    double nb0(vector<double> par);
    double nb1(vector<double> par);
    double nb2(vector<double> par);
    double sumGains(int st, double g2);
    double expectedGains(double g2);
    double expectedTime1_old(int rep, double l1, double eps1);
    double expectedTime1(int rep, double l1, double eps1);
    double expectedTime2(int rep, double g2, double l2, double eps2);
    void expectedTime(string file, double l1, double g2, double l2, double eps1, double eps2);
    vector<double> expectedTime2_times();
    vector<double> expectedTime2_aux(int rep, double g2, double l2, double eps2);
    void expectedTime2_store(string file, double g2, double l2, double eps2);

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
    vector<double> MLEpar;
    vector<double> MLEerr;
    double MLElik;
};
