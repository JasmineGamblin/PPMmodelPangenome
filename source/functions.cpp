/*
Functions to perform Maximum Likelihood inference under the PPM model.

This code is distributed under the GNU GPL license.

Author: Jasmine Gamblin
*/

#include "functions.hpp"
#include <cmath>
#include <queue>
#include <algorithm>
#include <numeric>
#include <stack>
#include <iostream>
#include <random>
#include <iomanip>
#include "oneapi/tbb.h"
#include <fstream>

using namespace std;



//initialization of the inference object
void Inference::addBl(RecTree* rt, vector<double>& bl)
{
    if (rt->getLeaf().size() == 0)
    {
        bl.push_back(rt->getBranchLengths().first);
        addBl(rt->getLeftChild(), bl);
        bl.push_back(rt->getBranchLengths().second);
        addBl(rt->getRightChild(), bl);
    }
}
bool Inference::hasAllLeaves(RecTree* rt, vector<string>& repLeaves)
{
    vector<string> treeLeaves;
    rt->getLeaves(treeLeaves);
    for (int i=0; i<repLeaves.size(); i++)
    {
        if (find(treeLeaves.begin(), treeLeaves.end(), repLeaves[i]) == treeLeaves.end()) {return false;}
    }
    return true;
}
int Inference::computeMRCA(RecTree* rt, PAmatrix& m, int i)
{
    vector<string> presentLeaves = m.extractPresentLeaves(i);

    // compute mrca of leaves containing the gene i
    int mrca;
    queue<RecTree*> q;
    q.push(rt);
    while (!q.empty())
    {
        if (hasAllLeaves(q.front(), presentLeaves)) {mrca = q.front()->getNodeID();}
        if (q.front()->getLeaf().size() == 0)
        {
            q.push(q.front()->getLeftChild());
            q.push(q.front()->getRightChild());
        }
        q.pop();
    }
    return mrca;
}
void Inference::depths(RecTree* rt, vector<double>& d)
{
    if (rt->getLeaf().size() == 0)
    {
        d[rt->getLeftChild()->getNodeID()] = d[rt->getNodeID()] + rt->getBranchLengths().first;
        d[rt->getRightChild()->getNodeID()] = d[rt->getNodeID()] + rt->getBranchLengths().second;
        depths(rt->getLeftChild(), d);
        depths(rt->getRightChild(), d);
    }
}
Inference::Inference(RecTree& rt, PAmatrix m)
{
    tree = &rt;

    mat = m;

    nRep = mat.getMatrix().size();
    nObs = accumulate(mat.getCounts().begin(), mat.getCounts().end(), 0);
    nNodes = tree->addNodeID();
    nLeaves = (nNodes + 1)/2;
    L = tree->treeLength();

    // add a row of zeros to compute proba to be lost
    mat.add0(); 

    // build tree likelihood subfunctions
    tLik = TreeLikFunctions(rt, mat);

    // store branchlengths
    branchLength = {0.};
    addBl(tree, branchLength);

    // store mrca
    mrca = vector<int>(nRep, 0);
    tbb::detail::d1::parallel_for(0, nRep, [&](int i)
    {
        mrca[i] = computeMRCA(tree, mat, i);
    });

    // store H and height
    height = vector<double>(nNodes);
    depths(tree, height);
    H = *max_element(height.begin(), height.end());
    transform(height.begin(), height.end(), height.begin(), bind1st(minus<double>(), H));

    // store order
    order = vector<int>(nNodes);   
    for (int i = 0 ; i < nNodes ; i++) {order[i] = i;}
    sort(order.begin(), order.end(), [&](const int& a, const int& b) {return (height[a] < height[b]);});

    // store frequency of genes
    freq = vector<int>(nRep, 0);
    meanGenes = 0.;
    tbb::detail::d1::parallel_for(0, nRep, [&](int i)
    {
        freq[i] = accumulate(mat.getMatrix()[i].begin(), mat.getMatrix()[i].end(), 0);
    });
    for (int i(0); i<nRep; i++)
    {
        meanGenes += freq[i]*mat.getCounts()[i];
    }
    meanGenes /= nLeaves;
    
    // store subtrees (IDs of subtrees roots)
    subtrees = vector<vector<int>>();
    vector<int> current = {0};
    for (int i = nNodes-1; i > 0; i--)
    {
        current.erase(remove(current.begin(), current.end(), order[i]), current.end());
        if (!tLik.getLikFunctions()[order[i]].isLeaf()) {
           current.push_back(tLik.getLikFunctions()[order[i]].getChildren().first);
           current.push_back(tLik.getLikFunctions()[order[i]].getChildren().second);
        }
        subtrees.push_back(current);
    }
    reverse(subtrees.begin(), subtrees.end());

    // store leavesID
    leavesID = vector<string>(nNodes); // leaf name is stored at the leaf ID ("" for internal nodes)
    tree->getLeavesID(leavesID);

    MLEN0 = 0.;
    MLEl0 = 0.;
    MLEi1 = 0.;
    MLEl1 = 0.;
    MLEi2 = 0.;
    MLEg2 = 0.;
    MLEl2 = 0.;
    MLEeps = {0., 0., 0.};
    MLElik = 0.;
}



// transition probabilities for Private genes
double Inference::pLost(double t, double l) {return log(1-exp(-l*t));}
double Inference::pRemain(double t, double l) {return -l*t;}



// likelihood functions for Private genes
double Inference::likSubtree(int nodeID, double l, double e)
{
    if (tLik.getLikFunctions()[nodeID].isLeaf())
    {
        return pRemain(branchLength[nodeID], l) + log(e);
    } else {
        return pRemain(branchLength[nodeID], l) + tLik.evaluateRootRep(nodeID, nRep, {0.,1.}, 0., l, 0., e, 1);
    }
}
double Inference::likUpper1(int gainNodeID, double l, double e)
{
    /* traverse the tree from root to gainNode and for each node
    strictly between the 2 calculate likSubtree */
    
    stack<int> nodes;
    int node = 0;
    pair<int, int> children = tLik.getLikFunctions()[node].getChildren();
    if (gainNodeID >= children.second) {node = children.second;}
    else {node = children.first;}
    while (node != gainNodeID)
    {
        nodes.push(node);
        children = tLik.getLikFunctions()[node].getChildren();
        if (gainNodeID >= children.second) {node = children.second;}
        else {node = children.first;}
    }
    
    int child;
    double sum_bl = branchLength[gainNodeID];
    double likSubtrees(0);
    double likU = log(1-exp(-l*branchLength[gainNodeID]));
    while (!nodes.empty())
    {
        children = tLik.getLikFunctions()[nodes.top()].getChildren();
        if (gainNodeID >= children.second) {child = children.first;}
        else {child = children.second;}
        likSubtrees += logSumExp(pLost(branchLength[child], l), this->likSubtree(child, l, e));
        likU = logSumExp(likU, log(1-exp(-l*branchLength[nodes.top()])) + likSubtrees - l*sum_bl);
        sum_bl += branchLength[nodes.top()];
        nodes.pop();
    }
    return likU;
}
double Inference::likRep1(int rep, double l, double e)
{
    // mrca is the root
    if (mrca[rep] == 0) {return -INFINITY;} 
    
    // mrca is a leaf
    else if (tLik.getLikFunctions()[mrca[rep]].isLeaf()) {return this->likUpper1(mrca[rep], l, e)+log(1-e);} 
    
    else
    {
        return tLik.evaluateRootRep(mrca[rep], rep, {0.,1.}, 0., l, 0., e, 1) + this->likUpper1(mrca[rep], l, e);
    }
}



// auxiliary functions for Mobile genes calculations
double Inference::p01(double t, double g, double l)
{
  if (g == 0.) {return 0.;}
  else {return g*(1-exp(-t*(g+l)))/(g+l);}
}
double Inference::p0leaf(int obs, double t, double g, double l, double e1, double e2)
{
  if (obs == 1)
  {
    return p01(t, g, l)*(1-e2) + (1-p01(t, g, l))*e1;
  } else {
    return p01(t, g, l)*e2 + (1-p01(t, g, l))*(1-e1);
  }
}



// likelihood function for Mobile genes
double Inference::likRep2(int rep, double g, double l, double e1, double e2)
{
    double sum,integral,w,f_aux,t,p1;
    int nSub,st;
    sum = 0;
    for (int rank = nLeaves-1; rank < nNodes-1; rank++) // begin at the last leaf because tree is ultrametric
    {
        integral = 0.;
        w = height[order[rank+1]] - height[order[rank]];
        if (w > 0.)
        {
            nSub = max(1,(int)ceil(w*100/H));
            for (int i = 0; i < nSub; i++)
            {
                f_aux = 0.;
                t = height[order[rank]]+(2*i+1)*w/(2*nSub);
                for (int j = 0; j < subtrees[rank].size(); j++)
                {
                    st = subtrees[rank][j];
                    if (tLik.getLikFunctions()[st].isLeaf())
                    {
                        f_aux += log(p0leaf(mat.getMatrix()[rep][mat.getColnames()[leavesID[st]]],
                        t-height[st], g, l, e1, e2));
                    } else {
                        p1 = p01(t - height[st], g, l);
                        f_aux += tLik.evaluateRootRep(st, rep, {1-p1,p1}, g, l, e1, e2, 2);
                    }
                }
                integral += exp(f_aux);
            }
            sum += integral*w/nSub;
        }
    }
    return log(sum);
}



// complete likelihood
double Inference::getPobs1(double l, double e)
{
    double sum(0);
    for (int i=1; i<branchLength.size(); i++)
    {
        if (tLik.getLikFunctions()[i].isLeaf())
        {
            sum += (1-exp(-l*branchLength[i]))*(1-e);
        }
        else
        {
            sum += (1-exp(-l*branchLength[i]))*(1-exp(tLik.evaluateRootRep(i, nRep, {0.,1.}, 0., l, 0., e, 1)));
        }
    }
    return sum/l;
}
double Inference::getPobs2(double g, double l, double e1, double e2)
{
    double sum,integral,w,f_aux,t,p1;
    int nSub,st;
    sum = 0;
    for (int rank = nLeaves-1; rank < nNodes-1; rank++) // begin at the last leaf because tree is ultrametric
    {
        integral = 0;
        w = height[order[rank+1]] - height[order[rank]];
        if (w > 0) // when the 2 nodes are not at the same height
        {
            nSub = max(1,(int)ceil(w*100/H));
            for (int i = 0; i < nSub; i++)
            {
                f_aux = 0;
                t = height[order[rank]]+(2*i+1)*w/(2*nSub);
                for (int j = 0; j < subtrees[rank].size(); j++)
                {
                    st = subtrees[rank][j];
                    if (tLik.getLikFunctions()[st].isLeaf())
                    {
                        f_aux += log(p0leaf(0, t-height[st], g, l, e1, e2));
                    } else {
                        p1 = p01(t - height[st], g, l);
                        f_aux += tLik.evaluateRootRep(st, nRep, {1-p1,p1}, g, l, e1, e2, 2);
                    }
                }
                integral += exp(f_aux);
            }
            sum += integral*w/nSub;
        }
    }  
    return H-sum;
}
double Inference::computeI2(double N0, double l0, double i1, double l1, double g2, double l2,
vector<double> eps)
{   
    double p0;
    if (l0 == 0.) {p0 = 1-pow(eps[0], nLeaves);}
    else {p0 = 1-exp(tLik.evaluateRootRep(0, nRep, {0,1}, 0., l0, 0., eps[0], 0));}
    double p1 = 1-exp(tLik.evaluateRootRep(0, nRep, {0,1}, 0., l1, 0., eps[1], 1));
    double A = this->getPobs1(l1, eps[1]);
    double B = this->getPobs2(g2, l2, eps[2], eps[2]);
    return (nObs-N0*p0-i1/l1*p1-i1*A)/B;
}
double Inference::likRepTot(int rep, double N0, double l0, double i1, double l1, double i2, double g2,
double l2, vector<double> eps)
{
    double lik0r = log(N0/nObs);
    if (l0 == 0. && (eps[0] > 0. || freq[rep] < nLeaves)) {
        lik0r += freq[rep]*log(1-eps[0]) + (nLeaves-freq[rep])*log(eps[0]);
    } else if (l0 > 0.) {
        lik0r += tLik.evaluateRootRep(0, rep, {0,1}, 0., l0, 0., eps[0], 0);
    }
    double lik1r = log(i1/(l1*nObs)) + tLik.evaluateRootRep(0, rep, {0,1}, 0., l1, 0., eps[1], 1);
    double lik1 = log(i1/(l1*nObs)) + this->likRep1(rep, l1, eps[1]);
    double lik2 = log(i2/nObs) + this->likRep2(rep, g2, l2, eps[2], eps[2]);
    return logSumExp(logSumExp(lik0r, lik1r), logSumExp(lik1, lik2));
}
double Inference::logLik(double N0, double l0, double i1, double l1, double g2, double l2, vector<double> eps)
{
    double i2 = this->computeI2(N0, l0, i1, l1, g2, l2, eps);
    if (i2 < 0.)
    {
        return -1e9*1. + i2;
    } else {
        vector<double> lik(nRep, 0.);
        tbb::detail::d1::parallel_for(0, nRep, [&](int i){
            lik[i] = this->likRepTot(i, N0, l0, i1, l1, i2, g2, l2, eps)*mat.getCounts()[i];
        });
        return accumulate(lik.begin(), lik.end(), 0.);
    }
}
bool Inference::anyAbove(vector<double> start, vector<double> upper)
{
    bool res = false;
    for (int i(0); i<start.size(); i++)
    {
        if (start[i]>upper[i]) {res = true;}
    }
    return res;
}
void Inference::optim7(int seed)
{
    const size_t nbPar = 7;
    
    mt19937 rng(seed);    // Random-number engine used (Mersenne-Twister in this case)
    uniform_int_distribution<int> uniN0(meanGenes/3,meanGenes);
    uniform_int_distribution<int> unii1(nObs/(L*5),nObs/L);
    exponential_distribution<> expl0(L);
    exponential_distribution<> expl1(L/50);
    exponential_distribution<> expl2(L/2500);
    exponential_distribution<> expe(100.);
    array<double,nbPar> start = {uniN0(rng)*1., expl0(rng), unii1(rng)*1., expl1(rng), expl1(rng),
    expl2(rng), expe(rng)};
    array<double,nbPar> lower = {0., 0., 0., 0., 0., 0., 0.};
    array<double,nbPar> upper = {meanGenes*2., 10/L, 3.*nObs/L, 500/L, 500/L, 25000/L, 0.1};
    while (this->computeI2(start[0], start[1], start[2], start[3], start[4], start[5],
    {start[6], start[6], start[6]})
    < 0. || anyAbove({start.begin(), start.end()}, {upper.begin(), upper.end()}))
    {
        start = {uniN0(rng)*1., expl0(rng), unii1(rng)*1., expl1(rng), expl1(rng), expl2(rng), expe(rng)};
    }

    nelder_mead_result<double,nbPar> res = nelder_mead<double,nbPar>([&](const array<double,nbPar>& par, int nbEval){
        cout << nbEval << "th iteration: ";
        tLik.clearLikValues();
        cout << par[0] << ' ' << par[1] << ' ' << par[2] << ' ' << par[3] << ' ' << par[4] << ' ' << par[5] << ' '
        << par[6] << ' ';
        double ll = this->logLik(par[0], par[1], par[2], par[3], par[4], par[5], {par[6], par[6], par[6]});
        cout << setprecision(9) << ll << endl;
        return -ll;
        }, start, lower, upper, 0.0001, 3000, 10);
    MLEN0 = res.xmin[0];
    MLEl0 = res.xmin[1];
    MLEi1 = res.xmin[2];
    MLEl1 = res.xmin[3];
    MLEg2 = res.xmin[4];
    MLEl2 = res.xmin[5];
    MLEeps = {res.xmin[6], res.xmin[6], res.xmin[6]};
    MLEi2 = computeI2(MLEN0, MLEl0, MLEi1, MLEl1, MLEg2, MLEl2, MLEeps);
    MLElik = -res.ymin;
    cout << "Found minimum: " << MLEN0 << ' ' << MLEl0 << ' ' << MLEi1 << ' ' << MLEl1 << ' '
    << MLEi2 << ' ' << MLEg2 << ' ' << MLEl2 << ' ' << MLEeps[0] << ' ' << MLEeps[1] << ' ' << MLEeps[2] << endl;
    cout << "upper limit:   " << upper[0] << ' ' << upper[1] << ' ' << upper[2] << ' ' << upper[3] << '\t' << upper[4]
    << ' ' << upper[5] << ' ' << upper[6] << '\t' << '\t' << endl;
    cout << "logLik value: " << MLElik << endl;
    cout << "Number of function evaluation: " << res.icount << endl;
    cout << "Message: " << res.message << endl;
}
void Inference::optim7(array<double,7>& start)
{
    const size_t nbPar = 7;
    
    array<double,nbPar> lower = {0., 0., 0., 0., 0., 0., 0.};
    array<double,nbPar> upper = {2.*meanGenes, 10/L, 3.*nObs/L, 500/L, 500/L, 25000/L, 0.1};

    nelder_mead_result<double,nbPar> res = nelder_mead<double,nbPar>([&](const array<double,nbPar>& par, int nbEval){
        cout << nbEval << "th iteration: ";
        tLik.clearLikValues();
        cout << par[0] << ' ' << par[1] << ' ' << par[2] << ' ' << par[3] << ' ' << par[4] << ' ' << par[5] << ' '
        << par[6] << ' ';
        double ll = this->logLik(par[0], par[1], par[2], par[3], par[4], par[5], {par[6], par[6], par[6]});
        cout << setprecision(9) << ll << endl;
        return -ll;
        }, start, lower, upper, 0.0001, 3000, 10);
    MLEN0 = res.xmin[0];
    MLEl0 = res.xmin[1];
    MLEi1 = res.xmin[2];
    MLEl1 = res.xmin[3];
    MLEg2 = res.xmin[4];
    MLEl2 = res.xmin[5];
    MLEeps = {res.xmin[6], res.xmin[6], res.xmin[6]};
    MLEi2 = computeI2(MLEN0, MLEl0, MLEi1, MLEl1, MLEg2, MLEl2, MLEeps);
    MLElik = -res.ymin;
    cout << "Found minimum: " << MLEN0 << ' ' << MLEl0 << ' ' << MLEi1 << ' ' << MLEl1 << ' '
    << MLEi2 << ' ' << MLEg2 << ' ' << MLEl2 << ' ' << MLEeps[0] << ' ' << MLEeps[1] << ' ' << MLEeps[2] << endl;
    cout << "upper limit:   " << upper[0] << ' ' << upper[1] << ' ' << upper[2] << ' ' << upper[3] << '\t' << upper[4]
    << ' ' << upper[5] << ' ' << upper[6] << '\t' << '\t' << endl;
    cout << "logLik value: " << MLElik << endl;
    cout << "Number of function evaluation: " << res.icount << endl;
    cout << "Message: " << res.message << endl;
}
void Inference::optim9(int seed)
{
    const size_t nbPar = 9;
    
    mt19937 rng(seed);    // Random-number engine used (Mersenne-Twister in this case)
    uniform_int_distribution<int> uniN0(meanGenes/3,meanGenes);
    uniform_int_distribution<int> unii1(nObs/(L*5),nObs/L);
    exponential_distribution<> expl0(L);
    exponential_distribution<> expl1(L/50);
    exponential_distribution<> expl2(L/2500);
    exponential_distribution<> expe(100.);
    array<double,nbPar> start = {uniN0(rng)*1., expl0(rng), unii1(rng)*1., expl1(rng), expl1(rng),
    expl2(rng), expe(rng), expe(rng), expe(rng)};
    array<double,nbPar> lower = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
    array<double,nbPar> upper = {meanGenes*2., 10/L, 3.*nObs/L, 500/L, 500/L, 25000/L, 0.1, 0.1, 0.1};
    while (this->computeI2(start[0], start[1], start[2], start[3], start[4], start[5],
    {start[6], start[7], start[8]})
    < 0. || anyAbove({start.begin(), start.end()}, {upper.begin(), upper.end()}))
    {
        start = {uniN0(rng)*1., expl0(rng), unii1(rng)*1., expl1(rng), expl1(rng), expl2(rng),
        expe(rng), expe(rng), expe(rng)};
    }

    nelder_mead_result<double,nbPar> res = nelder_mead<double,nbPar>([&](const array<double,nbPar>& par, int nbEval){
        cout << nbEval << "th iteration: ";
        tLik.clearLikValues();
        cout << par[0] << ' ' << par[1] << ' ' << par[2] << ' ' << par[3] << ' ' << par[4] << ' ' << par[5] << ' '
        << par[6] << ' ' << par[7] << ' ' << par[8] << ' ';
        double ll = this->logLik(par[0], par[1], par[2], par[3], par[4], par[5], {par[6], par[7], par[8]});
        cout << setprecision(9) << ll << endl;
        return -ll;
        }, start, lower, upper, 0.0001, 3000, 10);
    MLEN0 = res.xmin[0];
    MLEl0 = res.xmin[1];
    MLEi1 = res.xmin[2];
    MLEl1 = res.xmin[3];
    MLEg2 = res.xmin[4];
    MLEl2 = res.xmin[5];
    MLEeps = {res.xmin[6], res.xmin[7], res.xmin[8]};
    MLEi2 = computeI2(MLEN0, MLEl0, MLEi1, MLEl1, MLEg2, MLEl2, MLEeps);
    MLElik = -res.ymin;
    cout << "Found minimum: " << MLEN0 << ' ' << MLEl0 << ' ' << MLEi1 << ' ' << MLEl1 << ' '
    << MLEi2 << ' ' << MLEg2 << ' ' << MLEl2 << ' ' << MLEeps[0] << ' ' << MLEeps[1] << ' ' << MLEeps[2] << endl;
    cout << "upper limit:   " << upper[0] << ' ' << upper[1] << ' ' << upper[2] << ' ' << upper[3] << '\t' << upper[4]
    << ' ' << upper[5] << ' ' << upper[6] << ' ' << upper[7] << ' ' << upper[8] << endl;
    cout << "logLik value: " << MLElik << endl;
    cout << "Number of function evaluation: " << res.icount << endl;
    cout << "Message: " << res.message << endl;
}
void Inference::optim9(array<double,9>& start)
{
    const size_t nbPar = 9;
    
    array<double,nbPar> lower = {0., 0., 0., 0., 0., 0., 0., 0., 0.};
    array<double,nbPar> upper = {meanGenes*2., 10/L, 3.*nObs/L, 500/L, 500/L, 25000/L, 0.1, 0.1, 0.1};

    nelder_mead_result<double,nbPar> res = nelder_mead<double,nbPar>([&](const array<double,nbPar>& par, int nbEval){
        cout << nbEval << "th iteration: ";
        tLik.clearLikValues();
        cout << par[0] << ' ' << par[1] << ' ' << par[2] << ' ' << par[3] << ' ' << par[4] << ' ' << par[5] << ' '
        << par[6] << ' ' << par[7] << ' ' << par[8] << ' ';
        double ll = this->logLik(par[0], par[1], par[2], par[3], par[4], par[5], {par[6], par[7], par[8]});
        cout << setprecision(9) << ll << endl;
        return -ll;
        }, start, lower, upper, 0.0001, 3000, 10);
    MLEN0 = res.xmin[0];
    MLEl0 = res.xmin[1];
    MLEi1 = res.xmin[2];
    MLEl1 = res.xmin[3];
    MLEg2 = res.xmin[4];
    MLEl2 = res.xmin[5];
    MLEeps = {res.xmin[6], res.xmin[7], res.xmin[8]};
    MLEi2 = computeI2(MLEN0, MLEl0, MLEi1, MLEl1, MLEg2, MLEl2, MLEeps);
    MLElik = -res.ymin;
    cout << "Found minimum: " << MLEN0 << ' ' << MLEl0 << ' ' << MLEi1 << ' ' << MLEl1 << ' '
    << MLEi2 << ' ' << MLEg2 << ' ' << MLEl2 << ' ' << MLEeps[0] << ' ' << MLEeps[1] << ' ' << MLEeps[2] << endl;
    cout << "upper limit:   " << upper[0] << ' ' << upper[1] << ' ' << upper[2] << ' ' << upper[3] << '\t' << upper[4]
    << ' ' << upper[5] << ' ' << upper[6] << ' ' << upper[7] << ' ' << upper[8] << endl;
    cout << "logLik value: " << MLElik << endl;
    cout << "Number of function evaluation: " << res.icount << endl;
    cout << "Message: " << res.message << endl;
}



// post-optim
void Inference::writeParam(string seed, string file)
{
    fstream fileO;
    fileO.open(file, ios::app);
    fileO << seed << " " ;
    fileO << MLEN0 << " ";
    fileO << MLEl0 << " ";
    fileO << MLEi1 << " ";
    fileO << MLEl1 << " ";
    fileO << MLEi2 << " ";
    fileO << MLEg2 << " ";
    fileO << MLEl2 << " ";
    fileO << MLEeps[0] << " ";
    fileO << MLEeps[1] << " ";
    fileO << MLEeps[2] << " ";
    fileO << MLElik << endl;
    fileO.close();
}
void Inference::assignCat(string file)
{
    // compute most likely category for each pattern
    vector<double> cat(nRep, 0.);
    tbb::detail::d1::parallel_for(0, nRep, [&](int i)
    {
        // compute proba of the 3 categories
        vector<double> lik(3, 0.);
        lik[0] = log(MLEN0/nObs);
        if (MLEl0 == 0. && (MLEeps[0] > 0. || freq[i] < nLeaves)) {
            lik[0] += freq[i]*log(1-MLEeps[0]) + (nLeaves-freq[i])*log(MLEeps[0]);
        } else if (MLEl0 > 0.) {
            lik[0] += tLik.evaluateRootRep(0, i, {0,1}, 0., MLEl0, 0., MLEeps[0], 0);
        }
        lik[1] = logSumExp(log(MLEi1/(MLEl1*nObs)) + tLik.evaluateRootRep(0, i, {0,1}, 0., MLEl1,
        0., MLEeps[1], 1), log(MLEi1/(MLEl1*nObs)) + this->likRep1(i, MLEl1, MLEeps[1]));
        lik[2] = log(MLEi2/nObs) + this->likRep2(i, MLEg2, MLEl2, MLEeps[2], MLEeps[2]);
        // find most probable
        cat[i] = 0;
        if (lik[1]>lik[0]) {cat[i] = 1;}
        if (lik[2]>lik[cat[i]]) {cat[i] = 2;}
    });

    // store inferred categories
    fstream fileO;
    fileO.open(file, ios::out);
    for_each(cat.begin(), cat.end(), [&](const auto& elem)
    {
        fileO << elem << " ";
    });
}
double Inference::nb0(double N0, double l0, double eps)
{
    double p0;
    if (l0 == 0.) {p0 = 1-pow(eps, nLeaves);}
    else {p0 = 1-exp(tLik.evaluateRootRep(0, nRep, {0,1}, 0., l0, 0., eps, 0));}
    return N0*p0;
}
double Inference::nb1(double i1, double l1, double eps)
{
    double p1 = 1-exp(tLik.evaluateRootRep(0, nRep, {0,1}, 0., l1, 0., eps, 1));
    double A = this->getPobs1(l1, eps);
    return i1*p1/l1 + i1*A;
}
double Inference::nb2(double i2, double g2, double l2, double eps)
{
    double B = this->getPobs2(g2, l2, eps, eps);
    return i2*B;
}
void Inference::printPangenomeCompo()
{
    cout << "Expected pangenome composition:" << endl;
    cout << this->nb0(MLEN0, MLEl0, MLEeps[0]) << " Persistent genes" << endl;
    cout << this->nb1(MLEi1, MLEl1, MLEeps[1]) << " Private genes" << endl;
    cout << this->nb2(MLEi2, MLEg2, MLEl2, MLEeps[2]) << " Mobile genes" << endl;
}
