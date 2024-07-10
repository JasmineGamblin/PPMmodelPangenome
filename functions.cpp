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

    mat = m; // I assume that this matrix is restrained to tree leaves, sorted and compressed
    //(but I don't know if it's necessary)

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
    
    // store subtrees (in fact it's just IDs of subtrees roots)
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

    MLEpar = vector<double>();
    MLEerr = vector<double>();
    MLElik = 0.;
}



// transition probabilities for cat. 1 genes
double Inference::pLost(double t, double l) {return log(1-exp(-l*t));}
double Inference::pRemain(double t, double l) {return -l*t;}



// likelihood functions for cat. 1
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
    /* parcourir l'arbre de la racine jusqu'à gainNode et pour chaque noeud
    strictement compris entre les 2 calculer likSubtree */
    
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
    double p_sub;
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



// auxiliary functions of cat. 2 calculations
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



// likelihood function for cat. 2
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



// complete lik
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
double Inference::computeI2(const array<double,7>& par)
{
    // error parametrization
    array<double,4> error = {par[6], par[6], par[6], par[6]};
    // array<double,4> error = {par[5], par[5]+par[6], par[7], par[5]+par[7]};
    
    // double p0 = 1-pow(error[0], nLeaves);
    double p0 = 1-exp(tLik.evaluateRootRep(0, nRep, {0,1}, 0., par[1], 0., error[0], 0));
    double p1 = 1-exp(tLik.evaluateRootRep(0, nRep, {0,1}, 0., par[3], 0., error[1], 1));
    double A = this->getPobs1(par[3], error[1]);
    double B = this->getPobs2(par[4], par[5], error[2], error[3]);
    return (nObs-par[0]*p0-par[2]*p1/par[3]-par[2]*A)/B;
}
double Inference::computeI2_0(const array<double,2>& par)
{
    double p0 = 1-exp(tLik.evaluateRootRep(0, nRep, {0,1}, 0., par[0], 0., par[1], 0));
    return nObs/p0;
}
double Inference::computeI2_1(const array<double,2>& par)
{
    double p1 = 1-exp(tLik.evaluateRootRep(0, nRep, {0,1}, 0., par[0], 0., par[1], 1));
    double A = this->getPobs1(par[0], par[1]);
    return nObs/(p1/par[0]+A);
}
double Inference::computeI2_2(const array<double,3>& par)
{
    double B = this->getPobs2(par[0], par[1], par[2], par[2]);
    return nObs/B;
}
double Inference::likRepTot(int rep, vector<double>& par, array<double,4>& error)
{
    // double lik0r = log(par[0]/nObs) + freq[rep]*log(1-error[0]) + (nLeaves-freq[rep])*log(error[0]);
    // if (error[0] == 0. && freq[rep] == nLeaves) {lik0r = log(par[0]/nObs);}
    double lik0r = log(par[0]/nObs)
    + tLik.evaluateRootRep(0, rep, {0,1}, 0., par[1], 0., error[0], 0);
    double lik1r = log(par[2]/(par[3]*nObs))
    + tLik.evaluateRootRep(0, rep, {0,1}, 0., par[3], 0., error[1], 1);
    double lik1 = log(par[2]/(par[3]*nObs)) + this->likRep1(rep, par[3], error[1]);
    double lik2 = log(par[6]/nObs) + this->likRep2(rep, par[4], par[5], error[2], error[3]);
    return logSumExp(logSumExp(lik0r, lik1r), logSumExp(lik1, lik2));
}
double Inference::likRepTot_0(int rep, vector<double>& par, double error)
{
    return log(par[1]/nObs) + tLik.evaluateRootRep(0, rep, {0,1}, 0., par[0], 0., error, 0);
}
double Inference::likRepTot_1(int rep, vector<double>& par, double error)
{
    double lik1r = log(par[1]/(par[0]*nObs))
    + tLik.evaluateRootRep(0, rep, {0,1}, 0., par[0], 0., error, 1);
    double lik1 = log(par[1]/(par[0]*nObs)) + this->likRep1(rep, par[0], error);
    return logSumExp(lik1r, lik1);
}
double Inference::likRepTot_2(int rep, vector<double>& par, array<double,2>& error)
{
    return log(par[2]/nObs) + this->likRep2(rep, par[0], par[1], error[0], error[1]);
}
double Inference::verif(vector<double> par)
{
    double sum = 0.;
    double aux;
    fstream file;
    file.open("freq_cpp_I2_ultr2.txt", ios::out);
    
    // error parametrization
    array<double,4> error = {par[7], par[7], par[7], par[7]};
    // array<double,4> error = {par[6], par[6]+par[7], par[8], par[6]+par[9]};
    
    for (int rep = 0; rep<nRep; rep++)
    {
        // aux = exp(freq[rep]*log(1-error[0]) + (nLeaves-freq[rep])*log(error[0]))/(1-pow(error[0], nLeaves));
        
        // aux = exp(tLik.evaluateRootRep(0, rep, {0,1}, 0., par[1], 0., error[0], 0))
        // /(1-exp(tLik.evaluateRootRep(0, nRep, {0,1}, 0., par[1], 0., error[0], 0)));
        // aux = exp(tLik.evaluateRootRep(0, rep, {0,1}, 0., par[3], 0., error[1], 1))
        // /(1-exp(tLik.evaluateRootRep(0, nRep, {0,1}, 0., par[3], 0., error[1], 1)));
        // aux = exp(this->likRep1(rep, par[3], error[1]))/(this->getPobs1(par[3], error[1])*par[3]);
        aux = exp(this->likRep2(rep, par[4], par[5], error[2], error[3]))
        /this->getPobs2(par[4], par[5], error[2], error[3]);
        sum += aux;
        file << aux << endl;
    }
    file.close();
    return sum;
}
double Inference::logLik(const array<double,7>& par)
{
    vector<double> parAll(par.begin(), par.begin()+6);
    
    // error parametrization
    array<double,4> error = {par[6], par[6], par[6], par[6]};
    // array<double,4> error = {par[5], par[5]+par[6], par[7], par[5]+par[7]};
    
    parAll.push_back(this->computeI2(par));
    if (parAll[6] < 0.)
    {
        return -1e9*1. + parAll[6];
    } else {
        vector<double> lik(nRep, 0.);
        tbb::detail::d1::parallel_for(0, nRep,
        [&](int i){lik[i] = this->likRepTot(i, parAll, error)*mat.getCounts()[i];});
        return accumulate(lik.begin(), lik.end(), 0.);
    }
}
double Inference::logLik_0(const array<double,2>& par)
{
    vector<double> parAll(par.begin(), par.begin()+1);
    
    parAll.push_back(this->computeI2_0(par));
    if (parAll[1] < 0.)
    {
        return -1e9*1. + parAll[1];
    } else {
        vector<double> lik(nRep, 0.);
        tbb::detail::d1::parallel_for(0, nRep,
        [&](int i){lik[i] = this->likRepTot_0(i, parAll, par[1])*mat.getCounts()[i];});
        return accumulate(lik.begin(), lik.end(), 0.);
    }
}
double Inference::logLik_1(const array<double,2>& par)
{
    vector<double> parAll(par.begin(), par.begin()+1);

    parAll.push_back(this->computeI2_1(par));
    if (parAll[1] < 0.)
    {
        return -1e9*1. + parAll[1];
    } else {
        vector<double> lik(nRep, 0.);
        tbb::detail::d1::parallel_for(0, nRep,
        [&](int i){lik[i] = this->likRepTot_1(i, parAll, par[1])*mat.getCounts()[i];});
        return accumulate(lik.begin(), lik.end(), 0.);
    }
}
double Inference::logLik_2(const array<double,3>& par)
{
    vector<double> parAll(par.begin(), par.begin()+2);
    
    // error parametrization
    array<double,2> error = {par[2], par[2]};
    // array<double,4> error = {par[5], par[5]+par[6], par[7], par[5]+par[7]};
    
    parAll.push_back(this->computeI2_2(par));
    if (parAll[2] < 0.)
    {
        return -1e9*1. + parAll[2];
    } else {
        vector<double> lik(nRep, 0.);
        tbb::detail::d1::parallel_for(0, nRep,
        [&](int i){lik[i] = this->likRepTot_2(i, parAll, error)*mat.getCounts()[i];});
        return accumulate(lik.begin(), lik.end(), 0.);
    }
}
void Inference::optim(int seed)
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
    while (this->computeI2(start) < 0. || start[1]>upper[1] || start[3]>upper[3] || start[4]>upper[4]
    || start[5]>upper[5] || start[6]>upper[6])
    {
        start = {uniN0(rng)*1., expl0(rng), unii1(rng)*1., expl1(rng), expl1(rng), expl2(rng), expe(rng)};
    }

    nelder_mead_result<double,nbPar> res = nelder_mead<double,nbPar>([&](const array<double,nbPar>& par, int nbEval){
        cout << nbEval << "th iteration: ";
        tLik.clearLikValues();
        cout << par[0] << ' ' << par[1] << ' ' << par[2] << ' ' << par[3] << ' ' << par[4] << ' ' << par[5] << ' '
        << par[6] << ' ';
        double ll = this->logLik(par);
        cout << setprecision(9) << ll << endl;
        return -ll;
        }, start, lower, upper, 0.0001, 3000, 10);
    MLEpar = {res.xmin[0], res.xmin[1], res.xmin[2], res.xmin[3], res.xmin[4], res.xmin[5]};
    MLEpar.push_back(computeI2(res.xmin));
    MLEerr = {res.xmin[6]};
    MLElik = -res.ymin;
    cout << "Found minimum: " << MLEpar[0] << ' ' << MLEpar[1] << ' ' << MLEpar[2]
    << ' ' << MLEpar[3] << ' ' << MLEpar[4] << ' ' << MLEpar[5] << ' ' << MLEpar[6] << ' ';
    cout << MLEerr[0] << endl;
    cout << "upper limit:   " << upper[0] << ' ' << upper[1] << ' ' << upper[2] << ' ' << upper[3] << ' '
    << upper[4] << ' ' << upper[5] << ' ' << upper[6] << endl;
    cout << "logLik value: " << MLElik << endl;
    cout << "Number of function evaluation: " << res.icount << endl;
    cout << "Message: " << res.message << endl;
}
void Inference::optim_0(int seed)
{
    const size_t nbPar = 2;
    
    mt19937 rng(seed);    // Random-number engine used (Mersenne-Twister in this case)
    exponential_distribution<> expl0(L);
    exponential_distribution<> expe(100.);
    array<double,nbPar> start = {expl0(rng), expe(rng)};
    array<double,nbPar> lower = {0., 0.};
    array<double,nbPar> upper = {10/L, 0.1};
    while (this->computeI2_0(start) < 0. || start[1]>upper[1] || start[0]>upper[0])
    {
        start = {expl0(rng), expe(rng)};
    }

    nelder_mead_result<double,nbPar> res = nelder_mead<double,nbPar>([&](const array<double,nbPar>& par, int nbEval){
        cout << nbEval << "th iteration: ";
        tLik.clearLikValues();
        cout << par[0] << ' ' << par[1] << ' ';
        double ll = this->logLik_0(par);
        cout << setprecision(9) << ll << endl;
        return -ll;
        }, start, lower, upper, 0.0001, 3000, 10);
    MLEpar = {res.xmin[0]};
    MLEpar.push_back(computeI2_0(res.xmin));
    MLEerr = {res.xmin[1]};
    MLElik = -res.ymin;
    cout << "Found minimum: " << MLEpar[0] << ' ' << MLEpar[1] << ' ';
    cout << MLEerr[0] << endl;
    cout << "upper limit:   " << upper[0] << ' ' << upper[1] << endl;
    cout << "logLik value: " << MLElik << endl;
    cout << "Number of function evaluation: " << res.icount << endl;
    cout << "Message: " << res.message << endl;
}
void Inference::optim_1(int seed)
{
    const size_t nbPar = 2;
    
    mt19937 rng(seed);    // Random-number engine used (Mersenne-Twister in this case)
    exponential_distribution<> expl1(L/50);
    exponential_distribution<> expe(100.);
    array<double,nbPar> start = {expl1(rng), expe(rng)};
    array<double,nbPar> lower = {0., 0.};
    array<double,nbPar> upper = {500/L, 0.1};
    while (this->computeI2_1(start) < 0. || start[1]>upper[1] || start[0]>upper[0])
    {
        start = {expl1(rng), expe(rng)};
    }

    nelder_mead_result<double,nbPar> res = nelder_mead<double,nbPar>([&](const array<double,nbPar>& par, int nbEval){
        cout << nbEval << "th iteration: ";
        tLik.clearLikValues();
        cout << par[0] << ' ' << par[1] << ' ';
        double ll = this->logLik_1(par);
        cout << setprecision(9) << ll << endl;
        return -ll;
        }, start, lower, upper, 0.0001, 3000, 10);
    MLEpar = {res.xmin[0]};
    MLEpar.push_back(computeI2_1(res.xmin));
    MLEerr = {res.xmin[1]};
    MLElik = -res.ymin;
    cout << "Found minimum: " << MLEpar[0] << ' ' << MLEpar[1] << ' ';
    cout << MLEerr[0] << endl;
    cout << "upper limit:   " << upper[0] << ' ' << upper[1] << endl;
    cout << "logLik value: " << MLElik << endl;
    cout << "Number of function evaluation: " << res.icount << endl;
    cout << "Message: " << res.message << endl;
}
void Inference::optim_2(int seed)
{
    const size_t nbPar = 3;
    
    mt19937 rng(seed);    // Random-number engine used (Mersenne-Twister in this case)
    exponential_distribution<> expl1(L/50);
    exponential_distribution<> expl2(L/2500);
    exponential_distribution<> expe(100.);
    array<double,nbPar> start = {expl1(rng), expl2(rng), expe(rng)};
    array<double,nbPar> lower = {0., 0., 0.};
    array<double,nbPar> upper = {500/L, 25000/L, 0.1};
    while (this->computeI2_2(start) < 0. || start[1]>upper[1] || start[0]>upper[0] || start[2]>upper[2])
    {
        start = {expl1(rng), expl2(rng), expe(rng)};
    }

    nelder_mead_result<double,nbPar> res = nelder_mead<double,nbPar>([&](const array<double,nbPar>& par, int nbEval){
        cout << nbEval << "th iteration: ";
        tLik.clearLikValues();
        cout << par[0] << ' ' << par[1] << ' ' << par[2] << ' ';
        double ll = this->logLik_2(par);
        cout << setprecision(9) << ll << endl;
        return -ll;
        }, start, lower, upper, 0.0001, 3000, 10);
    MLEpar = {res.xmin[0], res.xmin[1]};
    MLEpar.push_back(computeI2_2(res.xmin));
    MLEerr = {res.xmin[2]};
    MLElik = -res.ymin;
    cout << "Found minimum: " << MLEpar[0] << ' ' << MLEpar[1] << ' ' << MLEpar[2] << ' ' ;
    cout << MLEerr[0] << endl;
    cout << "upper limit:   " << upper[0] << ' ' << upper[1] << ' ' << upper[2] << endl;
    cout << "logLik value: " << MLElik << endl;
    cout << "Number of function evaluation: " << res.icount << endl;
    cout << "Message: " << res.message << endl;
}
void Inference::writeParam(string seed, string file)
{
    fstream fileO;
    fileO.open(file, ios::app);
    fileO << seed << " " ;
    for_each(MLEpar.begin(), MLEpar.end(), [&](const auto& elem)
    {
        fileO << elem << " ";
    });
    for_each(MLEerr.begin(), MLEerr.end(), [&](const auto& elem)
    {
        fileO << elem << " ";
    });
    fileO << MLElik << endl;
    fileO.close();
}
// void Inference::writeParam(string seed, string file, double N0, double l0, double i1, double l1, double g2, double l2,
// double e)
// {
//     fstream fileO;
//     fileO.open(file, ios::app);
//     fileO << seed << " " ;
//     for_each(MLEpar.begin(), MLEpar.end(), [&](const auto& elem)
//     {
//         fileO << elem << " ";
//     });
//     for_each(MLEerr.begin(), MLEerr.end(), [&](const auto& elem)
//     {
//         fileO << elem << " ";
//     });
//     fileO << MLElik << " " << this->logLik({N0, l0, i1, l1, g2, l2, e}) << endl;
//     fileO.close();
// }
// void Inference::writeParam_2(string seed, string file, double g2, double l2, double e)
// {
//     fstream fileO;
//     fileO.open(file, ios::app);
//     fileO << seed << " " ;
//     for_each(MLEpar.begin(), MLEpar.end(), [&](const auto& elem)
//     {
//         fileO << elem << " ";
//     });
//     for_each(MLEerr.begin(), MLEerr.end(), [&](const auto& elem)
//     {
//         fileO << elem << " ";
//     });
//     fileO << MLElik << " " << this->logLik_2({g2, l2, e}) << endl;
//     fileO.close();
// }
void Inference::assignCat(string file)
{
    // error parametrization
    array<double,4> error = {MLEerr[0], MLEerr[0], MLEerr[0], MLEerr[0]};
    // array<double,4> error = {par[6], par[6]+par[7], par[8], par[6]+par[8]};
    
    // compute most probable category for each pattern
    vector<double> cat(nRep, 0.);
    tbb::detail::d1::parallel_for(0, nRep, [&](int i)
    {
        // compute proba of the 3 categories
        vector<double> lik(3, 0.);
        // lik[0] = log(MLEpar[0]/nObs) + freq[i]*log(1-error[0]) + (nLeaves-freq[i])*log(error[0]);
        lik[0] = log(MLEpar[0]/nObs) + tLik.evaluateRootRep(0, i, {0,1}, 0., MLEpar[1],
        0., error[0], 0);
        lik[1] = logSumExp(log(MLEpar[2]/(MLEpar[3]*nObs)) + tLik.evaluateRootRep(0, i, {0,1}, 0., MLEpar[3],
        0., error[1], 1), log(MLEpar[2]/(MLEpar[3]*nObs)) + this->likRep1(i, MLEpar[3], error[1]));
        lik[2] = log(MLEpar[6]/nObs) + this->likRep2(i, MLEpar[4], MLEpar[5], error[2], error[3]);
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
void Inference::computeProbaCat(string file, double N0, double l0, double i1, double l1, double i2, double g2, double l2,
double eps0, double eps1, double eps2)
{
    vector<double> row(3, 0.);
    vector<vector<double>> lik(nRep, row);
    tbb::detail::d1::parallel_for(0, nRep, [&](int i)
    {
        // compute proba of the 3 categories
        lik[i][0] = log(N0/nObs) + tLik.evaluateRootRep(0, i, {0,1}, 0., l0, 0., eps0, 0);
        lik[i][1] = logSumExp(log(i1/(l1*nObs)) + tLik.evaluateRootRep(0, i, {0,1}, 0., l1, 0., eps1, 1),
        log(i1/(l1*nObs)) + this->likRep1(i, l1, eps1));
        lik[i][2] = log(i2/nObs) + this->likRep2(i, g2, l2, eps2, eps2);
        double denominateur = logSumExp(lik[i][0], logSumExp(lik[i][1], lik[i][2]));
        lik[i][0] -= denominateur;
        lik[i][1] -= denominateur;
        lik[i][2] -= denominateur;
    });

    // store proba
    fstream fileO;
    fileO.open(file, ios::out);
    for_each(lik.begin(), lik.end(), [&](const auto& elem)
    {
        fileO << elem[0] << " " << elem[1] << " " << elem[2] << endl;
    });
}
void Inference::computeProbaPattern(string file, double N0, double l0, double i1, double l1, double i2, double g2, double l2,
double eps0, double eps1, double eps2)
{
    vector<double> row(3, 0.);
    vector<vector<double>> lik(nRep, row);
    tbb::detail::d1::parallel_for(0, nRep, [&](int i)
    {
        // compute proba of the 3 categories
        lik[i][0] = exp(log(N0/nObs) + tLik.evaluateRootRep(0, i, {0,1}, 0., l0, 0., eps0, 0));
        lik[i][1] = exp(logSumExp(log(i1/(l1*nObs)) + tLik.evaluateRootRep(0, i, {0,1}, 0., l1, 0., eps1, 1),
        log(i1/(l1*nObs)) + this->likRep1(i, l1, eps1)));
        lik[i][2] = exp(log(i2/nObs) + this->likRep2(i, g2, l2, eps2, eps2));
    });

    // store proba
    fstream fileO;
    fileO.open(file, ios::out);
    for_each(lik.begin(), lik.end(), [&](const auto& elem)
    {
        fileO << elem[0] + elem[1] + elem[2] << endl;
    });
}
double Inference::nb0(vector<double> par)
{
    double p0 = 1-exp(tLik.evaluateRootRep(0, nRep, {0,1}, 0., par[1], 0., par[2], 0));
    return par[0]*p0;
}
double Inference::nb1(vector<double> par) // warning: params are i1, l1, eps1
{
    double p1 = 1-exp(tLik.evaluateRootRep(0, nRep, {0,1}, 0., par[1], 0., par[2], 1));
    double A = this->getPobs1(par[1], par[2]);
    return par[0]*p1/par[1] + par[0]*A;
}
double Inference::nb2(vector<double> par)
{
    
    // double p2 = 1-exp(tLik.evaluateRootRep(0, nRep, {0,1}, par[1], par[2], par[3], par[3], 2));
    double B = this->getPobs2(par[1], par[2], par[3], par[3]);
    return B*par[0];
}
double Inference::sumGains(int st, double g2)
{
    LikFunction lf = tLik.getLikFunctions()[st];
    if (lf.isLeaf())
    {
        return -1;
    } else {
        return exp(g2*height[st]) + sumGains(lf.getChildren().first, g2) + sumGains(lf.getChildren().second, g2);
    }
}
double Inference::expectedGains(double g2)
{
    // compute expected number of gains
    double sum,sum2,w;
    int st;
    sum = 0;
    for (int rank = nLeaves-1; rank < nNodes-1; rank++) // begin at the last leaf because tree is ultrametric
    {
        w = height[order[rank+1]] - height[order[rank]];
        if (w > 0.)
        {
            sum2 = 0.;
            for (int j = 0; j < subtrees[rank].size(); j++)
            {
                st = subtrees[rank][j];
                sum2 += w + this->sumGains(st, g2)*(exp(-g2*height[order[rank]]) - exp(-g2*height[order[rank+1]]))/g2;
            }
            sum += sum2;
        }
    }
    return sum/H;
}
double Inference::expectedTime1_old(int rep, double l1, double eps1)
{
    int gainNodeID = mrca[rep];
    if (gainNodeID != 0) // if mrca is not the root
    {
        /* parcourir l'arbre de la racine jusqu'à gainNode et pour chaque noeud
        strictement compris entre les 2 calculer likSubtree */
        stack<int> nodes;
        int node = 0;
        nodes.push(node); // on ajoute la racine (contrairement à likUpper1)
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
        double p_sub,ta;
        double likSubtrees(0.);
        double tmrca = height[gainNodeID];
        double ta1 = height[nodes.top()];
        double likU = log(1/l1 - exp(-l1*(ta1-tmrca))*(ta1-tmrca+1/l1));
        while (nodes.top() != 0)
        {
            children = tLik.getLikFunctions()[nodes.top()].getChildren();
            ta = ta1;
            nodes.pop();
            ta1 = height[nodes.top()];
            if (gainNodeID >= children.second) {child = children.first;}
            else {child = children.second;}
            likSubtrees += logSumExp(pLost(branchLength[child], l1), this->likSubtree(child, l1, eps1));
            likU = logSumExp(likU, likSubtrees + log(exp(-l1*(ta-tmrca))*(ta-tmrca+1/l1)
            - exp(-l1*(ta1-tmrca))*(ta1-tmrca+1/l1)));
        }
        ;
        return (exp(tLik.evaluateRootRep(mrca[rep], rep, {0.,1.}, 0., l1, 0., eps1, 1))*exp(likU)
        + exp(tLik.evaluateRootRep(0, rep, {0,1}, 0, l1, 0, eps1, 1))*(H-tmrca+1/l1))/(l1*L+1);
    } else { // if mrca is the root
        return exp(tLik.evaluateRootRep(0, rep, {0,1}, 0, l1, 0, eps1, 1))*(1/l1)/(l1*L+1);
    }
}
double Inference::expectedTime1(int rep, double l1, double eps1)
{
    int gainNodeID = mrca[rep];
    if (gainNodeID != 0) // if mrca is not the root
    {
        /* parcourir l'arbre de la racine jusqu'à gainNode et pour chaque noeud
        strictement compris entre les 2 calculer likSubtree */
        stack<int> nodes;
        int node = 0;
        nodes.push(node); // on ajoute la racine (contrairement à likUpper1)
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
        double p_sub,ta,q_root;
        double likSubtrees(0.);
        double tmrca = height[gainNodeID];
        double ta1 = height[nodes.top()];
        double sumt = log(1/l1 - exp(-l1*(ta1-tmrca))*(ta1-tmrca+1/l1));
        double sum = log(1 - exp(-l1*(ta1-tmrca)));
        while (nodes.top() != 0)
        {
            children = tLik.getLikFunctions()[nodes.top()].getChildren();
            ta = ta1;
            nodes.pop();
            ta1 = height[nodes.top()];
            if (gainNodeID >= children.second) {child = children.first;}
            else {child = children.second;}
            likSubtrees += logSumExp(pLost(branchLength[child], l1), this->likSubtree(child, l1, eps1));
            sumt = logSumExp(sumt, likSubtrees + log(exp(-l1*(ta-tmrca))*(ta-tmrca+1/l1)
            - exp(-l1*(ta1-tmrca))*(ta1-tmrca+1/l1)));
            sum = logSumExp(sum, likSubtrees + log(exp(-l1*(ta-tmrca)) - exp(-l1*(ta1-tmrca))));
        }
        children = tLik.getLikFunctions()[nodes.top()].getChildren();
        if (gainNodeID >= children.second) {child = children.first;}
        else {child = children.second;}
        q_root = exp(likSubtrees + logSumExp(pLost(branchLength[child], l1), this->likSubtree(child, l1, eps1)));
        return (exp(sumt)+q_root*exp(-l1*(H-tmrca))*(H-tmrca+1/l1))/(exp(sum) + q_root*exp(-l1*(H-tmrca)));
    } else { // if mrca is the root
        return 1/l1;
    }
}
double Inference::expectedTime2(int rep, double g2, double l2, double eps2)
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
                t = height[order[rank]]+(2*i+1)*w/(2*nSub);
                f_aux = log(t); // only difference with likRep2 (f_aux = 0)
                for (int j = 0; j < subtrees[rank].size(); j++)
                {
                    st = subtrees[rank][j];
                    if (tLik.getLikFunctions()[st].isLeaf())
                    {
                        f_aux += log(p0leaf(mat.getMatrix()[rep][mat.getColnames()[leavesID[st]]],
                        t-height[st], g2, l2, eps2, eps2));
                    } else {
                        p1 = p01(t - height[st], g2, l2);
                        f_aux += tLik.evaluateRootRep(st, rep, {1-p1,p1}, g2, l2, eps2, eps2, 2);
                    }
                }
                integral += exp(f_aux);
            }
            sum += integral*w/nSub;
        }
    }
    return sum;
}
void Inference::expectedTime(string file, double l1, double g2, double l2, double eps1, double eps2)
{
    vector<double> row(2, 0.);
    vector<vector<double>> times(nRep, row);
    tbb::detail::d1::parallel_for(0, nRep, [&](int i)
    {
        // compute expected time of arrival for type 1 and 2 
        // times[i][0] = 0;
        times[i][0] = this->expectedTime1(i, l1, eps1);
        times[i][1] = 0;
        // times[i][1] = this->expectedTime2(i, g2, l2, eps2)/exp(this->likRep2(i, g2, l2, eps2, eps2));
    });

    // store proba
    fstream fileO;
    fileO.open(file, ios::out);
    for_each(times.begin(), times.end(), [&](const auto& elem)
    {
        fileO << elem[0] << endl;
    });
}
vector<double> Inference::expectedTime2_times()
{
    vector<double> times = {};
    double w;
    int nSub;
    for (int rank = nLeaves-1; rank < nNodes-1; rank++) // begin at the last leaf because tree is ultrametric
    {
        w = height[order[rank+1]] - height[order[rank]];
        if (w > 0.)
        {
            nSub = max(1,(int)ceil(w*100/H));
            for (int i = 0; i < nSub; i++) {times.push_back(height[order[rank]]+(2*i+1)*w/(2*nSub));}
        }
    }
    return times;
}
vector<double> Inference::expectedTime2_aux(int rep, double g2, double l2, double eps2)
{
    vector<double> f_val = {};
    double w,f_aux,t,p1;
    int nSub,st;
    for (int rank = nLeaves-1; rank < nNodes-1; rank++) // begin at the last leaf because tree is ultrametric
    {
        w = height[order[rank+1]] - height[order[rank]];
        if (w > 0.)
        {
            nSub = max(1,(int)ceil(w*100/H));
            for (int i = 0; i < nSub; i++)
            {
                t = height[order[rank]]+(2*i+1)*w/(2*nSub);
                f_aux = 0.;
                for (int j = 0; j < subtrees[rank].size(); j++)
                {
                    st = subtrees[rank][j];
                    if (tLik.getLikFunctions()[st].isLeaf())
                    {
                        f_aux += log(p0leaf(mat.getMatrix()[rep][mat.getColnames()[leavesID[st]]],
                        t-height[st], g2, l2, eps2, eps2));
                    } else {
                        p1 = p01(t - height[st], g2, l2);
                        f_aux += tLik.evaluateRootRep(st, rep, {1-p1,p1}, g2, l2, eps2, eps2, 2);
                    }
                }
                f_val.push_back(exp(f_aux));
            }
        }
    }
    return f_val;
}
void Inference::expectedTime2_store(string file, double g2, double l2, double eps2)
{
    vector<vector<double>> dist(nRep+1); // first line is time
    dist[0] = this->expectedTime2_times();
    dist[0].push_back(0);
    // tbb::detail::d1::parallel_for(0, nRep, [&](int i)
    // {
    //     // compute distribution of arrival time for type 2 
    //     dist[i+1] = this->expectedTime2_aux(i, g2, l2, eps2);
    //     dist[i+1].push_back(exp(this->likRep2(i, g2, l2, eps2, eps2))); // last column is denominator
    // });

    // store dist
    fstream fileO;
    fileO.open(file, ios::out);
    for (int i = 0; i < 1; i++)
    {
        for (int j = 0; j < dist[0].size(); j++)
        {
            fileO << dist[i][j] << " ";
        }
        fileO << endl;
    }
}