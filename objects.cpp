#include "objects.hpp"
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <numeric>

using namespace std;



double logSumExp(double a, double b)
{
    if (isinf(a))
    {
        if (isinf(b)) {return -INFINITY;}
        else {return b;}
    }
    else if (isinf(b)) {return a;}
    else if (b-a <= 100. && b-a >= -100.) {return a + log(1 + exp(b-a));}
    // else if (b-a > 700. && b-a <= 1400.) {return 700 + a + log(exp(-700) + exp(b-a-700));}
    // else if (b-a < -700. && b-a >= -1400.) {return -700 + a + log(exp(700) + exp(b-a+700));}
    // if b-a is too small or too big we neglect the smallest term
    // 100 should be enough, on one example we began to see a difference in the likelihood (precision 10-5) under 20
    else if (b-a < -100.) {return a;}
    else {return b;}
}



// RecTree class
int RecTree::firstPoint(string& s)
{
    int par = 0;
    int pos = 1;
    do
    {
        if (s[pos] == '(') {par ++;}
        else if (s[pos] == ')') {par --;}
        pos ++;
    } while (par > 0);
    while (s[pos] != ':') {pos ++;}
    return pos;
}
RecTree::RecTree(string text)
{
    if (text[0] != '(') // create leaf
    {
        leaf = text;
        branchLengths = {0.,0.};
        leftChild = NULL;
        rightChild = NULL;
    }
    else // recursive call to constructor
    {
        int first_point = firstPoint(text);
        int coma = first_point+1;
        while (text[coma] != ','){coma++;}
        int last_par = text.size()-1;
        while (text[last_par] != ')'){last_par--;}
        int second_point = last_par-1;
        while (text[second_point] != ':'){second_point--;}
        
        leaf = "";
        branchLengths = {stold(text.substr(first_point+1,coma-first_point-1)),
        stold(text.substr(second_point+1,last_par-second_point-1))};
        RecTree* lc = new RecTree(text.substr(1,first_point-1));
        RecTree* rc = new RecTree(text.substr(coma+1,second_point-coma-1));
        leftChild = lc;
        rightChild = rc;
    }
}
RecTree::~RecTree()
{
    delete leftChild;
    delete rightChild;
}
int RecTree::getNodeID()
{
    return nodeID;
}
string RecTree::getLeaf()
{
    return leaf;
}
pair<double, double> RecTree::getBranchLengths()
{
    return branchLengths;
}
RecTree* RecTree::getLeftChild()
{
    return leftChild;
}
RecTree* RecTree::getRightChild()
{
    return rightChild;
}
int RecTree::addNodeID(int ind)
{
    nodeID = ind;
    ind ++;
    if (leaf.size() == 0)
    {
        ind = leftChild->addNodeID(ind);
        ind = rightChild->addNodeID(ind);
    }
    return ind;
}
double RecTree::treeLength()
{
    if (leaf.size() != 0) {return 0.;}
    else
    {
        return branchLengths.first + leftChild->treeLength() + branchLengths.second + rightChild->treeLength();
    }
}
RecTree RecTree::readTree(string path)
{
    fstream file;
    file.open(path, ios::in);
    string line;
    getline(file, line);
    file.close();
    return RecTree(line);
}
void RecTree::getLeaves(vector<string>& v)
{
    if (leaf.size() > 0) {v.push_back(leaf);}
    else
    {
        leftChild->getLeaves(v);
        rightChild->getLeaves(v);
    }
}
void RecTree::getLeavesID(vector<string>& v)
{
    if (leaf.size() > 0) {v[nodeID] = leaf;}
    else
    {
        leftChild->getLeavesID(v);
        rightChild->getLeavesID(v);
    }
}



// PAmatrix class
bool PAmatrix::more(vector<int>& row1, vector<int>& row2)
{
    return accumulate(row1.begin(), row1.end(), 0) > accumulate(row2.begin(), row2.end(), 0);
}
PAmatrix::PAmatrix()
{
    matrix = vector<vector<int>>();
    colNames = map<string, int>();
    counts = vector<int>();
}
PAmatrix::PAmatrix(string path)
{
    matrix = vector<vector<int>>();
    colNames = map<string, int>();
    counts = vector<int>();
    
    int col = 0;
    fstream file;
    file.open(path, ios::in);
    string line, tok;
    
    // first line contains counts if reading a compressed matrix
    // a compressed matrix must be already sorted
    getline(file, line);
    if (line[0] == 'c' && line[1] == 'c' && line[2] == 'c')
    {
        stringstream ss(line);
        int row = 0;
        getline(ss, tok, ','); // first token is compressed signal "ccc"
        while(getline(ss, tok, ','))
        {
            counts.push_back(stoi(tok));
            row++;
        }
    }
    
    while (getline(file, line))
    {
        stringstream ss(line);
        
        // first token is colname
        getline(ss, tok, ','); 
        colNames[tok] = col;
        
        int row = 0;
        while(getline(ss, tok, ','))
        {
            if (col == 0) {matrix.push_back(vector<int>());}
            matrix[row].push_back(stoi(tok));
            row++;
        }
        col++;

    }
    file.close();

    // if matrix was not compressed
    if (counts.empty())
    {
        counts = vector<int>(matrix.size(),1);
    }
}
vector<vector<int>>& PAmatrix::getMatrix()
{
    return matrix;
}
map<string, int>& PAmatrix::getColnames()
{
    return colNames;
}
vector<int>& PAmatrix::getCounts()
{
    return counts;
}
vector<string> PAmatrix::extractPresentLeaves(int row)
{
    vector<string> l;
    map<string, int>::iterator it;
    for(it = colNames.begin(); it!=colNames.end(); ++it)
    {
        if (matrix[row][it->second] == 1)
        {
            l.push_back(it->first);
        }
    }
    return l;
}
vector<int> PAmatrix::extractCol(string colName)
{
    vector<int> col;
    int colID = colNames[colName];
    for (int i = 0; i<matrix.size(); i++)
    {
        col.push_back(matrix[i][colID]);
    }
    return col;
}
void PAmatrix::sortByFreq()
{
    sort(matrix.begin(), matrix.end(), more);
}
void PAmatrix::compress()
{
    int ind,i;
    ind = 0;
    while (ind < matrix.size())
    {
        i = ind+1;
        while (i<matrix.size())
        {
            if (matrix[ind] == matrix[i])
            {
                matrix.erase(matrix.begin()+i);
                counts[ind]++;
                counts.erase(counts.begin()+i);
            }
            else {i++;}
        }
        ind++;
    }   
}
void PAmatrix::add0()
{
    matrix.push_back(vector<int>());
    for (int i=0; i<matrix[0].size(); i++)
    {
        matrix[matrix.size()-1].push_back(0);
    }
}
void PAmatrix::restrainToLeaves(RecTree* t)
{   
    vector<string> leaves;
    t->getLeaves(leaves);

    map<string, int>::iterator it = colNames.begin();
    int nbRemoved(0);
    while (it != colNames.end())
    {
        colNames[it->first] -= nbRemoved;
        if (find(leaves.begin(), leaves.end(), it->first) == leaves.end())
        {
            // remove column
            for (int i = 0; i<matrix.size(); i++)
            {
                matrix[i].erase(matrix[i].begin()+it->second);
            }
            
            // remove colname
            colNames.erase(it++);
            nbRemoved++;
        }
        else
        {
            it++;
        }
        
    }

    // remove rows of zeros
    int i(0);
    while (i<matrix.size())
    {
        int sum = accumulate(matrix[i].begin(), matrix[i].end(), 0);
        if (sum == 0)
        {
            matrix.erase(matrix.begin()+i);
            counts.erase(counts.begin()+i);
        }
        else {i++;}
    }
}



// LikFunctionClass
LikFunction::LikFunction(vector<int> lv)
{
    leafValues = lv;
    branchLengths = pair<double, double>();
    children = pair<int, int>();
}
LikFunction::LikFunction(pair<double, double> bl, pair<int, int> c)
{
    leafValues = vector<int>();
    branchLengths = bl;
    children = c;
}
bool LikFunction::isLeaf()
{
    return !leafValues.empty();
}
vector<int>& LikFunction::getLeafValues()
{
    return leafValues;
}
int LikFunction::getLeafValue(int i)
{
    return leafValues[i];
}
pair<double, double>& LikFunction::getBranchLengths()
{
    return branchLengths;
}
pair<int, int> LikFunction::getChildren()
{
    return children;
}



// TreeLikFunctions
void TreeLikFunctions::add_func(RecTree* rt, PAmatrix& mat, vector<LikFunction>& lf)
{
    if (rt->getLeaf().size()>0) // leaf
    {
        LikFunction l(mat.extractCol(rt->getLeaf()));
        lf.push_back(l); 

    }
    else // recursive tree traversal
    {
        pair<double, double> branchLengths = rt->getBranchLengths();
        pair<int, int> children = {rt->getLeftChild()->getNodeID(),rt->getRightChild()->getNodeID()}; 
        LikFunction l(branchLengths, children);
        lf.push_back(l);
        add_func(rt->getLeftChild(), mat, lf);
        add_func(rt->getRightChild(), mat, lf);
    }
}
TreeLikFunctions::TreeLikFunctions()
{
    likFunctions = vector<LikFunction>();
    likValues = vector<vector<vector<double>>>();
}
TreeLikFunctions::TreeLikFunctions(RecTree& tree, PAmatrix& mat)
{
    vector<LikFunction> lf;
    add_func(&tree, mat, lf);
    likFunctions = lf;
    vector<double> row1 = {0.,0.,0.,0.,0.,0.}; // indices 0 and 1 for cat 0, 2 and 3 for cat 1, and 4 and 5 for cat 2
    vector<vector<double>> row2(mat.getMatrix().size(), row1);
    likValues = vector<vector<vector<double>>>(likFunctions.size(), row2);;
}
vector<LikFunction>& TreeLikFunctions::getLikFunctions()
{
    return likFunctions;
}
void TreeLikFunctions::evaluateNodeRep(int node, int rep, double g, double l, double e1, double e2, int cat)
{
    if (likValues[node][rep][2*cat] + likValues[node][rep][2*cat+1] == 0.)
    {
        if (likFunctions[node].isLeaf())
        {
            // lik values at the leaves
            likValues[node][rep][2*cat] = (likFunctions[node].getLeafValue(rep) == 0)? log(1-e1) : log(e1);
            likValues[node][rep][2*cat+1] = (likFunctions[node].getLeafValue(rep) == 1)? log(1-e2) : log(e2);
        }
        else
        {
            pair<int, int> children = likFunctions[node].getChildren();
            this->evaluateNodeRep(likFunctions[node].getChildren().first, rep, g, l, e1, e2, cat);
            this->evaluateNodeRep(likFunctions[node].getChildren().second, rep, g, l, e1, e2, cat);
            vector<double>& p_left = likValues[children.first][rep];
            vector<double>& p_right = likValues[children.second][rep];
            pair<double, double>& bl = likFunctions[node].getBranchLengths();
            double p0;
            double p1;
            if (g == 0.) // cat 0, 1 or 2
            {
                p0 = p_left[2*cat] + p_right[2*cat];
                p1 = logSumExp(log(1-exp(-l*bl.first))+p_left[2*cat], -l*bl.first+p_left[2*cat+1])
                + logSumExp(log(1-exp(-l*bl.second))+p_right[2*cat], -l*bl.second+p_right[2*cat+1]);
                likValues[node][rep][2*cat] = p0;
                likValues[node][rep][2*cat+1] = p1;
            } else { // cat 2
                double gl = g + l;
                p0 = logSumExp(log((l+g*exp(-gl*bl.first))/gl)+p_left[4], log((1-exp(-gl*bl.first))*g/gl)+p_left[5])
                + logSumExp(log((l+g*exp(-gl*bl.second))/gl)+p_right[4], log((1-exp(-gl*bl.second))*g/gl)+p_right[5]);
                p1 = logSumExp(log((1-exp(-gl*bl.first))*l/gl)+p_left[4], log((g+l*exp(-gl*bl.first))/gl)+p_left[5])
                + logSumExp(log((1-exp(-gl*bl.second))*l/gl)+p_right[4], log((g+l*exp(-gl*bl.second))/gl)+p_right[5]);
                likValues[node][rep][4] = p0;
                likValues[node][rep][5] = p1;
            } 
        }
    }
}
double TreeLikFunctions::evaluateRootRep(int node, int rep, pair<double, double> root, double g, double l,
double e1, double e2, int cat)
{
    this->evaluateNodeRep(node, rep, g, l, e1, e2, cat);
    return logSumExp(log(root.first)+likValues[node][rep][2*cat],
    log(root.second)+likValues[node][rep][2*cat+1]);
}
void TreeLikFunctions::clearLikValues()
{
    for (int i = 0; i < likValues.size(); i++)
    {
        for (int j = 0; j < likValues[i].size(); j++)
        {
            likValues[i][j][0] = 0.;
            likValues[i][j][1] = 0.;
            likValues[i][j][2] = 0.;
            likValues[i][j][3] = 0.;
            likValues[i][j][4] = 0.;
            likValues[i][j][5] = 0.;
        }
    } 
}
void TreeLikFunctions::clearLikValues0()
{
    for (int i = 0; i < likValues.size(); i++)
    {
        for (int j = 0; j < likValues[i].size(); j++)
        {
            likValues[i][j][0] = 0.;
            likValues[i][j][1] = 0.;
        }
    } 
}
void TreeLikFunctions::clearLikValues1()
{
    for (int i = 0; i < likValues.size(); i++)
    {
        for (int j = 0; j < likValues[i].size(); j++)
        {
            likValues[i][j][2] = 0.;
            likValues[i][j][3] = 0.;
        }
    } 
}
void TreeLikFunctions::clearLikValues2()
{
    for (int i = 0; i < likValues.size(); i++)
    {
        for (int j = 0; j < likValues[i].size(); j++)
        {
            likValues[i][j][4] = 0.;
            likValues[i][j][5] = 0.;
        }
    } 
}
 

