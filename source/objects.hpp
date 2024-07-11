#include <vector>
#include <string>
#include <map>
#include <cmath>

using namespace std;



// auxiliary function for log-sum-exp trick
double logSumExp(double a, double b);



// class for species tree
class RecTree
{
    public:

    static int firstPoint(string& s);
    RecTree(string text);
    ~RecTree();
    int getNodeID();
    string getLeaf();
    pair<double, double> getBranchLengths();
    RecTree* getLeftChild();
    RecTree* getRightChild();
    int addNodeID(int ind = 0);
    double treeLength();
    static RecTree readTree(string path);
    void getLeaves(vector<string>& v);
    void getLeavesID(vector<string>& v);


    private:
    
    int nodeID;
    string leaf;
    pair<double, double> branchLengths;
    RecTree* leftChild;
    RecTree* rightChild;
};



// class for presence/absence matrix
class PAmatrix
{
    public:

    static bool more(vector<int>& row1, vector<int>& row2);
    PAmatrix();
    PAmatrix(string path);
    vector<vector<int>>& getMatrix();
    map<string, int>& getColnames();
    vector<int>& getCounts();
    vector<string> extractPresentLeaves(int row);
    vector<int> extractCol(string colName);
    void sortByFreq();
    void compress();
    void add0();
    void restrainToLeaves(RecTree* t);

    private:

    vector<vector<int>> matrix; // columns are genomes and rows are genes
    // ie. matrix[i] is the repartition of gene i
    map<string, int> colNames;
    vector<int> counts;
};



// class for nodes of tree-like structure TreeLikFunctions
class LikFunction
{
    public:

    LikFunction(vector<int> leafValues);
    LikFunction(pair<double, double> branchLengths, pair<int, int> children);
    bool isLeaf();
    vector<int>& getLeafValues();
    int getLeafValue(int ind);
    pair<double, double>& getBranchLengths();
    pair<int, int> getChildren();

    private:

    vector<int> leafValues;
    pair<double, double> branchLengths;
    pair<int, int> children;
};



// class for storing likelihood calculations
class TreeLikFunctions
{
    public:

    static void add_func(RecTree* rt, PAmatrix& mat, vector<LikFunction>& lf);
    TreeLikFunctions();
    TreeLikFunctions(RecTree& tree, PAmatrix& mat);
    vector<LikFunction>& getLikFunctions();
    void evaluateNodeRep(int node, int rep, double gain, double loss, double error1, double error2, int cat);
    double evaluateRootRep(int node, int rep, pair<double, double> root, double gain, double loss,
    double error1, double error2, int cat);
    void clearLikValues();
    void clearLikValues0();
    void clearLikValues1();
    void clearLikValues2();

    private:

    vector<LikFunction> likFunctions;
    vector<vector<vector<double>>> likValues;
};

