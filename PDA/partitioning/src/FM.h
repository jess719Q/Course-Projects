#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <cmath>
#include <unordered_map>
#include <random>
#include <algorithm>
#include <chrono>

using namespace std;

struct Node;

struct cell{
    string name;
    int size;
    vector<int> nets;
    vector<int> nidx;
    int group;
    vector<Node*> pos;
    vector<int> gidx;
};

struct Node {
    cell* c;
    Node* front;
    Node* next;
};

struct record{
    cell* c;
    int cutsize;
    int fromG;
    int toG;
    int sizeDiff;
};


class FMEngine{
public:

    FMEngine(vector<cell*>& , vector<vector<cell*>>& , int, chrono::time_point<std::chrono::steady_clock>);
    void FiducciaMattheyses();

    vector<vector<string>> groupsAfter;
    vector<cell*> cellList;
    vector<vector<cell*>>& netList;
    int partitions=0, cutSize=0, totSize=0, maxP=0;
    chrono::time_point<std::chrono::steady_clock> startTime;

private:
    vector<int> groupSize;
    vector<vector<int>> netGroupNum;
    vector<vector<int>> bucketHead;
    vector<vector<vector< Node* >>> groupBucket;

    void MultiWayFM(vector<int>, int);
    void TwoWayInitFM(vector<int>, int);

    void removeFromBucket(cell*, int, int);
    void appendToBucket(cell*, int, int);
    void updateBucket(cell*, int, int, int);
    void InitializeGroupBucket();
    void Iniitalize2Group();
    void updateGain(cell*, int, int);
    int moveCell(vector<int>, vector<record>&, int, int);
    void move2anotherGroup(cell*, int, int);
    int moveCellforInitialize(int, int);
};