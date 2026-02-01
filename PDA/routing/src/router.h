#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>
#include <unordered_set>
#include <chrono>
#include <climits>
#include <queue>
#include <set>

using namespace std;

struct Segment{
    int x1;
    int y1;
    int x2;
    int y2;
};

struct Net{
    string name;
    int x1;
    int y1;
    int x2;
    int y2;

    int HPWL=0;
    int length=0;
    int overflow=0;

    vector<Segment> path;

    bool inQueue = false;
};

struct Gcell{
    int occH =0;
    int occV =0;
    float hisOccH =0;
    float hisOccV =0;
    vector<Net*> netPassH;
    vector<Net*> netPassV;
};

struct Node{
    int x;
    int y;
    float f;
};

class GlobalRouter{

public:
    GlobalRouter(int,int,int,int,vector<Net>&);

    void Routing(chrono::time_point<std::chrono::steady_clock>);

    int getWirelength(){ return wireLength; };


private:
    int C;
    int R;
    int N;
    int capH;
    int capV;
    vector<Net>& nets;

    int wireLength;

    float coffH;
    float coffV;
    float coffhis;
    float coffbbox;

    int iteration;
    int totOverflow;
    
    vector<Gcell> gcells;

    vector<int> parentMap;
    vector<float> costMap;
    vector<int> nodeTag;
    int currNetId;      

    struct  cmp{
        bool operator()(Node a, Node b) {
            return a.f >= b.f;
        }
    };

    inline void UpdateOccupy(Net*, Segment&, int);
    int AstarSearch(int, int, int, int, vector<Segment>&);

    void checkLimit();
};