#include "FM.h"
#include <fstream>
#include <omp.h>

int main(int argc, char* argv[]){

    auto start = chrono::steady_clock::now();  // start time

    if (argc != 4) {
        std::cerr << "Usage: " << argv[0] << " <input> <output> <number of partitions>\n";
        return 1;
    }

    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    string inputFile = argv[1];
    string outputFile = argv[2];
    int partitions = stoi(argv[3]);
    //===================================================================
    
    int NumCells, NumNets;
    unordered_map<string, cell> cellMap;
    vector<cell*> cellList;
    vector<vector<cell*>> netList;

    //===================================================================
    
    // read file
    
    ifstream file(inputFile);
    string word, cellName, netName;
    int value;

    file >> word >> NumCells;               // NumCells 12752
    cellList.reserve(NumCells);
    cellMap.reserve(NumCells);
    for (int i=0; i<NumCells; i++) {
        file >> word >> cellName >> value;                  //Cell C1 50
        cellMap[cellName]= {cellName,value,{},{},0,{},{}} ; // name, size, nets, nidx, group, pos, gidx
        cellMap[cellName].pos.resize(partitions, nullptr);
        cellMap[cellName].gidx.resize(partitions, -1);
        cellMap[cellName].pos.resize(partitions, nullptr);

        cellList.push_back(&cellMap[cellName]);
    }


    file >> word >> NumNets;                //NumNets 14111
    netList.resize(NumNets);
    for (int i=0; i<NumNets; i++) {
        file >> word >> netName >> value;   // Net N1 2
        for(int j=0; j<value; j++){
            file >> word >> cellName;       // Cell C1
            cell& cel = cellMap[cellName];
            netList[i].push_back(&cel);
            cel.nets.push_back(i);
            cel.nidx.push_back(j);
        }
    }
    file.close();

    //===================================================================

    int maxThreads = omp_get_max_threads();
    int numTrials=maxThreads>32 ? 32: maxThreads>0 ? maxThreads : 16;
    int bestCutSize=NumNets;
    vector<vector<string>> bestGroups;
    
    // Parallel partitioning with different initial conditions
    #pragma omp parallel for
    for (int t = 0; t < numTrials; t++) {
        try {
            // copy cellMap
            unordered_map<string, cell> localCellMap = cellMap;

            // build localCellList
            vector<cell*> localCellList;
            localCellList.reserve(cellList.size());
            for (auto* c : cellList){
                localCellList.push_back(&localCellMap[c->name]);
                for (int j = 0; j < partitions; ++j)
                    localCellMap[c->name].pos[j] = new Node();
            }

            // build localNetList
            vector<vector<cell*>> localNetList;
            localNetList.resize(netList.size());
            for (size_t n = 0; n < netList.size(); n++)
                for (cell* c : netList[n])
                    localNetList[n].push_back(&localCellMap[c->name]);
            
            // sort/suffle cellList
            if(t==0){
                std::sort(localCellList.begin(), localCellList.end(), [](const cell* a, const cell* b) {
                    return a->nets.size() < b->nets.size();
                });
            }
            else{
                random_device rd;
                mt19937 g(rd());
                shuffle(localCellList.begin(), localCellList.end(), g);
            }

            // run FM
            FMEngine fm(localCellList, localNetList, partitions, start);
            fm.FiducciaMattheyses();

            // store result (critical)
            #pragma omp critical
            {
                if (fm.cutSize < bestCutSize) {
                    bestCutSize = fm.cutSize;
                    bestGroups = fm.groupsAfter;
                }
            }
        }catch (const std::exception& e) {
            #pragma omp critical
            std::cerr << "Thread " << t << " crashed: " << e.what() << std::endl;
        }
    }

    //output=============================================================
    ofstream outfile(outputFile);
    if (!outfile) {
        cerr << "Error opening output file: " << outputFile << endl;
        return 1;
    }
    outfile<<"CutSize "<<bestCutSize<<"\n";
    for (int g = 0; g < partitions; ++g) {
        sort(bestGroups[g].begin(), bestGroups[g].end(), [](const string& a, const string& b) {
            if (a.length() != b.length()) return a.length() < b.length();
            return a < b;
        });
        outfile<<endl;
        outfile << "Group" << static_cast<char>(65+g) <<" " << bestGroups[g].size() << "\n";
        for (const auto& cellName : bestGroups[g]) {
            outfile << cellName << "\n";
        }
    }
    outfile.close();



    return 0;
}