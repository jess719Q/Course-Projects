#include "router.h"

#include <fstream>
#include <sstream>

int main(int argc, char* argv[]){

    auto startTime = chrono::steady_clock::now();   // start time

    ios::sync_with_stdio(false); // Faster C++ I/O
    cin.tie(nullptr);

    if (argc != 3) { // Check argument count
        std::cerr << "Usage: " << argv[0] << " <input> <output>\n";
        return 1;
    }
    string inputFile = argv[1];
    string outputFile= argv[2];

    ifstream fin(inputFile); // Open LEF file
    if (!fin) {
        cerr << "Cannot open file.\n";
        return 1;
    }

    int C, R, capH, capV, nNets;
    string tmp;
    vector<Net> nets;

    fin >> tmp >> C >> R;
    fin >> tmp >> capH >> capV;
    fin >> tmp >> nNets;

    nets.reserve(nNets);

    for(int i=0; i<nNets; i++){
        string name;
        int x1, y1, x2, y2;
        fin >> tmp >> name >> tmp;
        fin >> tmp >> tmp >> x1 >> y1;
        fin >> tmp >> tmp >> x2 >> y2;
        int hpwl = abs(x1 - x2) + abs(y1 - y2);
        nets.push_back({name, x1, y1, x2, y2, hpwl, 0, 0, {}, 0});
    }

    fin.close();

    //=======================================================

    GlobalRouter globalRouter(C, R, capH, capV, nets);
    globalRouter.Routing(startTime);

    //=======================================================

    ofstream fout(outputFile);
    if (!fout) {
        cerr << "Error opening output file: " << outputFile << endl;
        return 1;
    }

    fout << "Wirelength " << globalRouter.getWirelength() <<endl;
    
    for(auto& n : nets){
        fout<<"Net "<<n.name<<endl;
        for(auto& s : n.path){
            fout<<"Segment "<<s.x1<<" "<<s.y1<<" "<<s.x2<<" "<<s.y2<<endl;
        }
    }
    fout.close();

    // show total time
    // auto now = chrono::steady_clock::now();
    // auto elapsed = chrono::duration_cast<chrono::seconds>(now - startTime).count();
    // cout<<"time: "<<elapsed<<endl;
}