#include "router.h"
#include <iomanip>

GlobalRouter::GlobalRouter(int c, int r, int ch, int cv, vector<Net>& n)
    : C(c), R(r), capH(ch), capV(cv), nets(n) {
    
    N = C*R;
    gcells.resize(C*R, {0,0,0,0});

    parentMap.resize(C*R);
    costMap.resize(C*R);
    nodeTag.resize(C*R, 0);
    currNetId = 0;

    wireLength = 0;
}

void GlobalRouter::Routing(chrono::time_point<std::chrono::steady_clock> startTime){

    // Sort nets by HPWL (High to Low).
    sort(nets.begin(), nets.end(), [](Net a, Net b) {return a.HPWL > b.HPWL;});

    // Initialize penalty coefficients
    coffH=1;
    coffV=1;
    coffbbox=1;     // Bounding box size control
    coffhis=0.2;    // History cost increment base

    iteration=1;

    // initial routing =========================================================
    for(auto& net: nets){
        int length = AstarSearch(net.x1, net.y1, net.x2, net.y2, net.path);
        net.length = length;
        wireLength += length;

        // Mark the path usage on the grid
        for(auto& s : net.path) UpdateOccupy(&net, s,1);
    }

    // First Refinement, rip-up and reroute all nets once ========================
    for(auto& net: nets){
        
        for(auto& s : net.path) UpdateOccupy(&net, s, -1);

        net.path.clear();
        wireLength -= net.length;
        int length = AstarSearch(net.x1, net.y1, net.x2, net.y2, net.path);
        net.length = length;
        wireLength += length;

        for(auto& s : net.path) UpdateOccupy(&net, s, 1);
    }
    
    // Main loop to rip-up and reroute nets  ================================

    int lastNumOF=0;
    int stucked=0;
    int inc =1;

    for(int k=0; k<1000; k++){

        iteration++;
        vector<Net*> rerouteQueue;

        bool of = false; // overflow exist?
        int numOF = 0;   // total overflow amount

        for(auto& n: nets) {
            n.overflow=0;
            n.inQueue=false;
        }

        // If progress is stuck, increase penalty faster
        inc = (stucked>1)? 3:1;

        // Check every grid cell for overflow and exactly full
        for(auto& g: gcells){
            if(g.occH>capH){
                for(auto* n: g.netPassH) {
                    if(!n->inQueue){
                        rerouteQueue.push_back(n);
                        n->inQueue=true;
                    }
                    n->overflow++;
                }
                g.hisOccH +=inc;    // Increase history cost
                of = true;
                numOF++;
            }
            else if(g.occH==capH){
                for(auto* n: g.netPassH) {
                    if(!n->inQueue){
                        rerouteQueue.push_back(n);
                        n->inQueue=true;
                    }
                }
            }

            if(g.occV>capV){
                for(auto* n: g.netPassV) {
                    if(!n->inQueue){
                        rerouteQueue.push_back(n);
                        n->inQueue=true;
                    }
                    n->overflow++;
                }
                g.hisOccV +=inc;
                of = true;
                numOF++;
            }
            else if(g.occV==capV){
                for(auto* n: g.netPassV) {
                    if(!n->inQueue){
                        rerouteQueue.push_back(n);
                        n->inQueue=true;
                    }
                }
            }
        }

        // If no overflow, we are done!
        if(!of) break;

        // Stagnation Detection
        if(numOF==lastNumOF) stucked++;
        else{
            stucked=0;
            lastNumOF = numOF;
        }

        // Also reroute nets that are detouring too much
        for(auto& net: nets)
            if(net.length>net.HPWL*1.1){
                if(!net.inQueue){
                    rerouteQueue.push_back(&net);
                    net.inQueue=true;
                }
            }
        
        sort(rerouteQueue.begin(), rerouteQueue.end(), [](Net* a, Net* b) {
            if(a->length == b->length)
                return a->overflow > b->overflow;
            return a->length > b->length;
        });

        // Rip-up and Reroute
        for(auto* net: rerouteQueue){
            for(auto& s : net->path) UpdateOccupy(net, s, -1);

            net->path.clear();
            wireLength -= net->length;
            int length = AstarSearch(net->x1, net->y1, net->x2, net->y2, net->path);
            net->length = length;
            wireLength += length;

            for(auto& s : net->path)UpdateOccupy(net, s, 1);
        }

        // Adjust parameters for next iteration
        if(iteration>20){
            coffH*=1.08;
            coffV*=1.08;
        }
        if(coffbbox < 100) coffbbox *= 1.5;
        
        // check id is not overflow
        if(currNetId>2e9){
            for(auto& t: nodeTag) t=0;
            currNetId=0;
        }

        // Time limit check
        auto now = chrono::steady_clock::now();
        auto elapsed = chrono::duration_cast<chrono::seconds>(now - startTime).count();
        if(elapsed>400) break;
    }

    // Post-Refinement (Greedy) ======================================

    // Penalize overflow edges heavily so they become impassable.
    coffH =C*R;
    coffV =C*R;
    // Remove history bias for final length refinement.
    coffhis=0;

    int lastLength=0;

    for(int i=0; i<10; i++){
        vector<Net*> detourNets;
        // Only target nets that are not optimal (length > HPWL)
        for(auto& net: nets)
            if(net.length>net.HPWL) detourNets.push_back(&net);

        sort(detourNets.begin(), detourNets.end(), [](Net* a, Net* b) {return a->length > b->length;});

        for(auto* net: detourNets){
            for(auto& s : net->path) UpdateOccupy(net, s, -1);

            vector<Segment> newpath;
            int length = AstarSearch(net->x1, net->y1, net->x2, net->y2, newpath);

            if(length<net->length){
                wireLength -= net->length;
                net->path = newpath;
                net->length = length;
                wireLength += length;
            }

            for(auto& s : net->path) UpdateOccupy(net, s, 1);
        }

        // If no improvement, stop
        if(wireLength==lastLength) break;
        lastLength=wireLength;
    }

    // Final Report
    // checkLimit();

}

inline void GlobalRouter::UpdateOccupy(Net* n, Segment& s, int du){
    // du = 1 or -1

    auto removNet = [&](vector<Net*>& vec, Net* nn){
        for(size_t i = 0; i < vec.size(); ++i) {
            if(vec[i] == nn) {
                vec[i] = vec.back();
                vec.pop_back();
                return;
            }
        }
    };

    int x1=s.x1; int y1=s.y1;
    int x2=s.x2; int y2=s.y2;
    
    if(du>0){
        if     (x1 < x2) for(int x=x1; x<x2; x++){ gcells[y1*C+x].occH ++;  gcells[y1*C+x].netPassH.push_back(n); }
        else if(x1 > x2) for(int x=x2; x<x1; x++){ gcells[y1*C+x].occH ++;  gcells[y1*C+x].netPassH.push_back(n); }
        else if(y1 < y2) for(int y=y1; y<y2; y++){ gcells[y*C+x1].occV ++;  gcells[y*C+x1].netPassV.push_back(n); }
        else if(y1 > y2) for(int y=y2; y<y1; y++){ gcells[y*C+x1].occV ++;  gcells[y*C+x1].netPassV.push_back(n); }
    }
    else if(du<0){
        if     (x1 < x2) for(int x=x1; x<x2; x++){ gcells[y1*C+x].occH --;  removNet(gcells[y1*C+x].netPassH,n); }
        else if(x1 > x2) for(int x=x2; x<x1; x++){ gcells[y1*C+x].occH --;  removNet(gcells[y1*C+x].netPassH,n); }
        else if(y1 < y2) for(int y=y1; y<y2; y++){ gcells[y*C+x1].occV --;  removNet(gcells[y*C+x1].netPassV,n); }
        else if(y1 > y2) for(int y=y2; y<y1; y++){ gcells[y*C+x1].occV --;  removNet(gcells[y*C+x1].netPassV,n); }   
    }
    
}

int GlobalRouter::AstarSearch(int x1, int y1, int x2, int y2, vector<Segment>& path){

    if(x1==x2 && y1==y2) return 0;

    currNetId++; // Unique ID for this search session
    priority_queue<Node, vector<Node>, cmp> explored;

    // Process a neighbor node
    auto exploring = [&](int idx, int dx, int dy){
        int nIdx = idx + dy*C + dx; // Neighbor Index
        int nx = idx%C +dx;         // Neighbor X
        int ny = idx/C +dy;         // Neighbor Y

        // Check if node was visited in this session
        float gOld = (nodeTag[nIdx] == currNetId) ? costMap[nIdx] : INFINITY;
        float gNew =gOld+1;

        // --- Calculate Cost based on Congestion and History ---
        // Basic idea: Cost = Base + Congestion_Penalty + History_Penalty + Bend_Cost
        
        if(dx<0) {  // Moving Left
            float usage = (float)gcells[nIdx].occH / capH;
            float penalty = (gcells[ nIdx ].occH>=capH) ? coffH*usage*usage : 0.5*usage;
            penalty += coffhis*gcells[ nIdx ].hisOccH *usage*usage;

            gNew = costMap[idx] + penalty +1.5;
            if( parentMap[idx]!=-1 && parentMap[idx]/C != ny) gNew += 0.1;  // Add small cost for turning

        }
        else if(dx>0) { // Moving Right
            float usage = (float)gcells[idx].occH / capH;
            float penalty = (gcells[ idx ].occH>=capH) ? coffH*usage*usage : 0.5*usage;
            penalty += coffhis*gcells[ idx ].hisOccH *usage*usage;

            gNew = costMap[idx] + penalty +1.5;
            if( parentMap[idx]!=-1 && parentMap[idx]/C != ny) gNew += 0.1;
            
        }
        else if(dy<0) { // Moving Down
            float usage = (float)gcells[nIdx].occV / capV;
            float penalty = (gcells[ nIdx ].occV>=capV) ? coffV*usage*usage : 0.5*usage;
            penalty += coffhis*gcells[ nIdx ].hisOccV *usage*usage;

            gNew = costMap[idx] + penalty +1.5;
            if( parentMap[idx]!=-1 && parentMap[idx]%C != nx) gNew += 0.1;

        }
        else if(dy>0){ // Moving Up
            float usage = (float)gcells[idx].occV / capV;
            float penalty = (gcells[ idx ].occV>=capV) ? coffV*usage*usage : 0.5*usage;
            penalty += coffhis*gcells[ idx ].hisOccV *usage*usage;

            gNew = costMap[idx] + penalty +1.5;
            if( parentMap[idx]!=-1 && parentMap[idx]%C != nx) gNew += 0.1;

        }

        // If a cheaper path is found, update and push to queue
        if(gNew < gOld){
            float cost = gNew +  abs(x2-nx) + abs(y2-ny );
            explored.push({nx, ny, cost});
            parentMap[nIdx] = idx;
            costMap[nIdx] = gNew;
            nodeTag[nIdx] = currNetId;
        }
    };

    // Setup start node
    int sIdx = y1*C + x1;
    nodeTag[sIdx] = currNetId; 
    costMap[sIdx] = 0;
    parentMap[sIdx] = -1;
    explored.push({x1, y1, 0});

    bool found=true;

    // Define Bounding Box to limit search area
    float bb = (coffbbox-1);
    int bx1 = min(x1, x2) - abs(x1-x2)*bb;
    int by1 = min(y1, y2) - abs(y1-y2)*bb;
    int bx2 = max(x1, x2) + abs(x1-x2)*bb;
    int by2 = max(y1, y2) + abs(y1-y2)*bb;

    if(bx1<0) bx1=0;
    if(by1<0) by1=0;
    if(bx2>=C) bx2=C-1;
    if(by2>=R) by2=R-1;
    
    // Maze Search (A*)
    while(!explored.empty()){
        Node curr = explored.top(); explored.pop();
        int px = curr.x;
        int py = curr.y;
        int currIdx = py*C + px;

        if(px==x2 && py==y2) {
            found=true;
            break;
        }
        // Explore neighbors within bounding box
        if(px-1 >= bx1) exploring(currIdx, -1,  0);
        if(px+1 <= bx2) exploring(currIdx, +1,  0);
        if(py-1 >= by1) exploring(currIdx,  0, -1);
        if(py+1 <= by2) exploring(currIdx,  0, +1);

    }

    path.clear();

    if(found){
        // Backtrace to reconstruct the path from Target to Source
        int sx = x2;
        int sy = y2;

        int length=1; 
        bool isV = false; // true = vertical; false = horizontal
        int currIdx = parentMap[ sy*C + sx ];

        int tx = currIdx %C;
        int ty = currIdx /C;
        if(tx == sx) isV = true;

        int fromIdx;

        while(currIdx!=sIdx){
            fromIdx = parentMap[ currIdx ];
            if(fromIdx==-1) break; 

            length++;

            tx = fromIdx %C;
            ty = fromIdx /C;
            // Merge consecutive segments into one if direction is same
            if(isV){
                if(sx != tx){
                    path.push_back( {sx, sy, sx, ty} );
                    sy = ty;
                    isV = false;
                }
            }
            else{
                if(sy != ty){
                    path.push_back( {sx, sy, tx, sy} );
                    sx = tx;
                    isV = true;
                } 
            }
            currIdx = fromIdx;
        }

        path.push_back( {sx, sy, x1, y1} ); // Add final segment

        return length;
    }
    else return -1;
}

void GlobalRouter::checkLimit() {
    // --- Wirelength ---
    long long totalLength = 0;
    long long totalHPWL = 0;
    long long diff = 0;
    int perfectNets = 0;

    for (const auto& net : nets) {
        totalLength += net.length;
        totalHPWL += net.HPWL;
        diff += (net.length - net.HPWL);
        if (net.length == net.HPWL) perfectNets++;
    }

    double ratio = (double)totalLength / totalHPWL;

    // --- Overflow ---
    long long calcTotalOverflow = 0;
    int calcMaxOverflow = 0;
    int overflowEdges = 0;

    for(const auto& g : gcells){
        int ovH = g.occH - capH;
        if(ovH > 0) {
            calcTotalOverflow += ovH;
            if(ovH > calcMaxOverflow) calcMaxOverflow = ovH;
            overflowEdges++;
        }

        int ovV = g.occV - capV;
        if(ovV > 0) {
            calcTotalOverflow += ovV;
            if(ovV > calcMaxOverflow) calcMaxOverflow = ovV;
            overflowEdges++;
        }
    }

    cout << "========================================" << endl;
    cout << "           FINAL RESULT REPORT          " << endl;
    cout << "========================================" << endl;

    cout << "[Wirelength Metrics]" << endl;
    cout << "  Total Length : " << totalLength << endl;
    cout << "  Total HPWL   : " << totalHPWL << " (Lower Bound)" << endl;
    cout << "  Total Detour : " << diff << endl;
    cout << "  Detour Ratio : " << fixed << setprecision(5) << ratio << " (" << (ratio-1.0)*100 << "% extra)" << endl;
    cout << "  Perfect Nets : " << perfectNets << " / " << nets.size() 
         << " (" << (double)perfectNets/nets.size()*100 << "%)" << endl;
    cout << endl;

    cout << "[Congestion Metrics]" << endl;
    cout << "  Total Overflow : " << calcTotalOverflow << endl;
    cout << "  Max Overflow   : " << calcMaxOverflow << endl;
    cout << "  Ovfl. Edges    : " << overflowEdges << endl;

    cout << "========================================" << endl;
}