#include "detailedPM.h"

// Constructor for the FastDP placer class
FastDP::FastDP( unordered_map<string,vector<float>>& cs,
                unordered_map<string, Cell>& ce,
                vector<Row>& r,
                unordered_map<string, Component>& co,
                vector<Net>& n)
    // Initialize member references from main
    : CoreSite(cs), cells(ce), rows(r), componentMap(co), nets(n) {

        numRow=rows.size(); // Store number of rows
        // Sort rows primarily by y-coordinate, then x-coordinate
        sort(rows.begin(), rows.end(),
            [](Row a, Row b) {
                if (a.oy != b.oy)
                    return a.oy < b.oy;
                return a.ox < b.ox;
            });
        
        // Build a map from y-coordinate to row indices for fast lookup
        y2RowIdx.reserve(numRow);
        for(int i=0; i<numRow; i++){
            y2RowIdx[rows[i].oy].push_back(i);
        }

        // Create a flat vector of component pointers from the map
        components.reserve(componentMap.size());
        for(auto& c: componentMap){
            components.push_back(&c.second);
        }

        // placement vector to hold all rows
        placement.resize(numRow); 
        for(auto& comps: placement){
            comps.reserve(components.size()/numRow*2);
        }

        // Assign initial row index (ridx) and placement for each component
        for(auto* c:components){
            if(c->c->ctype==0){ // Movable components (ctype == 0)
                int r = findRowIndex(c->x, c->y); // Find initial row
                c->ridx = r;
                placement[ r ].push_back(c); // Add to placement vector
                rows[r].occupy += c->c->w / rows[r].dx; // Add to initial row occupancy
            }
            else{ // Fixed components (ctype != 0)

                int ridx;

                // Find the first row whose y-origin is >= component's y
                auto it = std::upper_bound(rows.begin(), rows.end(), c->y,
                                           [](int value, const Row& row){ return value < row.oy; });
                
                if(it == rows.begin()) {                      // Component is at or below the first row
                    if(c->y + c->c->h < rows[0].oy) continue; // Component is entirely below, skip
                    else ridx =0;                             // Belongs to first row
                }
                else 
                    ridx = (it - rows.begin() - 1); // Get index of the row just below

                // Check if component is entirely above the last row
                if(ridx==numRow-1){
                    if(c->y > rows[ridx].oy+rows[ridx].dy)
                        continue;
                }
                c->ridx = ridx; // Set initial row index

                int yMax=c->y + c->c->h; // Top edge of the component

                // Handle multi-row components (e.g., macros)
                while( yMax > rows[ridx].oy && ridx<numRow){
                    // Check if component is outside this row's x-bounds
                    if( c->x + c->c->w <= rows[ridx].ox || c->x >= rows[ridx].ox + rows[ridx].dx*rows[ridx].nw ) {
                        ridx++;
                        continue;
                    }

                    // Component overlaps with this row
                    placement[ ridx ].push_back(c);

                    // Calculate legalized x-span and occupancy this component covers *in this row*
                    int x1=legalizeX(c->x, rows[ridx].oy);
                    int x2=legalizeXback(c->x + c->c->w, rows[ridx].oy);

                    int doccupy = ( x2-x1 )/rows[ridx].dx;
                    rows[ridx].occupy += doccupy; // Add fixed component's occupancy

                    ridx++; // Check next row
                }
            }
        }

        // segment initialization
        for (int r=0; r<rows.size(); r++){
            auto& comps = placement[r];
            // Sort components by x to prepare for segment generation
            sort(comps.begin(), comps.end(), [](Component* a, Component* b){ return a->x < b->x; });

            
            auto& row = rows[r];
            int sidx = 0;       // Current segment index
            bool lastFix=false; // Flag if last component was fixed

            int rowEnd = row.ox+row.dx*row.nw; // Row's end x-coordinate

            // Start with one giant segment covering the whole row
            row.segments.push_back({row.ox, rowEnd, 0, 0});

            for(auto& c: comps){
                if(c->fixed) {
                    // less than row x origin
                    if(c->x + c->c->w < row.ox) continue;
                    // Handle case where multiple fixed cells are contiguous
                    if(lastFix && c->x <= row.segments[sidx].x_start){
                        int tXend=legalizeXback(c->x + c->c->w, row.oy);
                        if(tXend > row.segments[sidx].x_start){
                            // Shrink the start of the current segment
                            row.segments[sidx].x_start = tXend;
                            if(tXend > rowEnd){ // Check if segment is fully occluded
                                row.segments.pop_back(); // Remove invalid segment
                                sidx--;
                            }
                        }
                    }
                    else{
                        // Cut the current segment at the fixed cell's start
                        row.segments[sidx].x_end = legalizeX(c->x, row.oy);
                        row.segments[sidx].capacity = (row.segments[sidx].x_end - row.segments[sidx].x_start) / row.dx;
                        
                        // Create a new segment starting after this fixed cell
                        int tXend=legalizeXback(c->x + c->c->w, row.oy);
                        if(tXend<rowEnd){ // Only add if new segment is valid
                            row.segments.push_back( { tXend , rowEnd , 0 , 0} ); // (Bug: x_end is wrong)
                            sidx++;
                        }
                    }
                    lastFix=true;
                }
                else{ // Movable cell
                    // Add its occupancy to the current segment
                    row.segments[sidx].occupy += c->c->w / rows[r].dx;
                    lastFix=false;
                }
            }
            // Set capacity for the last segment
            row.segments[sidx].capacity = (row.segments[sidx].x_end-row.segments[sidx].x_start) / row.dx;
        }
}

void FastDP::Optimize(chrono::time_point<std::chrono::steady_clock> startTime){

    // Defines a structure to hold a complete placement solution state
    struct Solution{
        unsigned long long HPWL;
        unordered_map<string,Component> cMap;
        vector<Net> ns;
        vector<vector<Component*>> pl;
    };

    Solution bestSol;        // Holds the best solution found so far
    bestSol.HPWL=ULLONG_MAX; // Initialize best HPWL to "infinity"

    SingleSegmentClustering();

    // Sort components by net count to process high-degree nodes first
    sort(components.begin(), components.end(), [](Component* a, Component* b){ return a->net.size() > b->net.size(); });

    // Set Threshold based on problem size
    if(components.size()<300000)
        Threshold = -2000;
    else
        Threshold = -500;
    
    // Main optimization loop
    for(int i=0; i<100; i++){

        if(i>0){
            // Shuffle component processing order for randomness
            random_device rd;
            mt19937 g(rd());
            shuffle(components.begin(), components.end(), g);

            Threshold*=0.96; // cooling
        }

        // Perform placement improvement moves
        GlobalSwap();
        VerticalSwap();
        ReOrdering();

        // Calculate total HPWL for the current state
        unsigned long long currHPWL=0;
        for(const auto& n: nets)
            currHPWL += n.HPWL;
        // cout<<currHPWL<<endl;

        // Save the best solution found
        if(currHPWL < bestSol.HPWL){
            bestSol.HPWL=currHPWL;
            bestSol.cMap = componentMap;
            bestSol.ns = nets;           
            bestSol.pl = placement; 
        }

        // Check for timeout
        auto now = chrono::steady_clock::now();
        auto elapsed = chrono::duration_cast<chrono::seconds>(now - startTime).count();
        if (elapsed >= 250)
            break;
    }

    // Restore the best solution found
    if(bestSol.HPWL!=ULLONG_MAX){
        for(auto& c: componentMap){
            c.second.x = bestSol.cMap[c.first].x;
            c.second.y = bestSol.cMap[c.first].y;
            c.second.ridx = bestSol.cMap[c.first].ridx;
            if(!c.second.fixed)
                c.second.orient = rows[c.second.ridx].orient;
        }
        nets = bestSol.ns;
        placement = bestSol.pl;
    }

    SingleSegmentClustering();

    // Set final orientation for all components
    for(auto& c: components){
        if(!c->fixed)
            c->orient = rows[c->ridx].orient;
    }
}

// Find row index for a given (x,y) coordinate
int FastDP::findRowIndex(int x, int y){
    auto& idxList = y2RowIdx[y]; // Get candidate rows for this y
    if(idxList.size()==1) return idxList[0];

    if(y>=rows[numRow-1].oy) return numRow-1; // Above last row
    else if(y<=rows[0].oy) return 0;          // Below first row

    else{ // Check x-coordinate against rows with matching y
        for(int i : idxList){
            Row& row=rows[i];
            if( x>= row.ox && x < row.ox+row.dx*row.nw )
                return i;
        }
    }

    return -1; // Not found
}

// Find segment index for a given row and x-coordinate
int FastDP::findSegmentIndex(int r, int x) {
    vector<Segment>& segment = rows[r].segments;
    if (segment.size() ==1) return 0; // Optimization for single-segment rows
    else if (segment.size() <= 10) {
        // Small count -> linear scan
        for (int i = 0; i < segment.size(); ++i) {
            if (x >= segment[i].x_start && x < segment[i].x_end)
                return i;
        }
        return - 1; // Not in any segment
    } 
    else {
        // Large count -> binary search (lower_bound)
        auto it = std::lower_bound(segment.begin(), segment.end(), x,
            [](const Segment& a, int value) {
                return a.x_end < value; // Find first segment whose x_end is not less than x
            }
        );
        int idx = int(it - segment.begin());
        if (idx == 0) return 0; // Belongs to first segment
        if (idx >= (int)segment.size()) idx = segment.size() - 1; // Belongs to last segment
        return idx;
    }
}

// Check if a component's symmetry allows a given orientation
bool FastDP::ValidOrient(int sym, string orient){
    if(orient=="N") return true;
    else if(orient=="FS")
        if((sym & SYM_X) == SYM_X) return true;
    else if(orient=="FN")
        if((sym & SYM_Y) == SYM_Y) return true;
    else if(orient=="S")
        if ((sym & (SYM_X | SYM_Y)) == (SYM_X | SYM_Y)) return true;
    
    return false;
}

// Check if swapping two components is valid (capacity check)
bool FastDP::ValidSwap(Component* c1, Component* c2, int r1, int r2){
    int dOccupy1 = c1->c->w/rows[r1].dx - c2->c->w/rows[r1].dx;
    int dOccupy2 = c1->c->w/rows[r2].dx - c2->c->w/rows[r2].dx;

    // Check row-level capacity
    if(rows[r2].occupy + dOccupy2 > rows[r2].nw || rows[r1].occupy - dOccupy1 > rows[r1].nw) return false;
    
    // Find segments for both components
    Segment* s1 = &rows[r1].segments[ findSegmentIndex(r1,c1->x) ];
    Segment* s2 = &rows[r2].segments[ findSegmentIndex(r2,c2->x) ];

    // Check segment-level capacity
    if(s2->occupy + dOccupy2 > s2->capacity || s1->occupy - dOccupy1 > s1->capacity) return false;

    return true;
}

// Check if moving a component to a space is valid (capacity check)
bool FastDP::ValidSwap(Component* c1, int x2, int r1, int r2){
    if(findSegmentIndex(r2,x2)==-1) return false; // Target space is not in a valid segment

    int dOccupy1 = c1->c->w/rows[r1].dx;
    int dOccupy2 = c1->c->w/rows[r2].dx;

    // Check row-level capacity
    if(rows[r2].occupy + dOccupy2 > rows[r2].nw || rows[r1].occupy - dOccupy1 > rows[r1].nw) return false;
    
    // Find segments
    Segment* s1 = &rows[r1].segments[ findSegmentIndex(r1,c1->x) ];
    Segment* s2 = &rows[r2].segments[ findSegmentIndex(r2,x2) ];

    // Check segment-level capacity
    if(s2->occupy + dOccupy2 > s2->capacity || s1->occupy - dOccupy1 > s1->capacity) return false;

    return true;
}

// Update occupancy data after a cell swap
void FastDP::UpdateOccupy(Component* c1, Component* c2, int r1, int r2){

    int dOccupy1 = c1->c->w/rows[r1].dx - c2->c->w/rows[r1].dx;
    int dOccupy2 = c1->c->w/rows[r2].dx - c2->c->w/rows[r2].dx;
    
    rows[r1].occupy -= dOccupy1;
    rows[r2].occupy += dOccupy2;
    
    Segment* s1 = &rows[r1].segments[ findSegmentIndex(r1,c1->x) ];
    Segment* s2 = &rows[r2].segments[ findSegmentIndex(r2,c2->x) ];

    s1->occupy -= dOccupy1;
    s2->occupy += dOccupy2;

    return;
}

// Update occupancy data after a move-to-space
void FastDP::UpdateOccupy(Component* c1, int x2, int r1, int r2){
    
    int dOccupy1 = c1->c->w/rows[r1].dx;
    int dOccupy2 = c1->c->w/rows[r2].dx;
    
    rows[r1].occupy -= dOccupy1;
    rows[r2].occupy += dOccupy2;

    
    Segment* s1 = &rows[r1].segments[ findSegmentIndex(r1,c1->x) ];
    Segment* s2 = &rows[r2].segments[ findSegmentIndex(r2,x2) ];

    s1->occupy -= dOccupy1;
    s2->occupy += dOccupy2;

    return;
}

// Get component's start x-coordinate from placement vector
int FastDP::getX(int r, int i) {
    if (i < 0) return rows[r].ox; // Left boundary of row
    if (i >= placement[r].size()) return rows[r].ox + rows[r].dx * rows[r].nw; // Right boundary of row
    return placement[r][i]->x;
};

// Get component's end x-coordinate from placement vector
int FastDP:: getEndX(int r, int i) {
    if (i < 0) return rows[r].ox; // Left boundary of row
    if (i >= placement[r].size()) return rows[r].ox + rows[r].dx * rows[r].nw; // Right boundary of row
    return getX(r, i) + placement[r][i]->c->w;
};

// Legalize x by rounding down to the nearest site boundary
int FastDP::legalizeX(int x, int y){
    
    Row* row = &rows[ findRowIndex(x,y) ];
    if( (x-row->ox) % row->dx ==0)
        return x;
    else{
        int a=( (x-row->ox) / row->dx) * row->dx + row->ox;

        return a;
    }
}

// Legalize x by rounding up to the nearest site boundary
int FastDP::legalizeXback(int x, int y){
    
    Row* row = &rows[ findRowIndex(x,y) ];
    if( (x-row->ox) % row->dx ==0)
        return x;
    else{
        return (x - row->ox + row->dx - 1) / row->dx * row->dx + row->ox;
    }
}

// Find optimal x-range (median) for a component based on its nets
void FastDP::FindOptXs(Component* c2mov, int& x1, int& x2){
    vector<int> xList;

    xList.reserve(c2mov->net.size()*2);

    for(auto& net : c2mov->net){
        int nMinX, nMaxX, nMinY, nMaxY;
        if(net->pminX==-1){ // No fixed PINs
            nMinX=2147483647;
            nMaxX=0;
        }else{ // Has fixed PINs
            nMinX=net->pminX;
            nMaxX=net->pmaxX;
        }
        // Find BBox of other components in the net
        for(auto& c: net->c){
            if(c->x < nMinX) nMinX=c->x;
            if(c->x > nMaxX) nMaxX=c->x;
        }
        xList.push_back(nMinX);
        xList.push_back(nMaxX);
    }

    sort(xList.begin(),xList.end());

    // Set optimal range to the median
    x1=xList[xList.size()/2-1];
    x2=xList[xList.size()/2];
}

// Single Segment Clustering
void FastDP::SingleSegmentClustering(){

    #pragma omp parallel for schedule(dynamic)
    for(int ridx=0; ridx<numRow; ridx++){
        vector<Component*> &comps = placement[ridx];

        int cidx =0; // Component index

        for(auto& seg: rows[ridx].segments){ // Iterate through segments in this row

            // Find first movable component in this segment
            while(cidx<comps.size() && ( comps[cidx]->x < seg.x_start || comps[cidx]->fixed ) ) cidx++;
            if(cidx>=comps.size()) break; // No more components in row
            if(comps[cidx]->x > seg.x_end) continue; // No components in this segment

            int startIdx=cidx; // Start of a movable cluster
            int endIdx=-1;

            vector<int> optX1; // Store optimal x1 for each cell
            vector<int> optX2; // Store optimal x2 for each cell
            
            // Collect all movable components in this segment
            while( cidx<comps.size() ){
                if(comps[cidx]->x >= seg.x_end || comps[cidx]->fixed) break;

                int x1,x2;
                FindOptXs(comps[cidx],x1,x2); // Get optimal x-range
                x1=clamp(x1, seg.x_start, seg.x_end);
                x2=clamp(x2, seg.x_start, seg.x_end);
                optX1.push_back(x1);
                optX2.push_back(x2);

                cidx++;
            }
            cidx--;
            endIdx=cidx;

            // Perform clustering legalization
            vector<Cluster> clusters;
            for(int i=startIdx; i<=endIdx; i++){
                Component* c = comps[i];

                int x1=optX1[i-startIdx];
                int x2=optX2[i-startIdx];
                // Calculate cost for placing at each site
                vector<int> bounds( seg.capacity );
                for(int j=0; j<seg.capacity; j++){
                    int xj = seg.x_start + j*rows[ridx].dx;
                    if( xj >=x1 && xj<=x2) bounds[j] = 0;
                    else if( xj < x1 )     bounds[j] = xj-x1;
                    else                   bounds[j] = xj-x2;
                }
                
                // Add new cluster for this component
                clusters.push_back({i, legalizeXback(optX1[i-startIdx], rows[ridx].oy), 1, c->c->w, bounds});
                Cluster* curr = &clusters.back();
                Cluster* last;

                int idx = clusters.size()-2;
                
                // Merge with overlapping previous clusters
                while(idx>=0){
                    last = &clusters[idx];
                    if(curr->x >= last->x+last->w) break; // No overlap

                    
                    // Update bounds cost for merged cluster
                    for(int j=0; j<seg.capacity; j++){
                        if( j+last->w/rows[ridx].dx <seg.capacity){
                            last->bounds[j] += curr->bounds[j + last->w/rows[ridx].dx ];
                        }
                        else{
                            last->bounds[j] += seg.capacity*rows[ridx].dx*10; // High penalty
                        }
                    }

                    // Find new optimal position for the merged cluster
                    auto min_it = std::min_element(last->bounds.begin(), last->bounds.end(),
                        [](int a, int b){ return std::abs(a) < std::abs(b); });
                    int min_index = std::distance(last->bounds.begin(), min_it);


                    // Update merged cluster properties
                    last->x = legalizeX( seg.x_start + min_index*rows[ridx].dx , c->y );
                    last->e += curr->e;
                    last->w += curr->w;

                    curr=last;
                    clusters.pop_back(); // Remove the cluster that was merged

                    idx--;
                }
            }
            // Legalize cluster positions based on segment boundaries
            legalizeBoundary(clusters, seg.x_start, seg.x_end, rows[ridx].oy);


            // Apply new x-coordinates to components
            for(auto& clus : clusters){
                int num = clus.e;
                int currX = clus.x;
                for(int i=clus.idx; i<clus.idx+clus.e; i++){
                    comps[i]->x=currX;
                    currX += comps[i]->c->w;
                }
            }       
        }
    }

    // Re-calculate HPWL for all nets after legalization
    for(auto& net : nets){
        if(net.pminX==-1){
            net.minX=net.c[0]->x;
            net.maxX=net.c[0]->x;
        }else{
            net.minX=net.pminX;
            net.maxX=net.pmaxX;
        }

        for(auto& c: net.c){
            if(c->x < net.minX)
                net.minX=c->x;
            else if(c->x > net.maxX)
                net.maxX=c->x;
        }
        net.HPWL = (net.maxX-net.minX) + (net.maxY-net.minY);
    }

}

// Performs local window re-ordering (3-cell window)
void FastDP::ReOrdering(){

    int order[6][3]={{0,1,2},{0,2,1},{1,0,2},{1,2,0},{2,0,1},{2,1,0}};

    for (auto& comps : placement){ // Iterate all rows
        if(comps.size() < 3) continue;
        for(int i=0; i<comps.size()-2; i++){ // Sliding window of 3 cells

            if(comps[i]->fixed || comps[i+1]->fixed || comps[i+2]->fixed) continue;

            int minIdx =0;
            int minHPWL=0;

            // Collect affected nets
            vector<Net*> netChange;
            unordered_map<Net*, vector<Component*>> cellAffected;
            unordered_set<Net*> seen;

            Component* cTmp[3] = { comps[i], comps[i+1], comps[i+2] };
            
            for (auto* cn: cTmp)
                for (auto* n : cn->net){
                    if (seen.insert(n).second) // Returns pair<iterator,bool>
                        netChange.push_back(n);
                    cellAffected[n].push_back(cn);
                }

            // Check all 5 permutations (0 is original)
            for(int o=1; o<6; o++){

                Component* cNews[3];
                int xNews[3];

                // Get components in new order
                cNews[0]= comps[i+order[o][0]];
                cNews[1]= comps[i+order[o][1]];
                cNews[2]= comps[i+order[o][2]];

                // Calculate new x positions
                xNews[0]= comps[i]->x;
                xNews[2]= comps[i+2]->x + comps[i+2]->c->w - cNews[2]->c->w;
                xNews[1]= legalizeX((xNews[2] +xNews[0] +cNews[0]->c->w -cNews[1]->c->w)/2 , cNews[1]->y);

                int dHPWL=0;
                // Calculate dHPWL incrementally
                for(int nidx=0; nidx<netChange.size(); nidx++){
                    Net* net = netChange[nidx];
                    
                    int nMinX=net->minX;
                    int nMaxX=net->maxX;

                    bool affected = false;
                    
                    // Fast check: see if BBox boundaries are affected
                    for(auto* cn : cellAffected[net]){
                        for(int j=0; j<3; j++){
                            if(cn==cNews[j]){
                                if(xNews[j]<nMinX)
                                    nMinX=xNews[j];
                                else if(cn->x==net->minX){
                                    affected=true;
                                    break;
                                }
                                
                                if(xNews[j]>nMaxX)
                                    nMaxX=xNews[j];
                                else if(cn->x==net->maxX){
                                    affected=true;
                                    break;
                                }
                            }
                        }
                    }

                    if(affected){ // Slow check: recalculate BBox
                        if(net->pminX!=-1){
                            nMinX=net->pminX;
                            nMaxX=net->pmaxX;
                        }else{
                            nMinX=2147483647;
                            nMaxX=0;
                        }

                        for(auto* c: net->c){
                            
                            if(c==cNews[0]){
                                if(xNews[0] < nMinX) nMinX=xNews[0];
                                if(xNews[0] > nMaxX) nMaxX=xNews[0];
                            }
                            else if(c==cNews[1]){
                                if(xNews[1] < nMinX) nMinX=xNews[1];
                                if(xNews[1] > nMaxX) nMaxX=xNews[1];
                            }
                            else if(c==cNews[2]){
                                if(xNews[2] < nMinX) nMinX=xNews[2];
                                if(xNews[2] > nMaxX) nMaxX=xNews[2];
                            }
                            else{
                                if(c->x < nMinX) nMinX=c->x;
                                if(c->x > nMaxX) nMaxX=c->x;
                            }
                        }
                    }

                    dHPWL += (nMaxX-nMinX) - (net->maxX - net->minX); // Add delta xHPWL
                }               
                if(dHPWL<minHPWL){ // Check if this order is better
                    minIdx = o;
                    minHPWL  = dHPWL;
                }
            }

            if(minIdx==0) continue; // Original order was best

            // Apply the best permutation found
            Component* cNews[3];
            int xNews[3];

            cNews[0]= comps[i+order[minIdx][0]];
            cNews[1]= comps[i+order[minIdx][1]];
            cNews[2]= comps[i+order[minIdx][2]];

            xNews[0]= comps[i]->x;
            xNews[2]= comps[i+2]->x + comps[i+2]->c->w - cNews[2]->c->w;
            xNews[1]= legalizeX((xNews[2] +xNews[0] +cNews[0]->c->w -cNews[1]->c->w)/2 , cNews[1]->y );

            // Update component x-coords and placement vector
            cNews[0]->x = xNews[0];     cNews[1]->x = xNews[1];     cNews[2]->x = xNews[2];
            comps[i]    = cNews[0];     comps[i+1]  = cNews[1];     comps[i+2]  = cNews[2];

            // Update Net BBoxes and HPWL
            for(int nidx=0; nidx<netChange.size(); nidx++){
                Net* net=netChange[nidx];

                if(net->pminX!=-1){
                    net->minX=net->pminX;
                    net->maxX=net->pmaxX;
                }else{
                    net->minX=net->c[0]->x;
                    net->maxX=net->c[0]->x;
                }

                for(auto& c: net->c){
                    if(c->x < net->minX) net->minX=c->x;
                    else if(c->x > net->maxX) net->maxX=c->x;
                }

                net->HPWL = (net->maxX-net->minX) + (net->maxY-net->minY);
            }
        }
    }
}

// Find optimal region for a component based on its nets
void FastDP::FindOptRegion(Component* c2mov, int& x1, int& x2, int& y1, int& y2){
    vector<int> xList;
    vector<int> yList;

    xList.reserve(c2mov->net.size()*2);
    yList.reserve(c2mov->net.size()*2);

    for(auto& net : c2mov->net){
        bool changed=false;
        // Check if component is on the boundary of its net's BBox
        if(c2mov->x==net->minX || c2mov->x==net->maxX || c2mov->y==net->minY || c2mov->y==net->maxY)
            changed=true;

        if(changed){ // Recalculate BBox *without* this component
            int nMinX, nMaxX, nMinY, nMaxY;
            if(net->pminX==-1){
                nMinX=2147483647;
                nMaxX=0;
                nMinY=2147483647;
                nMaxY=0;
            }else{
                nMinX=net->pminX;
                nMaxX=net->pmaxX;
                nMinY=net->pminY;
                nMaxY=net->pmaxY;
            }
            for(auto& c: net->c){
                if(c!=c2mov){
                    if(c->x < nMinX) nMinX=c->x;
                    if(c->x > nMaxX) nMaxX=c->x;
                    if(c->y < nMinY) nMinY=c->y;
                    if(c->y > nMaxY) nMaxY=c->y;
                }
            }
            xList.push_back(nMinX);
            xList.push_back(nMaxX);
            yList.push_back(nMinY);
            yList.push_back(nMaxY);
        }
        else{ // Component is not on boundary, use existing BBox
            xList.push_back(net->minX);
            xList.push_back(net->maxX);
            yList.push_back(net->minY);
            yList.push_back(net->maxY);
        }
    }
    sort(xList.begin(),xList.end());
    sort(yList.begin(),yList.end());

    // Set optimal region to the median of BBox coordinates
    x1=xList[xList.size()/2-1];
    x2=xList[xList.size()/2];
    y1=yList[yList.size()/2-1];
    y2=yList[yList.size()/2];
}

// Calculate overlap penalty for swapping two cells
long long FastDP::SwapPenalty(Component* c1, Component* c2, int r, int idx){
    
    int wdiff = c1->c->w - c2->c->w; // Width difference

    // P1: Overlap with immediate neighbors
    long long P1 = wdiff - ( getX(r, idx+1) - getEndX(r, idx-1) - placement[r][idx]->c->w );

    // P2: Overlap including next-nearest neighbors
    long long P2 = P1;
    if(idx>0)
        P2 -= (getX(r, idx - 1) - getEndX(r, idx-2));

    if(idx<placement[r].size()-1)
        P2 -= (getX(r, idx+2) - getEndX(r, idx + 1));

    // Return weighted penalty
    if(P2>0) return P1*2 + P2*10;
    else if(P1>0) return P1*2;
    else     return 0; // No penalty if no overlap (P1 <= 0)
}

// Calculate overlap penalty for moving a cell to a space
long long FastDP::SwapPenalty(Component* c1, int r, int idx){
    
    int wdiff = c1->c->w;

    // P1: Overlap with immediate neighbors
    long long P1 = wdiff - ( getX(r, idx+1) - getEndX(r, idx));

    // P2: Overlap including next-nearest neighbors
    long long P2 = P1 - (getX(r, idx) - getEndX(r, idx-1)) - (getX(r, idx+2) - getEndX(r, idx + 1));

    // Return weighted penalty
    if(P2>0) return P1*2 + P2*10;
    else if(P1>0) return P1*2;
    else     return 0; // No penalty if no overlap (P1 <= 0)
}

// Incrementally calculate dHPWL for a list of moved components
long long FastDP::dHPWLafterSwap(vector<Component*>& cList, vector<int>& newX, vector<int>& newY, vector<NetInfo>& netInfos){
    
    vector<Net*> netChange; // List of nets affected
    unordered_map<Net*, vector<int>> idxAffected; // Map net -> index of cell in cList
    unordered_set<Net*> seen; // Set to track unique nets
    
    // Collect all unique nets affected by the components in cList
    for (int cidx=0; cidx<cList.size(); cidx++){
        for (auto* n : cList[cidx]->net){
            if (seen.insert(n).second) // Returns pair<iterator,bool>
                netChange.push_back(n);
            idxAffected[n].push_back(cidx);
        }
    }
    
    long long dHPWL=0; // Total change in HPWL
    for(auto* net: netChange){
        
        // Start with current net BBox
        int nMinX=net->minX;
        int nMaxX=net->maxX;
        int nMinY=net->minY;
        int nMaxY=net->maxY;

        bool affected = false; // Flag if a boundary component moved
        
        // Fast check: see if any moved component was on the BBox boundary
        for(int cidx : idxAffected[net]){

            Component* cn = cList[cidx];
            if(newX[cidx]<nMinX) nMinX=newX[cidx];
            else if(cn->x==net->minX){
                affected=true; break;
            }
            if(newX[cidx]>nMaxX) nMaxX=newX[cidx];
            else if(cn->x==net->maxX){
                affected=true; break;
            }

            if(newY[cidx]<nMinY) nMinY=newY[cidx];
            else if(cn->y==net->minY){
                affected=true; break;
            }
            if(newY[cidx]>nMaxY) nMaxY=newY[cidx];
            else if(cn->y==net->maxY){
                affected=true; break;
            }

        }

        if(affected){ // Slow check: BBox boundary moved, must recalculate
            // Reset BBox based on fixed PINs
            if(net->pminX!=-1){
                nMinX=net->pminX;
                nMaxX=net->pmaxX;
                nMinY=net->pminY;
                nMaxY=net->pmaxY;
            }else{ // No fixed PINs
                nMinX=2147483647;
                nMaxX=0;
                nMinY=2147483647;
                nMaxY=0;
            }

            // Recalculate BBox from all components in the net
            for(auto* c: net->c){
                bool ismov=false;
                // Check if this component is one of the moved components
                for(int cidx : idxAffected[net]){
                    if(c==cList[cidx]){
                        // Use its new position
                        if(newX[cidx] < nMinX) nMinX=newX[cidx];
                        if(newX[cidx] > nMaxX) nMaxX=newX[cidx];
                        if(newY[cidx] < nMinY) nMinY=newY[cidx];
                        if(newY[cidx] > nMaxY) nMaxY=newY[cidx];

                        ismov = true;
                    }
                }

                if(!ismov){ // Not a moved component, use its current position
                    if(c->x < nMinX) nMinX=c->x;
                    if(c->x > nMaxX) nMaxX=c->x;
                    if(c->y < nMinY) nMinY=c->y;
                    if(c->y > nMaxY) nMaxY=c->y;
                }
            }
        }

        int newHPWL = (nMaxX-nMinX + nMaxY-nMinY);
        dHPWL += net->HPWL - newHPWL;
        
        // Store new BBox info to be applied later
        netInfos.push_back( {net, nMinX, nMinY, nMaxX, nMaxY, newHPWL} );
    }

    return dHPWL;
}

// Legalize placement for all rows
void FastDP::legalizeRows(){

    #pragma omp parallel for schedule(dynamic) // Parallelize row legalization
    for(int ridx=0; ridx<numRow; ridx++){
        
        vector<Component*> comps = placement[ridx]; // Get components in this row
        vector<Cluster> clusters; // Vector to hold clusters of movable cells
        
        // Iterate through components in the row
        for(int ci=0; ci<comps.size(); ci++){
            Component* c = comps[ci];
            if(!c->fixed){ // Movable component
                // Create a new cluster for this component
                clusters.push_back({ci, c->x, 1, c->c->w});
                Cluster* curr = &clusters.back();

                int idx = clusters.size()-2;

                Cluster* last;
                
                // Merge with overlapping clusters to its left
                while(idx>=0){
                    last = &clusters[idx];
                    if(curr->x >= last->x+last->w) break; // No overlap

                    // Merge curr into last
                    last->x = legalizeX( (last->e*last->x + curr->e*(curr->x-last->w)) / (last->e + curr->e),c->y);
                    last->e += curr->e;
                    last->w += curr->w;

                    curr=last;
                    clusters.pop_back(); // Remove curr

                    idx--;
                }
            }
            else if(!clusters.empty()){ // Fixed component
                // Legalize the preceding cluster(s)
                legalizeBoundary(clusters,  getEndX(ridx, clusters[0].idx-1), c->x, rows[ridx].oy);

                // Apply legalized x-coordinates to components
                for(auto& clus : clusters){
                    int num = clus.e;
                    int currX = clus.x;
                    for(int i=clus.idx; i<clus.idx+clus.e; i++){
                        
                        comps[i]->x= legalizeXback(currX, c->y);
                        currX += comps[i]->c->w;
                    }
                }
                clusters.clear(); // Start new cluster
            }
        }

        if(!clusters.empty()){ // Legalize last cluster in the row
            legalizeBoundary(clusters,  getEndX(ridx, clusters[0].idx-1), getX(ridx, comps.size()), rows[ridx].oy);

            // Apply legalized x-coordinates
            for(auto& clus : clusters){
                int num = clus.e;
                int currX = clus.x;
                for(int i=clus.idx; i<clus.idx+clus.e; i++){
                    comps[i]->x=legalizeXback(currX, rows[ridx].oy);
                    currX += comps[i]->c->w;
                }
            }
        }
    }

    // Update HPWL for all nets after row legalization
    for(auto& net : nets){
        if(net.pminX==-1){
            net.minX=net.c[0]->x;
            net.maxX=net.c[0]->x;
        }else{
            net.minX=net.pminX;
            net.maxX=net.pmaxX;
        }

        for(auto& c: net.c){
            if(c->x < net.minX)
                net.minX=c->x;
            else if(c->x > net.maxX)
                net.maxX=c->x;
        }
        net.HPWL = (net.maxX-net.minX) + (net.maxY-net.minY);
    }
}

// Legalize a set of clusters within given start/end boundaries
void FastDP::legalizeBoundary(vector<Cluster>& clusters, int startX, int endX, int y){
    if(clusters.empty()) return;

    Cluster* last = &clusters.back(); // Get the rightmost cluster

    // Check for overflow at the right boundary (endX)
    if(last->x+last->w > endX){
        
        int idx = clusters.size()-2;
        Cluster* curr = last;
        curr->x = legalizeX(endX - curr->w, y); // Snap to boundary

        // Propagate the shift leftwards by merging
        while(idx>=0){
            last = &clusters[idx];
            if(last->x+last->w <= curr->x) break; // No overlap

            // Merge curr into last
            last->e += curr->e;
            last->w += curr->w;
            last->x = legalizeX(endX - last->w, y); // Re-snap merged cluster

            curr=last;
            clusters.pop_back(); // Remove curr

            idx--;
        }
    }

    Cluster* front = &clusters[0]; // Get the leftmost cluster

    // Check for overflow at the left boundary (startX)
    if(front->x < startX){
        front->x = legalizeXback(startX, y); // Snap to boundary

        // Propagate the shift rightwards by merging
        if(clusters.size() > 1){
            int merge_idx = 1; // Track next cluster to check
            while(merge_idx < clusters.size()){
                if(clusters[merge_idx].x < front->x + front->w){
                    // Merge cluster[merge_idx] into front (clusters[0])
                    front->e += clusters[merge_idx].e;
                    front->w += clusters[merge_idx].w;
                    merge_idx++; // Check next
                }
                else{
                    break; // No more overlap
                }
            }
            
            // If merging happened (merge_idx > 1), do deletion
            if(merge_idx > 1){
                // One-shot erase of all merged elements [1, merge_idx-1]
                clusters.erase(clusters.begin() + 1, clusters.begin() + merge_idx);
            }
        }
    }
}

void FastDP::GlobalSwap(){
    // Iterate through all movable components
    for(auto* cfrom : components){
        if(cfrom->fixed) continue; // Skip fixed components

        int x1, x2, y1, y2;
        // Find the optimal median region for this component based on its nets
        FindOptRegion(cfrom, x1, x2, y1, y2);

        // Skip if the component is already inside its optimal region
        if( x1 < cfrom->x && cfrom->x < x2 && y1 < cfrom->y && cfrom->y < y2)
            continue;

        int rfrom = cfrom->ridx; // Get the component's current row index
        
        // Find the component's index within its current row's placement vector
        auto& comps0 = placement[rfrom];
        auto it0 = std::find(comps0.begin(), comps0.end(), cfrom);
        int idxfrom = it0 - comps0.begin();

        // Find the range of target rows that fall within the optimal y-region
        // r1 = first row with oy >= y1
        auto it1 = lower_bound(
            rows.begin(), rows.end(), y1,
            [](const Row& a, int value){
                return a.oy < value;
            }
        );
        int r1 = it1 - rows.begin();

        // r2 = last row with oy <= y2
        auto it2 = upper_bound(
            rows.begin(), rows.end(), y2,
            [](int value, const Row& a){
                return value < a.oy; // upper_bound: first oy > y2
            }
        );
        int r2 = (it2 - rows.begin()) - 1;

        // Create a list of candidate row indices
        vector<int> rList;
        rList.reserve(r2-r1+1);
        for(int rr=r1; rr<=r2; rr++) rList.push_back(rr);

        // Shuffle the candidate rows to randomize the search
        random_device rd;
        mt19937 g(rd());
        shuffle(rList.begin(), rList.end(), g);

        bool found = false; // Flag to stop searching once a move is made
        for(int r : rList){ // Iterate through the shuffled target rows

            // Check if site type and orientation are compatible
            if(rows[r].site!=rows[rfrom].site || !ValidOrient(cfrom->c->sym, rows[r].orient)) continue;

            // Pick a random x-coordinate within the optimal x-region
            int randNum = rand()%(x2-x1 + 1) + x1;
            auto& comps = placement[r]; // Get components in the target row
            
            // Find the component at or just before the random x-coordinate
            auto it = lower_bound(comps.begin(), comps.end(), randNum,
                [](const Component* a, int value) {
                    return a->x < value;
                }
            );
            int idx = it-comps.begin() -1; // Start search from this index
            int step=0;

            bool goRight=true; // Spiral search flag

            // Spiral search (left and right) from the random starting index
            while(!found){
                step++;

                // --- Try to swap with another cell ---
                if(goRight) idx+=step;
                else idx-=step;
                
                if(idx>=comps.size()) break; // Out of bounds
                
                Component* cto = comps[idx]; // Get target component
                
                if(cto->x >x2 || cto->x<x1) break; // Target is outside optimal x-region
                if(cto->fixed) continue; // Can't swap with a fixed cell
                
                // Check orientation and capacity validity
                if( !ValidOrient(cto->c->sym, rows[rfrom].orient) || !ValidSwap(cfrom, cto, rfrom, r)) continue;

                // Prepare to calculate dHPWL
                vector<Component*> cList={cfrom, cto};
                vector<int> newX = {cto->x, cfrom->x}; 
                vector<int> newY = {cto->y, cfrom->y};

                vector<NetInfo> netInfos;
                long long dHPWL = dHPWLafterSwap(cList, newX, newY, netInfos);

                // acceptance check
                if(dHPWL>Threshold){

                    long long penalty=0;
                    // Calculate overlap penalty
                    if(cfrom->c->w > cto->c->w)
                        penalty = SwapPenalty(cfrom, cto, r, idx);
                    else
                        penalty = SwapPenalty(cto, cfrom, rfrom, idxfrom);

                    // Final check including penalty
                    if(dHPWL-penalty>Threshold){
                        // Commit the swap
                        for(auto& info : netInfos){ // Update net BBoxes
                            Net* n = info.net;
                            n->HPWL = info.HPWL;
                            n->minX = info.minX;
                            n->minY = info.minY;
                            n->maxX = info.maxX;
                            n->maxY = info.maxY;
                        }
                        UpdateOccupy(cfrom, cto, rfrom, r); // Update row/segment occupancy

                        // Swap coordinates and row indices
                        swap(cfrom->x, cto->x);
                        swap(cfrom->y, cto->y);

                        cfrom->ridx = r;
                        cto->ridx = rfrom;

                        // Update the placement vectors
                        placement[rfrom][idxfrom] = cto;
                        placement[r][idx] = cfrom;

                        found=true; // Mark as found
                        break; // Exit while loop
                    }
                }

                // --- Try to move into a space ---
                
                // Find the end x of the current cell (which is the start of a potential space)
                int spaceX = legalizeXback( getEndX(r,idx) , rows[r].oy);
                if(spaceX>x2) break; // Space is outside optimal x-region
                if(spaceX==getX(r,idx+1)) continue; // No space

                if(!ValidSwap(cfrom, spaceX, rfrom, r)) continue; // Check capacity

                // Prepare to calculate dHPWL
                vector<Component*> cList1={cfrom};
                vector<int> newX1 = {spaceX}; 
                vector<int> newY1 = {rows[r].oy};

                vector<NetInfo> netInfos1;
                dHPWL = dHPWLafterSwap(cList1, newX1, newY1, netInfos1);
                
                // SA acceptance check
                if(dHPWL>Threshold){

                    long long penalty=SwapPenalty(cfrom, r, idx); // Calculate overlap penalty
                    
                    if(dHPWL-penalty>Threshold){ // Final check
                        // Commit the move-to-space
                        for(auto& info : netInfos1){ // Update net BBoxes
                            Net* n = info.net;
                            n->HPWL = info.HPWL;
                            n->minX = info.minX;
                            n->minY = info.minY;
                            n->maxX = info.maxX;
                            n->maxY = info.maxY;
                        }

                        UpdateOccupy(cfrom, spaceX, rfrom, r); // Update occupancy

                        // Update component's coordinates and row
                        cfrom->x = spaceX;
                        cfrom->y = rows[r].oy;
                        cfrom->ridx = r;

                        if(r==rfrom && idx>=idxfrom){ // Handle intra-row move
                            placement[r].insert(placement[r].begin()+idx+1, cfrom);
                            placement[rfrom].erase(placement[rfrom].begin()+idxfrom);
                        }
                        else{ // Handle inter-row move
                            placement[rfrom].erase(placement[rfrom].begin()+idxfrom);
                            placement[r].insert(placement[r].begin()+idx+1, cfrom);
                        }
                        found=true; // Mark as found
                        break; // Exit while loop
                    }
                }

                // Alternate spiral search direction
                if(goRight) goRight=false;
                else goRight=true;
            }

            // --- Fallback: Try to place at the left edge of the optimal region ---
            if(!found && idx>0 && idx<comps.size()){
                if( comps[idx]->x < x1){ // If search ended to the left of optimal region

                    int spaceX = legalizeXback( getEndX(r,idx) , rows[r].oy);
                    if(spaceX<x1) spaceX=x1; // Snap to left edge of optimal region

                    if(spaceX!=getX(r,idx+1) && ValidSwap(cfrom, spaceX, rfrom, r)) {

                        vector<Component*> cList1={cfrom};
                        vector<int> newX1 = {spaceX}; 
                        vector<int> newY1 = {rows[r].oy};

                        vector<NetInfo> netInfos1;
                        long long dHPWL = dHPWLafterSwap(cList1, newX1, newY1, netInfos1);
                        
                        if(dHPWL>Threshold){
                            long long penalty=SwapPenalty(cfrom, r, idx);
                            
                            if(dHPWL-penalty>Threshold){
                                // Commit the move-to-space
                                for(auto& info : netInfos1){
                                    Net* n = info.net;
                                    n->HPWL = info.HPWL;
                                    n->minX = info.minX;
                                    n->minY = info.minY;
                                    n->maxX = info.maxX;
                                    n->maxY = info.maxY;
                                }
                                UpdateOccupy(cfrom, spaceX, rfrom, r);

                                cfrom->x = spaceX;
                                cfrom->y = rows[r].oy;
                                cfrom->ridx = r;

                                if(r==rfrom && idx>=idxfrom){
                                    placement[r].insert(placement[r].begin()+idx+1, cfrom);
                                    placement[rfrom].erase(placement[rfrom].begin()+idxfrom);
                                }
                                else{
                                    placement[rfrom].erase(placement[rfrom].begin()+idxfrom);
                                    placement[r].insert(placement[r].begin()+idx+1, cfrom);
                                }
                                found=true;
                                break; // Exit for loop
                            }
                        }
                    }
                }
            }
        }
    }
    legalizeRows(); // Legalize all rows after swaps
}

void FastDP::VerticalSwap(){
    // Iterate through all movable components
    for(auto* cfrom : components){
        if(cfrom->fixed) continue; // Skip fixed

        bool found = false;

        int x1, x2, y1, y2;
        // Find optimal region
        FindOptRegion(cfrom, x1, x2, y1, y2);

        // Skip if already in optimal y-range
        if(y1 < cfrom->y && cfrom->y < y2)
            continue;

        int rfrom = cfrom->ridx; // Get current row
        
        // Determine target row (one row up or down)
        int r = cfrom->y < y1 ? rfrom+1 : rfrom-1;
        if( r<0 || r>=rows.size() ) continue; // Out of bounds

        // If adjacent row is not compatible (site/orient)
        if(rows[r].site!=rows[rfrom].site || !ValidOrient(cfrom->c->sym, rows[r].orient)) {
            // Search a few rows further
            if(cfrom->y < y1){ // Search down
                for(int i=1; i<4; i++){
                    if(r+i>=rows.size()){
                        found=true; // Reached end, stop
                        break;
                    }
                    if(rows[r+i].site==rows[rfrom].site){
                        r=r+i; // Found compatible row
                        break;
                    }
                }
            }
            else{ // Search up
                for(int i=1; i<4; i++){
                    if(r-i<0){
                        found=true; // Reached end, stop
                        break;
                    }
                    if(rows[r-i].site==rows[rfrom].site){
                        r=r-i; // Found compatible row
                        break;
                    }
                }
            }
        }
        if(found) continue; // No compatible row found within 4 steps

        auto& comps0 = placement[rfrom]; // Source row components

        // Find index of component in source row
        auto it0 = std::find(comps0.begin(), comps0.end(), cfrom);
        int idxfrom = it0 - comps0.begin();
        
        // --- Find best move in the target row ---
        long long bestBenifit = 0;
        int bestIdx=-1;
        vector<NetInfo> bestNetInfos;
        bool isSpace=false;

        auto& comps = placement[r]; // Target row components
        if(comps.size()==0) continue; // Target row is empty

        // Find search window in target row based on component's x
        auto it = lower_bound(comps.begin(), comps.end(), cfrom->x,
            [](const Component* a, int value) {
                return a->x < value;
            }
        );
        int idx = it-comps.begin();
        int idxMax= idx<comps.size()-1 ? idx+1 : idx; // Search window of ~3 cells

        if(idx>0) idx--; // Widen window to the left

        idx--; // Start search from one left of window
        while(!found){
            idx++;
            if(idx>=idxMax) break; // End of search window
            Component* cto = comps[idx];
            if(cto->fixed) continue; // Skip fixed
            
            // Check orientation and capacity
            if( !ValidOrient(cto->c->sym, rows[rfrom].orient) || !ValidSwap(cfrom, cto, rfrom, r)) continue;

            // --- Check "swap cell" move ---
            vector<Component*> cList={cfrom, cto};
            vector<int> newX = {cto->x, cfrom->x}; 
            vector<int> newY = {cto->y, cfrom->y};

            vector<NetInfo> netInfos;
            long long dHPWL = dHPWLafterSwap(cList, newX, newY, netInfos);
            
            //  check
            if(dHPWL>0){
                long long penalty=0;
                if(cfrom->c->w > cto->c->w)
                    penalty = SwapPenalty(cfrom, cto, r, idx);
                else
                    penalty = SwapPenalty(cto, cfrom, rfrom, idxfrom);

                long long benefit=dHPWL-penalty;

                // Find the move with the best benefit
                if(benefit>0 && benefit>bestBenifit){
                    bestBenifit=benefit;
                    bestIdx=idx;
                    bestNetInfos=netInfos;
                    isSpace=false;
                }
            }

            // --- Check "move to space" move ---
            int spaceX =  legalizeXback( getEndX(r,idx) , rows[r].oy);
            if(spaceX>x2) break;
            if(spaceX==getX(r,idx+1)) continue; // No space

            if(!ValidSwap(cfrom, spaceX, rfrom, r)) continue;

            vector<Component*> cList1={cfrom};
            vector<int> newX1 = {spaceX}; 
            vector<int> newY1 = {rows[r].oy};

            vector<NetInfo> netInfos1;
            dHPWL = dHPWLafterSwap(cList1, newX1, newY1, netInfos1);
            
            if(dHPWL>0){
                long long penalty=SwapPenalty(cfrom, r, idx);
                
                long long benefit=dHPWL-penalty;
                if(benefit>0 && benefit>bestBenifit){
                    bestBenifit=benefit;
                    bestIdx=idx;
                    bestNetInfos=netInfos1;
                    isSpace=true;
                }
            }
        }
        // Commit the best move found in the window
        if(bestIdx!=-1){

            if(!isSpace){ // Best move was a cell swap
                Component* cto = comps[bestIdx];

                vector<Component*> cList={cfrom, cto};
                vector<int> newX = {cto->x, cfrom->x}; 
                vector<int> newY = {cto->y, cfrom->y};

                // Update net BBoxes
                for(auto& info : bestNetInfos){
                    Net* n = info.net;
                    n->HPWL = info.HPWL;
                    n->minX = info.minX;
                    n->minY = info.minY;
                    n->maxX = info.maxX;
                    n->maxY = info.maxY;
                }

                UpdateOccupy(cfrom, cto, rfrom, r); // Update occupancy

                // Swap coordinates and row indices
                swap(cfrom->x, cto->x);
                swap(cfrom->y, cto->y);

                cfrom->ridx = r;
                cto->ridx = rfrom;

                // Update placement vectors
                placement[rfrom][idxfrom] = cto;
                placement[r][bestIdx] = cfrom;

                found=true;
            }
            else{ // Best move was a move-to-space
                int spaceX = legalizeXback( getEndX(r,bestIdx) , rows[r].oy);

                vector<Component*> cList={cfrom};
                vector<int> newX = {spaceX}; 
                vector<int> newY = {rows[r].oy};

                // Update net BBoxes
                for(auto& info : bestNetInfos){
                    Net* n = info.net;
                    n->HPWL = info.HPWL;
                    n->minX = info.minX;
                    n->minY = info.minY;
                    n->maxX = info.maxX;
                    n->maxY = info.maxY;
                }
                UpdateOccupy(cfrom, spaceX, rfrom, r); // Update occupancy

                // Update component's coordinates and row
                cfrom->x = spaceX;
                cfrom->y = rows[r].oy;
                cfrom->ridx = r;

                placement[rfrom].erase(placement[rfrom].begin()+idxfrom);
                placement[r].insert(placement[r].begin()+bestIdx+1, cfrom);

                found=true;
            }
        }
    }
    legalizeRows(); // Legalize all rows after swaps
}