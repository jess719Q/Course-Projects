#include "FM.h"

FMEngine::FMEngine(vector<cell*>& cm, vector<vector<cell*>>& nl, int p, chrono::time_point<std::chrono::steady_clock> st)
    : cellList(cm), netList(nl), partitions(p), startTime(st) {
}

void FMEngine::MultiWayFM(vector<int> groups, int iter){
    // Run one pass of the FM algorithm for the given set of groups. 
    // 
    // Parameters:
    //   groups: the indices of groups allowed to exchange cells in this pass.
    //   iter  : the current stage index in the hierarchical multi-way partitioning.
    //           It represents how many groups have already been generated 
    //           (e.g., iter = 4 means we are refining the four-way partition stage).
    //           Note that 'iter' does NOT necessarily equal to groups.size();
    //           it only indicates the expansion depth of the current partitioning hierarchy.
    //           It is used to control the allowed group-size balance range.

    int minSize=totSize*pow(0.5,log2(iter))*pow(0.9,float(iter)/partitions);
    int maxSize=totSize*pow(0.5,log2(iter))*pow(1.1,float(iter)/partitions);

    // move cell
    int cutSize0 = cutSize;
    vector<record> movRecord;
    movRecord.reserve(cellList.size());
    while(true){
        int moveGain = moveCell(groups, movRecord, minSize, maxSize);
        cutSize -= moveGain-maxP;
        movRecord[movRecord.size() - 1].cutsize=cutSize;
        if(moveGain<0 || cutSize>cutSize0*10)
            break;
    }


    // find min cut size
    int minCutsize=cutSize0;
    int minIdx=-1;

    for(int i=0; i<movRecord.size(); i++){
        if(movRecord[i].cutsize<minCutsize){
            minCutsize=movRecord[i].cutsize;
            minIdx=i;
        }
        else if(movRecord[i].cutsize==minCutsize){
            if(movRecord[i].sizeDiff<movRecord[minIdx].sizeDiff){
                minCutsize=movRecord[i].cutsize;
                minIdx=i;
            }
        }
    }
    
    // recover cells from records
    for(int i=movRecord.size()-1; i>minIdx; i--){
        cell* re = movRecord[i].c;
        re->group=movRecord[i].fromG;
        groupSize[ movRecord[i].fromG ] += re->size;
        groupSize[ movRecord[i].toG ] -= re->size;
        for(int n: re->nets){
            netGroupNum[n][movRecord[i].fromG]++;
            netGroupNum[n][movRecord[i].toG]--;
        }
    }
    cutSize=minCutsize;
}

void FMEngine::TwoWayInitFM(vector<int> groups, int iter){
    // Initialize two groups by moving cells from the first group (g1)
    // to the second group (g2) one by one based on gain values.
    //
    // Parameters:
    //   groups: indices of the two groups involved in this initialization.
    //           The first group initially contains all cells, while
    //           the second group starts empty.
    //
    //   iter  : the current stage index in the hierarchical multi-way
    //           partitioning process. It indicates how many groups
    //           will be generated in total at this stage
    //           (e.g., iter = 4 corresponds to the four-way partition stage).
    //           It is used to control the allowed group-size balance range.

    int g0=groups[0], g1=groups[1];

    int minSize=totSize*pow(0.5,log2(iter))*pow(0.9, float(iter)/partitions);
    int maxSize=totSize*pow(0.5,log2(iter))*pow(1.1, float(iter)/partitions);

    
    // move cell
    while(true){
        int moveGain = moveCellforInitialize(g0,g1);
        cutSize -= moveGain-maxP;
        if(groupSize[g0]<=maxSize && groupSize[g1]>=minSize) break;
    }

}

void FMEngine::FiducciaMattheyses(){
    // Run the full Fiduccia–Mattheyses (FM) multi-level partitioning driver.
    //
    // Outputs:
    //   - Final cut size is stored in `cutSize`.
    //   - Final group assignment is stored in `groupsAfter`.
    //   - Also updates the following internal states:
    //       * `cellList`       – cell objects with updated group numbers
    //       * `groupSize`      – current cell count (or total size) per group
    //       * `netGroupNum`    – number of groups each net spans

    // initialize variables
    cutSize=0; 
    totSize=0; 
    maxP=0;

    groupSize.resize(partitions, 0);
    netGroupNum.resize(netList.size(),vector<int>(partitions,0));
    groupBucket.resize(partitions, vector<vector<Node*>>(partitions, vector<Node*>(2*maxP+1, nullptr)));
    bucketHead.resize(partitions, vector<int>(partitions, -1));

    for(cell* c : cellList) {
        totSize += c->size;
        if(c->nidx.size()>maxP)
            maxP=c->nidx.size();
    }

    groupSize[0]=totSize;

    for (int net=0; net<netList.size(); net++)
        netGroupNum[net][0]=netList[net].size();

    groupBucket.resize(partitions);
    for (int i=0; i<partitions; i++){
        groupBucket[i].resize(partitions);
        for(int j=0; j<partitions; j++)
            if(i!=j)
                groupBucket[i][j].resize(2 * maxP + 1, nullptr);
    }
    
    //===================================================================
    
    // initialize buckets

    for(cell* c: cellList){
        c->group=0;


        int gain=-c->nets.size();

        int gidx=gain+maxP;

        Node* tmp=c->pos[0];
        c->gidx[0]=-1;
        tmp->c=c;
        tmp->front=nullptr;
        tmp->next=nullptr;

        for(int j=1; j<partitions; j++){
            Node* tmp = c->pos[j];
            tmp->c = c;
            tmp->front=nullptr;
            c->gidx[j]=gidx;
            if(groupBucket[0][j][gidx]){
                tmp->next = groupBucket[0][j][gidx];
                groupBucket[0][j][gidx]->front=tmp;
            }
            else
                tmp->next=nullptr;
            groupBucket[0][j][gidx]=tmp;
        }
    }

    for(int g=2*maxP; g>=0; g--){
        if(groupBucket[0][1][g]){
            for(int i=1; i<partitions; i++)
                bucketHead[0][i]=g;
            break;
        }
    }
    //===================================================================
    
    // initial partitioning
    for(int p=1; p<partitions; p*=2){
        if(p>1)
            InitializeGroupBucket();

        for(int g=0; g<p; g++){
            vector<int> groups(2);
            groups[0]=g;
            groups[1]=g+p;
            TwoWayInitFM(groups,p*2);    


            int cutSizeLast=cutSize;
            while(true){
                InitializeGroupBucket();
                MultiWayFM(groups,p*2);

                if(cutSizeLast-cutSize<=netList.size()*0.0001)
                    break;
                cutSizeLast=cutSize;
            }
        }
    }
    
    //===================================================================
    
    // partitioning
    int cutSizeLast=cutSize;
    vector<int> groups(partitions);
    for(int i=0; i<partitions; i++) groups[i]=i;
    while(true){
        
        // random choose two groups to run partition
        if(partitions>2){
          random_device rd;
          mt19937 gen(rd());
          shuffle(groups.begin(), groups.end(), gen);
          vector<int> groups2(2);
          groups2[0]=groups[0];
          groups2[1]=groups[1];
          InitializeGroupBucket();
          MultiWayFM(groups2,partitions);
        }

        
        // partitioning all groups
        InitializeGroupBucket();
        MultiWayFM(groups,partitions);

        if(cutSize==cutSizeLast){
            break;
        }
        cutSizeLast=cutSize;
        
        // if run too long
        auto now = chrono::steady_clock::now();
        auto elapsed = chrono::duration_cast<chrono::seconds>(now - startTime).count();
        if (elapsed >= 50)
            break;
    }

    groupsAfter.resize(partitions);
    for(cell* c: cellList)
        groupsAfter[c->group].push_back(c->name);
}


void FMEngine::removeFromBucket(cell* cell2mov, int fromGroup, int toGroup){
    int gidx = cell2mov->gidx[toGroup];
    Node* n = cell2mov->pos[toGroup];
    if (!n) {
        return;
    }

    if(n->front) n->front->next=n->next;
    if(n->next) n->next->front=n->front;
    if(groupBucket[fromGroup][toGroup][gidx]==n)
        groupBucket[fromGroup][toGroup][gidx]=n->next;

    if(bucketHead[fromGroup][toGroup]==gidx && !groupBucket[fromGroup][toGroup][gidx]){
        int p=gidx;
        while(p>=0 && !groupBucket[fromGroup][toGroup][--p]);
        bucketHead[fromGroup][toGroup]=p;
    }
}

void FMEngine::appendToBucket(cell* cell2mov, int fromGroup, int toGroup){
    Node* n = cell2mov->pos[toGroup];
    int gidx=cell2mov->gidx[toGroup];
    n->front=nullptr;
    if(groupBucket[fromGroup][toGroup][gidx]){
        n->next = groupBucket[fromGroup][toGroup][gidx];
        groupBucket[fromGroup][toGroup][gidx]->front=n;
    }
    else
        n->next=nullptr;
    groupBucket[fromGroup][toGroup][gidx]=n;

    if(bucketHead[fromGroup][toGroup]<gidx)
        bucketHead[fromGroup][toGroup]=gidx;
    if(bucketHead[fromGroup][toGroup] == -1)
        bucketHead[fromGroup][toGroup] = gidx;
}

void FMEngine::updateBucket(cell* cell2mov, int fromGroup, int toGroup, int gchange){
    removeFromBucket(cell2mov, fromGroup, toGroup);
    cell2mov->gidx[toGroup]+=gchange;
    appendToBucket(cell2mov, fromGroup, toGroup);
    return;
}

void FMEngine::InitializeGroupBucket(){
    // Initialize the multi-way gain bucket structure for all group pairs.
    // It clears previous bucket contents and resets the bucket of each 
    // (fromGroup, toGroup) combination.

    for (auto& v1 : groupBucket)
        for (auto& v2 : v1)
            fill(v2.begin(), v2.end(), nullptr);
    

    for(cell* c: cellList){

        int gain=0;
        int selfGroup=c->group;

        for(int j=0; j<partitions; j++){
            if(j==selfGroup) {
                c->gidx[j]=-1;
                continue;
            }

            gain=0;
            for (int net: c->nets) {
                if (netGroupNum[net][selfGroup] == netList[net].size())
                    gain--;
                else if (netGroupNum[net][selfGroup]==1 && netGroupNum[net][j]+1==netList[net].size())
                    gain++;
            }
            
            int gidx=gain+maxP;

            Node* tmp = c->pos[j];
            tmp->front=nullptr;
            c->gidx[j]=gidx;
            if(groupBucket[selfGroup][j][gidx]){
                tmp->next = groupBucket[selfGroup][j][gidx];
                groupBucket[selfGroup][j][gidx]->front=tmp;
            }
            else
                tmp->next=nullptr;
            groupBucket[selfGroup][j][gidx]=tmp;
        }
    }

    for(int i=0; i<partitions; i++){
        for(int j=0; j<partitions; j++){
            if(i==j) continue;
            for(int g=2*maxP; g>=0; g--){
                if(groupBucket[i][j][g]){
                    bucketHead[i][j]=g;
                    break;
                }
            }
        }
    }
}

void FMEngine::updateGain(cell* cell2mov, int fromGroup, int toGroup){
    // Update the gain values of the affected cells connected through the same nets.
    
    for(int net: cell2mov->nets){
        netGroupNum[net][fromGroup]--;
        netGroupNum[net][toGroup]++;
        
        if(netGroupNum[net][toGroup]==1 && netGroupNum[net][fromGroup]+1==netList[net].size())
            for(auto* cel: netList[net])
                if(cel->gidx[toGroup]>=0 && cel!=cell2mov)
                    for(int g=0; g<partitions; g++)
                        if(g!=fromGroup)
                            updateBucket(cel, fromGroup, g, 1);

        if(netGroupNum[net][toGroup]==2 && netGroupNum[net][fromGroup]+2==netList[net].size())
            for(auto* cel: netList[net])
                if(cel->group==toGroup && cel->gidx[fromGroup]>=0 && cel!=cell2mov)
                    updateBucket(cel, toGroup, fromGroup, -1);


        if(netGroupNum[net][fromGroup]==0 && netGroupNum[net][toGroup]==netList[net].size())
            for(auto* cel: netList[net])
                if(cel->gidx[fromGroup]>=0 && cel!=cell2mov)
                    for(int g=0; g<partitions; g++)
                        if(g!=toGroup)
                            updateBucket(cel, toGroup, g, -1);

        if(netGroupNum[net][fromGroup]==1 && netGroupNum[net][toGroup]+1==netList[net].size())
            for(auto* cel: netList[net])
                if(cel->group==fromGroup && cel->gidx[toGroup]>=0 && cel!=cell2mov)
                    updateBucket(cel, fromGroup, toGroup, 1);


        if(netGroupNum[net][fromGroup]==0 && netGroupNum[net][toGroup]+1==netList[net].size())
            for(auto* cel: netList[net])
                if(cel->group!=fromGroup && cel->gidx[toGroup]>=0 && cel!=cell2mov)
                    updateBucket(cel, cel->group, toGroup, 1);

        if(netGroupNum[net][toGroup]==1 && netGroupNum[net][fromGroup]+2==netList[net].size())
            for(auto* cel: netList[net])
                if(cel->group!=toGroup && cel->gidx[fromGroup]>=0 && cel!=cell2mov)
                    updateBucket(cel, cel->group, fromGroup, -1);
                    
    }
}

int FMEngine::moveCell(vector<int> groups, vector<record>& movRecord, int minSize, int maxSize){
    int gidx=-1;
    int fromGroup=-1, toGroup=-1;
    bool canMove=false;
    cell* cell2mov = nullptr;

    vector<vector<int>> tempBucketHead = bucketHead;
    
    // find the cell to move
    while(!canMove){
        for(int i: groups){
            for(int j: groups){
                if(i==j) continue;
                else if(tempBucketHead[i][j]>gidx){
                    gidx=tempBucketHead[i][j];   fromGroup=i;    toGroup=j;
                }
                else if(tempBucketHead[i][j]==gidx){
                    if(groupSize[i]<groupSize[fromGroup] && groupSize[j]>groupSize[toGroup]){
                        continue;
                    }
                    else if(groupSize[i]>groupSize[fromGroup] && groupSize[j]<groupSize[toGroup]){
                        gidx=tempBucketHead[i][j];   fromGroup=i;    toGroup=j;
                    }
                    else if(groupSize[i]-groupSize[j] > groupSize[fromGroup]-groupSize[toGroup]){
                        gidx=tempBucketHead[i][j];   fromGroup=i;    toGroup=j;
                    }


                }
            }
        }

        if(gidx==-1) return -1; // no more cell can move
        Node* n = groupBucket[fromGroup][toGroup][gidx];  
        
        for(int t=0; t<2; t++) {
            if(groupSize[fromGroup]-n->c->size >= minSize && groupSize[toGroup]+n->c->size <= maxSize){
                canMove=true;    // found
                cell2mov=n->c;
                for(int k=0; k<partitions; k++)
                    if(k!=fromGroup){
                        if(cell2mov->gidx[k]>=0)
                            removeFromBucket(cell2mov, fromGroup, k);   // remove from buckets
                        cell2mov->gidx[k]=-1;
                    }
                
                break;
            }
            n=n->next;
            if(!n) break;
        }

        if(!canMove){
            tempBucketHead[fromGroup][toGroup]=-1;
            gidx=-1;
            fromGroup=-1;
            toGroup=-1;
        }
    }

    // move cell
    groupSize[fromGroup]-=cell2mov->size;
    groupSize[toGroup]+=cell2mov->size;
    cell2mov->group=toGroup;
    
    updateGain(cell2mov, fromGroup, toGroup);

    int maxG=0, minG=2147483647;
    for(int g: groups){
        if(groupSize[g]>maxG)
            maxG=groupSize[g];
        if(groupSize[g]<minG)
            minG=groupSize[g];
    }
    movRecord.push_back({cell2mov,0,fromGroup,toGroup, maxG-minG});
    return gidx;
}

void FMEngine::move2anotherGroup(cell* cell2mov, int fromGroup, int toGroup) {
    // Move the given cell from 'fromGroup' to 'toGroup'.
    // The cell is not locked after moving; it remains eligible for future moves.


    groupSize[fromGroup] -= cell2mov->size;
    groupSize[toGroup]   += cell2mov->size;
    cell2mov->group = toGroup;

    updateGain(cell2mov, fromGroup, toGroup);

    for (int g = 0; g < partitions; g++) {


        // remove old bucket node
        if(cell2mov->gidx[g]>=0){
            removeFromBucket(cell2mov, fromGroup, g);
            cell2mov->gidx[g] = -1;
        }

        if(g == toGroup) continue;


        // calculte new gain
        int gain = 0;
        for (int net: cell2mov->nets) {
            if (netGroupNum[net][toGroup] == netList[net].size())
                gain--;
            else if (netGroupNum[net][toGroup]==1 && netGroupNum[net][g]+1==netList[net].size())
                gain++;
        }

        int gidx = gain + maxP;
        
        // append to new bucket
        Node* tmp = cell2mov->pos[g];
        tmp->front = nullptr;
        if (groupBucket[toGroup][g][gidx]) {
            tmp->next = groupBucket[toGroup][g][gidx];
            groupBucket[toGroup][g][gidx]->front = tmp;
        } else {
            tmp->next = nullptr;
        }

        groupBucket[toGroup][g][gidx] = tmp;
        cell2mov->pos[g] = tmp;
        cell2mov->gidx[g] = gidx;
    }
}


int FMEngine::moveCellforInitialize(int fromGroup, int toGroup){

    int gidx=bucketHead[fromGroup][toGroup];
    cell* cell2mov = nullptr;
    Node* n = groupBucket[fromGroup][toGroup][gidx];
    
    cell2mov=n->c;
    move2anotherGroup(cell2mov, fromGroup, toGroup);

    return gidx;
}

