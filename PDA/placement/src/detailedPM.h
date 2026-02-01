#include <iostream>
#include <cmath>
#include <string>
#include <vector>
#include <unordered_map>
#include <algorithm>
#include <random>
#include <unordered_set>
#include <chrono>
#include <omp.h>
#include <climits>

using namespace std;


enum SymmetryFlag {
    SYM_X   = 1 << 0,  // 001
    SYM_Y   = 1 << 1,  // 010
    SYM_R90 = 1 << 2   // 100
};

// Represents a cell macro from LEF
struct Cell {
    int ctype = 0; // Cell type (e.g., 0 for CORE, 1 for BLOCK)
    int w = 0;     // Width
    int h = 0;     // Height
    int sym = 0;   // Symmetry flags
    vector<float>* site; // Pointer to site info
};

// Represents a placeable segment within a row
struct Segment {
    int x_start;    // segment start coordinate (inclusive)
    int x_end;      // segment end coordinate (exclusive)
    int capacity;   // Available capacity (in site units)
    int occupy;     // Occupied capacity (in site units)
};

// Represents a placement row
struct Row {
    vector<float>* site; // Pointer to site info
    int ox=0;           // Origin x
    int oy=0;           // Origin y
    string orient;      // Row orientation (e.g., N, FS)
    int nw=0;           // Number of sites (width)
    int nh=1;           // Number of sites (height)
    int dx=0;           // Site width
    int dy=0;           // Site height (row height)
    int occupy=0;       // Total occupied sites in row
    vector<Segment> segments; // Vector of placeable segments
};

// Represents a PIN location
struct Pin{
    int x=0;
    int y=0;
};

struct Net; // Forward declaration for Component

// Represents a component instance
struct Component {
    Cell* c;          // Pointer to the cell macro
    int x=0;          // x-coordinate
    int y=0;          // y-coordinate
    string orient;    // Component orientation
    vector<Net*> net; // Pointers to nets this component is part of
    bool fixed;       // Is the component fixed?

    string name;      // Instance name

    int ridx=-1;      // Row index
};

// Represents a net (connection)
struct Net {
    vector<Component*> c={}; // Pointers to components in this net
    vector<string> o={};     // Pin orientations (unused?)
    int minX=0;              // Bounding box min X
    int minY=0;              // Bounding box min Y
    int maxX=0;              // Bounding box max X
    int maxY=0;              // Bounding box max Y
    int pminX=-1;            // Bounding box min X (PINs only)
    int pminY=-1;            // Bounding box min Y (PINs only)
    int pmaxX=-1;            // Bounding box max X (PINs only)
    int pmaxY=-1;            // Bounding box max Y (PINs only)
    int HPWL=0;              // Half-Perimeter Wirelength

    string name;             // Net name
};

// Struct to hold temporary net info during swap calculations
struct NetInfo {
    Net* net;
    int minX=0;
    int minY=0;
    int maxX=0;
    int maxY=0;
    int HPWL=0;
};

// Struct for clustering during legalization
struct Cluster{
    int idx=0; // Start index in placement vector
    int x=0;   // Legalized x-coordinate
    int e=1;   // Number of elements in cluster
    int w=0;   // Total width of cluster
    vector<int> bounds; // (Unused?)
};

// Main placer class
class FastDP{

public:
    // Constructor
    FastDP(unordered_map<string,vector<float>>&,
           unordered_map<string, Cell>&,
           vector<Row>&,
           unordered_map<string, Component>&,
           vector<Net>&);

    // Main optimization function
    void Optimize(chrono::time_point<std::chrono::steady_clock>);

private:

    int numRow; // Total number of rows

    // References to data from main
    unordered_map<string,vector<float>>& CoreSite;
    unordered_map<string, Cell>& cells;
    unordered_map<string, Component>& componentMap;
    vector<Component*> components; // Vector of pointers to movable components
    vector<Net>& nets;
    vector<Row>& rows;
    
    // Internal data structures
    vector<vector<Component*>> placement;    // 2D vector [row_idx][component_ptr]
    unordered_map<int,vector<int>> y2RowIdx; // Map y-coordinate to row indices

    float Threshold = -2000;        // threshold parameter

    // Core optimization functions
    void SingleSegmentClustering(); // (Unused?)
    void GlobalSwap();              // Perform long-range swaps
    void VerticalSwap();            // Perform swaps with adjacent rows
    void ReOrdering();              // Re-order cells within rows

    // Utility functions for row/cell positions
    int getX(int, int);             // Get x of component at placement[r][idx]
    int getEndX(int, int);          // Get end x of component at placement[r][idx]
    int legalizeX(int,int);         // Legalize x to nearest site boundary (forward)
    int legalizeXback(int,int);     // Legalize x to nearest site boundary (backward)
    int findSegmentIndex(int, int); // Find segment index for a given row and x
    int findRowIndex(int, int);     // Find row index for a given x, y

    // Functions for calculating optimal regions and penalties
    void FindOptXs(Component*, int&, int&);                  // Find optimal x-range
    void FindOptRegion(Component*, int&, int&, int&, int&);  // Find optimal region (x,y)
    long long SwapPenalty(Component*, Component*, int, int); // Penalty for swapping cells
    long long SwapPenalty(Component*, int, int);             // Penalty for moving cell to space
    long long dHPWLafterSwap(vector<Component*>& , vector<int>& , vector<int>& , vector<NetInfo>&); // Calculate HPWL change

    // Legalization functions
    void legalizeRows();                                    // Legalize all rows
    void legalizeBoundary(vector<Cluster>&, int, int, int); // Legalize a cluster within boundaries

    // Validation functions
    bool ValidOrient(int, string);                          // Check if orientation is valid
    bool ValidSwap(Component*, Component*, int, int);       // Check if cell swap is valid (capacity)
    bool ValidSwap(Component*, int, int, int);              // Check if move-to-space is valid (capacity)
    void UpdateOccupy(Component*, Component*, int, int);    // Update occupancy after cell swap
    void UpdateOccupy(Component*, int, int, int);           // Update occupancy after move-to-space
};