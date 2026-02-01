#include "detailedPM.h"

#include <fstream>
#include <sstream>

void readLEF(const string& inputLEF, int& units,
             unordered_map<string,vector<float>>& CoreSite,
             unordered_map<string, Cell>& cells){
    // read LEF
    ifstream finLEF(inputLEF); // Open LEF file
    if (!finLEF) {
        cerr << "Cannot open file.\n";
        return;
    }

    string line;

    while (getline(finLEF, line)) { // Read line by line

        if (line.find("UNITS") != string::npos) { // Find UNITS section
            string tmp;

            getline(finLEF, line);
            stringstream ss(line);
            ss >> tmp >>  tmp >> units >> tmp;   //DATABASE MICRONS 2000 ;
            getline(finLEF, line);
        }
        else if (line.find("SITE") != string::npos) { // Find SITE section
            string siteName;
            float w,h;

            string tmp;
            stringstream ss(line);
            ss >> tmp >> siteName; // Get site name

            getline(finLEF, line);
            getline(finLEF, line);
            ss.clear(); ss.str(line);
            ss >> tmp >> w >> tmp >> h >> tmp;  //SIZE 0.1 BY 1.2 ;
            
            getline(finLEF, line);

            CoreSite[siteName]={w*units,h*units}; // Store site dimensions
        }
        else if (line.find("MACRO") != string::npos) { // Find MACRO section
            string macroName, ctype, site;
            float w,h;
            int sym = 0; // Symmetry flag

            string tmp;
            stringstream ss(line);
            ss >> tmp >> macroName; // Get macro name
            
            while(getline(finLEF, line)){ // Loop inside MACRO definition
                if (line.find("CLASS") != string::npos){
                    ss.clear(); ss.str(line);      
                    ss >> tmp >> ctype >> tmp;         //CLASS CORE ;
                }
                else if (line.find("SIZE") != string::npos){
                    ss.clear(); ss.str(line);
                    ss >> tmp >> w >> tmp >>h >> tmp;  //SIZE 1.400000 BY 1.200000 ;
                }
                else if (line.find("SYMMETRY") != string::npos) {
                    // Set symmetry bits
                    if (line.find("X") != string::npos)   sym |= SYM_X;
                    if (line.find("Y") != string::npos)   sym |= SYM_Y;
                    if (line.find("R90") != string::npos) sym |= SYM_R90;
                }
                else if (line.find("SITE") != string::npos){
                    ss.clear(); ss.str(line);
                    ss >> tmp >> site;  //SITE CoreSite ; 
                }

                if(line.find("END "+macroName) != string::npos) break; // End of MACRO
            }
            
            
            // Store cell info based on its class
            if(ctype=="CORE")
                cells[macroName]={0,(int)(w*units),(int)h*units,sym,&CoreSite[site]};
            else 
                cells[macroName]={1,(int)(w*units),(int)h*units,sym,&CoreSite[site]};
        }

    }
    finLEF.close();
}

void readDEF(const string& inputDEF,
             unordered_map<string, Cell>& cells,
             unordered_map<string,vector<float>>& CoreSite,
             vector<Row>& rows,
             unordered_map<string, Component>& components,
             vector<Net>& nets){
    
    ifstream finDEF(inputDEF); // Open DEF file
    if (!finDEF) {
        cerr << "Cannot open file.\n";
        return ;
    }

    unordered_map<string, Pin> pins; // Map to store PIN locations
    string line;

    while (getline(finDEF, line)){ // Read line by line
        if (line.find("ROW ") != string::npos) { // Find ROW section

            string site;
            int ox,oy;
            string orient;
            int nw,nh,dx,dy;
            string tmp;

            stringstream ss(line);
            //ROW CORE_ROW_0 CoreSite 2000 2000 FS DO 1464 BY 1 STEP 200 0 
            ss>>tmp>>tmp>>site>>ox>>oy>>orient>>tmp>>nw>>tmp>>nh>>tmp>>dx>>dy;
            rows.push_back({&CoreSite[site] , ox,oy,orient,nw,nh,dx,dy,0,{}}); // Add new row
        }

        else if (line.find("COMPONENTS ") != string::npos){ // Find COMPONENTS section
            stringstream ss(line);
            string tmp;
            int num;
            ss>>tmp>>num>>tmp; // Get number of components
            components.reserve(num); // Pre-allocate memory

            while (getline(finDEF, line)){
                if (line.find("END COMPONENTS") != string::npos) break; // End of COMPONENTS

                else if(line.find("-") != string::npos){ // Component entry
                    string tmp, name, cname, orient, mtype;
                    int x, y;

                    ss.clear(); ss.str(line);
                    ss>>tmp>>name>>cname; // Get component instance name and macro name

                    if(line.find("PLACED") != string::npos || line.find("FIXED") != string::npos){
                        // Parse PLACED or FIXED status
                        ss>>tmp;
                        while( tmp!="PLACED" && tmp!="FIXED" ) ss>>tmp;
                        mtype = tmp;
                        ss>>tmp>>x>>y>>tmp>>orient;
                        bool fixed = mtype=="FIXED" ? true:false;
                        components[name]={ &cells[cname],x,y,orient,{},fixed,name}; // Store component
                    }
                    else{
                        // Handle multi-line component definitions
                        while (getline(finDEF, line)){
                            if(line.find("PLACED") != string::npos || line.find("FIXED") != string::npos){
                                ss.clear(); ss.str(line);
                                ss>>tmp;
                                while( tmp!="PLACED" && tmp!="FIXED" ) ss>>tmp;
                                mtype = tmp;
                                ss>>tmp>>x>>y>>tmp>>orient>>tmp;
                                bool fixed = mtype=="FIXED" ? true:false;
                                components[name]={ &cells[cname],x,y,orient,{},fixed,name, -1};
                                if(y==0){
                                    cout<<"!!"<<line<<endl; // Debug print
                                }
                                break;
                            }
                        }
                    }


                }


            }
        }

        else if (line.find("PINS ") != string::npos){ // Find PINS section
            stringstream ss(line);
            string tmp;
            int num;
            ss>>tmp>>num>>tmp; // Get number of pins
            pins.reserve(num);

            while (getline(finDEF, line)){
                if (line.find("END PINS") != string::npos) break; // End of PINS

                string name;
                int x, y;

                ss.clear(); ss.str(line);
                ss>>tmp>>name; // Get pin name

                getline(finDEF, line); // Skip layer info

                getline(finDEF, line); // Get pin location
                ss.clear(); ss.str(line);
                ss>>tmp>>tmp>>tmp>>x>>y;

                pins[name]={x,y}; // Store pin
            }
        }

        else if (line.find("SPECIALNETS ") != string::npos) continue; // Skip SPECIALNETS
        
        else if (line.find("NETS ") != string::npos){ // Find NETS section
            stringstream ss(line);
            string tmp;
            int num;
            ss>>tmp>>num>>tmp; // Get number of nets
            nets.resize(num);

            int nidx=-1; // Net index

            while (getline(finDEF, line)){
                if (line.find("END NETS") != string::npos) break; // End of NETS
                else if (line.find("-") != string::npos) { // New net definition
                    nidx++;

                    string name;
                    ss.clear(); ss.str(line);
                    ss>>tmp>>name;
                    nets[nidx].name=name; // Store net name
                }
                
                else if (line.find("(") != string::npos){ // Net connections
                    ss.clear(); ss.str(line);

                    string token, name, orient;
                    while (ss >> token && token != ";") { // Parse all connections in line
                        if (token == "(") {
                            ss>>name>>orient;
                            if(name=="PIN") { // Connection is to a PIN
                                Pin p=pins[orient];
                                // Update net bounding box for PINs
                                if(nets[nidx].pminX==-1){ // First PIN in this net
                                    nets[nidx].pminX=p.x;
                                    nets[nidx].pmaxX=p.x;
                                    nets[nidx].pminY=p.y;
                                    nets[nidx].pmaxY=p.y;
                                }
                                else{ // Update existing BBox
                                    if     (p.x<nets[nidx].pminX) nets[nidx].pminX=p.x;
                                    else if(p.x>nets[nidx].pmaxX) nets[nidx].pmaxX=p.x;
                                    if     (p.y<nets[nidx].pminY) nets[nidx].pminY=p.y;
                                    else if(p.y>nets[nidx].pmaxY) nets[nidx].pmaxY=p.y;
                                }

                            }
                            else{ // Connection is to a component
                                nets[nidx].c.push_back(&components[name]);
                                nets[nidx].o.push_back(orient);
                                components[name].net.push_back(&nets[nidx]); // Cross-reference
                            }

                        }
                    }
                }
                else if (line.find(";") != string::npos){ // End of net entry
                    // Calculate initial HPWL
                    if(nets[nidx].pminX==-1){ // No PINs, use first component
                        nets[nidx].minX=nets[nidx].c[0]->x;
                        nets[nidx].minY=nets[nidx].c[0]->y;
                        nets[nidx].maxX=nets[nidx].c[0]->x;
                        nets[nidx].maxY=nets[nidx].c[0]->y;
                    }else{ // Has PINs, use PIN BBox
                        nets[nidx].minX=nets[nidx].pminX;
                        nets[nidx].minY=nets[nidx].pminY;
                        nets[nidx].maxX=nets[nidx].pmaxX;
                        nets[nidx].maxY=nets[nidx].pmaxY;
                    }

                    // Expand BBox with component locations
                    for(auto& c: nets[nidx].c){

                        if(c->x < nets[nidx].minX)
                            nets[nidx].minX=c->x;
                        else if(c->x > nets[nidx].maxX)
                            nets[nidx].maxX=c->x;
                        if(c->y < nets[nidx].minY)
                            nets[nidx].minY=c->y;
                        else if(c->y > nets[nidx].maxY)
                            nets[nidx].maxY=c->y;
                    }
                    // Calculate Half-Perimeter Wirelength (HPWL)
                    nets[nidx].HPWL = (nets[nidx].maxX-nets[nidx].minX) + (nets[nidx].maxY-nets[nidx].minY);
                }
            }
        }
    }

    finDEF.close();
}

void writeDEF(const string& inputDEF, const string& outputDEF,
              unordered_map<string, Component> components){

    ifstream finDEF(inputDEF); // Open original DEF
    if (!finDEF) {
        cerr << "Cannot open file.\n";
        return;
    }

    ofstream foutDEF(outputDEF); // Open output DEF
    if (!foutDEF) {
        cerr << "Cannot open file o.\n";
        return;
    }

    string line;

    while (getline(finDEF, line)){
        if (line.find("COMPONENTS") != string::npos){ // Find COMPONENTS section
            foutDEF << line << "\n"; // Write header

            string name, tmp;

            while (getline(finDEF, line)){ // Loop through components in input file
                if (line.find("END COMPONENTS") != string::npos) {
                    foutDEF << line << "\n";
                    break;
                }

                if(line.find("-") != string::npos){ // Component entry
                    stringstream ss(line);
                    ss>>tmp>>name>>tmp;

                    if(line.find("+") != string::npos){ // Check for multi-line properties

                        foutDEF<<"  - "<<name<<" "<<tmp; // Write component name

                        ss>>tmp;
                        while( tmp!="PLACED" && tmp!="FIXED" ) { // Copy properties
                            foutDEF<<" "<<tmp;
                            ss>>tmp;
                        }

                        Component c = components[name]; // Get updated component data
                        string mtype =  c.fixed?" FIXED":" PLACED";
                        // Write new PLACED/FIXED line with updated coordinates
                        foutDEF << mtype << " ( " << c.x <<" "<< c.y <<" ) "<< c.orient;
                        ss>>tmp>>tmp>>tmp>>tmp>>tmp; // Skip old coordinates
                        
                        while(ss>>tmp) foutDEF<<" "<<tmp; // Write remaining properties
                        foutDEF<<endl;
                    }
                    else{
                        foutDEF<<line<<endl; // Copy line as is
                    }
                }
                else if(line.find("PLACED") != string::npos || line.find("FIXED") != string::npos){
                    // Handle single-line component definition (part of multi-line)
                    stringstream ss(line);
                    ss>>tmp;
                    while( tmp!="PLACED" && tmp!="FIXED" ) {
                        foutDEF<<" "<<tmp;
                        ss>>tmp;
                    }
                    
                    Component c = components[name]; // Get updated component data
                    string mtype =  c.fixed?" FIXED":" PLACED";
                    // Write new PLACED/FIXED line
                    foutDEF <<mtype << " ( " << c.x <<" "<< c.y <<" ) "<< c.orient;
                    ss>>tmp>>tmp>>tmp>>tmp>>tmp; // Skip old coordinates

                    while(ss>>tmp) foutDEF<<" "<<tmp; // Write remaining properties
                    foutDEF<<endl;
                }
                else{
                    foutDEF<<line<<endl; // Copy line as is
                }
                

            }
        }
        else{
            foutDEF << line << "\n"; // Copy all other lines
        }
    }

    finDEF.close();
    foutDEF.close();
}

int main(int argc, char* argv[]) {

    auto startTime = chrono::steady_clock::now();   // start time

    ios::sync_with_stdio(false); // Faster C++ I/O
    cin.tie(nullptr);

    if (argc != 4) { // Check argument count
        std::cerr << "Usage: " << argv[0] << " <input> <output> <number of partitions>\n";
        return 1;
    }
    string inputLEF = argv[1];
    string inputDEF = argv[2];
    string outputDEF= argv[3];

    string line;
    int units = 0;
    // Data structures for design
    unordered_map<string,vector<float>> CoreSite;
    unordered_map<string, Cell> cells;
    vector<Row> rows;
    unordered_map<string, Component> componentMap;
    vector<Net> nets;

    // Parse input files
    readLEF(inputLEF,units,CoreSite,cells);
    readDEF(inputDEF,cells,CoreSite,rows,componentMap,nets);
    //=============================================================================

    // Create placer object and run optimization
    FastDP fastDP(CoreSite, cells, rows, componentMap, nets);
    fastDP.Optimize(startTime);

    //=============================================================================
    // Write output DEF file with new component locations
    writeDEF(inputDEF,outputDEF,componentMap);

    // show total time
    // auto now = chrono::steady_clock::now();
    // auto elapsed = chrono::duration_cast<chrono::seconds>(now - startTime).count();
    // cout<<"time: "<<elapsed<<endl;

    return 0;
}