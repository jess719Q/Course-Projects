#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <cstring>
#include <unordered_map>
#include <queue>
#include <bitset>
#include <algorithm>

#define PI 3.141592653589793
#define SQH 0.707106781186547  /* square root of 2 */
#define SWAP(a,b)  tempr=(a); (a) = (b); (b) = tempr

using namespace std;


bool readRawImage(string, vector<unsigned char>&);
void saveRawImage(string, const unsigned char*, int);
bool saveBitmap(string, string);
vector<unsigned char> RGB2YUV(const vector<unsigned char>&, int, int);
vector<unsigned char> YUV2RGB(const vector<unsigned char>& , int, int);
vector<int> quantDct2(vector<unsigned char>&, int , int, int, bool);
vector<unsigned char> iquantDct2(vector<int>&, int , int, int, bool);
string DCAC(vector<int>, int, int, bool);
vector<int> ACDCdecode(string, int, int, bool);