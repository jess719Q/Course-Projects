#include "myimage.h"
#include "HuffmanTable.h"

//Encode signed integers to binary (JPEG style)
string bincode(int val){
    string bin = "";
    if(val>=0){
        while (val > 0) {
            int bit = val%2;
            bin.push_back('0' + bit);
            val /= 2;
        }
    }
    else{
        val = abs(val);
        while (val > 0) {
            int bit = (val+1)%2;
            bin.push_back('0' + bit);
            val /= 2;
        }
    }

	reverse(bin.begin(), bin.end());
    return bin;
}
    
string DCAC(vector<int> img, int height, int width, bool gray=false) {
    string bitstream = "";

    // Encode luminance blocks
    for (int y = 0; y < height; y+=8) {
        for (int x = 0; x < width; x+=8) {
            int idx = y * width + x;

            // DC difference encoding (DPCM)
            int DIFF = (idx == 0) ? img[idx]
                     : (x == 0) ? img[idx] - img[idx - width * 8]
                                : img[idx] - img[idx - 8];

            bitstream += (DIFF == 0) ? luminanceDC[0]
                                     : luminanceDC[(int)log2(abs(DIFF)) + 1];
            bitstream += bincode(DIFF);

            // AC run-length and Huffman encoding
            int n0 = 0;
            for (int i = 1; i < 64; ++i) {
                idx = y * width + x + zigzagIndex[i][0] + zigzagIndex[i][1] * width;
                if (img[idx] == 0) {
                    n0++;
                } else {
                    while (n0 > 15) {
                        bitstream += luminanceAC[15 * 11]; // ZRL
                        n0 -= 15;
                    }
                    bitstream += luminanceAC[n0 * 11 + (int)log2(abs(img[idx])) + 1];
                    bitstream += bincode(img[idx]);
                    n0 = 0;
                }
            }
            bitstream += luminanceAC[0]; // End-of-block
        }
    }

    // Chrominance (U, V)
    if (!gray) {
        int framesize = height * width;
        for (int y = 0; y < height; y += 8) {
            for (int x = 0; x < width / 2; x += 8) {
                int idx = y * width/2 + x + framesize;

                int DIFF = (idx == 0) ? img[idx]
                         : (x == 0) ? img[idx] - img[idx - width * 4]
                                    : img[idx] - img[idx - 8];

                bitstream += (DIFF == 0) ? chrominanceDC[0]
                                         : chrominanceDC[(int)log2(abs(DIFF)) + 1];
                bitstream += bincode(DIFF);

                int n0 = 0;
                for (int i = 1; i < 64; ++i) {
                    idx = y * width/2 + x + framesize + zigzagIndex[i][0] + zigzagIndex[i][1] * width/2;
                    if (img[idx] == 0) {
                        n0++;
                    } else {
                        while (n0 > 15) {
                            bitstream += chrominanceAC[165]; // ZRL
                            n0 -= 15;
                        }
                        bitstream += chrominanceAC[n0 * 11 + (int)log2(abs(img[idx])) + 1];
                        bitstream += bincode(img[idx]);
                        n0 = 0;
                    }
                }
                bitstream += chrominanceAC[0]; // End-of-block
            }
        }
    }

    return bitstream;
}


struct Node {
    int symbol;
    Node* left;
    Node* right;
    Node() : symbol(-1), left(nullptr), right(nullptr) {}
};


// Insert a symbol and its binary code into the Huffman decoding tree
void insertCode(Node* root, const string& code, int symbol) {
    Node* node = root;
    for (char c : code) {
        if (c == '0') {
            if (!node->left) node->left = new Node();
            node = node->left;
        } else {
            if (!node->right) node->right = new Node();
            node = node->right;
        }
    }
    node->symbol = symbol;  // Set the symbol at the leaf node
}

Node* buildTree(const char** table, int size){

    vector<pair<string, int>> symbols;  // Pair of (symbol, code length)
    for (int i = 0; i < size; ++i) {
        symbols.emplace_back(table[i], i);
    }

    sort(symbols.begin(), symbols.end(), [](auto& a, auto& b) {
        if (a.second == b.second) return a.first < b.first;
        return a.second < b.second;
    });

    // Step 3: Build Huffman decoding tree from the canonical codebook
    Node* root = new Node();
    for (auto& kv : symbols) {
        insertCode(root, kv.first, kv.second);
    }

    return root;
}


//Decode JPEG binary string to integer
int bindecode(string bin){
    int val=0;

    if(bin[0]=='0'){
        for(char c : bin){
            val = val<<1;
            if(c=='0') val++;
        }
        val = -val;
    }
    else{
        for(char c : bin){
            val = val<<1;
            if(c=='1') val++;
        }
    }

    return val;
}

vector<int> ACDCdecode(string bitstream, int height, int width, bool gray=false) {
    int framesize = height * width;

    // Build Huffman decoding trees for luminance and chrominance DC/AC
    Node* luDC = buildTree(luminanceDC, 12);
    Node* chDC = buildTree(chrominanceDC, 12);
    Node* luAC = buildTree(luminanceAC, 176);
    Node* chAC = buildTree(chrominanceAC, 176);

    // Vector to store decoded DPCM coefficients (in zigzag order)
    vector<int> decoded;
    Node* node = luDC; // Start with luminance DC tree

    int mod = 0; // 0: expecting DC symbol, 1: expecting AC symbol or zeros
    int typ = 0; // 0: luminance, 1: chrominance
    string num = ""; // Holds the binary number for DIFF or AC coefficient
    int numlen = 0;  // Number of bits remaining to read for DIFF or AC value

    // Main decoding loop over bitstream
    for (char bit : bitstream) {

        // If we are currently reading the actual value (after Huffman length code)
        if (numlen > 0) {
            num += bit;
            numlen--;
            if (numlen <= 0) {
                decoded.push_back(bindecode(num)); // Convert binary string to int
                num = "";
            }
            continue;
        }

        // Huffman tree traversal
        node = (bit == '0') ? node->left : node->right;

        // Safety check for invalid traversal
        if (!node) {
            cerr << "Decoding error: invalid traversal in Huffman tree.\n";
            break;
        }

        // Leaf node found — process symbol
        if (node->symbol != -1) {
            if (mod == 0) { // Handling DC component
                numlen = node->symbol;
                if (numlen == 0) decoded.push_back(0); // Zero difference
                mod = 1; // Next we will read AC components
                node = (typ == 0) ? luAC : chAC; // Switch to AC tree
            } else if (mod == 1) { // Handling AC component
                int idx = node->symbol;

                if (idx == 0) { // End of Block (EOB)
                    while (decoded.size() % 64 != 0) decoded.push_back(0); // Fill rest of block with zeros
                    mod = 0; // Next block will begin with DC again
                    if (decoded.size() >= framesize) typ = 1; // Switch to chrominance after luminance
                    node = (typ == 0) ? luDC : chDC; // Switch to correct DC tree
                    continue;
                } else {
                    int n0 = idx / 11; // Run of zeros
                    while (n0-- > 0) decoded.push_back(0);
                    numlen = idx % 11; // Bits to read for the non-zero value
                    node = (typ == 0) ? luAC : chAC; // Remain in AC tree
                }
            }
        }
    }

    // Reconstruct image from decoded coefficients
    vector<int> imgOut;
    if (gray) imgOut.resize(framesize);
    else imgOut.resize(framesize * 3 / 2); // Account for chroma subsampling

    // Decode luminance blocks
    for (int y = 0; y < height; y += 8) {
        for (int x = 0; x < width; x += 8) {
            int idx = y * width + x;
            int idx0 = (y / 8 * width / 8 + x / 8) * 64;

            // Reconstruct DC coefficient
            if (idx == 0)
                imgOut[idx] = decoded[idx0];
            else if (x == 0)
                imgOut[idx] = imgOut[idx - width * 8] + decoded[idx0];
            else
                imgOut[idx] = imgOut[idx - 8] + decoded[idx0];

            // Place remaining AC coefficients in zigzag order
            for (int i = 0; i < 8; ++i) {
                for (int j = 0; j < 8; ++j) {
                    int k = i * 8 + j;
                    if (k == 0) continue; // DC already handled
                    int idx2 = y * width + x + zigzagIndex[k][0] + zigzagIndex[k][1] * width;
                    imgOut[idx2] = decoded[idx0 + k];
                }
            }
        }
    }

    // Decode chrominance blocks if not grayscale
    if (!gray) {
        for (int y = 0; y < height; y += 8) {
            for (int x = 0; x < width / 2; x += 8) {
                int idx = y * width / 2 + x + framesize;
                int idx0 = (y / 8 * width / 16 + x / 8) * 64 + framesize;

                // Reconstruct DC coefficient
                if (idx == 0)
                    imgOut[idx] = decoded[idx0];
                else if (x == 0)
                    imgOut[idx] = imgOut[idx - width * 4] + decoded[idx0];
                else
                    imgOut[idx] = imgOut[idx - 8] + decoded[idx0];

                // Place AC coefficients
                for (int i = 0; i < 8; ++i) {
                    for (int j = 0; j < 8; ++j) {
                        int k = i * 8 + j;
                        if (k == 0) continue;
                        int idx2 = y * width / 2 + framesize + x + zigzagIndex[k][0] + zigzagIndex[k][1] * width / 2;
                        imgOut[idx2] = decoded[idx0 + k];
                    }
                }
            }
        }
    }

    return imgOut;
}
