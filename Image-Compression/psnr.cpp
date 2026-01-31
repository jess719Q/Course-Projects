#include <stdio.h>
#include <math.h>
#include <iostream>
#include <cstring>
#include <vector>
#include <fstream>

using namespace std;

bool readRawImage(string filename, vector<unsigned char>& img) {
    ifstream file(filename, ios::binary);
    if (!file) {
        cerr << "Failed to open " << filename << endl;
        return false;
    }

    file.read(reinterpret_cast<char*>(&img[0]), img.size());
    if (!file) {
        cerr << "Error reading file or file too short." << endl;
        return false;
    }

    file.close();
    return true;
}

float calculatePSNR(vector<unsigned char> original, vector<unsigned char> reconstructed, int size, float maxPixelValue) {
    float mse = 0.0f;
    for (int i = 0; i < size; ++i)
        mse += pow(static_cast<float> (original[i] - reconstructed[i]), 2);
    mse /= size;
    if (mse == 0) return INFINITY;  // no difference
    return 10.0f * log10((maxPixelValue * maxPixelValue) / mse);
}

int main(int argc, char* argv[]){

    string originFile, compressFile;
    bool grayscale=false;

    // Parse command line arguments
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-a") == 0) {
            originFile = argv[++i];
        }
        else if (strcmp(argv[i], "-b") == 0) {
            compressFile = argv[++i];
        }
        else if (strcmp(argv[i], "-c") == 0) {
            if(strcmp(argv[++i], "gray") == 0) grayscale=true;
        }
    }
    
    int framesize = grayscale? 512*512 : 512*512*3;
    
    vector<unsigned char> originImage(framesize);
    if (!readRawImage(originFile, originImage)) return 1;
    vector<unsigned char> compressImage(framesize);
    if (!readRawImage(compressFile, compressImage)) return 1;

    float psnr = calculatePSNR(originImage,compressImage,framesize,256);

    cout<<"psnr="<<psnr<<endl;
}