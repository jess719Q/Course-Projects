#include "myimage.h"

void saveRawImage(string filename, const unsigned char* img, int size) {
    ofstream ofs(filename, ios::binary);
    if (!ofs) {
        cerr << "Error opening file for writing: " << filename << endl;
        return;
    }

    ofs.write(reinterpret_cast<const char*>(img), size);
    ofs.close();

    cout << "Saved raw image to: " << filename << endl;
}

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

bool saveBitmap(string bitstream, string outputFile){
    ofstream fout(outputFile, ios::binary);
    if (!fout) {
        cerr << "Cannot open output file.\n";
        return false;
    }

    // Pad to 8 bits and write bytes to file
    while (bitstream.size() % 8 != 0) bitstream += '0';
    for (size_t i = 0; i < bitstream.size(); i += 8) {
        bitset<8> byte(bitstream.substr(i, 8));
        unsigned char b = static_cast<unsigned char>(byte.to_ulong());
        fout.write(reinterpret_cast<char*>(&b), 1);
    }

    fout.close();

    return true;
}