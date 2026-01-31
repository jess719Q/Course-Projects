#include "src/myimage.h"

int main(int argc, char* argv[]) {
    string inputFile, outputFile;
    int QF = 50;              // Default Quality Factor
    bool grayscale = false;   // Whether to use grayscale mode

    // Parse command line arguments
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-o") == 0) {
            outputFile = argv[++i];  // Get output file name
        } else if (strcmp(argv[i], "-qf") == 0) {
            QF = atoi(argv[++i]);    // Get quality factor
        } else if (strcmp(argv[i], "-c") == 0) {
            if (strcmp(argv[++i], "gray") == 0) grayscale = true; // Set grayscale flag
        } else if (argv[i][0] != '-') {
            inputFile = argv[i];     // First non-flag argument is input file
        }
    }

    // Convert quality factor to quantization scale
    if (QF < 50)
        QF = 5000 / QF;
    else
        QF = 200 - 2 * QF;

    const int width = 512;    // Image width
    const int height = 512;   // Image height
    const int framesize = width * height; // Grayscale image size in bytes

    // Open the input file
    ifstream fin(inputFile, ios::binary);
    if (!fin) {
        cerr << "Failed to open input file.\n";
        return 1;
    }

    // Read the encoded bitstream from the file
    vector<unsigned char> encodedData((istreambuf_iterator<char>(fin)), {});
    string bitstream = "";
    for (unsigned char byte : encodedData) {
        bitset<8> bits(byte);
        bitstream += bits.to_string();  // Convert each byte to an 8-bit binary string
    }

    // Decode the bitstream into a coefficient array (DC + AC)
    vector<int> decoded = ACDCdecode(bitstream, height, width, grayscale);

    // Perform inverse quantization and inverse DCT to reconstruct the image
    vector<unsigned char> image = iquantDct2(decoded, QF, height, width, grayscale);

    // Save the reconstructed image
    if (grayscale)
        saveRawImage(outputFile, &image[0], framesize); // Save grayscale image
    else {
        vector<unsigned char> imageRGB = YUV2RGB(image, width, height); // Convert YUV to RGB
        saveRawImage(outputFile, &imageRGB[0], framesize * 3); // Save color image
    }

    return 0;
}