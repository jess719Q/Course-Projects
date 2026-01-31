#include "src/myimage.h"

int main(int argc, char* argv[]) {
    string inputFile, outputFile;
    int QF = 50;               // Default Quality Factor
    bool grayscale = false;   // Flag for grayscale mode

    // Parse command line arguments
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-o") == 0) {
            outputFile = argv[++i];  // Output file path
        } else if (strcmp(argv[i], "-qf") == 0) {
            QF = atoi(argv[++i]);    // Read Quality Factor
        } else if (strcmp(argv[i], "-c") == 0) {
            if (strcmp(argv[++i], "gray") == 0) grayscale = true; // Enable grayscale mode
        } else if (argv[i][0] != '-') {
            inputFile = argv[i];     // The first non-option argument is the input file path
        }
    }

    // Convert Quality Factor to quantization scale (JPEG standard approximation)
    if (QF < 50) 
        QF = 5000 / QF;
    else      
        QF = 200 - 2 * QF;

    const int width = 512;     // Fixed image width
    const int height = 512;    // Fixed image height
    const int framesize = width * height;

    vector<int> imageDCT;      // Vector to hold quantized DCT coefficients

    if (grayscale) {
        // Grayscale mode
        vector<unsigned char> imageGRAY(framesize);  // Allocate grayscale image buffer
        if (!readRawImage(inputFile, imageGRAY))     // Read raw grayscale image
            return 1;
        imageDCT = quantDct2(imageGRAY, QF, height, width, grayscale); // Perform DCT and quantization
    } else {
        // Color mode
        vector<unsigned char> imageRGB(framesize * 3);  // Allocate RGB image buffer
        if (!readRawImage(inputFile, imageRGB))         // Read raw RGB image
            return 1;
        vector<unsigned char> imageYUV = RGB2YUV(imageRGB, width, height); // Convert RGB to YUV
        imageDCT = quantDct2(imageYUV, QF, height, width, grayscale);      // Perform DCT and quantization on YUV
    }

    // Encode the DCT coefficients into a bitstream (DC + AC encoding)
    string bitstream = DCAC(imageDCT, height, width, grayscale);

    // Save the bitstream to a file (e.g., a custom bitmap or binary format)
    if (saveBitmap(bitstream, outputFile))
        cout << "Compressed bitstream saved to " << outputFile << endl;

    return 0;
}
