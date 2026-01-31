#include "myimage.h"

vector<unsigned char> RGB2YUV(const vector<unsigned char>& imageRGB, int width = 0, int height = 0) {

    vector<unsigned char> imageY, imageU, imageV;

    imageY.resize(width * height);
    imageU.resize(width * height / 4);
    imageV.resize(width * height / 4);

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            int idx = y * width + x;
            int r = imageRGB[3 * idx];
            int g = imageRGB[3 * idx + 1];
            int b = imageRGB[3 * idx + 2];

            int Y = r * 0.299000 + g * 0.587000 + b * 0.114000;
            int U = r *-0.168736 + g *-0.331264 + b * 0.500000 + 128;
            int V = r * 0.500000 + g *-0.418688 + b *-0.081312 + 128;


            imageY[idx] = static_cast<unsigned char>(Y);


            if (y % 2 == 0 && x % 2 == 0) {
                int groupIdx = (y / 2) * (width / 2) + (x / 2);
                imageU[groupIdx] = static_cast<unsigned char>(U);
                imageV[groupIdx] = static_cast<unsigned char>(V);
            }
        }
    }

    vector<unsigned char> imageYUV;
    imageYUV.reserve(width * height + width * height / 2);
    imageYUV.insert(imageYUV.end(), imageY.begin(), imageY.end());
    imageYUV.insert(imageYUV.end(), imageU.begin(), imageU.end());
    imageYUV.insert(imageYUV.end(), imageV.begin(), imageV.end());

    return imageYUV;
}

vector<unsigned char> YUV2RGB(const vector<unsigned char>& imageYUV, int width = 0, int height = 0) {

    const int frameSize = width * height;
    const int chromaWidth = width / 2;
    const int chromaHeight = height / 2;

    const unsigned char* Y_plane = imageYUV.data();
    const unsigned char* U_plane = Y_plane + frameSize;
    const unsigned char* V_plane = U_plane + (frameSize / 4);

    vector<unsigned char> imageRGB;
    imageRGB.resize(frameSize * 3);

    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            int idx = y * width + x;
            int chromaIdx = (y / 2) * chromaWidth + (x / 2);

            int Y = Y_plane[idx];
            int U = U_plane[chromaIdx];
            int V = V_plane[chromaIdx];

            // Convert YUV to RGB using float math and clamp the result
            float c = Y;
            float d = U - 128;
            float e = V - 128;

            int R = clamp(static_cast<int>(c + 1.402f * e + 0.5f), 0, 255);
            int G = clamp(static_cast<int>(c - 0.344136f * d - 0.714136f * e + 0.5f), 0, 255);
            int B = clamp(static_cast<int>(c + 1.772f * d + 0.5f), 0, 255);

            imageRGB[3 * idx    ] = static_cast<unsigned char>(R);
            imageRGB[3 * idx + 1] = static_cast<unsigned char>(G);
            imageRGB[3 * idx + 2] = static_cast<unsigned char>(B);
        }
    }

    return imageRGB;
}