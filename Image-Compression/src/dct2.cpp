/* 2-D DCT (type II) program */
/* This file contains two subprograms. The first one is "dct2(x,n)",
which performs the forward 2-D DCT, and the second one is "idct2(x,n)", 
which performs the inverse 2-D DCT.  The program, dct2 (or idct2),
will replace the input x (2-D square array [0..n-1][0..n-1]) by 
its discrete cosine transform (or inverse discrete cosine transform).
The array size is n*n where n must be an integer power of 2. */

#include "myimage.h"

void dct1(float *x, int n);
void idct1(float *x, int n);

void dct2(float **x, int n)
{
  int i,j;
  float *y;

  y = (float *) calloc (n,sizeof(float));
  if (y == NULL) {
   printf("allocation failure\n");
   exit(1);
  }
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) 
      y[j] = x[j][i];
    dct1(y,n);
    for (j=0;j<n;j++) 
      x[j][i] = y[j];
  }   /* end of loop i */

  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) 
      y[j] = x[i][j];
    dct1(y,n);
    for (j=0;j<n;j++) 
      x[i][j] = y[j];
  }   /* end of loop i */

  free(y);
}

/* ----------------------------------------------- */

void idct2(float **x, int n)
{
  int i,j;
  float *y;

  y = (float *) calloc (n,sizeof(float));
  if (y == NULL) {
   printf("allocation failure\n");
   exit(1);
  }
  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) 
      y[j] = x[j][i];
    idct1(y,n);
    for (j=0;j<n;j++) 
      x[j][i] = y[j];
  }   /* end of loop i */

  for (i=0;i<n;i++) {
    for (j=0;j<n;j++) 
      y[j] = x[i][j];
    idct1(y,n);
    for (j=0;j<n;j++) 
      x[i][j] = y[j];
  }   /* end of loop i */

  free(y);
}

/* ----------------------------------------------- */

const int luminanceQuantMatrix[8][8] = {
    {16, 11, 10, 16,  24,  40,  51,  61},
    {12, 12, 14, 19,  26,  58,  60,  55},
    {14, 13, 16, 24,  40,  57,  69,  56},
    {14, 17, 22, 29,  51,  87,  80,  62},
    {18, 22, 37, 56,  68, 109, 103,  77},
    {24, 35, 55, 64,  81, 104, 113,  92},
    {49, 64, 78, 87, 103, 121, 120, 101},
    {72, 92, 95, 98, 112, 100, 103,  99}
};

const int chrominanceQuantMatrix[8][8] = {
    {17, 18, 24, 47, 99, 99, 99, 99},
    {18, 21, 26, 66, 99, 99, 99, 99},
    {24, 26, 56, 99, 99, 99, 99, 99},
    {47, 66, 99, 99, 99, 99, 99, 99},
    {99, 99, 99, 99, 99, 99, 99, 99},
    {99, 99, 99, 99, 99, 99, 99, 99},
    {99, 99, 99, 99, 99, 99, 99, 99},
    {99, 99, 99, 99, 99, 99, 99, 99}
};

vector<int> quantDct2(vector<unsigned char>& img, int QF, int height, int width, bool gray = false) {
    int framesize = width * height;

    vector<int> imgOut;
    if (gray)
        imgOut.resize(framesize);       // Only luminance (Y) component
    else
        imgOut.resize(framesize * 3 / 2); // Y + subsampled U + V components (YUV 4:2:0 layout)

    // Allocate an 8x8 block for DCT input
    float** arr = (float**)malloc(8 * sizeof(float*));
    for (int i = 0; i < 8; ++i)
        arr[i] = (float*)malloc(8 * sizeof(float));

    // Process the luminance (Y) channel in 8x8 blocks
    for (int y = 0; y < height; y += 8) {
        for (int x = 0; x < width; x += 8) {
            int idx = y * width + x;

            // Copy and level-shift the 8x8 block from the image
            for (int i = 0; i < 8; ++i)
                for (int j = 0; j < 8; ++j)
                    arr[i][j] = static_cast<float>(img[idx + i * width + j] - 128);  // Center range to [-128,127]

            dct2(arr, 8); // Perform 2D DCT on the block

            // Quantize the DCT coefficients using the luminance quantization matrix
            for (int i = 0; i < 8; ++i)
                for (int j = 0; j < 8; ++j)
                    imgOut[idx + i * width + j] =
                        static_cast<int>(round(arr[i][j] / luminanceQuantMatrix[i][j] * QF / 100.0));
        }
    }

    if (!gray) {
        // Process chrominance (U and V) channels in 8x8 blocks
        for (int y = 0; y < height; y += 8) {
            for (int x = 0; x < width / 2; x += 8) {
                int idx = y * (width / 2) + x + framesize;  // Start index for U/V blocks

                // Copy and level-shift the 8x8 block
                for (int i = 0; i < 8; ++i)
                    for (int j = 0; j < 8; ++j)
                        arr[i][j] = static_cast<float>(img[idx + i * (width / 2) + j] - 128);

                dct2(arr, 8); // Perform DCT

                // Quantize using the chrominance quantization matrix
                for (int i = 0; i < 8; ++i)
                    for (int j = 0; j < 8; ++j)
                        imgOut[idx + i * (width / 2) + j] =
                            static_cast<int>(round(arr[i][j] / chrominanceQuantMatrix[i][j] * QF / 100.0));
            }
        }
    }

    // Free dynamically allocated memory
    for (int i = 0; i < 8; ++i)
        free(arr[i]);
    free(arr);

    return imgOut;
}


/* ----------------------------------------------- */

vector<unsigned char> iquantDct2(vector<int>& img, int QF, int height, int width, bool gray = false) {
    int framesize = width * height;

    vector<unsigned char> imgOut;
    if (gray)
        imgOut.resize(framesize);        // For grayscale: only Y channel
    else
        imgOut.resize(framesize * 3 / 2); // For color: Y + subsampled U and V (YUV 4:2:0)

    // Allocate an 8x8 block for inverse DCT
    float** arr = (float**)malloc(8 * sizeof(float*));
    for (int i = 0; i < 8; ++i)
        arr[i] = (float*)malloc(8 * sizeof(float));

    // Process the luminance (Y) channel in 8x8 blocks
    for (int y = 0; y < height; y += 8) {
        for (int x = 0; x < width; x += 8) {
            int idx = y * width + x;

            // Dequantize each coefficient by multiplying with quantization matrix
            for (int i = 0; i < 8; ++i)
                for (int j = 0; j < 8; ++j)
                    arr[i][j] = static_cast<float>(
                        img[idx + i * width + j] * luminanceQuantMatrix[i][j] * 100.0 / QF);

            idct2(arr, 8); // Apply 2D inverse DCT

            // Add 128 to shift back from [-128,127] to [0,255]
            for (int i = 0; i < 8; ++i)
                for (int j = 0; j < 8; ++j)
                    imgOut[idx + i * width + j] = static_cast<unsigned char>(round(arr[i][j] + 128));
        }
    }

    if (!gray) {
        // Process chrominance (U and V) channels in 8x8 blocks
        for (int y = 0; y < height; y += 8) {
            for (int x = 0; x < width / 2; x += 8) {
                int idx = y * (width / 2) + x + framesize;

                // Dequantize using chrominance quantization matrix
                for (int i = 0; i < 8; ++i)
                    for (int j = 0; j < 8; ++j)
                        arr[i][j] = static_cast<float>(
                            img[idx + i * (width / 2) + j] * chrominanceQuantMatrix[i][j] * 100.0 / QF);

                idct2(arr, 8); // Apply 2D inverse DCT

                // Re-shift to [0,255] range
                for (int i = 0; i < 8; ++i)
                    for (int j = 0; j < 8; ++j)
                        imgOut[idx + i * (width / 2) + j] = static_cast<unsigned char>(round(arr[i][j] + 128));
            }
        }
    }

    // Free dynamically allocated memory
    for (int i = 0; i < 8; ++i)
        free(arr[i]);
    free(arr);

    return imgOut;
}
