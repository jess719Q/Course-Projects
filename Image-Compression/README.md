# JPEG-like Image Compression

This project implements a **still image transform coding** system, a prototype of JPEG compression. You can use any programming language you are familiar with. The project supports both **gray level** and **color images**.

## Implementation Details

- 8x8 DCT is provided (TMN version optimized for H.263 video coding).
- Apply quantization and coding to compress the images.
- Quantization tables can be adjusted using the **Quality Factor (QF)**.
- Compressed images can be recovered to `.raw` format for viewing.
- Calculate **PSNR** between the original and compressed images.

## Compilation

Compile the programs using `g++` (C++17 standard):

### Compile

```
g++ ./encode.cpp ./src/*.cpp -o encode.exe -std=c++17
g++ ./decode.cpp ./src/*.cpp -o decode.exe -std=c++17
g++ ./psnr.cpp -o psnr.exe
```

### Usage

Encode image:
```
./encode.exe image.raw -o imgJPG.bmp -qf QF (-c gray)
```

Decode image:
```
./decode.exe imgJPG.bmp -o imgBack.raw -qf QF (-c gray)
```

Calculate PSNR:
```
./psnr.exe -a image.raw -b imgBack.raw (-c gray)
```


### Parameters explanation

- `image.raw` original image file

- `imgJPG.bmp` compressed image file

- `imgBack.raw` : recovered image file

- `QF` Quality Factor (1–100)

- `-c gray` optional flag for gray level images

## Results

The PSNR results show the quality of compressed images at different QFs. Higher QF → better image quality.

Below is an example of compressed images:

<img width="1059" height="1318" alt="image" src="https://github.com/user-attachments/assets/cc4ac3aa-7f20-459b-a957-bac211361ee8" />
