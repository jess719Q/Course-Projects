# Image Stitching

This script performs automatic stitching of multiple images in a folder to create a panorama using feature matching, homography, and blending. It builds on the `stitch.py` module for keypoint detection, image warping, and blending.

## Features

- Automatically reads all images from a folder (`.png` or `.jpg`).
- Detects keypoints and matches between images.
- Computes transformation matrices (homography) to align images.
- Supports linear and pyramid blending.
- Optional cylindrical or spherical warping.
- Can display intermediate stitching process with matplotlib.

## Usage

```
python stitch_folder.py --folder <input_folder> --output <output_file> \
                        [--focal FOCAL_LENGTH] [--count MIN_MATCHES] \
                        [--blender {linear,pyramid}] [--warp {none,cylindrical,spherical}] \
                        [--show]
```

### Arguments

* `--folder` (required): Path to the folder containing input images.
* `--output` (required): Output filename for the stitched panorama.
* `--focal` (optional, default=1000): Focal length used for cylindrical/spherical warping.
* `--count` (optional, default=50): Minimum number of keypoint matches to consider a valid link between images.
* `--blender` (optional, default=`linear`): Blending method. Options:
    * `linear` – simple linear blending
    * `pyramid` – Laplacian pyramid blending
* `--warp` (optional, default=`none`): Warp method. Options:
    * `none` – no warp
    * `cylindrical` – cylindrical projection
    * `spherical` – spherical projection
* `--show` (optional): Show intermediate stitching results.

## Results
<img width="1711" height="835" alt="image" src="https://github.com/user-attachments/assets/73043533-1905-47b0-bd72-eb83f2466ebf" />
