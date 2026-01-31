import os
import argparse
import numpy as np
import cv2
from matplotlib import pyplot as plt
from stitch import *  # 你原本的 stitch.py

def getimgs(folder):
    # List image files
    file_names = [f for f in os.listdir(folder) if f.lower().endswith((".png", ".jpg", ".jpeg"))]
    file_names.sort()
    imgs = []
    kps = []
    for file_name in file_names:
        file_path = os.path.join(folder, file_name)
        img = cv2.imread(file_path)
        imgs.append(img)
        kps.append(keypoints(img))
    return imgs, kps

def findImgLink(kps, count):
    n = len(kps)
    img_link = [[] for _ in range(n)]
    goods = [[[] for _ in range(n)] for __ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            good = Stich.findMatch(kps[i], kps[j])
            if len(good) > count:
                img_link[i].append(j)
                img_link[j].append(i)
                goods[i][j] = good
    return img_link, goods

def FolderStiching(folder, output, f=1000, count=50, blender=None, warpper=None, showProcess=False):
    imgs, kps = getimgs(folder)
    img_link, goods = findImgLink(kps, count)
    print(f"Image links: {img_link}")
    # Choose the most connected image as base
    mid = max(range(len(img_link)), key=lambda i: len(img_link[i]))
    done = [False] * len(img_link)
    done[mid] = True

    if warpper:
        result_img = warpper(imgs[mid], f)
    else:
        result_img = cv2.cvtColor(imgs[mid], cv2.COLOR_BGR2BGRA)

    while sum(done) < len(img_link):
        for i, links in enumerate(img_link):
            if done[i]:
                for j in links:
                    if not done[j]:
                        next_idx = j
                        break
                if not done[next_idx]:
                    break

        done[next_idx] = True

        last_img = result_img
        next_img = warpper(imgs[next_idx], f) if warpper else cv2.cvtColor(imgs[next_idx], cv2.COLOR_BGR2BGRA)
        result_img = Stich.wrap2img(last_img, next_img, blender, count)

        if showProcess:
            print(f"Processed {sum(done)}/{len(done)} images")
            plt.imshow(result_img)
            plt.show()

    cv2.imwrite(output, result_img)
    print(f"Result saved to {output}")
    return result_img

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Batch image stitching")
    parser.add_argument("--folder", required=True, help="Input folder containing images")
    parser.add_argument("--output", required=True, help="Output file name")
    parser.add_argument("--focal", type=float, default=1000, help="Focal length")
    parser.add_argument("--count", type=int, default=50, help="Minimum keypoint matches")
    parser.add_argument("--blender", choices=["linear", "pyramid"], default="linear", help="Blending method")
    parser.add_argument("--warp", choices=["none", "cylindrical", "spherical"], default="none", help="Warp method")
    parser.add_argument("--show", action="store_true", help="Show intermediate results")

    args = parser.parse_args()

    # Map blender choice to function
    blender_func = Blender.linear if args.blender == "linear" else Blender.pyramid

    # Map warp choice to function
    warp_func = None
    if args.warp == "cylindrical":
        warp_func = Stich.cylindricalWarp
    elif args.warp == "spherical":
        warp_func = Stich.sphericalWarp

    FolderStiching(args.folder, args.output, f=args.focal, count=args.count,
                   blender=blender_func, warpper=warp_func, showProcess=args.show)
