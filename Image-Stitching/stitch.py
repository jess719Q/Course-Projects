import numpy as np
import cv2
from matplotlib import pyplot as plt
import os

class keypoints:
    def __init__(self, img):
        sift = cv2.SIFT_create()
        self.kp, self.des = sift.detectAndCompute(img,None)

class Blender:
    
    def __init__(self, image1, image2):
        self.img1 = image1
        self.img2 = image2
        
    def linear(self):
        
        and_img = cv2.bitwise_and(self.img1, self.img2)
        _, overlap = cv2.threshold(and_img[:,:,3], 1, 255, cv2.THRESH_BINARY)
        overlap = overlap.astype(np.uint8)
        h, w = overlap.shape[:2]

        left_edge = np.argmax(overlap[h//2]== 255)
        if self.img1[h//2][left_edge-10][3]!=0:
            left_img = self.img1
            right_img = self.img2
        else:
            left_img = self.img2
            right_img = self.img1
        
        weighting = -np.ones((h, w, 3))
        
        for i in range(h):
            l=0
            for j in range(w):
                if(overlap[i][j]==255):
                    weighting[i][j][:]=l
                    l+=1
            if l>0:
                weighting[i] /= l
        
        overlap_img = right_img*np.concatenate([weighting,np.ones((h,w,1))], axis=2) \
                    + left_img*np.concatenate([(1-weighting),np.ones((h,w,1))], axis=2)
        overlap_img = cv2.bitwise_and(overlap_img,overlap_img,mask=overlap).astype(np.uint8)
        result = cv2.bitwise_or(left_img,right_img,mask=cv2.bitwise_not(overlap))
        result = cv2.bitwise_or(result,overlap_img)
        return result
    
    def pyramid(self):
        and_img = cv2.bitwise_and(self.img1, self.img2)
        _, overlap0 = cv2.threshold(and_img[:,:,3], 1, 255, cv2.THRESH_BINARY)
        overlap0 = overlap0.astype(np.uint8)
        h, w = overlap0.shape[:2]
        
        new_w = (w // 64) * 64 if w % 64 == 0 else ((w // 64) + 1) * 64
        new_h = (h // 64) * 64 if h % 64 == 0 else ((h // 64) + 1) * 64
        
        left_edge = np.argmax(overlap0[h//2]== 255)
        if self.img1[h//2][left_edge-10][3]!=0:
            left_img = cv2.resize(self.img1, (new_w, new_h))
            right_img = cv2.resize(self.img2, (new_w, new_h))
        else:
            left_img = cv2.resize(self.img2, (new_w, new_h))
            right_img = cv2.resize(self.img1, (new_w, new_h))
            
        and_img = cv2.bitwise_and(left_img, right_img)
        _, overlap = cv2.threshold(and_img[:,:,3], 1, 255, cv2.THRESH_BINARY)
        overlap = overlap.astype(np.uint8)
        
        
        # Find the indices of 1s
        indices_of_ones = np.where(overlap[h//2] == 255)[0]
        # Get the middle ones
        middle = indices_of_ones[len(indices_of_ones) // 2]
        
        # generate Gaussian pyramid for A
        G = left_img.copy()
        gpA = [G]
        for i in range(6):
            G = cv2.pyrDown(G)
            gpA.append(G)
        
        # generate Gaussian pyramid for B
        G = right_img.copy()
        gpB = [G]
        for i in range(6):
            G = cv2.pyrDown(G)
            gpB.append(G)
        
        # generate Laplacian Pyramid for A
        lpA = [gpA[5]]
        for i in range(5,0,-1):
            GE = cv2.pyrUp(gpA[i])
            L = cv2.subtract(gpA[i-1],GE)
            lpA.append(L)
        
        # generate Laplacian Pyramid for B
        lpB = [gpB[5]]
        for i in range(5,0,-1):
            GE = cv2.pyrUp(gpB[i])
            L = cv2.subtract(gpB[i-1],GE)
            lpB.append(L)
        
        # Now add left and right halves of images in each level
        LS = []
        a=5
        for la,lb in zip(lpA,lpB):
            mid = middle//2**a
            ls = np.hstack((la[:,0:mid], lb[:,mid:]))
            LS.append(ls)
            
            a-=1
        
        # now reconstruct
        ls_ = LS[0]
        for i in range(1,6):
            ls_ = cv2.pyrUp(ls_)
            ls_ = cv2.add(ls_, LS[i])
            
        ls_[:,:,3] = np.where(ls_[:,:,3] > 100, 255, 0)
        ls_[ls_[:,:,3] < 100, 0:3] = 0
        
        weighting = np.zeros((new_h, new_w,3))
        weighting0 = np.zeros((new_h, new_w,3))
        weighting1 = np.zeros((new_h, new_w,3))
        for i in range(new_h):
            l=0
            r=0
            for j in range(new_w):
                if(overlap[i][j]==255):
                    if(j<middle-new_w//100):
                        weighting0[i][j][:]=l
                        weighting[i][j][:]=l
                        l+=1
                    elif(j>=middle+new_w//100):
                        weighting1[i][j][:]=r
                        weighting[i][j][:]=r
                        r-=1
                    elif(j>=middle-new_w//10 and j<middle+new_w//10):
                        weighting[i][j][:]=1
            if(l>0):
                weighting0[i][:middle-new_w//100] = 1-weighting0[i][:middle-new_w//100]/l
                weighting[i][:middle-new_w//100] /=l
            if(r<0):
                weighting1[i][middle+new_w//100:] /=r
                weighting[i][middle+new_w//100:] = 1-weighting[i][middle+new_w//100:]/r

        overlap_img = ls_*np.concatenate([weighting,np.ones((new_h, new_w,1))], axis=2)\
                    + left_img*np.concatenate([(weighting0),np.ones((new_h, new_w,1))], axis=2)\
                    + right_img*np.concatenate([(weighting1),np.ones((new_h, new_w,1))], axis=2)
        overlap_img = cv2.bitwise_and(overlap_img,overlap_img,mask=overlap).astype(np.uint8)
        result0 = cv2.bitwise_or(left_img,right_img,mask=cv2.bitwise_not(overlap))
        result = cv2.bitwise_or(result0,overlap_img)
        result[:,:,3] = np.where(result[:,:,3] > 1, 255, 0)
            
        return result
        
        

class Stich:
    def __init__(self, count):
        self.MIN_MATCH_COUNT = count

    def gaussian(x, mu, sig):
        return (
            1.0 / (np.sqrt(2.0 * np.pi) * sig) * np.exp(-np.power((x - mu) / sig, 2.0) / 2)
        )

    def imgPyr(h, w):
        arr = np.zeros((h,w))
        for i in range(h):
            for j in range(w):
                # arr[i][j] = min(j,w-j-1)
                arr[i][j] = min(i,j,h-i-1,w-j-1)
        return arr

    def findKeypoint(img1, img2):
        # SIFT detector
        sift = cv2.SIFT_create()
        # find the keypoints and descriptors with SIFT
        kp1, des1 = sift.detectAndCompute(img1,None)
        kp2, des2 = sift.detectAndCompute(img2,None)
        flann = cv2.DescriptorMatcher_create(cv2.DescriptorMatcher_FLANNBASED)
        matches = flann.knnMatch(des1,des2,k=2)

        # store all the good matches as per Lowe's ratio test.
        good = []
        for m,n in matches:                 
            if m.distance < 0.75*n.distance:
                good.append(m)
        return good, kp1, kp2

    def findMatch(kp1,kp2):
        flann = cv2.DescriptorMatcher_create(cv2.DescriptorMatcher_FLANNBASED)
        matches = flann.knnMatch(kp1.des,kp2.des,k=2)

        # store all the good matches as per Lowe's ratio test.
        good = []
        for m,n in matches:                 
            if m.distance < 0.75*n.distance:
                good.append(m)
        return good

    def findTransMatrix(kp1, kp2, good, match_count=50):
        if len(good)>match_count:
            pts1 = np.float32([ kp1[m.queryIdx].pt for m in good ]).reshape(-1,1,2)
            pts2 = np.float32([ kp2[m.trainIdx].pt for m in good ]).reshape(-1,1,2)
            H21, mask = cv2.findHomography(pts2, pts1, cv2.RANSAC,5.0)
            
        else:
            print( f"Not enough matches are found - {len(good)}/{match_count}")
            H21 = None
        return H21

    def warp_new(img1, img2, H, blender=None):
        h1, w1 = img1.shape[:2]
        h2, w2 = img2.shape[:2]

        corners1 = np.float32([[0, 0], [0, h1], [w1, h1], [w1, 0]]).reshape(-1, 1, 2)
        corners2 = np.float32([[0, 0], [0, h2], [w2, h2], [w2, 0]]).reshape(-1, 1, 2)
        
        warped_corners2 = cv2.perspectiveTransform(corners2, H)

        corners = np.concatenate((corners1, warped_corners2), axis=0)
        
        [xmin, ymin] = np.int32(corners.min(axis=0).ravel() - 0.5)
        [xmax, ymax] = np.int32(corners.max(axis=0).ravel() + 0.5)
        wn, hn = (xmax - xmin, ymax - ymin)

        t = [-xmin, -ymin]
        Ht = np.array([[1, 0, t[0]], [0, 1, t[1]], [0, 0, 1]])
        
        warped_img2 = cv2.warpPerspective(img2, Ht @ H, (wn, hn))
        
        warped_img2[(warped_img2[:,:,3]<255) & (warped_img2[:,:,3]>0)]=0
        
        warped_img1 = np.zeros((hn, wn,4),dtype="uint8")
        warped_img1[t[1]:h1 + t[1], t[0]:w1 + t[0]] = img1
        
        if blender == None:
            _, binary_alpha = cv2.threshold(warped_img2[:, :, 3], 1, 255, cv2.THRESH_BINARY_INV)
            binary_alpha = binary_alpha.astype(np.uint8)
            
            warped_img1 = cv2.bitwise_and(warped_img1,warped_img1, mask=binary_alpha)
            out = cv2.bitwise_or(warped_img1, warped_img2)
            
        else:
            b=Blender(warped_img1, warped_img2)
            out = blender(b)

        
        return out

    def wrap2img(img1, img2, blender=None, count=50):
        img1_grap = cv2.cvtColor(img1, cv2.COLOR_BGRA2GRAY)
        img2_grap = cv2.cvtColor(img2, cv2.COLOR_BGRA2GRAY)
        
        good, kp1, kp2 = Stich.findKeypoint(img1_grap, img2_grap)
        H = Stich.findTransMatrix(kp1, kp2, good, count)
        if(H is None):
            return None
        result_img = Stich.warp_new(img1, img2, H, blender)
        
        return result_img

    def cylindricalWarp(img, f):
        '''This function returns the cylindrical warp for a given image and intrinsics matrix K'''
        h_,w_ = img.shape[:2]
        
        K = np.array([[f, 0, w_/2],
                    [0, f, h_/2],
                    [0, 0, 1]],dtype="float32")
        
        # pixel coordinates
        y_i, x_i = np.indices((h_,w_))
        X = np.stack([x_i,y_i,np.ones_like(x_i)],axis=-1).reshape(h_*w_,3) # to homog
        Kinv = np.linalg.inv(K) 
        X = Kinv.dot(X.T).T # normalized coords
        # calculate cylindrical coords (sin\theta, h, cos\theta)
        A = np.stack([np.sin(X[:,0]),X[:,1],np.cos(X[:,0])],axis=-1).reshape(w_*h_,3)
        B = K.dot(A.T).T # project back to image-pixels plane
        # back from homog coords
        B = B[:,:-1] / B[:,[-1]]
        # make sure warp coords only within image bounds
        B[(B[:,0] < 0) | (B[:,0] >= w_) | (B[:,1] < 0) | (B[:,1] >= h_)] = -1
        B = B.reshape(h_,w_,-1)
        
        # img_rgba = cv2.cvtColor(img,cv2.COLOR_BGR2BGRA) # for transparent borders...
        # warp the image according to cylindrical coords
        img_warp = cv2.remap(cv2.cvtColor(img,cv2.COLOR_BGR2BGRA), B[:,:,0].astype(np.float32), 
                             B[:,:,1].astype(np.float32), cv2.INTER_AREA, borderMode=cv2.BORDER_TRANSPARENT)
        return img_warp
    
    def sphericalWarp(img, f):
        '''This function returns the cylindrical warp for a given image and intrinsics matrix K'''
        h_,w_ = img.shape[:2]
        
        K = np.array([[f, 0, w_/2],
                    [0, f, h_/2],
                    [0, 0, 1]],dtype="float32")
        
        # pixel coordinates
        y_i, x_i = np.indices((h_,w_))
        X = np.stack([x_i,y_i,np.ones_like(x_i)],axis=-1).reshape(h_*w_,3) # to homog
        Kinv = np.linalg.inv(K) 
        X = Kinv.dot(X.T).T # normalized coords
        # calculate spherical coords
        A = np.stack([np.sin(X[:,0])*np.cos(X[:,1]),np.sin(X[:,1]),np.cos(X[:,0])*np.cos(X[:,1])],axis=-1).reshape(w_*h_,3)
        
        B = K.dot(A.T).T # project back to image-pixels plane
        # back from homog coords
        B = B[:,:-1] / B[:,[-1]]
        # make sure warp coords only within image bounds
        B[(B[:,0] < 0) | (B[:,0] >= w_) | (B[:,1] < 0) | (B[:,1] >= h_)] = -1
        B = B.reshape(h_,w_,-1)
        
        # img_rgba = cv2.cvtColor(img,cv2.COLOR_BGR2BGRA) # for transparent borders...
        # warp the image according to cylindrical coords
        img_warp = cv2.remap(cv2.cvtColor(img,cv2.COLOR_BGR2BGRA), B[:,:,0].astype(np.float32), 
                             B[:,:,1].astype(np.float32), cv2.INTER_AREA, borderMode=cv2.BORDER_TRANSPARENT)
        return img_warp