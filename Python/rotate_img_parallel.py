#!/usr/bin/env python
from glob import glob
from PIL import Image
from time import time

ImgLst = glob('/Users/PLS/Downloads/Bab2/*')
def rotateImage(ImgNm):
    try:
        img = Image.open(ImgNm)
        width, height = img.size
        if (width > height):
            out = img.transpose(Image.ROTATE_90)
            out.save(ImgNm)
            out.close()
        img.close()
    except IOError:
        pass

if __name__ == '__main__':
    from  multiprocessing import Pool
    startTime = time()
    pool = Pool(processes=4)
    pool.map(rotateImage, ImgLst)

    endTime = time()
    print "%d images"%len(ImgLst)
    print "execute time: %.2fs"%(endTime-startTime)
