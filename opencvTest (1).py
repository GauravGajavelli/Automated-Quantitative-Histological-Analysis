import numpy as np
import cv2
from matplotlib import pyplot as plt
import sys
import os
import glob

shamArea = 0
intraArea = 0
periArea = 0
vehicleArea = 0
numSham = 0
numIntra = 0
numPeri = 0
numVehicle = 0
currentSham = 0
currentIntra = 0
currentPeri = 0
currentVehicle = 0

def getBrainArea(filepath, showImg):
    img = cv2.imread(filepath,0)
    img2 = cv2.imread(filepath)
    # img = cv2.medianBlur(img,7)
    ret,thresh1 = cv2.threshold(img,208,255,cv2.THRESH_BINARY)
    titles = ['OG','Basic Thresholding']
    images = [img2,thresh1]
    pixelCount = 0
    for i in range(thresh1.shape[0]):
        for j in range(thresh1.shape[1]):
            pixel = thresh1.item(i,j)
            if pixel == 0:
                pixelCount = pixelCount+1
    if (showImg == 1):
        # print(pixelCount)
        for i in range(2):
            plt.subplot(2,1,i+1),plt.imshow(images[i],'gray')
            plt.title(titles[i])
            plt.xticks([]),plt.yticks([])
        plt.show()
    return pixelCount


# path = 'F:\Documents 2\GUAVVY\2020 REsearch\Images\Lesion'
def pywalker(path):
    global shamArea
    global intraArea
    global periArea
    global vehicleArea

    global currentSham
    global currentIntra
    global currentPeri
    global currentVehicle

    global numSham
    global numIntra
    global numPeri
    global numVehicle
    
    for root, dirs, files in os.walk(path): 
        for dir_ in dirs:
            dirPath = root + '/' + dir_
            # print(dirPath)
            for root2, dirs2, files2 in os.walk(dirPath):
                # print('in a new ' + dir_[7:] + ' dir of size: ' + str(len(next(os.walk(dirPath))[2])))
                numFiles = len(glob.glob1(dirPath,"*.jpg"))
                # print('in a new ' + dir_[7:] + ' dir of size: ' + str(numFiles))
                if numFiles == 22:
                    for file_ in files2:
                        if (file_[-3:] == 'jpg'):
                            print(dirPath+'/'+file_)
                            sum = numSham + numIntra + numVehicle + numPeri
                            show = 0
                            if sum % 20 == 0:
                                show = 1
                            dirType = dir_[7:]
                            if(dirType == 'sham'):
                                shamArea = shamArea + getBrainArea(dirPath+'/'+file_,show)
                                numSham = numSham + 1
                                currentSham = currentSham + getBrainArea(dirPath+'/'+file_,show)
                            elif (dirType == 'intra'):
                                intraArea = intraArea + getBrainArea(dirPath+'/'+file_,show)
                                numIntra = numIntra + 1
                                currentIntra = currentIntra + getBrainArea(dirPath+'/'+file_,show)
                            elif (dirType == 'peri'):
                                periArea = periArea + getBrainArea(dirPath+'/'+file_,show)
                                numPeri = numPeri + 1
                                currentPeri = currentPeri + getBrainArea(dirPath+'/'+file_,show)
                            elif (dirType == 'vehicle'):
                                vehicleArea = vehicleArea + getBrainArea(dirPath+'/'+file_,show)
                                numVehicle = numVehicle + 1
                                currentVehicle = currentVehicle + getBrainArea(dirPath+'/'+file_,show)
                    print('currentSham: ' + str(float(currentSham)))
                    print('currentIntra: ' + str(float(currentIntra)))
                    print('currentPeri: ' + str(float(currentPeri)))
                    print('currentVehicle: ' + str(float(currentVehicle)))
                    print('-----------------------------------------------')
                    currentSham = 0
                    currentIntra = 0
                    currentPeri = 0
                    currentVehicle = 0

if __name__ == '__main__':
    pywalker('F:/Documents 2/GUAVVY/2020 REsearch/Images/Lesion')
    print('sham: ' + str((float(shamArea)/float(numSham))))
    print('intra: ' + str((float(intraArea)/float(numIntra))))
    print('peri: ' + str((float(periArea)/float(numPeri))))
    print('vehicle: ' + str((float(vehicleArea)/float(numVehicle))))


