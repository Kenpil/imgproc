#@title Default title text
from IPython.display import clear_output
clear_output()
import scipy
import numpy
from scipy.io import loadmat
depthMapOri = scipy.io.loadmat('/content/drive/My Drive/develop/Sea-thru/depthMap5.mat')
BcRemovedImOri = scipy.io.loadmat('/content/drive/My Drive/develop/Sea-thru/BcRemovedM.mat')

import cupy as cp
import numpy as np
height = 1827
width = 2737
depthMap = cp.array(depthMapOri['depthMap'])
BcRemovedIm = cp.array(BcRemovedImOri['BcRemovedIm'])
p = 1/5000
maxDepth = 50
maxDepthThre = 0.2
acM = cp.zeros(height*width*3).reshape(height,width,3).astype('double')
repetition = 5000
print('Go!!')

NeM = cp.zeros(height*width).reshape(height,width).astype('f')
index1 = cp.where((abs(depthMap[1:height-1,1:width-1]-depthMap[1:height-1,0:width-2]) < maxDepthThre) & (depthMap[1:height-1,0:width-2] < maxDepth))
index2 = cp.where((abs(depthMap[1:height-1,1:width-1]-depthMap[1:height-1,2:width]) < maxDepthThre) & (depthMap[1:height-1,2:width] < maxDepth))
index3 = cp.where((abs(depthMap[1:height-1,1:width-1]-depthMap[0:height-2,1:width-1]) < maxDepthThre) & (depthMap[0:height-2,1:width-1] < maxDepth))
index4 = cp.where((abs(depthMap[1:height-1,1:width-1]-depthMap[2:height,1:width-1]) < maxDepthThre) & (depthMap[2:height,1:width-1] < maxDepth))
NeM[index1[0]+1,index1[1]+1] += 1
NeM[index2[0]+1,index2[1]+1] += 1
NeM[index3[0]+1,index3[1]+1] += 1
NeM[index4[0]+1,index4[1]+1] += 1
index = cp.where(NeM > 0)

for i in range (repetition):
    acMtmp = cp.zeros(height*width*3).reshape(height,width,3).astype('double')
    acMtmp[index1[0]+1,index1[1]+1,:] += acM[index1[0]-1+1,index1[1]+1,:]
    acMtmp[index2[0]+1,index2[1]+1,:] += acM[index2[0]+1+1,index2[1]+1,:]
    acMtmp[index3[0]+1,index3[1]+1,:] += acM[index3[0]+1,index3[1]-1+1,:]
    acMtmp[index4[0]+1,index4[1]+1,:] += acM[index4[0]+1,index4[1]+1+1,:]

    for j in range(3):
        acMtmp[index[0],index[1],j] /= NeM[index[0],index[1]]
        
    acM = p*BcRemovedIm + (1-p)*acMtmp
    if (i+1)%100 == 0:
        print(i+1)
        #fileString = '/content/drive/My Drive/develop/Sea-thru/acM' +str(i) + '1000.mat'
        #print(fileString)
        #scipy.io.savemat(fileString, {'acM':acM})
acMnp = cp.asnumpy(acM)
scipy.io.savemat('/content/drive/My Drive/develop/Sea-thru/acM5000.mat', {'acM':acMnp})
