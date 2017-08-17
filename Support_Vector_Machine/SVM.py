import csv
import os
import sys

import numpy as np
from sklearn.model_selection import train_test_split
from sklearn import svm
import pylab as pl
from sklearn.decomposition import PCA

os.getcwd()
os.chdir('/Users/GSU/Desktop')

from numpy import genfromtxt
my_data = genfromtxt('SQV_SVM.csv', delimiter=',')        

Data = np.delete(my_data, 0, 1)

Target = my_data[:,[0]]

X_train, X_test, y_train, y_test = train_test_split(Data, Target, test_size=0.4, random_state=0)

Data.shape, Target.shape # (5387, 210), (5387, 1) where 1 is resistance 100 and above, 0 is resistance below 1

X_train.shape, y_train.shape # (3232, 210), (3232, 1)

X_test.shape, y_test.shape # (2155, 210), (2155, 1)

clf = svm.SVC(kernel='linear', C=1).fit(X_train, y_train)

clf.score(X_test, y_test) # 0.99721577726218102

clf.score(X_test, y_test) # 0.90545808966861596 where 1 is resistance 100 and above, 0 is resistance below 1 and rest are classified as 3

clf.score(X_test, y_test) # SQV 0.97807017543859653 where resistant fold < 3.0 were classified as non- resistant (susceptible), denoted as 0; while those with the relative resistant fold ≥ 3.0 were classified as resistant, denoted as 1


# SQV 0.97807017543859653, TPV 0.97522063815342841, LPV 0.987548828125, 

pca = PCA(n_components=2).fit(X_train)
pca_2d = pca.transform(X_train)
svmClassifier_2d =   svm.LinearSVC(random_state=111).fit(   pca_2d, y_train)


for i in range(0, pca_2d.shape[0]):
    if y_train[i] == 0:
        c1 = pl.scatter(pca_2d[i,0],pca_2d[i,1],c='r',    s=50,marker='+')
    elif y_train[i] == 1:
        c2 = pl.scatter(pca_2d[i,0],pca_2d[i,1],c='g',    s=50,marker='o')

pl.legend([c1, c2], ['Non-Resistant', 'Resistant'])
x_min, x_max = pca_2d[:, 0].min() - 1,   pca_2d[:,0].max() + 1
y_min, y_max = pca_2d[:, 1].min() - 1,   pca_2d[:, 1].max() + 1
xx, yy = np.meshgrid(np.arange(x_min, x_max, .01),   np.arange(y_min, y_max, .01))
Z = svmClassifier_2d.predict(np.c_[xx.ravel(),  yy.ravel()])
Z = Z.reshape(xx.shape)
pl.contour(xx, yy, Z)
pl.title('Support Vector Machine Decision Surface for ATV')
#pl.axis('off')
pl.show()

     
   
