clear;
close all;
clc;



%% preprocess
sequenceFileName = 'C:/Users/spawar2/Desktop/NFV_All_Trial.txt';
width = 100;

outputFileName = sequenceFileName;
outputFileName = strrep(outputFileName, '.txt', '_preprocessed.txt');

expend = 1;
preprocessSequenceTextData1(sequenceFileName, width, outputFileName, expend);


