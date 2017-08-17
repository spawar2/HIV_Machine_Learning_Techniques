clear;
close all;
clc;

templatePDBFileName = pdbread('http://files.rcsb.org/download/3OXC.pdb');

%% compute 210 vectors for a data set, ATV
datePrefix = 'C:/Users/spawar2/Desktop/TPV';
numOfDecoy = 7363; % read the *_preprocessd.txt file and see how many lines are there, chose a number <= numOfLines


wildType = computeAverageDistanceVector(templatePDBFileName);
resistantLevels = load([datePrefix, '_all_preprocessed_resistant.txt']);
tic;
SQVVectors = getVector210FromDecoy(templatePDBFileName, [datePrefix, '_all_preprocessed.txt'], numOfDecoy); 
toc
save([datePrefix, '_all_preprocessed_vectors.mat']);



