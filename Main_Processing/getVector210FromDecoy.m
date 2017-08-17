function averageDistanceVectors = getVector210FromDecoy(wildTypePDBFileName, decoySequenceList, numOfDecoy)

wildTypePDBFileName = pdbread('3OXC.pdb');

[caXYZList, resNames] = getCAXYZListFromPdbFile(wildTypePDBFileName);


numOfCA = 198;

%% only take the first 198, coz the decoy sets only has 198
if size(caXYZList, 2) ~= numOfCA
    error('Error: 198');
end

if numel(resNames) ~= numOfCA
    error('Error: need 198 residue names');
end

x = caXYZList(1, :)';
y = caXYZList(2, :)';
z = caXYZList(3, :)';

dt = DelaunayTri(x, y, z);

% tetramesh(dt.Triangulation, dt.X, 'FaceAlpha', 0.5, 'edgecolor', 'interp');
% tetramesh(dt.Triangulation, dt.X, 'FaceAlpha', 0.5);
% cameratoolbar;

%% Decoys only has sequence, not structure. The sequences are stored in a text file.
fid = fopen(decoySequenceList);
decoys = textscan(fid, '%s', 'delimiter', '\n');
fclose(fid);

decoys = decoys{1}; % now that decoys is a cell-array of strings

n = numel(decoys);

if n < numOfDecoy
    error('Error: number of decoy is less than required (numOfDecoy)');
end

n = numOfDecoy;

averageDistanceVectors = zeros(n, 210);

matlabpool;
parfor it = 1:n
% for it = 1:n
    it
    decoyEesNames = [decoys{it}, decoys{it}];
    
    averageDistanceVectors(it, :) = computeAverageDistanceVectorFromDelauNayAndResidueNames(dt, decoyEesNames);
end
matlabpool close;