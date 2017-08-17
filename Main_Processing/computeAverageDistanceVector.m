function averageDistanceVector = computeAverageDistanceVector(pdbFileName)

pdbFileName = 'C:/Users/spawar2/Desktop/3OXC.pdb';

[caXYZList, resNames] = getCAXYZListFromPdbFile(pdbFileName);
x = caXYZList(1, :)';
y = caXYZList(2, :)';
z = caXYZList(3, :)';

% figure;
% plot3(x, y, z, '.');
% axis equal;
% cameratoolbar;


%% do Delaunay
dt = DelaunayTri(x, y, z);

% figure;
% tetramesh(dt.Triangulation, dt.X, 'FaceAlpha', 0.05);
% cameratoolbar;


% AAList = {'Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Glu', 'Gln', 'Gly', 'His', 'Ile', 'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val'};

% AAListThreeLetter = {'Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Glu', 'Gln', 'Gly', 'His', 'Ile', 'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val'};
AAListOneLetter = {'A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'};


distance = zeros(20, 20);
occurance = zeros(20, 20);

numTetra = size(dt.Triangulation, 1);
for iTetra = 1:numTetra
    thisTetra = dt.Triangulation(iTetra, :);
    
    for it = 1:4
        for itt = (it + 1):4
            firstAtomID = thisTetra(it);
            firstAtomResName = resNames(firstAtomID);
            firstAtomResID = getResidueID(firstAtomResName, AAListOneLetter);
            firstAtomCoord = dt.X(firstAtomID, :);
            
            secondAtomID = thisTetra(itt);
            secondAtomResName = resNames(secondAtomID);
            secondAtomResID = getResidueID(secondAtomResName, AAListOneLetter);
            secondAtomCoord = dt.X(secondAtomID, :);
            
            dist = sqrt(sum((firstAtomCoord - secondAtomCoord).^2));
            
            
            distance(min([firstAtomResID, secondAtomResID]), max([firstAtomResID, secondAtomResID])) = distance(min([firstAtomResID, secondAtomResID]), max([firstAtomResID, secondAtomResID])) + dist;
            occurance(min([firstAtomResID, secondAtomResID]), max([firstAtomResID, secondAtomResID])) = occurance(min([firstAtomResID, secondAtomResID]), max([firstAtomResID, secondAtomResID])) + 1;
        end
    end
end

averageDistance = distance./(occurance + eps);
% figure;
% imshow(averageDistance, []);
% colormap jet;

ii = 0;
averageDistanceVector = zeros(1, 210);
for it = 1:20
    for itt = it:20
        ii = ii + 1;
        averageDistanceVector(ii) = averageDistance(it, itt);
    end
end
