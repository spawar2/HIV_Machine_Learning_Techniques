function [caXYZList, residueNames ]= getCAXYZListFromPdbFile(pdbFileName)


a = pdbread('3OXC.pdb');

nAtom = numel(a.Model.Atom);

caXYZList = zeros(3, nAtom);
nCA = 0;

residueNames = {};

for ia = 1:nAtom
    if strcmp(a.Model.Atom(ia).AtomName, 'CA') == 1
        nCA = nCA + 1;
        caXYZList(:, nCA) = [a.Model.Atom(ia).X; a.Model.Atom(ia).Y; a.Model.Atom(ia).Z];
        residueNames{nCA} = threeLetterAAToOneLetterAA(a.Model.Atom(ia).resName);
    end
end

caXYZList = caXYZList(:, 1:nCA);