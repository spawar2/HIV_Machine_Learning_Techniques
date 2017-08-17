function it = getResidueID(secondAtomResName, AAList)

numStr = numel(AAList);

for it = 1:numStr
    if strcmpi(secondAtomResName, AAList{it})
        return;
    end
end