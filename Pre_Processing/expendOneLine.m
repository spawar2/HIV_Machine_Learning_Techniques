function expendedRows = expendOneLine(wildTypeTxt, thisRow, limitOfNumOfVariations)

nc = size(thisRow, 2);

%% case 1. see a ".", or "~", or "#", or "*", means various things, just REMOVE that line
for ic = 1:nc
    str = char(thisRow{ic}); % now str is a string
    
    if ~isempty(strfind(str, '.')) || ~isempty(strfind(str, '*')) || ~isempty(strfind(str, '~')) || ~isempty(strfind(str, '#'))
        expendedRows = [];
        return;
    end
end



%% case 2. process -, which means "just use the Amino acid at the top row"
numOfVariations = 1;
positionOfVar = [];
for ic = 1:nc
    str = char(thisRow{ic}); % now str is a string
    numOfVariations = numOfVariations*numel(thisRow{ic});
    
    if (numel(thisRow{ic}) >= 2)
        positionOfVar = [positionOfVar, [ic; numel(thisRow{ic})]];
    end
    
    if strcmp(str, '-')
        %             txt(ir,ic) = {'ASDF'}; % need to enclose string by {}
        thisRow(ic) = wildTypeTxt(ic);
    end
end

%% case 3, if no duplicates, then just return this row
% if numOfVariations == 1
%     expendedRows = thisRow;
% else
%     expendedRows = thisRow;
% end

numOfVariations = min([numOfVariations, limitOfNumOfVariations]);

expendedRows = repmat(thisRow, [numOfVariations, 1]);

for iv = 1:numOfVariations
    indexToChoose = getPosition1(positionOfVar, iv);
    for ic = 1:size(positionOfVar, 2)
        thisStr = expendedRows{iv, positionOfVar(1, ic)};
        expendedRows(iv, positionOfVar(1, ic)) = {thisStr(indexToChoose(ic))};
    end
end