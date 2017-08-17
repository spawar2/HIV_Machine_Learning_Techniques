function preprocessSequenceTextData1(sequenceFileName, width, outputTxtFileName, expand)


if nargin <= 3
    expand = 0;
end

%% read sequences from txt file to a cell matrix
% sequenceFileName = './data/xlsx20111107/APV.txt';
% width = 101; % the number of cells in a line. It should be 100, but for some reasons there is a tailing empty string which task a cell. So 101

fid = fopen(sequenceFileName);
% txt = textscan(fid, '%s', 'delimiter', '\t');
txt = textscan(fid, '%s', 'delimiter', ',');
fclose(fid);

txt = txt{1};

height = numel(txt)/width;
txt = reshape(txt, [width, height])';

% txt = txt(:, 1:(end-1)); % echo the above width = 101, we need to remove the last non-sense column. If width = 100, dont need this.

resistantValueStr = txt(:, 1);
resistantValue = zeros(size(resistantValueStr));
for it = 1:numel(resistantValue)
    resistantValue(it) = str2double(resistantValueStr{it});
end

txt = txt(:, 2:end); % The first column is resistant values, take off


if expand == 1
    
    %% process the strings, has these cases:
    % case 1. see a "-", means we want to copy from top line
    % case 2. see a ".", or "~", or "#", or "*", means various things, just REMOVE that line
    % case 3. see a cell with more than one character, means there are more than one choices for that position, expend them.
    
    [nr, nc] = size(txt);
    
    limitOfNumOfVariationsOfOneLine = 256;
    
    allTxt = cell(nr*limitOfNumOfVariationsOfOneLine, nc);
    allResistant = zeros(nr*limitOfNumOfVariationsOfOneLine, 1);
    wildTypeTxt = txt(1, :);
    currentPassEndPosition = 1;
    for ir = 2:nr
        disp(['have expended to ', int2str(currentPassEndPosition-1), ' rows, now working on row ', int2str(ir), ' of ', int2str(nr)]);
        expendedThisLine = expendOneLine(wildTypeTxt, txt(ir, :), limitOfNumOfVariationsOfOneLine);
        if ~isempty(expendedThisLine)
            allResistant(currentPassEndPosition:(currentPassEndPosition + size(expendedThisLine, 1) - 1)) = ones(size(expendedThisLine, 1), 1)*resistantValue(ir);
            allTxt(currentPassEndPosition:(currentPassEndPosition + size(expendedThisLine, 1) - 1), :) = expendedThisLine;
            
            currentPassEndPosition = currentPassEndPosition + size(expendedThisLine, 1);
        end
    end
    
    allResistant = allResistant(1:(currentPassEndPosition-1));
    allTxt = allTxt(1:(currentPassEndPosition-1), :);
    
    
    %% output to file
    
    fid = fopen(outputTxtFileName, 'w');
    sequences = [];
    irr = 1;
    for ir = 1:size(allTxt, 1)
        thisSequence = [];
        
        for ic = 1:size(allTxt, 2)
            
            str = allTxt(ir, ic); % str is a cell of string
            str = char(str{1}); % now str is a string
            
            thisSequence = [thisSequence, str];
        end
        
        fprintf(fid, '%s\n', thisSequence);
    end
    
    fclose(fid);
    
    % output resistant values
    outputTxtFileName = strrep(outputTxtFileName, '.txt', '_resistant.txt');
    dlmwrite(outputTxtFileName, allResistant, 'delimiter', ' ');
else
    % expend = 0
    
    [nr, nc] = size(txt);
    
    
    for ir =  1:nr
        for ic = 1:nc
            
            str = txt(ir, ic); % str is a cell of string
            str = char(str{1}); % now str is a string
            
            %% process -, which means "just use the Amino acid at the top row"
            if strcmp(str, '-')
                %             txt(ir,ic) = {'ASDF'}; % need to enclose string by {}
                txt(ir,ic) = txt(1, ic);
            end
            
        end
    end
    
    
    %% There are positions with two letters, indicating changing to either or them. Skip them for now.
    
    fid = fopen(outputTxtFileName, 'w');
    
    
    outputResistantFileName = strrep(outputTxtFileName, '.txt', '_resistant.txt');
    fid1 = fopen(outputResistantFileName, 'w');
    
    
    
    sequences = [];
    irr = 1;
    for ir = 1:nr
        thisSequence = [];
        for ic = 1:nc
            
            str = txt(ir, ic); % str is a cell of string
            str = char(str{1}); % now str is a string
            
            thisSequence = [thisSequence, str];
            
        end
        
        thisSequence = strrep(thisSequence, '.', ''); % '.' means deletion. skip them also.
        
        if length(thisSequence) == nc
            sequences{irr} = thisSequence;
            irr = irr + 1;
            
            fprintf(fid, '%s\n', thisSequence);
            fprintf(fid1, '%f ', resistantValue(ir));
        end
    end
    
    
    fclose(fid);
    
    fprintf(fid1, '\n');
    fclose(fid1);
end
