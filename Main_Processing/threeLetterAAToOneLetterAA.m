function oneLetterAA = threeLetterAAToOneLetterAA(name)

AAListThreeLetter = {'Ala', 'Arg', 'Asn', 'Asp', 'Cys', 'Glu', 'Gln', 'Gly', 'His', 'Ile', 'Leu', 'Lys', 'Met', 'Phe', 'Pro', 'Ser', 'Thr', 'Trp', 'Tyr', 'Val'};

AAListOneLetter = {'A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V'};

numStr = numel(AAListThreeLetter);

for it = 1:numStr
    if strcmpi(name, AAListThreeLetter{it})
        oneLetterAA = AAListOneLetter{it};
        return;
    end
end