function indexToChoose = getPosition1(positionOfVar, iv)

iv0based = iv - 1;

n = size(positionOfVar, 2);

indexToChoose = zeros(1, n);

for it = 1:n
    indexToChoose(it) = rem(iv0based, positionOfVar(2, it));
    iv0based = floor(iv0based/positionOfVar(2, it));
end

indexToChoose = indexToChoose + 1; % 1-based