function iNew = repack2DColorMap(I)

S = size(I);
newSize = S(1) * S(2);
iNew = NaN(newSize,3);
itemp = 1;
for i = 1:S(1)
    for j = 1:S(2)
        iNew(itemp,:) = I(i,j,:);
        itemp = itemp + 1;
    end
end

