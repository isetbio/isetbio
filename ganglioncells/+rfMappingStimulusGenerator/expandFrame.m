function theFullFrame = expandFrame(theFrame, expansionFactor)

    rowsNum = size(theFrame,1);
    colsNum = size(theFrame,2);
    theFullFrame = zeros(rowsNum * expansionFactor, colsNum * expansionFactor);
    yy = 1:expansionFactor;
    xx = 1:expansionFactor;
    for row = 1:rowsNum
        for col = 1:colsNum
            theFullFrame((row-1)*expansionFactor+yy, (col-1)*expansionFactor+xx) = theFrame(row,col);
        end
    end
end
