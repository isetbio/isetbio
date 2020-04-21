function rfValues = lookUpValues(rfPositions, tabulatedPositions, tabulatedValue, neighborsNum)
    neighborsNum = 1;
    [D, I] = pdist2(tabulatedPositions, rfPositions, 'euclidean', 'Smallest', neighborsNum);

    if (neighborsNum > 1)
        totalD = sum(D,1);
        b = zeros(neighborsNum, size(D,2));
        for k = 1:neighborsNum
            b(k,:) = (totalD - D(k,:)) ./ totalD;
        end
        rfValues = sum(b.*tabulatedValue(I),1);
    else
        rfValues = tabulatedValue(I);
    end
    
    rfValues = rfValues';
end