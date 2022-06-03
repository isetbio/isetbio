function overlapCoeff = overlap(weights1, weights2)
    commonArea = sum(min(weights1(:), weights2(:)));
    jointArea = sum(max(weights2(:), weights2(:)));
    %jointArea = 0.5*(sum(weights1(:))+sum(weights2(:)));
    overlapCoeff = commonArea/jointArea;
end
