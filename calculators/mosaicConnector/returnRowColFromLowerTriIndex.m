function [rows,cols] = returnRowColFromLowerTriIndex(indices, N)
    totalElementsNum = (N*(N-1))/2;
    reversedIndices = totalElementsNum-indices;
    k = floor( (sqrt(1+8*reversedIndices)-1)/2);
    j = reversedIndices - k.*(k+1)/2;
    rows = N-j;
    cols = N-(k+1);
end
