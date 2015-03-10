function res = rgcIsCovarianceMatrix(R)
    res = 1;
    res = res && all(size(R) == [2 2]);
    res = res && (abs(R(1,2)-R(2,1)) < 1e-6);
    res = res && (det(R) > 0);
    res = res && (trace(R) > 0);
end