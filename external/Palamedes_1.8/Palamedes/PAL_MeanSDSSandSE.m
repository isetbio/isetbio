%
%PAL_MeanSDSSandSE   Descriptive statistics
%
%   syntax: [Mean SD SS SE] = PAL_MeanSDSSandSE(x)
%
%   Given vector or matrix 'x', [Mean SD SS SE] = PAL_MeanSDSSandSE(x) 
%   returns the mean, standard deviation (using n - 1, i.e., the estimate 
%   of standard deviation of a population from which sample 'x' was drawn), 
%   the sum of squared deviations from mean, and estimate of the standard 
%   error of the mean (again assuming 'x' to be sample from population).
%
%   In case 'x' is a matrix, function will return row vectors with
%   statistics calculated across columns of 'x'.
%
%   Example: [Mean SD SS SE] = PAL_MeanSDSSandSE([1 1 2; 1 3 4]) returns:
%
%   Mean =  1     2     3
%
%   SD =    0    1.4142    1.4142
%
%   SS =    0     2     2
%
%   SE =    0     1     1
%
%Introduced: Palamedes version 1.0.0 (NP)
%Modified: Palamedes version 1.6.3 (see History.m)

function [Mean, SD, SS, SE] = PAL_MeanSDSSandSE(x)

if size(x,1) == 1
    x = x';
end

SS = zeros(1,size(x,2));
SD = zeros(1,size(x,2));
SE = zeros(1,size(x,2));

Mean = sum(x,1)/size(x,1);
for i = 1:size(x,2)
    SS(i) = (sum((x(:,i) - Mean(i)).^2));
    SD(i) = (SS(i)./(size(x,1)-1)).^.5;
    SE(i) = sqrt((SS(i)./(size(x,1)-1)))/sqrt(size(x,1));
end

if isempty(x)
    Mean = [];
    SS = [];
    SD = [];
    SE = [];
end