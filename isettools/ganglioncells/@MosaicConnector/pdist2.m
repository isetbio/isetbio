function [D,idx] = pdist2(A, B, varargin)
% Calculate pairwise distance between two sets of observations
%
% Syntax:
%   [D,idx] = pdist2(A, B, '', varargin)
%
% Description:
% This is an alternate to the pdist2 native MATLAB function 
% which exists in the Statistics and Machine Learning Toolbox
% Based on the fact that: ||u-v||^2 = ||u||^2 + ||v||^2 - 2*u.v
%
% Inputs:
%    A: [Mobservations x kDimensions] matrix
%    B: [Nobservations x kDimensions] matrix
%
% Outputs:
%    D: [Mobservations x Nobservations] matrix of distance from all 
%        M points in A to all Npoints in B
%    idx: indices of points in B that are closest to points in A
%         (non-empty if 'smallest' is non-empty)
%
% Optional key/value pairs
%   'smallest' : integer.  Return the smallest number points in A that are
%                          closest to each point in B
%   
% Example usage:
%{
    B = [3 3];
    A = [2 2; 2 3; 4 4; 5 5; 3 2.2];
    [D, idx] = pdist2(A,B, '', 'smallest', 1);
%} 

    assert(size(A,2) == 2, 'A matrix must be an N x 2 matrix');
    assert(size(B,2) == 2, 'B matrix must be an M x 2 matrix');
    
    p = inputParser;
    p.addOptional('method', '', @(x)(isempty(x)||(ischar(x))));
    p.addParameter('smallest', [], @(x)(isempty(x) || (isnumeric(x))));
    p.parse(varargin{:});
    smallest = p.Results.smallest;
    method = p.Results.method;

    % Compute all pairwise distances
    D = sqrt( bsxfun(@plus,sum(A.^2,2),sum(B.^2,2)') - 2*(A*B') );

    % Return smallest distances if so desired
    if (~isempty(smallest)) && (smallest > 0)
        dimension = 1;
        [D,idx] = sort(D, dimension, 'ascend');
        smallest = min([smallest size(D,1)]);
        D = D(1:smallest,:);
        idx = idx(1:smallest,:);
    else
        idx = [];
    end
    
end