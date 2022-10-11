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
%
%{
    A = [2 2; 2 3; 4 4; 5 5; 3 2.2];
    B = [3 3];
    [D, idx] = pdist2(A,B, 'smallest', 1);

    A = [2 2; 2 3; 4 4; 5 5; 3 2.2];
    B = [];
    [D, idx] = pdist2(A,B, 'fromPosition', 'maxAbsPosition');

%} 

    assert(size(A,2) == 2, 'A matrix must be an N x 2 matrix');
    assert((isempty(B))||(size(B,2) == 2), 'B matrix must be either empty or an M x 2 matrix');
    
    p = inputParser;
    p.addParameter('smallest', [], @(x)(isempty(x) || (isnumeric(x))));
    p.addParameter('fromPosition', '', @(x) ((ischar(x) && (ismember(x,{'maxAbsPosition'}))) || ...
                                             (isnumeric(x)&&(numel(x) == 2)&&(size(x,2) == 2)))  );
    p.parse(varargin{:});
    smallest = p.Results.smallest;
    fromPosition = p.Results.fromPosition;

    % Compute all pairwise distances
    AA = sum(A.^2,2);
    if (isempty(B))
        if (isempty(fromPosition))
            error('MosaicConnector.pdist2:: if B is empty, the ''fromPosition'' optional argument must be set to a valid value');
        else
            if (ischar(fromPosition))
                switch (fromPosition)
                    case 'maxAbsPosition'
                        [~,idx] = max(AA);
                        B = A(idx,:);
                    otherwise
                        error('MosaicConnector.pdist2:: if B is empty, the ''fromPosition'' optional argument must be set to a valid value. ''%s'' is not a valid value',  fromPosition);
                end
            else
                B = fromPosition;
            end
        end

    end

    BB = sum(B.^2,2);
    D2 = bsxfun(@plus,AA,BB') - 2*(A*B');

    % Return smallest distances if so desired
    if (~isempty(smallest)) && (smallest > 0)
        dimension = 1;
        [D2,idx] = sort(D2, dimension, 'ascend');
        smallest = min([smallest size(D2,1)]);
        D2 = D2(1:smallest,:);
        idx = idx(1:smallest,:);
    else
        idx = [];
    end
    
    D = sqrt(D2);
end