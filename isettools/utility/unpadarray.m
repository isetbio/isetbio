function outArray = unpadarray(inArray, unpadSize)
% Inversion of padarray
%
% Syntax:
%   outArray = unpadarray(inArray, unpadSize)
%
% Description:
%    This function is an inversion of the function padarray. Unpads evenly
%    from both side of the array. 
%
% Inputs:
%    inArray   - Input array to unpad
%    unpadSize - Can be a single integer (rows) or a 1x2 vector (for rows &
%                columns) - the number by which you unpad both sides
%
% Outputs:
%    outArray  - the unpadded array
%
% Notes:
%    * [Note: XXX - Inverts padarray. Not much tested yet. Work on it.]
%    * [Note: XXX - Roughly an inverse. Tested quite thoroughly. Hmmm.]
%

% History:
%    11/20/17  jnm  Formatting, example & heading text

% Examples:
%{
    j = padarray([1 2; 3 4], [2,2]);
    unpadarray(j, [2,2])
%}

if notDefined('inArray'), error('Input array required'); end
if notDefined('unpadSize'), error('unpad size required'); end

if length(unpadSize) < 2, unpadSize(2) = 0; end

r = size(inArray, 1);
c = size(inArray, 2);

rows = (unpadSize(1) + 1):(r - unpadSize(1));
cols = (unpadSize(2) + 1):(c - unpadSize(2));

if ndims(inArray) == 3
    outArray = inArray(rows, cols, :);
else
    outArray = inArray(rows, cols);
end

end

