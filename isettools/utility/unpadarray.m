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
%    Examples in code.
%
% Inputs:
%    inArray   - Input array to unpad
%    unpadSize - Can be a single integer (rows) or a 1x2 vector (for rows &
%                columns) - the number by which you unpad both sides
%
% Outputs:
%    outArray  - the unpadded array
%
% Optional key/value pairs:
%    None.
%
% Notes:
%
% See also: padarray

% History:
%    11/20/17  jnm  Formatting, example & heading text
%    12/26/17   BW  more tests in example
%    01/17/18  jnm  Formatting update to match Wiki.

% Examples:
%{
    test = rand(128, 128);
    j = padarray(test, [2, 2]);
    test2 = unpadarray(j, [2, 2]);
    isequal(test, test2)

    pd = [3, 4];
    j = padarray(test, pd);
    test2 = unpadarray(j, pd);
    isequal(test, test2)

    pd = [9, 9];
    j = padarray(test, pd);
    test2 = unpadarray(j, pd);
    isequal(test, test2)
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

