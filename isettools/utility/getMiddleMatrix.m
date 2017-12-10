function middleM = getMiddleMatrix(m, sz)
% Extract values near middle of a matrix.
%
% Syntax:
%   middleM = getMiddleMatrix(m, sz)
%
% Description:
%    Data values from the middle of a matrix are returned. The total number
%    of extracted pixels is 1 + round(sz / 2) * 2.  This is awkward  for
%    small numbers of sz, but OK for bigger numbers.
%
% Inputs:
%    m       - The matrix to extract values from
%    sz      - The size of the middle to extract
%
% Outputs:
%    middleM - The requested middle data
% 
% Notes:
%    * [Note: HJ/BW - Changed from round to floor on March 19, 2015. We are
%      worried, but it had no effect on validationfastall and the routine
%      is only called in a few places where this should be OK. Delete this
%      comment after a few months.]
%    * [Note: JNM - 1) Is it safe to delete the prior note? 2) The change
%      to floor seems to only be in one of two instances? Is this half-
%      change accidental or intentional?]
%

% History:
%    xx/xx/03        Copyright ImagEval Consultants, LLC, 2003.
%    03/19/15  HJ/BW Changed from round to floor. No effect on
%                    validationfastall. Should be ok?
%    12/01/17  jnm   Formatting
%

% Examples:
%{
    m = reshape([1:(9 * 9)], 9, 9);
    foo = getMiddleMatrix(m, 3)
%}
%{
    m = reshape([1:(9 * 9 * 3)], 9, 9, 3);
    foo = getMiddleMatrix(m, 5)
%}

% Check that the variables have been provided
if notDefined('m'), error('Matrix must be provided.'); end
if notDefined('sz'), error('Size must be provided.'); end

sz = floor(sz / 2);

center = round(size(m) / 2);
rMin = max(1, center(1) - sz);
rMax = min(size(m, 1), center(1) + sz);
cMin = max(1, center(2) - sz);
cMax = min(size(m, 2), center(2) + sz);

r = rMin : rMax;
c = cMin : cMax;
% w = size(m, 3);
middleM = m(r, c, :);

end
