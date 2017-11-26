function j = wvfZernikeMNToOSAIndex(n,m)
% Convert from Zernike 2-index standard index to OSA single-zernike index
%
% Syntax:
%   j = wvfZernikeMNToOSAIndex(n,m)
%
% Description:
%    Convert from the Zernike 2 index standard indexing to the OSA single
%    Zernike index (starting at j = 0)
% 
%    Uses equation 4 from the OSA numbering document
%
% Inputs:
%    n - The radial order
%    m - The angular frequency
%
% Outputs:
%    j - Can be a vector, in which case n and m would also be vectors.
%
% Notes:
%    * Validation code s postpended to wvfZernikeNMToOSAIndex
%    * [Note: JNM - Which is the correct name then?]
%
% See Also:
%    wvfZernikeNMToOSAIndex, zernfun
%

% History:
%    07/29/12  dhb  Wrote it.
%    11/09/17  jnm  Formatted

% Examples:
%{
    % NMT naming
    osa_ind = wvfZernikeNMToOSAIndex(0,0)
%}
%{
    % MNT naming
    osa_ind = wvfZernikeMNToOSAIndex(4,6)
%}

% Get j
j = (n .* (n + 2) + m) / 2;
