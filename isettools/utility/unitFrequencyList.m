function fList = unitFrequencyList(N)
% Calculate a vector of normalized frequencies for an N-vector
%
% Syntax:
%   fList = unitFrequencyList(N)
%
% Description:
%    The list is calculated given an input that that is N-samples long.
%    This routine handles the case when N is even and odd.
%
%    The range that comes back, say when N = 100 or N = 101, the DC term is
%    at location 51. The general rule is that the center point is at
%    floor(n/2)+1.
% 
%    The main purpose is to get the zero (DC) term into the proper position
%    where Matlab expects it. Then, once we know the maximum frequency
%    (Nyquist) in our measurement, we can multiply the returned list here
%    times that maximum frequency.
%
%    Examples in the code.
%
% Inputs:
%    N     - Number of samples to make a list for.
%
% Outputs:
%    fList - Frequency List
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    11/20/17  jnm  Formatting
%    12/12/17  dhb  Expand on comments, make example work.
%    12/21/17  dhb  Slicker finding of center, floor(N/2)+1.  Works
%                   for odd and even, gives same answer as what was here.
%    01/17/18  jnm  Formatting update to match Wiki.

% Examples:
%{
    nyquistFrequency = 10;
    fList = unitFrequencyList(50);  
    fList = unitFrequencyList(51);
    dataFrequencies = fList * nyquistFrequency
%}

if notDefined('N'), N = 6; end

% Figure out which entry represents DC
mid = floor(N / 2) + 1;

% Here is a list of frequencies
c = 1:N;

% Subtract so that the one at the mid position is zero
c = c - c(mid);

% Normalize so that the largest value is +/-1
fList = c / max(abs(c(:)));

end