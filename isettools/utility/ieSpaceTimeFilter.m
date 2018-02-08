function resp = ieSpaceTimeFilter(sig, kernel, varargin)
% Apply a filter kernel to the space-time signal
%
% Syntax:
%   resp = ieSpaceTimeFilter(sig, kernel, [varargin])
%
% Description:
%    We use this to calculate either just the spatial pooling or
%    spatial-temporal pooling for a signal. If the kernel is 2D, say the
%    spatial receptive field, the the calculation is space-only. If the
%    kernel is 3D, then the calculation is space-time.
%
%    For example, motion analysis can be a kernel that is slanted in
%    space-time. Pure spatial is a kernel that is 2D.
%
% Inputs:
%    sig      - Signal
%    kernel   - Two or three-dimensional field indicating whether the
%               calculation is in space-only, or space-time.
%
% Outputs:
%    resp     - the requested filtered signal
%
% Optional key/value pairs
%    pad      - The pad parameter for the convolution call. Options are
%               'full', 'same' or 'valid'
%
% Notes:
%    * TODO: Think about when we can use the separable case.
%

% History:
%    xx/xx/16  JRG/HJ/BW  Isetbio Team, 2016
%    11/21/17  jnm  Formatting

% Examples:
%{
    sig = randn(128, 128, 64);
	kernel = ones(3, 3);
	resp = ieSpaceTimeFilter(sig, kernel, 'pad', 'same');
	var(resp(:))
    % Should be close to 9 (numel(kernel)), except for border
%}

%% Parse parameters
p = inputParser;
p.addRequired('sig', @isnumeric);
p.addRequired('kernel', @isnumeric);

% Parameter/Value options
vFunc = @(x)(ismember(x, {'full', 'same', 'valid'}));
p.addParameter('pad', 'same', vFunc);
p.parse(sig, kernel, varargin{:});
pad = p.Results.pad;

%% Apply convn typically using the pad = 'same' flag 
resp = convn(sig, kernel, pad);

end
