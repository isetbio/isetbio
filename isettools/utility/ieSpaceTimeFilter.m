function resp = ieSpaceTimeFilter(sig,kernel,varargin)
% Apply a kernel to the space-time signal
%
%    resp = 
%
% We use this to calculate either just the spatial pooling or
% spatial-temporal pooling for a signal.  
%
% For example, motion analysis can be a kernel that is slanted in
% space-time. Pure spatial is a kernel that is 2D.
%
% TODO: Think about when we can use the separable case.
%
%  sig = randn(128,128,64);
%  kernel = ones(3,3);
%  resp = ieSpaceTimeFilter(sig,kernel,'pad','same');
%  var(resp(:)) % Should be close to 9, except for border
%
% JRG/HJ/BW Isetbio Team, 2016

%% Parse parameters
p = inputParser;
p.addRequired('sig',@isnumeric);
p.addRequired('kernel',@isnumeric);

% Parameter/Value options
p.addParameter('pad','same',@ischar);

p.parse(sig,kernel,varargin{:});
sig = p.Results.sig;
kernel = p.Results.kernel;
pad = p.Results.pad;

%% Apply convn using the 'same' flag 

resp = convn(sig,kernel,pad);

end
