function imageSubsample = ieImageSubsample(I, stride, varargin)
% Subsample an image with specific stride
% 
% This can take as input a 2D or 3D matrix, and the subsampling is only
% applied along the first two dimensions.
%
% This function is too simple to be an independent function. We might
% consider remove it in the future. (HJ)
% 
% JRG/HJ/BW, ISETBIO Team, 2016

%% Parse parameters
p = inputParser;
p.addRequired('image',@isnumeric);
p.addRequired('stride',@isnumeric);

p.parse(I,stride,varargin{:});

%%  Set up the sub sampled image for output
imageSubsample = I(1:stride:size(I, 1), 1:stride:size(I, 2), :);

end