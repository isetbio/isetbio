function imageSubsample = ieImageSubsample(image, stride, varargin)
% Usually takes the output of ieSpatioTemporalFilter, and subsamples its
% output. 
% 
% This can take as input a 2D or 3D matrix, and the subsampling is only
% applied along the first two dimensions.
% 
% JRG/HJ/BW 4/2016 (c) isetbio team


%% Parse parameters
p = inputParser;
p.addRequired('image',@isnumeric);
p.addRequired('stride',@isnumeric);

p.parse(image,stride,varargin{:});
image = p.Results.image;
stride = p.Results.stride;

%%  Set up the sub sampled image for output

sz = size(image);
rSamps = 1:stride:sz(1);
cSamps = 1:stride:sz(2);
imageSubsample = zeros(length(rSamps),length(cSamps),sz(3));

%
for k = 1:sz(3)
    imageSubsample(:,:,k) = image(rSamps,cSamps,k);
end

end
