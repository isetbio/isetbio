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

%%

imageFrames= size(image,3);

switch imageFrames
    case 1        
        imageSubsample = image(1:stride:end,1:stride:end);        
    otherwise        
        for k = 1:imageFrames        
            imageSubsample(:,:,k) = image(1:stride:end,1:stride:end,k);           
        end
end

