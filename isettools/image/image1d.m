function [img, params] = image1d(pattern,varargin)
% Make an image from a 1D pattern.
%
%    img = image1d(pattern,varargin)
%
% Parameters
%  pattern:   1D pattern
%  varargin:  param/val
%     rows
%     rgb
%     
% Create a 2D image from a 1D pattern.
%
% Examples
%   pattern = ones(1,129)*0.5;  pattern(65) = 1;
%   img = image1d(pattern); 
%   vcNewGraphWin; imshow(img)
%
%   x = (-63:64)/128; f = 2;
%   pattern = 0.5*cos(2*pi*f*x) + 0.5;
%   [img,p] = image1d(pattern,'rows',256); imshow(img)
%   
%   [img,p] = image1d(pattern,'rgb',[.2 .2 1]); imshow(img)
%
% BW ISETBIO Team, Copyright 2015

% Is there a matlab parameter parser we should be using?

% Init parameters
params.rgb = [1 1 1];
params.rows = length(pattern);
params.pad  = [0 0];   % Extend the pattern

for ii=1:2:length(varargin)
    thisP = ieParamFormat(varargin{ii});
    switch thisP
        case 'rgb'
            params.rgb = varargin{ii+1};
        case 'rows'
            params.rows = varargin{ii+1};
        otherwise
            error('Unknown parameter %s\n',thisP);
    end
end

%% Make the image

% 1D image
img = repmat(pattern(:)',[params.rows,1]);

% Extend the pattern here?
% Look at params.pad

% Color it in RGB space.  
% Perhaps we should color the ccontrast, and not the mean?
% We could pull out the mean, multiply the contrast, and then add them
% together.
[img, r, c] = RGB2XWFormat(img);
img = img*params.rgb(:)';
img = XW2RGBFormat(img,r,c);

end
