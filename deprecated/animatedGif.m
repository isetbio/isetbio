function fName = animatedGif(pImg, fName)
% Save an animated gif from a 3D volume of gray scale images
%
%
% In the future, if we want to save RGB image sequences we need to apply
%
%    [M  c_map]= rgb2ind(RGB,256);  % If we have an rgb, do this
%
% to the RGB image before saving it in the gif, below
% 
% See also:  p_Microsoft2013XX
%
% BW,HJ PDCSOFT Team 2013

warning('%s: deprecated, use sensorVisualize', mfilename);

if notDefined('pImg'), error('3D volume of data needed'); end
if notDefined('fName'), fName = fullfile(pwd,'tmp.gif'); end

% Look forever
loops=65535;

% Gray scale
c_map = gray(256);

% Default frame delays
delay1 = .1;
delay2 = .1;

for ii=1:size(pImg,3)  
    M = uint8(pImg(:,:,ii));
    if ii==1
        imwrite(M,c_map,fName,'gif','LoopCount',loops,'DelayTime',delay1)
    end
    imwrite(M,c_map,fName,'gif','WriteMode','append','DelayTime',delay2)
end

fprintf('Finished creating animated gif in file %s\n',fName);

end