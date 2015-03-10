function vignetting = pixelVignetting(ISA,OI)
%Compute pixel vignetting across the image sensor array
%
%   vignetting = pixelVignetting(ISA,OI)
%
% PURPOSE:	
%   Calculate the correction factor for the irradiance of off-axis pixels
%   due to pixel vignetting.  The routine calculates the vignetting at 100
%   samples from the image center to the periphery along the X-axis.  THese
%   values are interpolated and then circularly symmetric 
%   function across the rest of the array.  Calculating every pixel takes
%   a long time, and this function is invariably smooth.
%
% Copyright ImagEval Consultants, LLC, 2003.
OPTICS = oiGet(OI,'optics');
PIXEL = sensorGet(ISA,'pixel');
% Setting up local variables
f = opticsGet(OPTICS,'focallength'); 	    % Focal Length [m]
[overlap,numGrid,numGridSpot] = pvFullOverlap(OPTICS,PIXEL);
% fullOverlapForPixelVignetting seems to be doing the right thing (PC)
% mesh(overlap);
[pixelCentersX,pixelCentersY] = sensorPixelCoord(ISA,'upper-right');
% Subsample to speed things up.  Perhaps we should put the sub-sampling
% option in sensorPixelCoord?
pixelCentersX = pixelCentersX(1:4:end,1:4:end);
pixelCentersY = pixelCentersY(1:4:end,1:4:end);
urCol = length(pixelCentersX); urRow = length(pixelCentersY);
[pixelCentersMeshX,pixelCentersMeshY] = meshgrid(pixelCentersX,pixelCentersY);
pixelCenters = [pixelCentersMeshX(:) pixelCentersMeshY(:)];
% 1st order
%%%%%%%%%%%%%
diodeAngle(:,1) = atan(pixelCentersMeshX(:)/f); % [rad]
diodeAngle(:,2) = atan(pixelCentersMeshY(:)/f); % [rad]
% This routines seems to work. We plot the angles for the main diagonal pixels
% plot(diodeAngle(find(abs(pixelCentersMeshX(:)) == abs(pixelCentersMeshY(:)))))
% 2nd order (We'll do this later)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% n = 1.5
% h = pixelGet(PIXEL,'depth');                % Distance from surface to diode [m]
% f = 180 - ( (f-h) * tan(diodeAngle(:,1)) + h * tan(asin(sin(diodeAngle(:,1)/n)));
% Compute the ratio of the collected to the incident irradiance.
% It appears to me that this function could be significantly faster if we
% find a way to avoid the looping, and if we compute just along, say, the
% center and then copy the data in a circularly symmetric function (which
% this is) or we just make it a 1D plot, which is just as informative.
wbar = waitbar(0,'Computing sensor vignetting...');
nPixels = length(pixelCentersMeshX(:));
    
for ii=1:nPixels      
   ratio(ii) = pvReduction(overlap,numGrid,numGridSpot,pixelCenters(ii,:),diodeAngle(ii,:),OPTICS,PIXEL);
   waitbar(ii/nPixels,wbar);
end
close(wbar);
nRows = sensorGet(ISA,'rows'); 
nCols = sensorGet(ISA,'cols');
% Don't ask me why.
ratio = rot90(reshape(ratio,urRow,urCol),1);
vignetting = upperQuad2FullMatrix(ratio,nRows,nCols);
end
