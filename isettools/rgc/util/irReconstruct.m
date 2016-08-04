function [stimulusReconstruction, params] = irReconstruct(innerRetina, varargin)
% Reconstructs the stimulus according to the inner retina representation by
% retinal ganglion cell mosaics.
% 
% As of now, this executes only a dumb linear reconstruction where the STRF
% of each cell is added to the stimulus when a spike occurs.
% 
% Inputs: an inner retina object and the method of reconstruction.
% 
% Outputs: a movie showing the reconstruced stimulus from the RGC
% representation.
% 
% Examples:
%   irReconstruct(innerRetina);
%   irReconstruct(innerRetina,'model','linear');
% 
% 2/2016 JRG (c) isetbio team
%
%% Parse input
p = inputParser;
p.addRequired('innerRetina');
p.addParameter('model', 'linear', @ischar);
p.addParameter('tuningWoff', 1, @isnumeric);
p.addParameter('percentDead', 0, @isnumeric);
p.parse(innerRetina, varargin{:});

innerRetina = p.Results.innerRetina;
model = p.Results.model;
tuningWoff = p.Results.tuningWoff;
percentDead = p.Results.percentDead;

%% 


switch model
case{'linear'}

% Initialize figure
% % vcNewGraphWin([],'upperleftbig');
% figure; set(gcf,'position',[160 60 1070 740]);
% hold on;

% TODO: Loop over all mosaics
nX = 0; nY = 0;
cellTypeInd = 1;

[nX,nY,~] = size(innerRetina.mosaic{cellTypeInd}.responseLinear);
nFrames = length(innerRetina.mosaic{cellTypeInd}.responseLinear(1,1,:));
% nX = nX + nXi;
% nY = nY = nYi;

% Find max positions
% allpos = vertcat(innerRetina.mosaic{cellTypeInd}.cellLocation{:});
% maxx = max(allpos(:,1))/1; maxy = max(allpos(:,2));

% metersPerPixel = 

% Find spatial RF size
[nPixX, nPixY] = size(innerRetina.mosaic{cellTypeInd}.sRFcenter{1,1});
nFramesRF = length(innerRetina.mosaic{cellTypeInd}.tCenter{1});


for cellTypeInd = 1:length(innerRetina.mosaic)
    centerCoords = innerRetina.mosaic{cellTypeInd}.cellLocation{1,1};
    mincoordmosaic(cellTypeInd) = abs(ceil(centerCoords(1) - (nPixY/2)));

end
mincoord = max(mincoordmosaic);

stimulusReconstruction = zeros(nPixX*nX +mincoord, nPixY*nY +mincoord, nFrames + nFramesRF);

for cellTypeInd = 1:length(innerRetina.mosaic)
    
    if cellTypeInd == 2 || cellTypeInd == 4
        tuningWeight = 1;%0.01;
    else 
        tuningWeight = tuningWoff;
    end
    
[nY,nX,~] = size(innerRetina.mosaic{cellTypeInd}.responseLinear);

deadIndicesAll = randperm(nX*nY);
numberDead = percentDead*(nX*nY);
deadIndices = deadIndicesAll(1:numberDead);

% [deadRow, deadCol] = ind2sub([nX nY], deadIndices);

maxx = 0; maxy = 0;
cellCtr = 0;
% Loop through each cell and plot spikes over time
for xc = 1:nX
    for yc = 1:nY
        
        cellCtr = cellCtr+1;
        if ~any(cellCtr == deadIndices)
        
        [nPixX, nPixY] = size(innerRetina.mosaic{cellTypeInd}.sRFcenter{1,1});
        nFramesRF = length(innerRetina.mosaic{cellTypeInd}.tCenter{1});
        
        % Build the STRF of the cell
        sRF = innerRetina.mosaic{cellTypeInd}.sRFcenter{xc,yc} - innerRetina.mosaic{cellTypeInd}.sRFsurround{xc,yc};
        tRF(1,1,:) = innerRetina.mosaic{cellTypeInd}.tCenter{1};
        strf = repmat(sRF,[1 1 nFramesRF]).*repmat(tRF, [nPixX, nPixY, 1]);
        
        % Get the appropriate spike data
        spPlot=innerRetina.mosaic{cellTypeInd}.responseSpikes{xc,yc,1,1};
        % spPlot=(median(horzcat(innerRetina.mosaic{3}.spikeResponse{xc,yc,:,2})'));
        if length(spPlot)>0
        % Add the STRF to the stimulus reconstruction for each spike
        for iFrame = 1:length(spPlot)
            
            % ycoords = (yc-1)*nPixY + 1 : yc*nPixY;
            % xcoords = (xc-1)*nPixX + 1 : xc*nPixX;
            
            centerCoords = innerRetina.mosaic{cellTypeInd}.cellLocation{xc,yc};            
            centerCoordsOut(xc,yc,:) = centerCoords;
            ycoords1 = mincoord + (ceil(centerCoords(1) - (nPixY/2)) : floor(centerCoords(1) + (nPixY/2))); 
            xcoords1 = mincoord + (ceil(centerCoords(2) - (nPixX/2)) : floor(centerCoords(2) + (nPixX/2))); 
            
            tcoords = ceil(1*spPlot(iFrame)) : ceil(1*spPlot(iFrame))+nFramesRF-1;
            xcgoodind = xcoords1>0; ycgoodind = ycoords1>0;
            xcoords = xcoords1(xcgoodind); ycoords = ycoords1(ycgoodind);
            stimulusReconstruction(xcoords, ycoords, tcoords) = ...
                stimulusReconstruction(xcoords, ycoords, tcoords) + ...
                tuningWeight*strf(xcgoodind,ycgoodind,:);
        end%iFrame
        maxx = max([maxx xcoords]); maxy = max([maxy ycoords]);
        end
        end
    end%nX
end%nY
end

case{'otherwise'}
    error('Model does not exist');
end

maxR = max(stimulusReconstruction(:));
minR = min(stimulusReconstruction(:));

params.maxx = maxx;
params.maxy = maxy;
params.minR = minR;
params.maxR = maxR;

% % Play the movie
% for iFrame = 1:size(stimulusReconstruction,3)
%     imagesc(stimulusReconstruction(1:maxx,1:maxy,iFrame));
%     colormap gray
%     caxis([minR maxR]);
%     pause(0.1);
% end

