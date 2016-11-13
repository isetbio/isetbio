function meanRate = coneMeanIsomerizations(cMosaic,varargin)
% Calculate the mean photon rate (R*/sec) for the 3 cone types in a mosaic
% 
% Input
%   cMosaic  - Cone Mosaic object (required)
%
% Return
%   meanRate - Three vector of mean rates for the (L,M,S) cones (R*/sec)
% 
% 11/2016 JRG (c) isetbio team

%% We use the absorptions in the cMosaic
%
p = inputParser; 
p.addRequired('cMosaic', @(x) isa(x, 'coneMosaic'));
p.parse(cMosaic,varargin{:});
cMosaic = p.Results.cMosaic;

coneType = cMosaic.pattern;
pRate    = cMosaic.absorptions;

%% Locations of each cone type

lConeIndices = find(coneType == 2);
mConeIndices = find(coneType == 3);
sConeIndices = find(coneType == 4);

% Reshape from 3D (x,y,t) to space x nCones
pRateXW = RGB2XWFormat(pRate);

% Get the individual cones
lConeAbsorptions = pRateXW(lConeIndices,:); %#ok<FNDSB>
mConeAbsorptions = pRateXW(mConeIndices,:); %#ok<FNDSB>
sConeAbsorptions = pRateXW(sConeIndices,:); %#ok<FNDSB>

%% Compute means for the given integration time

if ~isempty(lConeAbsorptions), lMean = mean(lConeAbsorptions(:));  
else                           lMean = 0;
end;

if ~isempty(mConeAbsorptions), mMean = mean(mConeAbsorptions(:)); 
else                           mMean = 0;
end;

if ~isempty(sConeAbsorptions), sMean = mean(sConeAbsorptions(:)); 
else                           sMean = 0;
end

%% Correct so returned units are absorptions per second

meanRate = [lMean, mMean, sMean]/cMosaic.integrationTime;

end
