function meanRate = coneMeanIsomerizations(cMosaic,varargin)
% Calculate the mean photon rate (R*/sec) for the 3 cone types in a mosaic
% 
% Input
%   cMosaic  - Cone Mosaic object (required)
%
% Parameters
%   perSample - Normally the returned rate is mean per sec.  But setting
%              'perSample',true  makes the mean rate per temporal sample bins
%
% Return
%   meanRate - Three vector of mean rates for the (L,M,S) cones (R*/sec)
% 
% 11/2016 JRG (c) isetbio team

%% Validate input parameters

p = inputParser; 
p.addRequired('cMosaic', @(x) isa(x, 'coneMosaic'));
p.addParameter('perSample',false,@islogical);

p.parse(cMosaic,varargin{:});
cMosaic    = p.Results.cMosaic;
perSample  = p.Results.perSample;

% Default
lMean = 0; mMean = 0; sMean = 0;
meanRate = [lMean, mMean, sMean];

%% Locations of each cone type

coneType = cMosaic.pattern;
pRate    = cMosaic.absorptions;  % Absorptions per sample
if isempty(pRate), return; end   % Return 0 when no absorptions

%% Compute
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

%% Correct so returned units are absorptions per second or per time bin

if perSample
    % Return mean number of absorptions per temporal sample
    meanRate = [lMean, mMean, sMean];    
else
    % Convert to R*/sec.  This is the default
    meanRate = [lMean, mMean, sMean]/cMosaic.integrationTime;
end

end
