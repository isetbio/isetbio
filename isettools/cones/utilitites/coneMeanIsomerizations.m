function meanRate = coneMeanIsomerizations(cMosaic,varargin)
%%coneMeanIsomerizations  Calculate the spatial mean photon rate (R*/sec) for the 3 cone types in a mosaic
% 
% Syntax:
%    meanRate = coneMeanIsomerizations(cMosaic);
%
% Description:
%    Calculate the spatial mean photon rate (R*/sec by default) for the 3
%    cone types in a mosaic.
%
%    [DHB NOTE: I can't figure out just from the code here whether this can
%    act on a time sequence or works just one frame at at time. There is a
%    key comment that says "Reshape from 3D (x,y,t) to space x nCones" that
%    is very confusing, because by the description time on the left has
%    turned into cones on the right. This comment needs to be
%    expanded/fixed and then this description needs to make clear whether
%    or not a time sequence is involved in the input and output, both for
%    the normal case and for when the absorptions are passed in as a
%    parameter via the keyword absorptionsInXWFormat.]
%    
% Input:
%    cMosaic          coneMosaic object
%
% Output:
%    meanRate         Mean absorption rate in R*/sec.
%
% Optional key/value pairs:
%    perSample                   Normally the returned rate is mean per sec.  Setting
%                                'perSample' to true  makes the mean rate per temporal
%                                sample bins (default false).
%
%    absorptionssInXWFormat      If empty (default), works on absorptions in cMosaic.absorptions.
%                                If this is passed, it acts on what is passed as this parameter,
%                                which is taken to be the absorptions in XW format. See RGB2XWFormat
%                                for a description of XW format.
%              
% Return
%   meanRate - Three vector of mean rates for the (L,M,S) cones (R*/sec)

% 11/2016   JRG   (c) Isetbio team
% 08/06/17  dhb   Comment cleaning pass.  Added notes where I could not figure it out.

%% Validate and parse input parameters
p = inputParser; 
p.addRequired('cMosaic', @(x) isa(x, 'coneMosaic'));
p.addParameter('perSample',false,@islogical);
p.addParameter('absorptionsInXWFormat', [], @isnumeric);
p.parse(cMosaic,varargin{:});
cMosaic    = p.Results.cMosaic;
perSample  = p.Results.perSample;

%% Default values to return for null case, namely 0 all the way.
lMean = 0; mMean = 0; sMean = 0;
meanRate = [lMean, mMean, sMean];

%% Locations of each cone type
coneType = cMosaic.pattern;

%% Compute
if (~isempty(p.Results.absorptionsInXWFormat))
    % [DHB NOTE: A commment here about the format would be great.]
    pRateXW = p.Results.absorptionsInXWFormat;
    nonNullConeIndices = find(cMosaic.pattern > 1);
    nonNullConeTypes = coneType(nonNullConeIndices);
    lConeIndices = find(nonNullConeTypes == 2);
    mConeIndices = find(nonNullConeTypes == 3);
    sConeIndices = find(nonNullConeTypes == 4);
else
    % Reshape from 3D (x,y,t) to space x nCones
    % [DHB NOTE: This comment does not parse for me - how does time turn
    % into cones?]
    pRate = cMosaic.absorptions;             % Absorptions per sample
    if isempty(pRate), return; end           % Return 0 when no absorptions
    pRateXW = RGB2XWFormat(pRate);
    lConeIndices = find(coneType == 2);
    mConeIndices = find(coneType == 3);
    sConeIndices = find(coneType == 4);
end

%% Get the individual cones
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
if (perSample)
    % Return mean number of absorptions per temporal sample
    meanRate = [lMean, mMean, sMean];    
else
    % Convert to R*/sec.  This is the default
    meanRate = [lMean, mMean, sMean]/cMosaic.integrationTime;
end

end
