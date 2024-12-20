function meanRate = coneMeanIsomerizations(cMosaicRect, varargin)
% Calculate spatial mean photon rate(R*/sec) for the 3 cone types in mosaic
% 
% Syntax:
%   meanRate = coneMeanIsomerizations(cMosaicRect);
%
% Description:
%    Calculate the spatial mean photon rate (R*/sec by default) for the 3
%    cone types in a mosaic.
%
%    Examples are contained in the code. To access, type 'edit
%    coneMeanIsomerizations.m' into the Command Window.
%
% Inputs:
%    cMosaicRect              - coneMosaicRect object
%
% Outputs:
%    meanRate                 - Three vector absorption rates for the L, M, 
%                               and S cones in R*/sec.
%
% Optional key/value pairs:
%    'perSample'              - Boolean. Normally the returned rate is mean
%                               per sec. Setting 'perSample' to true  makes
%                               the mean rate per temporal sample bins
%                               (default false)
%
%    'absorptionssInXWFormat' - Matrix. If empty (default), works on
%                               absorptions in cMosaic.absorptions. If this
%                               is passed, it acts on what is passed as
%                               this parameter, which is taken to be the
%                               absorptions in XW format. See RGB2XWFormat
%                               for a description of XW format.
%
% Notes:
%    * [NOTE: DHB - I can't figure out just from the code here whether this
%      can act on a time sequence or works just one frame at at time. There
%      is a key comment that says "Reshape from 3D (x, y, t) to space x
%      nCones" that is very confusing, because by the description time on
%      the left has turned into cones on the right. This comment needs to
%      be expanded/ fixed and then this description needs to make clear
%      whether or not a time sequence is involved in the input and output,
%      both for the normal case and for when the absorptions are passed in
%      as a parameter via the keyword absorptionsInXWFormat.]
%
%    * [NOTE: DHB - There is a comment, "Compute means for given
%      integration time", in the code below. I am not sure what this means.
%      I think it means, that we are computing the spatial mean for the
%      numbers in the mosaic, which in turn correspond to a particular
%      integration time. But the comment as written is confusing to me,
%      because I am tempted to interpret it as meaning that an integration
%      time could be passed to this routine, which I don't think it can.
%      Not changing comment because I'm not sure I fully understand what
%      the code is doing.]
%
% See Also:
%   RGB2XWFormat
%

% History:
%    11/xx/16  jrg  (c) Isetbio team
%    08/06/17  dhb  Comment cleaning pass.
%                   Added notes where I could not figure it out.
%    10/26/17  dhb  Reviewed jm changes, added new note, formatted if then
%                   else statements in a way I like better. Accepted some
%                   suggestions from Code Analyzer to remove obsolete
%                   warning supression (I'm running 2017a), and removed
%                   some stray semi-colons after some "end" statements.
%    10/26/17  dhb  Make example work, and add comment to example.
%    02/09/18  jnm  Formatting

% Examples:
%{
   % Create default scene, oi, mosaic and get the mean LMS isomerizations
   % from the mosaic.
   scene = sceneCreate;
   oi = oiCreate('human');
   oi = oiCompute(oi,scene,'pad value','mean');
   cMosaicRect = coneMosaicRect;
   cMosaicRect.compute(oi);
   tmp = coneMeanIsomerizations(cMosaicRect);
%}

%% Validate and parse input parameters
p = inputParser; 
p.addRequired('cMosaic', @(x) isa(x, 'coneMosaicRect'));
p.addParameter('perSample', false, @islogical);
p.addParameter('absorptionsInXWFormat', [], @isnumeric);
p.parse(cMosaicRect, varargin{:});
cMosaicRect = p.Results.cMosaic;
perSample = p.Results.perSample;

%% Default values to return for null case, namely 0 all the way.
lMean = 0;
mMean = 0;
sMean = 0;
meanRate = [lMean, mMean, sMean];

%% Locations of each cone type
coneType = cMosaicRect.pattern;

%% Compute
if (~isempty(p.Results.absorptionsInXWFormat))
    % The rows are spatial position.  The columns are the time points.  So
    % XW is a misnomer here.  It is really XT.
    % The action performed here should really be a 'get' method.
    pRateXW = p.Results.absorptionsInXWFormat;
    nonNullConeIndices = find(cMosaicRect.pattern > 1);
    nonNullConeTypes = coneType(nonNullConeIndices);
    lConeIndices = find(nonNullConeTypes == 2);
    mConeIndices = find(nonNullConeTypes == 3);
    sConeIndices = find(nonNullConeTypes == 4);
else
    % [NOTE: DHB - This comment does not parse for me - how does time turn
    % into cones?]
    % Response to DHB:  This appears to convert (x, y, t) to (space x time)
    % Then we pull out the spatial locations for each cone type separately
    % and average them.
    pRate = cMosaicRect.absorptions;     % Absorptions per sample
    if isempty(pRate), return; end   % Return 0 when no absorptions
    pRateXW = RGB2XWFormat(pRate);
    lConeIndices = find(coneType == 2);
    mConeIndices = find(coneType == 3);
    sConeIndices = find(coneType == 4);
end

%% Get the individual cones
lConeAbsorptions = pRateXW(lConeIndices, :); 
mConeAbsorptions = pRateXW(mConeIndices, :); 
sConeAbsorptions = pRateXW(sConeIndices, :); 

%% Compute means for the given integration time
if ~isempty(lConeAbsorptions)
    lMean = mean(lConeAbsorptions(:));  
else
    lMean = 0;
end
if ~isempty(mConeAbsorptions)
    mMean = mean(mConeAbsorptions(:)); 
else
    mMean = 0;
end
if ~isempty(sConeAbsorptions)
    sMean = mean(sConeAbsorptions(:)); 
else
    sMean = 0;
end

%% Correct so returned units are absorptions per second or per time bin
if (perSample)
    % Return mean number of absorptions per temporal sample
    meanRate = [lMean, mMean, sMean];    
else
    % Convert to R*/sec. This is the default
    meanRate = [lMean, mMean, sMean] / cMosaicRect.integrationTime;
end

end
