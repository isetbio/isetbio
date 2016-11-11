function current = osCompute(obj, pRate, coneType, varargin)
% function current = osCompute(obj, pRate, coneType, varargin)
% Compute the response of the outer segments using the linear model 
%
%    current = osCompute(obj, pRate, coneType, varargin)
%
% We use the linear model for cases in which there is a uniform background,
% as we typically find in psychophysical experiments.  When the images are
% more complex (e.g., natural scenes) you should use the osBioPhys model,
% not the linear model.
%
% The impulse response function of the linear model calculated here matches
% the linear osBioPhys model given that the observer is in a steady-state
% on the mean background (specified by the three cone pRate values).
%
% The basic idea is this.  When we have a mean field with a specific photon
% absorption rate, we get the base current returned from the osBioPhys()
% model.  We then get the temporal impulse response function that we expect
% when we add 1 photon to the mean level.  The predicted current, then, is
%
%   (absorptions - meanAbsorptions) * impulseResponse + (base current)
% 
% This converts isomerizations (R*) to outer segment current (pA). If the
% noiseFlag is set to true, this method adds noise to the current output
% signal. See Angueyra and Rieke (2013, Nature Neuroscience) for details.
%
%
% Inputs: 
%   obj       - osLinear class object
%   pRate     - R*/sec (x,y,t)
%   coneType  - coneMosaic.pattern 
%
% Optional input (key-value pairs)
%   append   - logical, compute new or to append to existing. When append
%              is true, the computed current is appended to the existing
%              photocurrent data in the object. The returned current value
%              only contains photocurrent for input pRate.
% 
% Outputs:
%   current  - outer segment current in pA
% 
% JRG/HJ/BW, ISETBIO TEAM, 2016

% We are hoping to switch to the arguments: [obj, cMosaic, varargin]
% When we do, this should be the parsing
%
% p = inputParser; 
% p.KeepUnmatched = true;
% p.addRequired('obj', @(x) isa(x, 'osLinear'));
% p.addRequired('cMosaic', @(x) isa(x, 'coneMosaic'));
% % To remove and write a script to check model compatibility with Juan's stuff
% % p.addParameter('linearized', true, @islogical);
% 
% p.parse(obj,cMosaic,varargin{:});
% 
% pRate = cMosaic.absorptions/cMosaic.integrationTime;
% coneType = cMosaic.pattern;


% parse inputs
p = inputParser; 
p.KeepUnmatched = true;
p.addRequired('obj', @(x) isa(x, 'osLinear'));
p.addRequired('pRate', @isnumeric);
p.addRequired('coneType', @ismatrix);  % Comes from coneMosaic parent

p.parse(obj, pRate, coneType, varargin{:});
coneType   = p.Results.coneType;
pRate      = p.Results.pRate;

nHistFrames = 0;


%% This is the place where we get the linear filters given the mean rate

% meanRate units are isomerizatons per section (R*/sec)
meanRate = coneMeanIsomerizations('pRate',pRate,'coneType',coneType);

% These convert a single photon increment on mean to a photocurrent IRF
[lmsFilters, meanCur] = obj.generateBioPhysFilters(meanRate, varargin{:});
% vcNewGraphWin; plot(obj.lmsConeFilter)

%% We need to ask whether we can't put this into the units

% This converts the mean absorption rate into the base current
% maxCur = 0.01*20.5^3; % Dunn & Rieke (2007, Nature)
% meanCur = maxCur * (1 - 1./(1 + 45000./meanRate));

%%  The predicted photocurrent is
%

% First entry is trial.  We are showing only the first trial here.
[absorptions, r, c] = RGB2XWFormat(pRate);

current = zeros(r*c, size(pRate, 3));

% convolve the filters with the isomerization data
for ii = 2 : 4  % loop for LMS, cone type 1 is black / blank
    % pull out the linear filter for current cone type.
    filter = lmsFilters(:, ii-1);
    
    % locate cones with specific type and convolve with temporal filter
    index = find(coneType==ii);
    if ~isempty(index)
        % The mean absorptions produces the meanCur
        % The differences from the mean produce the impulse response
        % So the formula convolves the differences from the mean with the
        % filter and then adds the mean background current
        %
        %  conv(absorptions - meanAbsorptions,lmsFilters) + baseCurrent
        
        % We are nervous about the size.
        dAbsorptions = absorptions(index,:) - meanRate(ii);
        current(index,:)  = conv2(dAbsorptions,lmsFilters(:,ii-1)','same') + meanCur(ii-1);
        
        % curData = conv2(absorptions(index, :), filter') - meanCur(ii-1);
        % current(index, :) = curData(:, nHistFrames+1+(1:size(pRate, 3)));
    end
end

% record only recent history in obj
newStart = max(size(pRate, 3) - length(filter) + 1, 1);
pRate = pRate(:, :, newStart:end);

% reshape the output signal matrix.
current = XW2RGBFormat(current, r, c);

% Add noise
% The osAddNoise function expects and input to be isomerization rate.
% This is handled properly because the params has the time sampling
% rate included.
if osGet(obj,'noiseFlag') == 1
    disp('Current noise added.')
    params.sampTime = obj.timeStep;
    current = osAddNoise(current, params);
else
    disp('No current noise added.')
end


obj.coneCurrentSignal = current;

end