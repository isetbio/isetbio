function current = osCompute(obj, cMosaic, varargin)
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


%% parse inputs
p = inputParser; 
p.KeepUnmatched = true;
p.addRequired('obj', @(x) isa(x, 'outerSegment'));
p.addRequired('cMosaic', @(x) isa(x, 'coneMosaic'));

p.parse(obj,cMosaic);

% p.addRequired('pRate', @isnumeric);
% p.addRequired('coneType', @ismatrix);  % Comes from coneMosaic parent

coneType   = cMosaic.pattern;
meanRate   = coneMeanIsomerizations(cMosaic,'perSample',true);  % R*/sample
tSamples   = size(cMosaic.absorptions,3);

%% This is the place where we get the linear filters given the mean rate

% These convert a single photon increment on mean to a photocurrent impulse
% response function
[lmsFilters, meanCur] = obj.linearFilters(cMosaic);
% obj.plot('current filters','meancurrent',meanCur)

%%  The predicted photocurrent is

% Convert (x,y,t) to (space,t)
[absorptions, r, c] = RGB2XWFormat(cMosaic.absorptions);   % Per sample
% vcNewGraphWin; plot(absorptions(100,:));

% We will store the current here
current = zeros(r*c, tSamples);

% convolve the filters with the isomerization data
for ii = 2 : 4  % loop for LMS, cone type 1 is black / blank
    
    % locate cones with specific type and convolve with temporal filter
    index = find(coneType==ii);
    if ~isempty(index)
        % The mean absorptions produces a mean current that was returned
        % above when we calculated the lmsFilters (meanCur).
        % 
        % We calculate the total photocurrent by convolving the difference
        % of the absorption rate from the mean absorption rate and the LM
        % or S filter.  We then add in the mean background current
        %
        %  conv(absorptions - meanAbsorptions,lmsFilters) + meanCur
        
        dAbsorptions = absorptions(index,:) - meanRate(ii-1);
        % The difference should be distributed around 0
        %
        % vcNewGraphWin; hist(dAbsorptions(:));
        % mean(dAbsorptions(:))
        %
        % vcNewGraphWin; plot(dAbsorptions(10,:))
        
        % Convolve and then add in the mean
        tmpCurrent = conv2(dAbsorptions,lmsFilters(:,ii-1)') + meanCur(ii-1);
        % vcNewGraphWin; plot(tmpCurrent(10,:))
        
        % Store it
        current(index,:) = tmpCurrent(:,1:tSamples);
        
    end
end

% record only recent history in obj
% newStart = max(size(pRate, 3) - length(filter) + 1, 1);
% pRate = pRate(:, :, newStart:end);

% Reshape the current back into (x,y,t) format
current = XW2RGBFormat(current, r, c);

% Add noise or not.

% The osAddNoise function expects the input to be current and it needs to
% know the time sampling.  Hence, we send that in, too.
if osGet(obj,'noiseFlag') == 1
    disp('Current noise added.')
    current = osAddNoise(current, 'sampTime',obj.timeStep);
else
    disp('No current noise added.')
end


obj.coneCurrentSignal = current;

end