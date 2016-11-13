function current = osCompute(obj, cMosaic, varargin)
% Compute the  outer segment photocurrents using the linear model
%
%    current = osCompute(obj, cMosaic, varargin)
%
% Inputs: 
%   obj       - osLinear class object
%   cMosaic   - parent object of the outersegment
%
% Output:
%   current  - outer segment photocurrent current in pA
% 
% We use  osLinear.osCompute (linear model) for experiments in which there
% is a uniform background, as we typically find in psychophysical
% experiments. When the images are more complex (e.g., natural scenes), use
% the osBioPhys model, not the linear model.
%
% The impulse response function of the linear model calculated here matches
% the small signal from the osBioPhys model, assuming the observer is in a
% steady-state on the mean background.
%
% The impulse response is calculated like this. We calculate the (base
% current) returned by the osBioPhys() model to a uniform stimulus.  We
% then calculate the response when we add 1 photon to the mean level for
% one sampling bin. This difference between these two signals is the
% photocurrent impulseResponse to a photon.
%
% To compute the current from a general stimulus we calculate the
% differences from the mean, convolve with the impulse response, and then
% add in the base current.
%
%   (absorptions - meanAbsorptions) * impulseResponse + (base current)
% 
% This converts isomerizations (R*) to outer segment current (pA).
%
% If the os.noiseFlag is true, this method adds noise to the current output
% signal. 
%
% See Angueyra and Rieke (2013, Nature Neuroscience) for details.  And see
% osAddNoise() for specification and validation of the noise model
%
% JRG/HJ/BW, ISETBIO TEAM, 2016

%% parse inputs
p = inputParser; 
p.KeepUnmatched = true;
p.addRequired('obj', @(x) isa(x, 'outerSegment'));
p.addRequired('cMosaic', @(x) isa(x, 'coneMosaic'));

p.parse(obj,cMosaic);

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

% Reshape the current back into (x,y,t) format
current = XW2RGBFormat(current, r, c);

% Noise anyone?
if osGet(obj,'noiseFlag') == 1
    % The osAddNoise function expects the input to be current and it needs to
    % know the time sampling.
    disp('Current noise added.')
    current = osAddNoise(current, 'sampTime',obj.timeStep);
else
    disp('No current noise added.')
end

obj.coneCurrentSignal = current;

end