function [current, interpFilters] = osCompute(obj, cMosaic, varargin)
% Linear model computing outer segment photocurrent from isomerizations (R*) 
%
%    [current, interpFilters] = osCompute(obj, cMosaic, varargin)
%
% We use  osLinear.osCompute (linear model) for experiments in which there
% is a uniform background, as we often find in psychophysical experiments.
% When the images are more complex (e.g., natural scenes), use the
% osBioPhys model, not the linear model.
%
% Inputs: 
%   obj       - osLinear class object
%   cMosaic   - parent object of the outerSegment
%
% Output:
%   current       - outer segment photocurrent current in pA
%   interpFilters - Interpolated impulse response functions (to integration
%                   time samples)
%
% The linear model impulse response function is the small signal of the
% osBioPhys model. The impulse response depends on the on mean
% isomerization rate.
%
% We calculate the osBioPhys current response to
%
%   * a constant abosrption rate
%   * a constant with 1 photon added to one sampling bin. 
%
% The difference between these two signals is the photocurrent
% impulseResponse to a photon.
%
% The current predicted to an arbitrary stimulus is the current from the
% mean isomerization rate plus the current from small deviations around the
% the mean isomerization rate.
%
%  (mean current) + convolve((absorptions - meanAbsorptions),impulseResponse) 
%
% If the os.noiseFlag is true, the method adds noise to the current output
% signal. See osAddNoise() for specification and validation of the noise
% model. See Angueyra and Rieke (2013, Nature Neuroscience) for details.
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

%% Get the linear filters for the mean rate

% These convert a single photon increment on mean to a photocurrent impulse
% response function
[lmsFilters, meanCur] = obj.linearFilters(cMosaic);

% obj.plot('current filters','meancurrent',meanCur)

%% Interpolate the stored lmsFilters to the time base of the absorptions

absTimeAxis   = cMosaic.timeAxis;
osTimeAxis    = obj.timeAxis;
interpFilters = zeros(cMosaic.tSamples,3);
for ii=1:3
    % Interpolation assumes that we are accounting for the time sample bin
    % width elsewhere.  Also, we extrapolate the filters with zeros to make
    % sure that they extend all the way through the absorption time axis.
    % See the notes in s_matlabConv2.m for an explanation of why.
    interpFilters(:,ii) = interp1(osTimeAxis(:),lmsFilters(:,ii),absTimeAxis(:),'linear',0);
end

%%  The predicted photocurrent is

% Convert (x,y,t) to (space,t)
[absorptions, r, c] = RGB2XWFormat(cMosaic.absorptions);   % Per sample
% vcNewGraphWin; plot(absorptions(100,:));

% We will store the current here
current = zeros(r*c, tSamples);

% convolve the filters with the isomerization data
for ii = 2 : 4  % loop for LMS, cone type 1 is black / blank, so we skip
    
    % locate cones with specific type and convolve with temporal filter
    index = find(coneType==ii);
    if ~isempty(index)
        % The mean absorptions produces a mean current. This current level
        % was returned above when we calculated the lmsFilters (meanCur).
        % 
        % Here we calculate the time-varying photocurrent by  
        %   * Convolving the difference of the absorption rate from the
        %   mean absorption rate with the L, M or S filter
        %   * Adding in the mean background current
        %
        %  conv(absorptions - meanAbsorptions,lmsFilters) + meanCur
        
        % dAbsorptions is [nCones by nTime]
        dAbsorptions = absorptions(index,:) - meanRate(ii-1);
        % The difference should be distributed around 0
        %
        %   vcNewGraphWin; hist(dAbsorptions(:));
        %   mean(dAbsorptions(:))
        
        % Convolve and  add in the mean.  The general conv2 produces a new
        % time series that is longer than nTimes.  We only want the
        % convolution up to the final absorption.  Not really sure if we
        % want 'same' here, or we want circonv, or ... (BW).
        % tmpCurrent = conv2(dAbsorptions,interpFilters(:,ii-1)','same') + meanCur(ii-1);
        tmpCurrent = conv2(interpFilters(:,ii-1)',dAbsorptions) + meanCur(ii-1);
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

end