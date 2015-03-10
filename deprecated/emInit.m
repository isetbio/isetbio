function [sensor, params] = emInit(emType, sensor, params)
% Init eye movement parameters in the sensor structure
%
%   [sensor, params] = emInit(emType, sensor, params)
%
% emType:  Eye movement type {fixation brownian, fixation gaussian}
% sensor:  The sensor
% params:  Depends on type
%   fixation brownian:
%      center, sdx, sdy, center, nSamples, randSeed (optional).  
%      The sd units are in deg of visual angle.
%
% General Process:
%   1. Check eyeMoveType and set random seed
%   2. Generate moving postion x and y
%   3. Generate frames per position according to distribution and nSamples
%   4. Generate linear eye movement for testing
%
%
% Output Parameter:
%   sensor       - sensor with eye movement related parameters set
%
% Example:
%    sensor = sensorCreate('human'); clear params
%    params.center = 0;
%    params.Sigma = 0.02; 
%    params.nSamples = 50;
%    params.fov = 1;  
%    sensor = emInit('fixation brownian', sensor, params)
%
% (HJ) Copyright PDCSOFT TEAM 2013

%% Check inputs and Init
warning('This function is deprecated. Use emCreate, eyemoveInit instead');
if notDefined('emType'), error('eye movement type required'); end
if notDefined('sensor'), error('sensor required.'); end
if notDefined('params'), error('parameters required.'); end

% Initialize random number generation
if ~isfield(params, 'center'), params.center = [0 0]; end
if ~isfield(params, 'Sigma')
    error('Covariance matrix for eye movement missing');
end
if ~isfield(params, 'nSamples'), params.nSamples = 1000; end
if ~isfield(params, 'fov'), error('Field of view in params required'); end

if isfield(params,'randSeed') % Set up random seed
    rng(randSeed);
else % Store the newly generated random seed
    params.randSeed = rng('shuffle');
end

if ~isfield(params, 'samplesPerTrial')
    params.samplesPerTrial = params.nSamples;
end

%% Initialize eye movements
% Each case builds the (x,y) and count variables for every position
emType = ieParamFormat(emType);
switch emType
    case {'fixationbrownian'}
        % The eye wanders around the center. The positions are random in a
        % disk around the center. the distances here are in deg of visual
        % angle. 
        % Eye movement is simulated by Brownian motion with increment of
        % Gaussian N(0, Sigma). The initial value is given by params.center
        
        % Compute sensor fov and size
        sz  = sensorGet(sensor, 'size');
        
        % Generate gaussian move
        pos  = ieMvnrnd(params.center, params.Sigma, params.nSamples);
        
        % Adding up to form brownian motion
        % It's a little tricky here. If using brownian motion, we could
        % goes to a very fall distance with high probability. So, we should
        % make it bounce back to center if it gets too large
        for cPos = 1 : params.samplesPerTrial : params.nSamples
            indx = cPos:min(cPos + param.samplesPerTrial, params.nSamples);
            pos(indx,:) = cumsum(pos(indx,:) - ...
                repmat(params.center, [length(indx),1]))+...
                repmat(params.center, [length(indx),1]); % in degree
        end   
        % Now we need to set the bounce back. When human eye gets too far
        % away from the fixation point, it will jump back to the fixation
        % point
        boundV = min(params.fov / 2, 0.6);
        indx = find(pos > boundV | pos < -boundV, 1, 'first');
        while ~isempty(indx)
            if pos(indx) > boundV
                pos(indx:end) = pos(indx:end) - pos(indx) + ...
                    mod(pos(indx), boundV);
            else
                pos(indx:end) = pos(indx:end) - pos(indx) + ...
                    mod(pos(indx), -boundV);
            end
            indx = find(pos > boundV | pos < -boundV, 1, 'first');
        end
        
        % For efficiency, we round the calculations
        % that are centered less than 1 detector's width.
        pos  = round((pos/params.fov).*repmat(sz,[params.nSamples,1])).*...
            params.fov ./ repmat(sz, [params.nSamples, 1]);
        
        % Group the same positions
        [pos,~,ic]  = unique(pos,'rows');
        
        % Compute frame per position
        f    = hist(ic,unique(ic));      % frames per position
        f(1) = f(1) + params.nSamples - sum(f); % make sure sum(f)=nSamples
    case 'fixationgaussian'
        sz  = sensorGet(sensor, 'size');
        
        % Generate gaussian move
        pos  = ieMvnrnd(params.center, params.Sigma, params.nSamples);
        
        % For efficiency, we round the calculations
        % that are centered less than 1 detector's width.
        pos  = round((pos/params.fov).*repmat(sz,[params.nSamples,1])).*...
            params.fov ./ repmat(sz, [params.nSamples, 1]);
        
        % Group the same positions
        [pos,~,ic]  = unique(pos,'rows');
        
        % Compute frame per position
        f    = hist(ic,unique(ic));      % frames per position
        f(1) = f(1) + params.nSamples - sum(f); % make sure sum(f)=nSamples
    otherwise
        error('Unknown emType %s\n',emType);
end

% Set sensor movement positions.
sensor = sensorSet(sensor,'movement positions', pos);
sensor = sensorSet(sensor,'frames per position',f);

end