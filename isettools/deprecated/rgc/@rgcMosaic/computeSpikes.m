function [rgcM, nTrialsSpikeResponseM] = computeSpikes(rgcM, varargin)
% Generate RGC spikes from the linear mosaic response
%
% Syntax:
%   [rgcM, nTrialsSpikeResponseM] = computeSpikes(rgcM, [varargin])
%
% Description:
%    Calculate the spikes from the RGC linear response for rgcGLM and
%    rgcLNP mosaics.
%
% Inputs:
%    rgcM          - Object. A rgc mosaic object.
%
% Outputs:
%    rgcM          - Object. The modified rgcM object.
%    nTrialsSpikeResponseM
%                  - Matrix. The spike response matrix for each trial.
%
% Optional key/value pairs:
%    coupling      - Boolean. Whether to couple. Default true.
%    nTrialsLinear - Numeric. The number of linear trials. Default [].
%
% See Also:
%   s_vaRGC.m in WL/WLVernierAcuity
%

% History:
%    XX/XX/15  BW   (c) isetbio team, 2015
%    06/18/19  JNM  Documentation pass

%% Parse
p = inputParser;
p.addRequired('rgcM', @(x)(isa(x, 'rgcMosaic')));  % Inner retina object
p.addParameter('coupling', true, @islogical);
p.addParameter('nTrialsLinear', [], @isnumeric);
p.parse(rgcM, varargin{:});

coupling = p.Results.coupling;
nTrialsLinear = p.Results.nTrialsLinear;

%% Required for Pillow code
% The refresh rate refers to the subsampling of the linear response. Each
% sample of the linear response is broken into N bins, where N is the
% refresh rate. This is used in simGLM.m and simGLMcpl.m from Pillow.
% BW - This should be fixed.

global RefreshRate
RefreshRate = 10;

% For every IR, this could be a vector in the future
% nRepeats = rgcM.get('number trials');

nRepeats = 1;  % Not implemented, but will be used
tt = 1;        % Will be used when we have repeats

%% Loop on the mosaics in the inner retina
if ~isempty(nTrialsLinear), nTrials = size(nTrialsLinear, 1);
else, nTrials = 1;
end

for iTrial = 1:nTrials
    % nTrialsSpikeResponse = cell(length(rgcM.mosaic), 1);

    % We will get rid of this and just use rgcM throughout.
    mosaic = rgcM;
    if ~isempty(nTrialsLinear)
        responseLinear = squeeze(nTrialsLinear(iTrial, :, :, :));
    else
        responseLinear = rgcM.get('response linear');
    end

    nSamples = size(responseLinear, 3);
    nCells = mosaic.get('mosaic samples');

    % This is a vector of times that each cell spiked
    spikeTimes = cell([nCells, nRepeats]);

    % Temporal sample of the voltage response
    respVolts = zeros(nCells(1), nCells(2), ...
        RefreshRate * nSamples, nRepeats);

    glminput = RGB2XWFormat(responseLinear);

    if coupling  % Only for GLM case
        % Call the Pillow code to generate spikes for the whole mosaic
        % using the coupled GLM
        glmprs = setGLMprs(mosaic, 'coupling', coupling);

        if ieSessionGet('wait bar')
            wbar = waitbar(0, sprintf('Trial 1  of %d', nRepeats));
        end

        % Run Pillow code
        if ieSessionGet('wait bar')
            wbar = waitbar(0, sprintf('Coupled GLM (Trial %d of %d)', ...
                tt, nRepeats));
        end
        [responseSpikesVec, Vmem] = simGLMcpl(glmprs, glminput');

        % Put the data in the outputs
        cellCtr = 0;
        for xc = 1:nCells(2)
            for yc = 1:nCells(1)
                cellCtr = cellCtr+1;
                % Vector of times when the cell spiked
                spikeTimes{yc, xc, tt} = responseSpikesVec{1, cellCtr};
                respVolts(yc, xc, :, tt) = Vmem(:, cellCtr);
            end
        end
        if ieSessionGet('wait bar')
            waitbar(tt / nRepeats, wbar, ...
                sprintf('Finished trial %d of %d', tt, nRepeats));
        end

        if ieSessionGet('wait bar'), delete(wbar); end
    else
        % No coupling, much faster

        % Run without the slow coupling component by looping on the
        % simGLM, not simGLMcp
        glmprs = setGLMprs(mosaic, 'coupling', coupling);
        cellCtr = 0;   % Reset the cell counter

        for xc = 1:nCells(2)
            for yc = 1:nCells(1)
                cellCtr = cellCtr + 1;
                % Pull out linear response of a cell
                thisCellInput = glminput(cellCtr, :);

                % Run Pillow code
                [responseSpikesVec, Vmem] = simGLM(glmprs, thisCellInput');

                % Vector of times when the cell spiked
                spikeTimes{yc, xc, tt} = responseSpikesVec;

                % Voltages as a function of time
                respVolts(yc, xc, :, tt) = Vmem;
            end
        end
    end

    % Set mosaic spike times
    rgcM.set('responseSpikes', spikeTimes);

    % The nonlinear voltage; only set in the GLM model
    if isa(rgcM, 'rgcGLM'), rgcM.set('response voltage', respVolts); end

    if ~isempty(nTrialsLinear)
        if iTrial == 1 % && ii == 1
            nTrialsSpikeResponse = ...
                zeros([nTrials, size(rgcM.responseSpikes)]);
        end
        spikesTemp = rgcM.get('spikes');
        nTrialsSpikeResponseM(iTrial, :, :, 1:size(spikesTemp, 3)) = ...
            spikesTemp;
    else
        nTrialsSpikeResponseM = [];
    end

end

end
