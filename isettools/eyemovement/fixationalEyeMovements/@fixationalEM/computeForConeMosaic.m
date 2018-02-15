function computeForConeMosaic(obj, coneMosaic, eyeMovementsPerTrial, varargin)
    p = inputParser;
    p.addRequired('coneMosaic', @(x)(isa(x, 'coneMosaic')));
    p.addRequired('eyeMovementsPerTrial', @isscalar);
    p.addParameter('nTrials',1,@isscalar);
    p.addParameter('computeVelocity', false, @islogical);
    p.addParameter('rSeed', [], @isscalar);
    p.addParameter('useParfor', false, @islogical);
    p.parse(coneMosaic, eyeMovementsPerTrial, varargin{:});
    
    % Set optional parameters based on input
    obj.randomSeed = p.Results.rSeed;
    nTrials = p.Results.nTrials;
    computeVelocitySignal = p.Results.computeVelocity;
    useParfor = p.Results.useParfor;
    
    % Extract sampleDuration from coneMosaic's integration time
    sampleDurationSeconds = coneMosaic.integrationTime;
    
    % Compute emDurationSeconds
    emDurationSeconds = eyeMovementsPerTrial * sampleDurationSeconds;
    
    % Generate fixational eye movements (obj.emPathsArcMin)
    compute(obj, emDurationSeconds, sampleDurationSeconds, nTrials, ...
        computeVelocitySignal, 'useParfor', useParfor);
    
    % Subsample to cone positions
    conePatternSampleMicrons = coneMosaic.patternSampleSize(1)*1e6;
    emPosMicronsNotRounded = (obj.emPosArcMin/60)*coneMosaic.micronsPerDegree;
    % Return the paths in units of cone mosaic pattern size
    obj.emPos = sign(emPosMicronsNotRounded) .* round(abs(emPosMicronsNotRounded)/conePatternSampleMicrons);
    % Also return the path in units of microns
    obj.emPosMicrons = obj.emPos * conePatternSampleMicrons;
    obj.emPosArcMin = obj.emPosMicrons / coneMosaic.micronsPerDegree * 60;
    
    if (computeVelocitySignal)
        iTrial = 1;
        velocityArcMinTrial1 = obj.computeVelocity(squeeze(obj.emPosArcMin(iTrial,:,:)));
        allTrialsVelocityArcMin = zeros(nTrials, length(velocityArcMinTrial1), 1);
        allTrialsVelocityArcMin(iTrial,:) = velocityArcMinTrial1;
    
        allTrialsEMPosArcMin = obj.emPosArcMin;
        if (useParfor)
            parfor iTrial = 2:nTrials
                allTrialsVelocityArcMin(iTrial,:) = fixationalEM.computeVelocity(squeeze(allTrialsEMPosArcMin(iTrial,:,:)));
            end
        else
            for iTrial = 2:nTrials
                allTrialsVelocityArcMin(iTrial,:) = obj.computeVelocity(squeeze(allTrialsEMPosArcMin(iTrial,:,:)));
            end
        end
        obj.velocityArcMin = allTrialsVelocityArcMin;
    end
end