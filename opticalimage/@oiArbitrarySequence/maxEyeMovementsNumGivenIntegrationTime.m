function maxEyeMovementsNum = maxEyeMovementsNumGivenIntegrationTime(...
    obj, integrationTime, varargin)
% Generate eye movement sequence for all oi's
%
% Syntax:
%   maxEyeMovementsNum = maxEyeMovementsNumGivenIntegrationTime(...
%       obj, integrationTime)
%
% Description:
%    Method to compute the maximum number of eye movements given an
%    oiSequence and a mosaic's intergration time. The optional parameter
%    stimulusSamplingInterval' is only used when the oiSequence has a
%    length of 1, in which case we do not have information regrading the
%    duration of that single frame.
%
% Inputs:
%    obj                - Object. An OI sequence object.
%    integrationTime    - Numeric. The integration time.
%    varargin           - (Optional). VARIES. Additional parameters in
%                         key/value fashion as needed.
%
% Outputs:
%    maxEyeMovementsNum - Numeric. The maximum number of eye movements.
%
% Optional key/value pairs:
%    None required.
%

%% Parse arguments
p = inputParser;
p.addParameter('stimulusSamplingInterval', [], @isnumeric);
p.parse(varargin{:});

% Generate eye movement sequence for all oi's
if numel(obj.timeAxis) == 1
    if ~isempty(p.Results.stimulusSamplingInterval)
        stimulusSamplingInterval = p.Results.stimulusSamplingInterval;
    else
        % No information about what the stimulus sampling interval is, so
        % arbitrarily set it to the integrationTime
        stimulusSamplingInterval = integrationTime;
    end
else
    stimulusSamplingInterval = obj.timeAxis(2) - obj.timeAxis(1);
end

eyeMovementsNumPerOpticalImage = stimulusSamplingInterval / ...
    integrationTime;
maxEyeMovementsNum = ...
    round(eyeMovementsNumPerOpticalImage * obj.length);

if maxEyeMovementsNum < 1
    error(['Less than 1 eye movement per oiSequence !!! \nStimulus' ...
        ' sampling interval:%g ms Cone mosaic integration time: ' ...
        '%g ms\n'], 1000 * stimulusSamplingInterval, ...
        1000 * integrationTime);
else
    % fprintf(['Optical image sequence contains %2.0f eye ' ...
    %    'movements (%2.2f eye movements/oi)\n'], ...
    %    maxEyeMovementsNum, eyeMovementsNumPerOpticalImage);
end
end