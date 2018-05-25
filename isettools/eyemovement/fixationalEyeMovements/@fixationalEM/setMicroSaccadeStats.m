function setMicroSaccadeStats(obj)
% Set the microSaccade statistics
%
% Syntax:
%   obj.setMicroSaccadeStats()
%   setMicroSaccadeStats(obj)
%
% Description:
%    Set the stats for the microSaccade. These include both referenced
%    numbers and arbitrary values.
%
%    The following are arbitrary values:
%       MicroSaccade minimum duration (assigned value: 2)
%       MicroSaccade target jitter arcMin (assigned value: 0.3)
%       MicroSaccade direction jitter degrees (assigned value: 15)
%
%    The following are known values:
%       Inter-saccade interval distribution(s)
%
%    The following are based off of references:
%       MicroSaccade amplitude distribution(s) - (Martinez 2009)
%       MicroSaccade speed distribution(s) - (Martinez 2008)
%
% Inputs:
%    obj - Object. The fixationalEM object.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    None.
%
% History:
%    01/03/18  NPC  ISETBIO Team, 2018
%    05/15/18  jnm  Formatting
%    05/24/18  NPC  Comments

% Set the inter-saccde interval distribution
obj.microSaccadeMeanIntervalSeconds = 0.45;
obj.microSaccadeIntervalGammaShapeParameter = 5;

% Set the microsaccade amplitude distribution
% (consistent with Fig 1 in Martinez 2009)
obj.microSaccadeMeanAmplitudeArcMin = 8;
obj.microSaccadeAmplitudeGammaShapeParameter = 4;

% 39 deg/second (See Martinez 2008, Table in Fig 3D)
obj.microSaccadeMeanSpeedDegsPerSecond = 39;
obj.microSaccadeStDevSpeedDegsPerSecond = 2;

% Minimum duration of a microsaccade - Arbitrary value
obj.microSaccadeMinDurationMilliSeconds = 2;

% Sigma of microsaccade target jitter - Arbitrary value
obj.microSaccadeTargetJitterArcMin = 0.3;

% Sigma of corrective microsaccade direction jitter - Arbitrary value
obj.microSaccadeDirectionJitterDegs = 15;

end