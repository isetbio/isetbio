function [irf, timeAxis] = currentIRF(meanIso,ecc, timeSamples)
% Return the photocurrent impulse response function
%
% Under development (see BW)
% 
% Synopsis
%   [irf, timeAxis] = currentIRF(meanIso,ecc,timeSamples)
%
% Inputs
%   meanIso     -  Mean isomerizations per temporal bin
%   ecc         - (Optional) Visual field eccentricity (default = 0)
%   timeSamples - (Optional) Temporal samples in seconds.
%
% Outputs
%   irf - Impulse response function, one for each mean ISO level
%   timeAxis - Sample times (secs)
%
% Description
%   Calls the biophysical model of the Rieke photocurrent.  Returns the
%   small signal impulse response function as a function of the
%   isomerization rate and the eccentricity.
%
%   Many things to check, such as why the variable depends on the rate per
%   bin rather than R*/sec.
%
% See also
%   @osLinear, t_osLinearize
%

% Examples:
%{
  mn = [10,50,100];
  [irf, timeAxis] = currentIRF(mn,15);
  ieNewGraphWin;
  plot(timeAxis,irf); xlabel('Time (s)');
  grid on; ylabel('Current (pA)');
  legend(strsplit(num2str(mn)))
%}
%{
  mn = 50;
  timeSamples = (0:0.002:0.4);
  irf = currentIRF(mn,15,timeSamples);
  ieNewGraphWin;
  plot(timeSamples,irf,'-o'); xlabel('Time (s)');
  grid on; ylabel('Current (pA)');
  legend(strsplit(num2str(mn)))
%}

%% Dummy coneMosaic object

% This enables us to access the outersegment methods.
cMosaic = coneMosaic('os', osLinear, 'pattern', [2 2 2]);

if ~exist('ecc','var'), ecc = 0; end

%% Allocate space for the impulse responses

% A litte ugly, we get an impulse response and find its size.
nTimeSamples = size(cMosaic.os.linearFilters(cMosaic), 1);

% The IRFs in the columns
irf = zeros(nTimeSamples, length(meanIso));

%%  Loop on different background rates and and compute
os = cMosaic.os;
for ii = 1:length(meanIso)
    % Set mean background for start of current output
    cMosaic.absorptions = repmat(meanIso(ii), 1, 3);

    % Compute outer segment currents for the fovea. Do this by pulling out
    % the first column of return, corresponding to the L cones.
    tmp = os.linearFilters(cMosaic, 'eccentricity', ecc);
    irf(:, ii) = tmp(:, 1);

end

timeAxis = cMosaic.os.timeAxis;

if exist('timeSamples','var')
    % Interpolate from current time axis to the provided time samples
    irf = interp1(timeAxis(:),irf,timeSamples);
end
    
end
