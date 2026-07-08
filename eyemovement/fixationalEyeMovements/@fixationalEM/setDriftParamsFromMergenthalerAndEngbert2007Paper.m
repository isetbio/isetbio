function setDriftParamsFromMergenthalerAndEngbert2007Paper(obj)
% Use the aforementioned paper as a reference to set drift parameters
%
% Syntax:
%   setDriftParamsFromMergenthalerAndEngbert2007Paper(obj)
%   obj.setDriftParamsFromMergenthalerAndEngbert2007Paper()
%
% Description:
%    Use the aforementioned paper as a reference to set drift parameters
%
% Inputs:
%    obj - Object. A fixationalEM object.
%
% Outputs:
%    None.
%
% Optional key/value pairs:
%    None.
%
% Notes:
%    * All references are from "Mergenthaler & Engbert", 2007
%

% Control gain (gamma)
obj.controlGamma = 0.25;

% Control noise, ksi, with zero MEAN and sigma STD
ksi = struct('mean', 0, 'sigma', 0.075);
obj.controlNoiseMean = ksi.mean;
obj.controlNoiseStd = ksi.sigma;

% Position noise, eta with zero MEAN and rho STD
eta = struct('mean', 0, 'rho', 0.35);
obj.positionNoiseMean = eta.mean;
obj.positionNoiseStd = eta.rho;

%Feedback gain (lambda)
obj.feedbackGain = 0.15;

% Feedback steepness (epsilon)
obj.feedbackSteepness = 1.1;

% Feedback delay for x- and y-pos
obj.feedbackXposDelaySeconds = 70 / 1000;  % tau_hor
obj.feedbackYposDelaySeconds = 40 / 1000;  % tau_vert

end