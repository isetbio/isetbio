function k = buildTemporalImpulseResponse(samplingTime)
%% buildTemporalImpulseResponse
% builds the temporal filter function for cells in an rgcMosaic.
% 
%     ir.tCenter{rgbInd,1} = buildTemporalImpulseResponse(samplingTime);
% 
% Create a default (temporal) stimulus filter using the gamma PDF.
% 
% Inputs: the sampling time from the outersegment object.
% 
% Output: a temporal impulse reponse;
% 
% This function incorporates code by Pillow available at 
%       http://pillowlab.princeton.edu/code_GLM.html
% under the GNU General Public License.
% 
% Example: 
%     obj.tCenter{rgbInd,1} = multFactor*buildTemporalImpulseResponse(samplingTime);
% 
% 09/2015 JRG 

%% Build temporal impulse response
% From code by J. Pillow:
% 
% DTsim = .01; % Bin size for simulating model & computing likelihood.
% nkt = 20;  % Number of time bins in filter;
DTsim = samplingTime; % Bin size for simulating model & computing likelihood.
filterLength = 0.2;
nkt = ceil(filterLength/samplingTime);  % Number of time bins in filter;
timeAxis = samplingTime:samplingTime:filterLength;
tk = [0:nkt-1]';
b1 = nkt/32; b2 = nkt/16;
k1 = 1/(gamma(6)*b1)*(tk/b1).^5 .* exp(-tk/b1);  % Gamma pdfn
k2 = 1/(gamma(6)*b2)*(tk/b2).^5 .* exp(-tk/b2);  % Gamma pdf
k = (k1-k2./1.5);
k = 1.2*(k./max(k));
