function [val,timeSamples] = meanIsomerizations(obj,excitations)
% Isomerizations (R*) per second
%
% Inputs
%  obj - cMosaic 
%
% Returns
%   val - mean isomerizations per sec (R*/sec) in L,M,S cones
%   timeSamples - Times in sec
%
% See also

% excitations: nInstances by nTimepoints by cone

% Mean across instances, time and across all the cones of that type
tmp = squeeze(mean(excitations(:,:,obj.lConeIndices),1));
val(1) = mean(tmp(:));

tmp = squeeze(mean(excitations(:,:,obj.mConeIndices),1));
val(2) = mean(tmp(:));

tmp = squeeze(mean(excitations(:,:,obj.sConeIndices),1));
val(3) = mean(tmp(:));

% Per second
val = val/obj.integrationTime;

timeSamples = (1:size(excitations,2))*obj.integrationTime;

end