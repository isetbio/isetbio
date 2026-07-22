function [LMS, factorsLMS] = MacBoynToLMS(lsY,T_cones,T_lum) 
% [LMS, factorsLMS] = MacBoynToLMS(lsY,T_cones,T_lum)
%
% Convert Macleod-Boynton chromaticity together with luminance to
% cone coordinates.
%
% This is designed to be used with the modern form of LMSTOMacBoyn,
% namely
%   lsY = LMSToMacBoyn(LMS2,T_cones,T_Y,1);
% as this form returns both ls and the corresponding luminance.  Passing
% both cone fundamentals and luminance sensitivity (T_cones and T_Y) keeps
% the routine general.  The scaling of s is also consistent when the calls
% are done in this way.
%
% See usage example in LMSToMacBoyn.

% History
%   06/xx/25 fh   Wrote it
%   06/16/25 dhb  Adjusted LMSToMacBoyn and this to mesh together more
%                 smoothly.
    
    % Magic call into LMSToMacBoyn to get the scaling factors
    % we need to invert.  See LMSToMacBoyn.
    [~,factorsLMS] = LMSToMacBoyn([],T_cones,T_lum);

    % Once we have the scaling factors, the inversion is pretty simple,
    % because we have the luminance returned by LMSToMacBoyn as well as the
    % ls coordinates.
    LMS = zeros(size(lsY));
    LMS(1,:) = lsY(1,:).*lsY(3,:)/factorsLMS(1);
    LMS(2,:) = (1-lsY(1,:)).*lsY(3,:)/factorsLMS(2);
    LMS(3,:) = lsY(2,:).*lsY(3,:)/factorsLMS(3);
   
end