function lsf = WestLSFMinutes(distance)
% lsf = WestLSFMinutes(distance)
%
% Compute Westheimer's LSF function as a function
% of passed distance.  Dadius passed in minutes of arc.
%
% 9/4/97  dhb  Wrote it.

lsf = 0.47*exp(-3.3*(distance.^2)) + ...
       0.53*exp(-0.93*(abs(distance)));
