function psf = WestPSFMinutes(radius)
% psf = WestPSFMinutes(radius)
%
% Compute Westheimer's PSF function as a function
% of passed radius.  Radius passed in minutes of arc.
%
% 7/11/94		dhb		Added comments, changed name.

psf = 0.952*exp(-2.59*(radius.^1.36)) + ...
       0.048*exp(-2.43*(radius.^1.74));
