function ieFontInit(fig)
% Initialize the font size based on the preference setting for ISET
%
%   ieFontInit(fig)
%
% There is a default font size in Matlab when a window opens.  For ISETBIO,
% we store a user-defined preference for a font size that modulates the
% Matlab default size.  This size appears different across platforms.
%
% ISETBIO font size preferences are managed using the preference mechanism
% in Matlab. The ISET preferences are obtained using getpref('ISET'). The
% font size is a field, fontDelta, in the ISET preference structure.
%
% This routine is called when an ISET window is opened. We assume that the
% default Matlab font size is in place, and we apply the preferred change.
%
% The preferred font size change is stored across sessions and windows.
%
% Copyright ImagEval Consultants, LLC, 2003.

try
    fSize = getpref('ISETBIO','fontSize');
catch
    % First time, it may not be there
    fSize = 12;
    setpref('ISETBIO','fontSize',fSize);
end

% fprintf('Setting font size to %d\n',fSize);

% No idea why I have to run this twice.  Must figure out ...
ieFontSizeSet(fig,fSize);
ieFontSizeSet(fig,fSize);

end
