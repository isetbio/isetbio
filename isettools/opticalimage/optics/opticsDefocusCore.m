function [otf, sampleSFmm] = opticsDefocusCore(optics,  sampleSF, D)
%Compute the optical transfer function for dioptric power D0 and defocus D 
%
%   [otf, sampleSFmm] = opticsDefocusCore(optics,sampleSF,D)
%
% A computation for a defocused OTF .
%
% Inputs
%  optics:    Optics structure
%  sampleSF:  Spatial frequencies in cycles/deg
%  D:         Defocus in diopters for each wavelength
%
%Returns
% otf:  Optical transfer function (actually, this is the MTF, just a set of
%       scale factors.  We assume there is no frequency-dependent phase
%       shift.
% sampleSFmm:  Sample spatial frequency in cyc/millimeters
%
% See also: humanCore (this routine derived from that), opticsDefocusedMTF,
%           defocusMTF, s_opticsDefocus
%
% Example:
%   REWRITE:
%   sampleSF = 1:60;
%   D = zeros(size(wave));             %No wavelength dependent defocus
%   D = 0.1*[1:length(wave)]/length(wave);   % A little.
%   otf = opticsDefocusCore(XXX);
%   vcNewGraphWin; mesh(sampleSF,wave,otf)
%
% Copyright ImagEval Consultants, LLC, 2005.

% Check parameters
if notDefined('optics'), error('Require optics'); end
if notDefined('sampleSF'), error('Spatial frequency needed'); end
if notDefined('D'), error('Defocus vector needed'); end

p  = opticsGet(optics, 'pupil radius', 'm');
D0 = opticsGet(optics, 'diopters', 'm');

% Converts the defocus in diopters to the Hopkins w20 parameter for a given
% pupil radius in meters, defocus (D, diopters), and dioptric power (D0).
% The explanation for this formula is in Marimont and Wandell, Appendix C:
% Converting from w20 to Defocus in diopters
w20 = (p^2/2)*(D0.*D)./(D0+D);
% plot(wave,w20); 

% Re-write so we can get sampleSF in cycles/mm directly without these two
% steps.
c = opticsGet(optics,'deg per dist', 'm');
% 1/(atand(1) * (1/D0));  %  deg per meter (rad/meter)

% The units are: 
% cycles/meter = (cycles/deg) * (deg/meter) 
cSF = sampleSF * c;

% The formulae in the opticsDefocusedMTF appears to have a problem handling
% the SF=0 value. At SF=0 it uses the classic diffraction formula and at
% other values it uses the formula from Marimont and Wandell.  THat formula
% seems to be off by a small scale factor that we don't understand.
% By doing this, we get a smooth OTF that doesn't have a value of 1 at DC.
% In the code we force a DC value of 1 - but that is a hack.  We should
% figure out what's up with the formula.  Also remember that the paper by
% Subbaro claims that the scale factor on the formula is wrong.  So there's
% two of out there.
% If we have SF=0, we replace it with a very small number.
ii = (cSF == 0);
cSF(ii) = min(cSF(~ii))*1e-12;

% Note: When D0 = 60, as for human, the number is:
%    c = 3434.07;
% This logic is repeated in the humanCore routine

lambda = opticsGet(optics,'wave','m');     % Wavelength in meters
s      = zeros(length(lambda),length(sampleSF));
alpha  = zeros(size(s));
otf    = zeros(size(s));

for ii = 1:length(lambda)
    
    % We should probably convert sampleSF to sampleSFmm above, and then get
    % rid of the 'c' parameter.
    %
    % Appendix B from Marimont and Wandell
    % Compute the reduced spatial frequency (0,2)
    %            m *        (m/m) *    cy/m  - Dimensionless in the end
    s(ii, :) = (lambda(ii) /(D0*p)) * cSF; 
    
    % Methods:  Marimont and Wandell
    % Related to the defocus specified by w20, which in turn depends on p,
    % D and D0.
    alpha(ii,:) = (4*pi./(lambda(ii))).*w20(ii).*abs(s(ii,:));
    
    % We put the vector of sample SF into this array.
    % Then we interpolate to the full 2D array outside of this loop.
    otf(ii,:) = opticsDefocusedMTF(s(ii,:),abs(alpha(ii,:)));
    % plot(otf(ii,:))
end

% Some parameters are out of range and yield complex values. We force them
% to 0 here. This is done implicitly for the diffraction limited case by
% using the incoherent cutoff frequency.

% Complex values to zero.
l = (angle(otf) ~= 0);
otf(l) = 0;

% Convert to cyc/mm, which is used in opticsGet/Set
% We convert (1/D0)*1000 to make it in millimeters, rather than meters.
degPerMillimeter = (1/(tand(1) * (1/D0)*1000));
%             cyc/deg * deg/mm
sampleSFmm = sampleSF * degPerMillimeter;

end