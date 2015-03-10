function [OTF2D, fSupport,inCutoffFreq] = dlMTF(oi,fSupport,wave,units)
% Gateway routine that assembles parameters to compute diffraction limited
% OTF
% 
% [OTF2D, fSupport,inCutoffFreq] = dlMTF(oi,[fSupport],[wave],[units])
%
% Compute the diffraction limited 2D OTF (OTF2D) at each wavelength for the
% optics within the optical image. 
%
% The diffraction limited OTF only depends on the optics.  But for some
% units and conditions we need to know properties of the optical image to
% perform the calculation.  If you know these parameters, we make it
% possible to call this routine using dlMTF(optics, ...).  This is possible
% because we check at the beginning of the routine to see whether the first
% argument is of type optical image or of type optics.
%
% The frequency support(fSupport) and the incoherent cutoff frequency as a
% function of wavelength (in nm) can also be calculated and returned.  The
% units for the frequency support, cycles/{meters,millimeters,microns}, can
% be specified (units).
% 
% The formulae are described in dlCore.m
%
% Examples:
% If you send in only one argument, it must be the optical image structure.
% In this first example, We calculate all of the OTFs (one for each
% wavelength). 
%
%   OTF2D = dlMTF(oi);
%
% To plot the result, we center DC using fftshift
%
%  figure(1); 
%  mesh(fSupport(1,:,1),fSupport(:,1,2),fftshift(OTF2D(:,:,1))); 
%  colorbar; xlabel('cyc/mm'); ylabel('cyc/mm');
%
% If you have pre-computed the arguments, then you can use the
% optics structure, as in this example
%
%  oi = vcGetObject('oi'); wavelength = 400; unit = 'mm';
%  fSupport = oiGet(oi,'fSupport',unit);
%  OTF2D = dlMTF(optics,fSupport,wavelength,unit);
%
%  figure(1); mesh(fSupport(1,:,1),fSupport(:,1,2),fftshift(OTF2D));
%  colorbar; xlabel('cyc/mm'); ylabel('cyc/mm');
%
% But, if you want units in cycles/deg, you must sent in the optical image.
%
%  See also
%    oitCalculateOTF, dlCore.m
%
% Copyright ImagEval Consultants, LLC, 2003.

if notDefined('oi'), error('Optics or optical image required.'); end

% The user can send in the OI or OPTICS.  We only need OPTICS.
if strcmpi(opticsGet(oi,'type'),'opticalimage')
    optics = oiGet(oi,'optics');
else
    if nargin < 4
        error('All dlMTF arguments must be sent if optics is first');
    end
    optics = oi;
    clear oi
end

if notDefined('wave'),  wave = oiGet(oi,'wavelength'); end
if notDefined('units'), units = 'cyclesPerDegree'; end
if notDefined('fSupport'), fSupport = oiGet(oi, 'fSupport', units); end

apertureDiameter = opticsGet(optics, 'aperturediameter');
fpDistance = opticsGet(optics, 'focalPlaneDistance');

fx = fSupport(:,:,1);
fy = fSupport(:,:,2);

%  Distance of each frequency pair from the origin (rho,theta)
%  (dc = [0,0]).  The frequency support is in cycles/deg. 
rho = sqrt(fx.^2 + fy.^2);

% Wavelength is stored in nanometers.  This converts it to meters, the same
% units as the apertureDimaeter. 
wave = wave *  1e-9;

% This formula assumes that a lens that is free of spherical aberration and
% coma.  When the source is far away, fpDistance is the focal length and
% the ratio (apertureDiameter / fpDistance) is the f#.
%
%  see discussion in dlCore.m  
%
% http://spie.org/x34304.xml  - Cutoff frequency
% http://spie.org/x34468.xml  - Airy disk
%
inCutoffFreq = (apertureDiameter / fpDistance) ./ wave;  

switch lower(units)
    case {'cyclesperdegree','cycperdeg'}
        % cycle/meter * meter/deg -> cycles/deg
        % Used for human calculations?
        inCutoffFreq = inCutoffFreq * ...
                    oiGet(oi, 'distancePerDegree', 'meters');
        
    case {'meters','m','millimeters','mm','microns','um'}
        inCutoffFreq = inCutoffFreq/ieUnitScaleFactor(units);
        
    otherwise
        error('Unknown units.');
end

% Now, both rho and the cutoff frequency are in cycles/degree.
OTF2D = dlCore(rho,inCutoffFreq);

end