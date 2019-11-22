function optics = opticsSet(optics, parm, val, varargin)
% Set optics structure parameters
%
% Syntax:
%   optics = opticsSet(optics, paramName, val, [varargin])
%
% Description:
%    The optics structure contains the basic optics parameters used to
%    control image formation. The parameters define parameters used in the
%    diffraction-limited or shift-invariant optics models. See opticsGet
%    for further information about these models.
%
%    The optics structure is normally part of the optical image and can be
%    retrieved using
%
%        oi = vcGetObject('OI');
%        optics = oiGet(oi, 'optics');
%
%    Often, we get and set optics properties from the oiSet/Get commands,
%    such as
%
%         fnumber = oiGet(oi, 'optics fnumber')
%         oi = oiSet(oi, 'optics fnumber', 2.8);
%
%    Those calls act via opticsGet and opticsSet
%
%    To set the aperture you must change either the focal length or the
%    f# = fL/aperture, so aperture = fL/f#
%
% Inputs:
%    optics   - Struct. An optics structure
%    parm     - String. The parameter you wish to change. The options and
%               their value types include:
%        Optics model
%           {'model'}         - String. One of the following:
%                               'diffractionLimited' or 'ShiftInvariant'
%        Diffraction limited optics specifications.
%           {'name'}          - String. This optics name
%           {'type'}          - String. Always 'optics'
%           {'fnumber'}       - Numeric. f# is focal length/aperture value
%                               is dimensionless.
%           {'focallength'}   - Numeric. The focal distance in meters for
%                               image at infinity
%           {'transmittance'} - Matrix. Wavelength transmittance  ([0, 1])
%        OTF Information for shift-invariant optics model
%           {'otfdata'}       - Array. Used to store custom data.
%                               Row x Col x Wave
%           {'otffx'}         - Vector. Frequency samples across otfdata
%                               cols (in cyc/mm)
%           {'otffy'}         - Vector. Frequency samples down otfdata rows
%                               (in cyc/mm)
%           {'otfwave'}       - Vector. The otf wavelengths
%        Lens
%           {'lens'}                - Object. For human optics, we store a
%                                     human lens object. The public
%                                     properties are the lens name, wave,
%                                     and density. Other values
%                                     (transmittance, absorptance,
%                                     absorbance) are derived.
%           {'transmittance wave'}  - Vector. For diffraction limited
%                                     optics, we allow the transmittance to
%                                     be aribtrary and we store the
%                                     wavelength samples in
%                                     optics.transmittance.wave.
%           {'transmittance scale'} - Vector. For diffraction limited
%                                     optics, we store the scale value at
%                                     each wavelength
%                                     optics.transmittance.scale.
%        Relative illumination data
%           {'relillummethod'}  - Boolean. A poorly-named offAxis flag.
%                                 (See below).
%           {'off axis method'} - String. Options for this are: 'Skip' to
%                                 turn off or 'cos4th'.
%           {'cos4thdata'}      - Matrix. The cached cos4th data.
%    val      - VARIES. The value to assign to the parameter, following the
%               form specified by the parameter examples above.
%    varargin - (Optional) VARIES. Additional arguments that may be
%               required. Some examples include units.
%
% Outputs:
%    optics   - Struct. The modified optics structure.
%
% Optional key/value pairs:
%    **Needs to be filled out**
%
% See Also:
%    oiGet oiSet opticsGet
%

% History:
%    xx/xx/05       Copyright ImagEval Consultants, LLC, 2005.
%    03/13/18  jnm  Formatting
%    07/25/18  dhb  Fix broken example (lens -> Lens).
%    06/26/19  JNM  Documentation update and minor formatting fixes.

% Examples:
%{
    oi = oiCreate('diffraction limited');
    oiGet(oi, 'optics fnumber')
    oi = oiSet(oi, 'optics fnumber', 8);
    oiGet(oi, 'optics fnumber')
%}
%{
    % Lens transmittance plot for human case
    l = Lens;
    vcNewGraphWin;
    hold on;
    plot(l.wave, l.transmittance);
    l.density = 0.3;
    oi = oiCreate;
    oiSet(oi, 'optics lens', l);
    l = oiGet(oi, 'lens');
    plot(l.wave, l.transmittance);
%}
%{
    oi = oiCreate('diffraction limited');
    oi = oiSet(oi, 'optics transmittance wave', 400:10:700);
    oi = oiSet(oi, 'optics transmittance scale', ones(31, 1));
%}

%%
if ~exist('optics', 'var') || isempty(optics)
  error('No optics specified.');
end
if ~exist('parm', 'var') || isempty(parm)
    error('No parameter specified.');
end
if ~exist('val', 'var'), error('No value.'); end

%%
parm = ieParamFormat(parm);
switch parm
    case 'name'
        optics.name = val;

    case 'type'
        % Should always be 'optics'
        if ~strcmp(val, 'optics')
            warning('Non standard optics type setting');
        end
        optics.type = val;

    case {'model', 'opticsmodel'}
        % Set the optics model type
        %
        % diffractionlimited and shiftinvariant are the legitimate
        % options. We allow 'raytrace' this is mainly in service of
        % TL's comment below.
        %
        % TL: I put back the ray trace model for now, since if we switch to
        % shiftinvariant it will try to find the lens info in the oi.optics
        % structure, which is non-existent for ray-tracing. Because of this
        % it will throw errors when you try to run oiWindow. This requires
        % more thought about how to handle this...
        %
        % BW: If we need a special case for sceneeye, let's make that.
        % ray trace is used for a different meaning in ISETCAM.  We
        % should use a different name for the scene eye calculation.

        % Remove white space and force lower case
        val = ieParamFormat(val);
        valid = {'diffractionlimited', 'shiftinvariant', ...
            'iset3d', 'raytrace'};

        if validatestring(val, valid)
            optics.model = ieParamFormat(val);
        else
            error('Invalid model %s\n', val);
        end

    case {'fnumber', 'f#'}
        optics.fNumber = val;

    case {'focallength', 'flength'}
        optics.focalLength = val;

    case {'spectrum'}
        % Spectrum structure
        warning('optics spectrum set, line 92')
        optics.spectrum = val;

    case {'wavelength', 'wave'}
        % This appears to be unnecessary. In fact, the whole
        % optics.spectrum slot may be unnecessary.
        %
        % We used to change the OTF at the same time. But this is not
        % necessary because the OTF structure has a wave term, and when we
        % need the OTF(w) we simply interpolate it from the OTF structure.
        %
        % I am not sure we ever use this particular wave for SI data, where
        % we use the OTF.wave and the oi.spectrum.wave. This one is kind
        % of caught in the middle. ISETBIO doesn't use rt, so the other
        % thing to check is whether it is used for diffraction.

        % Set new wavelength
        warning('optics spectrum wave set, line 110')
        optics.spectrum.wave = val(:);

    case {'lens'}
        % New lens object. This should replace transmittance.
        optics.lens = val;

    case {'transmittancewave'}
        % For human optics, there is a lens object but no
        % transmittance slot.  For diffraction limited optics,
        % however, there is a transmittance slot you can set.

        % The transmittance wave must be set before the scale is set.
        if isfield(optics, 'transmittance')
            optics.transmittance.wave = val;
        else
            fprintf('Set transmittance only for diffraction');
            fprintf('Set lens density for human case')
            fprintf('Optics name: %s', optics.name);
        end

    case {'transmittancescale'}
        % The wave must be set before you set the transmittance scale
        if isfield(optics, 'transmittance')
            if numel(val) ~= numel(optics.transmittance.wave)
                error(strcat('The transmittance scale does not match ', ...
                    'the number of wave slots'));
            end
            optics.transmittance.scale = val;
        else
            fprintf('Set transmittance only for diffraction');
            fprintf('Set lens density for human case')
            fprintf('Optics name: %s', optics.name);
        end

    % ---- Relative illumination calculations
    case {'offaxis', 'offaxismethod', 'relillummethod', 'cos4thflag'}
        % Flag determining whether you use the cos4th method
        % Bad naming because of history.
        optics.offaxis = val;

    case {'cos4thfunction', 'cos4thmethod'}
        % We only have cos4th offaxis implemented, and this probably is all
        % we will need.
        optics.cos4th.function = val;

    case {'cos4th', 'cos4thdata', 'cos4thvalue'}
        % Numerical values. Should change field to data from value.
        optics.cos4th.value = val;

    % ---- OTF information for shift-invariant calculations
    case {'otffunction', 'otfmethod'}
        % This should probably not be here.
        % We should probably only be using the 'model' option
        % But it is used, so we need to carefully debug
        % Current choices are 'dlmtf' ... - MP, BW
        optics.OTF.function = val;

    case {'otf', 'otfdata'}
        % Fraction of amplitude transmitted
        optics.OTF.OTF = val;

    case {'otffx'}
        % Units are cyc/mm
        %   - frequency samples across col of otfdata
        %   - these seem to have 0 at entry floor(n / 2) + 1, while the
        %     actual otf is not shifted this way.
        optics.OTF.fx = val;

    case {'otffy'}
        % Units are cyc/mm
        %   - frequency samples down rows of otfdata
        optics.OTF.fy = val;

    case {'otfwave'}
        % - otf wavelengths (nm)
        optics.OTF.wave = val;

    otherwise
        error('Unknown parameter')
end

end