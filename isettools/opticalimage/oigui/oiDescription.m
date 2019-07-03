function txt = oiDescription(oi)
% Generate optical image text description for window
%
% Syntax:
%   txt = oiDescription(oi)
%
% Description:
%    Writes the OI description for the box on the upper right of the
%    window. Manages cases of diffraction-limited, shift-invariant, and ray
%    trace separately.
%
% Inputs:
%    oi  - Struct. An optical image structure.
%
% Outputs:
%    txt - String. The generated optical image's description.
%
% Optional key/value pairs:
%    None.
%

% History:
%    xx/xx/03       Copyright ImagEval Consultants, LLC, 2003.
%    03/19/18  jnm  Formatting
%    07/01/19  JNM  Formatting update

txt = sprintf('\nOptical image\n');

if isempty(oi)
    txt = 'No image'; return;
else
    sz = oiGet(oi, 'size');
    if isempty(oiGet(oi, 'photons'))
        txt = addText(txt, sprintf('  No image\n'));
    else
        str = sprintf('  Size:       [%.0f, %.0f] samples\n', ...
            sz(1), sz(2));
        txt = addText(txt, str);

        u = round(log10(oiGet(oi, 'height', 'm')));
        if (u >= 0 )
            str = sprintf('  Hgt, wdth: [%.2f, %.2f] m\n', ...
                oiGet(oi, 'height', 'm'), oiGet(oi, 'width', 'm'));
        elseif (u >= -3)
            str = sprintf('  Hgt, wdth: [%.2f, %.2f] mm\n', ...
                oiGet(oi, 'height', 'mm'), oiGet(oi, 'width', 'mm'));
        else
            str = sprintf('  Hgt, wdth: [%.2f, %.2f] um\n', ...
                oiGet(oi, 'height', 'um'), oiGet(oi, 'width', 'um'));
        end
        txt = addText(txt, str);

        u = round(log10(oiGet(oi, 'sampleSize')));
        if (u >= 0 )
            str = sprintf('  Sample:  %.2f  m\n', ...
                oiGet(oi, 'sampleSize', 'm'));
        elseif (u >= -3)
            str = sprintf('  Sample:  %.2f mm\n', ...
                oiGet(oi, 'sampleSize', 'mm'));
        else
            str = sprintf('  Sample:  %.2f um\n', ...
                oiGet(oi, 'sampleSize', 'um'));
        end
        txt = addText(txt, str);

        wave = oiGet(oi, 'wave');
        spacing = oiGet(oi, 'binwidth');
        str = sprintf('  Wave:     %.0d:%.0d:%.0d nm\n', ...
            min(wave(:)), spacing, max(wave(:)));
        txt = addText(txt, str);

        meanIll = oiGet(oi, 'meanilluminance');
        if ~isempty(meanIll)
            txt = addText(txt, sprintf('  Illum:       %.1f lux\n', ...
                meanIll));
        end

        fov = oiGet(oi, 'hfov');
        txt = addText(txt, sprintf('  FOV (wide):    %.1f deg\n', fov));
    end
end

% Write out the optics parameters, either DL, SI or RT
optics = oiGet(oi, 'optics');
opticsModel = opticsGet(optics, 'model');

switch lower(opticsModel)
    case {'diffractionlimited', 'dlmtf'}
        txt = [txt, sprintf('Optics (DL)\n')];
        txt = [txt, sprintf('  Mag:  %.2e\n', ...
            opticsGet(optics, 'magnification'))];
        diameter = opticsGet(optics, 'aperturediameter', 'mm');
        txt = [txt, sprintf('  Diameter:  %.2f mm\n', diameter)];
    case 'shiftinvariant'
        txt = [txt, sprintf('Optics (SI)\n')];
        % See above
        % txt = [txt, sprintf('  Mag:  %.2e\n', ...
        %     opticsGet(optics, 'magnification'))];
        diameter = opticsGet(optics, 'aperture diameter', 'mm');
        txt = [txt, sprintf('  Diameter:  %.2f mm\n', diameter)];
        if checkfields(oi, 'optics', 'lens')
            d = oiGet(oi, 'lens density');
            txt = [txt, sprintf('  Lens density:  %.2f \n', d)];
        end
    case 'iset3d'
        % Like shift invariant for now
        txt = [txt, sprintf('Optics (iset3d)\n')];
        diameter = opticsGet(optics,'aperture diameter','mm');
        txt = [txt, sprintf('  Diameter:  %.2f mm\n',diameter)];
        if checkfields(oi,'optics','lens')
            try
                d = oiGet(oi,'lens density');
                txt = [txt, sprintf('  Lens density:  %.2f \n',d)];
            catch
                warning('No lens density set in this optical image.');
            end
        end
    otherwise
        error('Unknown optics model %s. ', opticsModel);
end

return;