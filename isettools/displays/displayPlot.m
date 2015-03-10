function [uData, g] = displayPlot(d,param,varargin)
% Gateway routine for plotting display structure
%
%   [uData, g] = displayPlot(d,param,varargin)
%
% Example:
%  d = displayCreate('CRT-Dell');
%  displayPlot(d,'gamut');
%  
%  displayPlot(d,'spd')
%
%  d = displayCreate('crt');
%  displayPlot(d,'gamma table')
%
% (c) Imageval Consulting, 2013

if notDefined('d'), error('Display required'); end

% format parameter - lower case and no space
param = ieParamFormat(param);

switch param
    case 'spd' % Plot spectral power distribution of the display
        spd = displayGet(d,'spd primaries');
        wave = displayGet(d,'wave');
        g = vcNewGraphWin;
        cOrder = {'r','g','b','k','y'}; % color order
        hold on
        for ii=1:size(spd,2)
            plot(wave,spd(:,ii),cOrder{ii});
        end
        
        xlabel('Wavelength (nm)');ylabel('Energy (watts/sr/m2/nm)');
        grid on; uData.wave = wave; uData.spd = spd;
        set(g,'userdata',uData);
        
    case {'gammatable','gamma'} % Plot display gamma curve
        gTable = displayGet(d,'gamma table');
        g = vcNewGraphWin; plot(gTable);
        xlabel('DAC'); ylabel('Linear');
        grid on
        
        uData = gTable;
        set(g,'userdata',uData);
    case {'inversegamma', 'invertgamma'} % Plot invert gamma table
        invG = displayGet(d, 'inverse gamma table');
        steps = linspace(0,1,length(invG));
        g = vcNewGraphWin; plot(steps, invG);
        xlabel('Linear'); ylabel('DAC');
        grid on;
        
        uData = invG;
        set(g,'userdata',uData);
        
    case 'gamut'  % Plot color gamut in chromaticity (xy) space
        spd = displayGet(d, 'spd primaries');
        wave = displayGet(d, 'wave');
        XYZ = ieXYZFromEnergy(spd', wave);
        xy = chromaticity(XYZ);
        
        % eliminate black primary
        indx = (sum(xy, 2) < 0.1);
        xy(indx, :) = [];
        
        g = chromaticityPlot(xy, 'gray', 256);
        indx = convhull(xy(:,1), xy(:,2));
        xy = xy(indx, :);
        l = line(xy(:,1),xy(:,2));
        set(l,'color',[.6 .6 .6]);
        
        % Store data in figure
        uData.xy = xy;
        set(g,'userdata',uData);
        
    case 'gamut3d'  % Plot display gamut in Lab space (3D)
        % check number of primaries in display
        if displayGet(d, 'n primaries') > 3
            warning('Display has more than 3 primaries');
            disp('Only first 3 primaries are used');
        end
        % create new graph
        vcNewGraphWin; ha = gca;
        
        % samples dac values
        nSamp  = 30;
        gTable = displayGet(d, 'gamma table');
        dac = linspace(0, length(gTable)-1, nSamp + 1);
        dac = round(dac(2:end));
        
        % We transform dac to nSamp^3 x 1 x3
        dacM = meshgrid(dac, dac, dac);
        dac = dacM(:);
        dacM = permute(dacM, [2 3 1]);
        dac = [dac dacM(:)];
        dacM = permute(dacM, [2 3 1]);
        dac = [dac dacM(:)];
        dac = dac(:, [2 1 3]);
        dac = reshape(dac, [nSamp^3 1 3]);
        
        % convert to display rgb with gamma table
        rgb = ieLUTDigital(dac, gTable);
        rgb = squeeze(rgb);
        rgb = padarray(rgb, [0 displayGet(d,'nprimaries') - 3], 0, 'post');
        
        % convert rgb to Lab
        XYZ = rgb * displayGet(d, 'rgb2xyz');
        Lab = ieXYZ2LAB(XYZ, displayGet(d, 'white xyz'));
        
        % generate delaunayTriangle
        delTri  = delaunayTriangulation(Lab);
        [ch, ~] = convexHull(delTri);
        
        % plot
        ptx = delTri.Points;
        trisurf(ch, ptx(:,2), ptx(:,3), ptx(:,1), ...
                'parent', ha, 'Tag', 'gamutPatch');
        axis image;
        
        % assign color to the patches
        hPat = findobj('Tag', 'gamutPatch', '-and', 'parent', ha);
        
        for ii = 1 : length(hPat)            
            % set color for the vertices
            set(hPat(ii), ...
                'FaceColor', 'interp', ...
                'FaceVertexCData', rgb(:, 1:3), ...
                'EdgeColor', 'interp', ...
                'LineWidth', 1, ...
                'facealpha', 0.5);
        end
        
        % title
        title('Color Gamut of display in Lab space');
        xlabel('a'); ylabel('b'); zlabel('L');
        
    otherwise
        error('Unknown parameter %s\n',param);
end

end

    