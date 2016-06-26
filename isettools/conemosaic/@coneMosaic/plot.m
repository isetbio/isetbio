function hf = plot(obj, type, varargin)
% Plot function for coneMosaic
%    hf = coneMosaic.plot(type, varargin)
%
% Inputs:
%   type - string, type of plot
%
% Optional input (key-val pairs in varargin):
%   'newFig' - whether to create a new figure
%
% Outputs:
%   hf    - figure handle
%
% type can be chosen from
%   'cone mosaic'
%   'cone fundamentals'
%   'macular transmittance'
%   'macular absorptance'
%   'spectral qe'
%   'absorptions'
%   'current'
%
% Example:
%   cm = coneMosaic;
%   cm.plot('macular transmittance')
%
% HJ/BW, ISETBIO TEAM, 2016

% parse input
p = inputParser;
p.addRequired('type', @isstr);
p.addParameter('newFig', true, @islogical);

p.parse(type, varargin{:});
newFig = p.Results.newFig;

% plot
if newFig, hf = vcNewGraphWin; else hf = []; end

switch ieParamFormat(type)
    case 'conemosaic'
        [~,~,~,coneMosaicImage] = conePlot(obj.coneLocs*1e6, obj.pattern);
        imagesc(coneMosaicImage); axis off;
    case 'conefundamentals'
        plot(obj.wave, obj.cone.absorptance, 'LineWidth', 2); grid on;
        xlabel('Wavelength (nm)'); ylabel('Cone Quanta Efficiency');
    case 'maculartransmittance'
        plot(obj.wave, obj.macular.transmittance, 'LineWidth', 2);
        xlabel('Wavelength (nm)'); ylabel('Transmittance'); grid on;
    case 'macularabsorptance'
        plot(obj.wave, obj.macular.absorptance, 'LineWidth', 2);
        xlabel('Wavelength (nm)'); ylabel('Absorptance'); grid on;
    case 'spectralqe'
        plot(obj.wave, obj.qe, 'LineWidth', 2); grid on;
        xlabel('Wavelength (nm)'); ylabel('Quanta Efficiency');
    case 'absorptions'
        if isempty(obj.absorptions)
            error('no absorption data computed');
        end
        % [~,~,~,absImage] = conePlot(obj.coneLocs*1e6,obj.pattern,[], ...
        %    [], [], [], obj.absorptions(:, :, 1));
        % imagesc(absImage);
        imagesc(obj.absorptions(:, :, 1)); axis off;
    case 'current'
        if isempty(obj.current)
            error('no current data computed');
        end
        imagesc(obj.current(:, :, 1)); axis off;
    otherwise
        error('unsupported plot type');
end

end