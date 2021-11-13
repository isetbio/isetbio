function figHdl = plot(cm,plotType,varargin)
% Interface method for plotting cMosaic parameters
%
% Synopsis
%    cm = cMosaic;
%    cm.plot(varargin);
%
% Inputs
%
% Optional key/val plot
%
% Output
%   figHdl
%
% See also
%   cMosaic.visualize

%% Parse input parameters
varargin = ieParamFormat(varargin);

p = inputParser;
p.addRequired('cm');
p.addRequired('plotType',@ischar);
p.addParameter('lens',[]);
p.parse(cm,plotType,varargin{:});

if isequal(ieParamFormat(plotType),'help')
    % Return and display a cell array of possible plot types
    % Maybe make this a table?
    figHdl = {'pigment absorptance','spectral qe'};
    return;
end

%%

figHdl = ieNewGraphWin;

switch ieParamFormat(plotType)
    case 'pigmentquantalefficiency'
        plot(cm.wave,cm.pigment.quantalEfficiency,'LineWidth',2);
        xlabel('Wavelength (nm)'); ylabel('Quantal efficiency'); grid on;
        
    case 'spectralqe'
        % The spectral QE includes the lens, macular and photo pigments.
        %
        % The lens is not part of the cMosaic, which makes this plot a
        % challenge. In principle, it should come from the OI or be passed
        % in. We handle various cases here.
        if isempty(p.Results.lens)
            thisLens = Lens('wave',cm.wave);            
        else
            % The user passed the lens in. So use it.
            thisLens = p.Results.lens;
            thisLens.wave = cm.wave;
        end
        
        % Match the lens to the cMosaic wave
        lensT      = thisLens.transmittance;        
        % ieNewGraphWin; plot(cm.wave,lensT);
        
        % The cMosaic qe is derived from macular and pigment information.
        % But it now seems it is not updated when I change the pigment
        % absorbance.
        plot(cm.wave,diag(lensT) * cm.pigment.quantalEfficiency,'LineWidth',2);
        xlabel('Wavelength (nm)'); ylabel('Spectral quantal efficiency'); grid on; 
        
    otherwise
        error('Unknown plotType %s\n',plotType);
end

end

