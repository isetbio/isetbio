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
        % The lens transmittance should come from the OI or be passed in.
        % We handle various other cases here.
        if isempty(p.Results.lens)
            oi = vcGetObject('oi');
            if isempty(oi), thisLens = Lens('wave',cm.wave);
            else,           thisLens = oiGet(oi,'lens');
            end
        else
            thisLens = p.Results.lens;
        end
        
        lensT      = thisLens.transmittance;
        macularT   = cm.macular.transmittance;
        spectralQE = diag(lensT.*macularT)*cm.pigment.quantalEfficiency;
        plot(cm.wave,spectralQE,'LineWidth',2);
        xlabel('Wavelength (nm)'); ylabel('Spectral quantal efficiency'); grid on;
    case 'qe'
        thisLens = Lens('wave',cm.wave);
        lensT      = thisLens.transmittance;

        plot(cm.wave,diag(lensT)*cm.qe,'LineWidth',2);
        xlabel('Wavelength (nm)'); ylabel('Spectral quantal efficiency'); grid on;
    otherwise
        error('Unknown plotType %s\n',plotType);
end

end

