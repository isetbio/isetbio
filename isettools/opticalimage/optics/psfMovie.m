function psfMovie(oi,figNum,delay)
%Show a movie of the pointspread functions 
%
%   psfMovie([oi],[figNum=1],[delay = 0.2])
%
% Show a movie of the psf as a function of wavelength
%
% Copyright ImagEval, LLC, 2005

%Examples:
%{
 oi = oiCreate('shift invariant');
 psfMovie(oiGet(oi,'optics'));
%}
%{
 oi = oiCreate('diffraction limited');
 psfMovie(oiGet(oi,'optics'));
%}

%%
if notDefined('oi'), oi = vcGetObject('oi'); end
if notDefined('figNum'), figNum = vcSelectFigure('GRAPHWIN'); end
if notDefined('delay'), delay = 0.2; end

figure(figNum)
set(figNum,'name','PSF Movie');

%%
optics      = oiGet(oi,'optics');
opticsModel = oiGet(oi,'optics model');

switch lower(opticsModel)
    case 'diffractionlimited'
        % This needs to be checked.  The psf doesn't seem to be changing
        % correctly with wavelength. 
        %
        % The OTF is computed on the fly for the diffraction limited case.
        %
        nSamp = 16;
        freqOverSample = 4;
        fSupport = oiGet(oi,'fsupport','um')*freqOverSample;
        samp = (-nSamp:(nSamp-1));
        %         [X,Y] = meshgrid(samp,samp);
        %         deltaSpace = 1/(2*max(fSupport(:)));
        %         sSupport(:,:,2) = Y*deltaSpace;
        %         sSupport(:,:,1) = X*deltaSpace;

        deltaSpace = 1/(2*max(fSupport(:)));
        x = samp*deltaSpace;
        y = x;
        wave = opticsGet(optics,'wave');

        vcNewGraphWin;
        for ii=1:length(wave)
            psf = opticsGet(optics,'diffraction limited psf data',...
                wave(ii),'um',nSamp,freqOverSample);
            %{
            imagesc(y,x,psf(:,:));
            xlabel('Position (um)');
            ylabel('Position (um)');
            grid on; axis image
            title(sprintf('Wave %.0f nm',wave(ii)));
            pause(delay);
                %}
        end
        
        
    case 'shiftinvariant'
        % We get the psf data all at once in this case
        psf  = opticsGet(optics,'shift invariant psf data');
        support = opticsGet(optics,'psf support','um');
        y = support{1}(:); x = support{2}(:);
        wave = opticsGet(optics,'wavelength');
        w = size(psf,3);

        for ii=1:w
            imagesc(y,x,psf(:,:,ii));
            xlabel('Position (um)');
            ylabel('Position (um)');
            grid on; axis image
            title(sprintf('Wave %.0f nm',wave(ii)));
            pause(delay);
        end
        %{
    case 'raytrace'
        name = opticsGet(optics,'rtname');
        figNum = vcNewGraphWin;
        set(figNum,'name',sprintf('%s: PSF movie',name));
        colormap(gray(256));

        wave   = opticsGet(optics,'rt psf wavelength');
        imgHgt = opticsGet(optics,'rt psf field height','um');
        psf    = opticsGet(optics,'rt psf data');
        c = opticsGet(optics,'rt psf support col','um');
        r = opticsGet(optics,'rt psf support col','um');

        % Should we plot them on a single image and move them, or centered
        % like this?
        gColor = [.5 .5 0];
        for jj=1:length(wave)
            for ii=1:length(imgHgt)
                imagesc(r + imgHgt(ii),c + imgHgt(ii),squeeze(psf(:,:,ii,jj)));
                set(gca,'yticklabel',[]); xlabel('Position (um)');
                set(gca,'xcolor',gColor,'ycolor',gColor);
                grid on; axis image
                title(sprintf('Wave %.0f nm\nField height %.2f um',wave(jj),imgHgt(ii)));
                pause(delay);
            end
        end
        %}
    otherwise
        error('Unknown model %s\n',opticsModel);
end

end