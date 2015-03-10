function ISA = sensorReadFilter(filterType,ISA,fname)
%Read color filter data of various sorts
%
%   isa = sensorReadFilter(filterType,isa,fname)
%
% Read and set the sensor wavelength properties.  The files are assumed to
% have variables with specific names defining the filter properties.  The
% specific variables and formats for each type are described below.
%
% Data defining color filters, the photodector spectral quantum efficiency,
% infrared cutoff filter, and the color filter array, can all be read by
% this program.
%
% Resetting any of these functions changes the currently the properties of
% the image sensor array.
%
% The wavelength sampling of the data matches the current specification in
% ISA.  To make the data match the properties of the optical image, make
% sure that
%
% Data formats:
%
%    FILL THIS IN
%
% Example:
%   [val,ISA] = vcGetSelectedObject('ISA');
%   ISA = sensorReadFilter('infrared',ISA);
%   ISA = sensorClearData(ISA);
%   vcReplaceObject(ISA,val);

if notDefined('ISA'), [~, ISA] = vcGetSelectedObject('ISA'); end
if notDefined('filterType'), filterType = 'cfa'; end
if notDefined('fname'), fname = []; end

wave = sensorGet(ISA,'wave');

switch lower(filterType)
    case 'cfa'
        if isempty(fname)
            fname = vcSelectDataFile(['sensor',filesep,'CFA']);
            if isempty(fname),  return;
            else
                if exist(fname,'file'),tmp = load(fname);
                else                    error('No such file %s\n',fname);
                end
            end
        end
        if checkfields(tmp,'spectrum')
            sWave = sensorGet(ISA,'wave');
            cWave = tmp.spectrum.wave;
            if isequal(sWave,cWave)
                if checkfields(tmp,'color'), ISA.color = tmp.color;
                else warning('No color structure in file.');end
                if checkfields(tmp,'cfa'), ISA.cfa = tmp.cfa;
                else warning('No cfa structure in file.'); end
            else
                error('CFA wavelength information does not match current sensor.');
            end
        else
            error('No spectrum (wavelength) information in CFA file');
        end

    case 'pdspectralqe'
        data = ieReadColorFilter(wave,fname);
        if isempty(data), return; end
        pixel = sensorGet(ISA,'pixel');
        pixel = pixelSet(pixel,'pdspectralqe',data);
        ISA = sensorSet(ISA,'pixel',pixel);
        % ISA = sensorSet(ISA,'pdspectralqe',data);

    case {'infrared','irfilter'}
        data = ieReadColorFilter(wave,fname);
        if isempty(data), return; end
        ISA = sensorSet(ISA,'infrared',data);

    case {'colorfilters','colorfilter'}
        [data,newFilterNames] = ieReadColorFilter(wave,fname);
        if isempty(data), return; end
        pattern = sensorGet(ISA,'pattern');
        nCols = size(data,2);
        pattern = min(nCols,pattern);
        ISA = sensorSet(ISA,'filterspectra',data);
        ISA = sensorSet(ISA,'filternames',newFilterNames);
        ISA = sensorSet(ISA,'pattern',pattern);

    otherwise
        error('Unknown color type.');
end

return;
