function keyValues = wvfKeySynonyms(keyValues)
% Convert cell array of key value pairs to canonical form
%
% Syntax:
%   keyValues = wvfKeySynonyms(keyValues)
%
% Description:
%   Check each string in the odd entries of the passed string cell
%   array, and convert each to the canonical form understood by the
%   wvf code.
%
%   Typical usage would be to pass varargin through this, resulting
%   in a cell array of strings that could be passed to the input parser,
%   where the input parser understood only one of each set of synonyms.
%
% Inputs:
%   keyValues - a cell array of strings.  Odd entries are converted.
%
% Outputs:
%   KeyValues - a cell array of strings, after conversion.
%
% See also:
%  ieParamFormat, wvfCreate, wvfGet, wvfSet

% History:
%   12/07/17  dhb  Wrote.

%% Loop over keys and convert each according to synoyms
for kk = 1:2:length(keyValues)
    switch (keyValues{kk})
        case {'name'}
            keyValues{kk} = 'name';
        case {'type'}
            keyValues{kk} = 'type';
        case {'measuredpupilsize', 'measuredpupil', 'measuredpupilmm', ...
                'measuredpupildiameter'}
            keyValues{kk} = 'measuredpupil';
        case {'measuredwave', 'measuredwl', 'measuredwavelength'}
            keyValues{kk} = 'measuredwl';
        case {'measuredopticalaxis', 'measuredopticalaxisdeg'}
            keyValues{kk} = 'measuredopticalaxis';
        case {'measuredobserveraccommodation', ...
                'measuredobserveraccommodationdiopters'}
            keyValues{kk} = 'measuredobserveraccommodation';
        case {'measuredobserverfocuscorrection', ...
                'measuredobserverfocuscorrectiondiopters'}
            keyValues{kk} = 'measuredobserverfocuscorrection';
        case {'zcoeffs', 'zcoeff', 'zcoef'}
            keyValues{kk} = 'zcoeffs';
        case {'sampleintervaldomain'}
            keyValues{kk} = 'sampleintervaldomain';
        case {'numberspatialsamples', 'spatialsamples', 'npixels', ...
                'fieldsizepixels'}
            keyValues{kk} = 'spatialsamples';
        case {'refpupilplanesize', 'refpupilplanesizemm', 'fieldsizemm'}
            keyValues{kk} = 'refpupilplanesize';
        case {'refpupilplanesampleinterval', 'fieldsamplesize', ...
                'refpupilplanesampleintervalmm', 'fieldsamplesizemmperpixel'}
            keyValues{kk} = 'refpupilplanesampleinterval';
        case {'refpsfsampleinterval' 'refpsfarcminpersample', ...
                'refpsfarcminperpixel'}
            keyValues{kk} = 'refpsfsampleinterval';           
        case {'calcpupilsize', 'calcpupildiameter', 'calculatedpupil', ...
                'calculatedpupildiameter'}    
            keyValues{kk} = 'calcpupilsize';       
        case {'calcopticalaxis'}
            keyValues{kk} = 'calcopticalaxis'; 
        case {'calcobserveraccommodation'}
            keyValues{kk} = 'calcobserveraccommodation';
     case {'calcobserverfocuscorrection', 'defocusdiopters'}
            keyValues{kk} = 'calcobserverfocuscorrection';
     case {'calcwave', 'calcwavelengths', 'wavelengths', 'wavelength', ...
                'wls', 'wave'}
            keyValues{kk} = 'calcwavelengths';         
        case {'calcconepsfinfo'}
            keyValues{kk} = 'calcconepsfinfo';         
        case {'sceparams', 'stilescrawford'}
            keyValues{kk} = 'sceparams';     
    end
end