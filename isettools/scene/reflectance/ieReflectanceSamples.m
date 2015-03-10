function [reflectances, sSamples] = ieReflectanceSamples(sFiles,sSamples,wave,replacement)
% Return a sample of reflectances 
%
%   [reflectances, sSamples] = ieReflectanceSamples(sFiles,sSamples,[wave],[sampling])
%
% The reflectances are drawn from the cell array of names in sFiles{} and
% returned in the columns of the matrix reflectances.
%
% sFiles:   Cell array of file names with reflectance spectra
% sSamples: This can either be
%      - Vector indicating how many surfaces to sample from each file
%      - A cell array of specifying the list of samples from each file
% wave:     Wavelength Samples
% sampling: 'r' means with replacement (default)
%
% See also: macbethChartCreate, sceneReflectanceChart
%
%Example:
%  sFiles = cell(1,2);
%  sFiles{1} = fullfile(isetRootPath,'data','surfaces','reflectances','MunsellSamples_Vhrel.mat');
%  sFiles{2} = fullfile(isetRootPath,'data','surfaces','reflectances','DupontPaintChip_Vhrel.mat');
%  sSamples = [12,12]*5; 
%  [reflectances, sSamples] = ieReflectanceSamples(sFiles,sSamples);
%  plot(wave,reflectances); grid on
%
% Copyright ImagEval Consultants, LLC, 2010.

if notDefined('sFiles'), error('Surface files required'); 
else                     nFiles = length(sFiles);
end

if notDefined('sSamples'), sSamples = zeros(1,nFiles);
elseif length(sSamples) ~= nFiles
    error('Mis-match between number of files and sample numbers');
end
if notDefined('wave'), wave = 400:10:700; end
if notDefined('replacement'), replacement = 'r'; end % With replacement

% sSamples might be a vector, indicating the number of samples, or a cell
% array specifying which samples.
if iscell(sSamples)
    nSamples = 0;
    for ii=1:nFiles
        nSamples = length(sSamples{ii}) + nSamples;
    end
else
    nSamples = sum(sSamples);
end

% Read the surface reflectance data. At the moment, we allow duplicates,
% though we should probably change this.
last = 0;
sampleList = cell(1,nFiles);
reflectances = zeros(length(wave),nSamples);

for ii=1:nFiles

    allSurfaces = ieReadSpectra(sFiles{ii},wave);
    nRef = size(allSurfaces,2);
    
    % Generate the random list of surfaces.  They are sampled with
    % replacement.
    if ~iscell(sSamples)
        if strncmp(replacement, 'r', 1)  % With replacement
            % randi doesn't exist in 2008 Matlab.
            if exist('randi', 'builtin')
                sampleList{ii} = randi(nRef,[1 sSamples(ii)]);
            else
                sampleList{ii} = ceil(rand([1 sSamples(ii)])*nRef);
            end
        else  % Without replacement
            if sSamples(ii) > nRef
                error('Not enough samples in %s\n',sFiles{ii});
            else
                sampleList{ii} = randperm(nRef, sSamples(ii));
            end
        end
        % fprintf('Choosing %d of %d samples\n',sSamples(ii),nRef);
        % fprintf('Unique samples:   %d\n',length(unique(sampleList{ii})));
    else
        sampleList{ii} = sSamples{ii};
    end
    
    this = last + 1;
    last = this + (length(sampleList{ii}) - 1);
    reflectances(:,this:last) = allSurfaces(:,sampleList{ii});
end

if max(reflectances(:)) > 1, error('Bad reflectance data'); end

% With this variable and the the filenames, you can run the same function
% and get the same reflectances.
sSamples = sampleList;

end