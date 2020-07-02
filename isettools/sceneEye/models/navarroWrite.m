function navarroWrite(thisR)
% Write out the navarro lens and support IoR files in the output directory
%
% Synopsis
%
% Description
%
% Input
%
% Optional key/val pai8rs
%
% Outputs
%
% See also
%

%%

% Writes out the navarro.dat file in the lens directory of the output
fullfilename = navarroLensCreate(filename,varargin);

%% Now write out the IoR files in the spds/lens directory?

% Or create an spds folder inside of the lens directory?
%
wave = (400:10:800); % nm
accom = thisR.get('object distance','m');

% Our convention (hard coded in writeNavarroLensFile) is always:
% ior1 --> cornea
% ior2 --> aqueuous
% ior3 --> lens
% ior4 --> vitreous
% We assume the eye is accommodated to the object distance
[cornea, aqueuous, lens, vitreous] = navarroRefractiveIndices(wave, 1/accom);
ior = [cornea(:),aqueuous(:), lens(:), vitreous(:)];
%% Write the spectrum to file.

iorNames = {'ior1.spd','ior2.spd','ior3.spd','ior4.spd'};

nSamples = numel(wave);
for ii=1:4
    filename = fullfile(thisR.get('lensfile directory'),iorNames{ii});
    fid = fopen(filename, 'w');
    for jj = 1:nSamples
        fprintf(fid, '%d %f\n', wave(jj), ior(jj,ii));
    end
    fclose(fid);
end

% Set the recipe slots correctly
thisR.camera.ior1.value = fullfile(workingFolder, iorNames{1});
thisR.camera.ior2.value = fullfile(workingFolder, iorNames{2});
thisR.camera.ior3.value = fullfile(workingFolder, iorNames{3});
thisR.camera.ior4.value = fullfile(workingFolder, iorNames{4});

thisR.camera.lensfile.value = fullfile(workingFolder, lensFile);
thisR.camera.lensfile.type = 'string';

end
