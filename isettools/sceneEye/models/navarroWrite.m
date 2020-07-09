function navarroWrite(thisR)
% Write out the navarro lens and support IoR files in the output directory
%
% Synopsis
%    navarroWrite(thisR);
%
% Description
%   The navarro.dat file and associated index of refraction files
%   (iorX.spd) are written into the lens rendering directory. The navarro
%   model accounts for accommodation in the ior and lens files.  We use the
%   'object distance' slot to define the accommodation.  Until the
%   accommodation is less than 0.5 m, the impact of that factor is very
%   small on the IOR.
%
% Input
%  thisR:  The rendering recipe.  It should include an object distance.
%
% Optional key/val pairs
%   N/A
%
% Outputs
%   N/A
%
% See also
%   navarroLensCreate, navarroRefractiveIndices

% Examples:
%{
  myScene = sceneEye('slantedBar');
  thisR = myScene.recipe;
  navarroWrite(thisR);
%}

%% Until we have another idea, we use this wave default as per TL
wave = (400:10:800); % nm

%% Writes out the navarro.dat file in the lens directory of the output
lensFile = fullfile(thisR.get('lens dir output'),'navarro.dat');
accom    = (1 / thisR.get('focal distance','m'));

navarroLensCreate(lensFile,'accommodation',accom);  % Diopters

% The file should be there, so no warning should come from this.
thisR.set('lens file',lensFile);

%% Now write out the IoR files

% Our convention (which was hard coded in writeNavarroLensFile) is
% ior1 --> cornea
% ior2 --> aqueuous
% ior3 --> lens
% ior4 --> vitreous
iorNames = {'ior1.spd','ior2.spd','ior3.spd','ior4.spd'};

% We assume the eye is accommodated to the object distance.  There is only
% a very very small impact of accommodation until the object is very close
% (less than 0.5 m).
[cornea, aqueuous, lens, vitreous] = navarroRefractiveIndices(wave, accom);
ior = [cornea(:),aqueuous(:), lens(:), vitreous(:)];

% We will put these files next to the lens file (navarro.dat).
nSamples = numel(wave);
for ii=1:4
    filename = fullfile(thisR.get('lens dir output'),iorNames{ii});
    fid = fopen(filename, 'w');
    for jj = 1:nSamples
        fprintf(fid, '%d %f\n', wave(jj), ior(jj,ii));
    end
    fclose(fid);
    
    % Update the recipe with the ior files
    [~,str,~] = fileparts(filename);
    thisR.set(str,filename);
end

fprintf('Wrote Navarro lens information with accommodation %f\n',accom);

end
