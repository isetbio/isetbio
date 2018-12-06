function [] = writeArizonaLensFile(A,filename)
% Write out the Arizona lens file for a given accomodation
%
% Syntax:
%   writeArizonaLensFile(A,filename)
%
% Description:
%    Write out the Arizona lens file for a given accomodation
%
% Inputs:
%    A        - The accommodation, in diopters
%    filename - The filename to write to
%
% Outputs:
%    None.
%

% All following equations are from :
% https://photonengr.com/wp-content/uploads/kbasefiles/ArizonaEyeModel.pdf
anteriorRadius = 12-(0.4*A);
posteriorRadius = -5.22+(0.2*A); 
aqueousThickness = 2.97-(0.04*A) ;
lensThickness = 3.77+(0.04*A) ;
anteriorAsph = -7.52+(1.29*A);
posteriorAsph = -1.35-(0.43*A);

% Columns are: radiusX, radiusY, thickness, materialIndex, semiDiameter,
% conicConstantX, and conicConstantY
corneaA = [-7.8   -7.8   0.55  1  4.820  -0.25   -0.25];
corneaP = [-6.5   -6.5   2.97  2  4.341  -0.25   -0.25];
pupil =   [ 0      0     0     2  2       0       0];
lensA =   [-12    -12    3.77  3  3.750  -7.52   -7.52];
lensP =   [ 5.22   5.22  0     4  3.750  -1.35   -1.35];

% flip these because of sign conventions
lensA(1:2) = [anteriorRadius anteriorRadius] * -1;
lensP(1:2) = [posteriorRadius posteriorRadius] * -1;

corneaP(3) = aqueousThickness;
lensA(3) = lensThickness;
lensA(6:7) = [anteriorAsph anteriorAsph];
lensP(6:7) = [posteriorAsph posteriorAsph];

%% Build matrix
lensMatrix = [corneaA; corneaP; pupil; lensA; lensP];

focalLength = 1 / (60.0 + A) * 10 ^ 3; % mm
fid = fopen(filename,'w');

str = sprintf('# Focal length (mm) \n');
fprintf(fid,'%s',str);
str = sprintf('%.3f\n',focalLength);
fprintf(fid,'%s',str);
str = sprintf(['# radiusX radiusY thickness materialIndex semiDiameter' ...
    ' conicConstantX conicConstantY\n']);
fprintf(fid,'%s',str);
for ii=1:size(lensMatrix,1)
    fprintf(fid,'%f\t%f\t%f\t%f\t%f\t%f\t%f\n', ...
        lensMatrix(ii,1), lensMatrix(ii, 2), lensMatrix(ii,3), ...
        lensMatrix(ii,4), lensMatrix(ii,5), lensMatrix(ii,6), ...
        lensMatrix(ii,7));
end
fclose(fid);

end
