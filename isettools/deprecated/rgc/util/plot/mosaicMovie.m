function psthMovie = mosaicMovie(irPSTH, ir, params)
% Use cell PSTH to create an RGC mosaic movie.
%
% Syntax:
%   psthMovie = mosaicMovie(irPSTH, ir, params)
%
% Description:
%    Makes a movie of the RGC mosaic according to each cell's PSTH. This is
%    useful for visualizing the spatial as
%
% Inputs:
%    irPSTH    - Cell. The cell array of inner retina PSTH objects.
%    ir        - Object. The inner retina object.
%    params    - Struct. A parameters structure.
%
% Outputs:
%    psthMovie - Matrix. A 3D matrix representing the PSTH movie.
%
% Optional key/value pairs:
%    None.
%

for i = 1:length(ir.mosaic{1}.cellLocation)
    loc(i, :) = params.inputScale .* ir.mosaic{1}.cellLocation{i}; 
end
mloc = max(loc);

rfRad = 1 * ir.mosaic{1}.rfDiaMagnitude / 2;
psthMovie = zeros(2 * ceil(mloc(1)), ...
    2 * ceil(mloc(2)), length((irPSTH{1})'));
for i = 1:length(ir.mosaic{1}.cellLocation)
    loc(i, :) = params.inputScale .* ir.mosaic{1}.cellLocation{i}; 
    % psthMovie(round(2 * loc(i, 1)), round(2 * loc(i, 2)), :) = irPSTH{i};
    xSelStart = round(2 * loc(i, 1) - rfRad / 2);
    xSelEnd = round(2 * loc(i, 1) + rfRad / 2);
    ySelStart = round(2 * loc(i, 2) - rfRad / 2);
    ySelEnd = round(2 * loc(i, 2) + rfRad / 2);
    psthMovieCell(1, 1, :) = irPSTH{i};
    pRep = repmat( psthMovieCell(1, 1, :), ...
        [length(xSelStart:xSelEnd), length(ySelStart:ySelEnd), 1]);
    psthMovie(xSelStart:xSelEnd, ySelStart:ySelEnd, :) = pRep;
end
% ieMovie(psthMovie(3:end - 2, :, 200:end));
