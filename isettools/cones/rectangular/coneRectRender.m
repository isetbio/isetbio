function axesHandle = coneRectRender(cRect,varargin)
% Render the coneRectMosaic cone array
%
% Synopsis
%   axesHandle = coneRectRender(cRect,varargin);
%
% Brief description
%   Based on the visualize method in cMosaic by NC.
%
% Input
%  cRect - coneMosaicRect object
%
% Optional Key/val pairs
%  axes handle    - Where to draw the image
%  aperture shape - Circle by default
%  unit           - Spatial unit for axes ('um' by default)
%  edge color - 
%  line width
%  face alpha
%  edge alpha
%
% Output
%   axesHandle - Where the 
%
% See also
%   visualize (renderPatchArray)
%

% Examples:
%{

 cRect = coneMosaicRect('center',[0 0]);
 cRect.patternSampleSize*1e6
 cRect.width*1e6
 coneRectRender(cRect);
%}
%{
 % Notice that the same row/col size spans a larger distance on the
 % retina because the cone apertures are larger
 cRect = coneMosaicRect('center',[1.0 0.0]*1e-3);
 cRect.patternSampleSize*1e6
 cRect.width*1e6
 coneRectRender(cRect);

%}
%% Parse
varargin = ieParamFormat(varargin);

p = inputParser;
p.addRequired('cRect',@(x)(isa(x,'coneMosaicRect')));

p.addParameter('axeshandle',[],@(x)(isa(x,'matlab.graphics.axis.Axes')));
p.addParameter('apertureshape',[],@isstruct);
p.addParameter('edgecolor',[0.1, 0.1, 0.1],@isvector);
p.addParameter('linewidth',1,@isnumeric);
p.addParameter('facealpha',1,@isnumeric);
p.addParameter('edgealpha',1,@isnumeric);
p.addParameter('unit','um',@ischar);

p.parse(cRect,varargin{:});

edgeColor = p.Results.edgecolor;
lineWidth = p.Results.linewidth;
faceAlpha = p.Results.facealpha;
edgeAlpha = p.Results.edgealpha;
axesHandle = p.Results.axeshandle;

if isempty(axesHandle)
    ieNewGraphWin;
    axesHandle = gca;
end

if isempty(p.Results.apertureshape)
    deltaAngle = 45;
    iTheta = (0:deltaAngle:360) / 180 * pi;
    apertureShape.x = cos(iTheta);
    apertureShape.y = sin(iTheta);
else
    apertureShape = p.Results.apertureshape;
end

locs = cRect.coneLocs * ieUnitScaleFactor(p.Results.unit);

% Pattern sample size is the aperture of each cone (square).  In this
% case, they are all the same.  
diameter = cRect.patternSampleSize(1)*ieUnitScaleFactor(p.Results.unit);  

%% Render

for ii=1:4
    faceColors = ii; 
    theseIndices = (cRect.pattern == ii);

    theseLocs = locs(theseIndices,:);  % Meters
    if numel(theseLocs) > 0

        apertureRadii = ones(size(theseLocs,1),1)*diameter/2;
        conesNum = numel(apertureRadii);
        verticesPerCone = numel(apertureShape.x);

        verticesList = zeros(verticesPerCone * conesNum, 2);
        facesList = [];

        if (numel(faceColors) == 1)
            colors = repmat(faceColors, [verticesPerCone*conesNum 1]);
        else
            colors = [];
        end

        for coneIndex = 1:conesNum
            idx = (coneIndex - 1) * verticesPerCone + (1:verticesPerCone);
            verticesList(idx, 1) = apertureShape.x*apertureRadii(coneIndex) + theseLocs(coneIndex,1);
            verticesList(idx, 2) = apertureShape.y*apertureRadii(coneIndex) + theseLocs(coneIndex,2);
            if ((numel(faceColors) == conesNum) && (conesNum > 1))
                colors = cat(1, colors, repmat(faceColors(coneIndex), [verticesPerCone 1]));
            end
            facesList = cat(1, facesList, idx);
        end


        S.Vertices = verticesList;
        S.Faces = facesList;
        S.FaceVertexCData = colors;
        S.FaceColor = 'flat';
        S.EdgeColor = edgeColor;
        S.FaceAlpha = faceAlpha;
        S.EdgeAlpha = edgeAlpha;
        S.LineWidth = lineWidth;
        patch(S, 'Parent', axesHandle);

    end

end

% There is something off because I am not accounting for k
% correctly. (BW)
cMap = [...
    1.0000    0.1000    0.5000;
    0.1000    1.0000    0.5000;
    0.1000    0.1000    1.0000];

colormap(axesHandle,cMap);
axis image;

end
