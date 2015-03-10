function plotSensorVignetting;
%
%   plotSensorVignetting
%
% Author: ImagEval
% Purpose:
%    Plot the sensor vignetting function.  
%

disp('Obsolete')
evalin('Caller',mfilename)
return;


[valISA,ISA] = vcGetSelectedObject('ISA');
[valOI,OI] = vcGetSelectedObject('OPTICALIMAGE');

vignetting = pvVignetting(ISA,OI);

row = sensorGet(ISA,'rows');
col = sensorGet(ISA,'cols');
meshX = linspace(1,col,size(vignetting,1));
meshY = linspace(1,row,size(vignetting,2));

figNum = vcSelectFigure('GRAPHWIN');
PlotSetUpWindow(figNum);
colormap('default');
mesh(meshY,meshX,vignetting);
view([30,2]);

xlabel('Sensor col'); ylabel('Sensor row'); zlabel('Tranmission')

ttl = title('Light loss due to 3D pixel geometry (vignetting)');
set(ttl,'Fontsize',10);
txt = [];

pixel = sensorGet(ISA,'pixel'); pixelDepth = pixelGet(pixel,'depth');
newText = sprintf('Pixel Depth:  %.01f um\n',pixelDepth*10^6);
txt = addText(txt,newText);

optics = sceneGet(OI,'optics'); fnumber = opticsGet(optics,'fnumber');
newText = sprintf('Optics f#:  %.02f',fnumber);
txt = addText(txt,newText);

%mx = max(vignetting(:)); mn = min(vignetting(:)); z = mx - (mx-mn)*.3;
%mx = max(row); mn = min(row); y = mn + (mx-mn)*.2;
%mx = max(col); mn = min(col); x = mn + (mx-mn)*.2;

t = text(0.8,0.8,0.8,txt)
set(t,'Background','w','Fontsize',8);

udata = vignetting;
set(gcf,'userdata',udata);

ISA.data.vignetting = vignetting;
vcReplaceObject(ISA,valISA);

return;
