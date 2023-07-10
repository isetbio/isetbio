%% s_fazMask
%
%{
June 5th email

My primary goal is to apply a mask on the cone mosaic such that those
cones covered by the mask do not participate in the computation of the
excitation or photon absorption.  

I attempted to bootstrap already programmed parameters in cMosaic.m to
accomplishing such a 2D mask. However, these parameters affect all the
cones within the mosaic.  
      1. eccVaryingConeAperture
      2. eccVaryingOuterSegmentLength
I then found a function within the scene folder labeled scenceCrop.m
which I thought I could crop specific regions from an unchanged
cMosaic which would mimic a masked cone mosaic.   

Questions: 
1. Is there a parameter in cMosiac that can accomplish a 2D mask?
      Is there a function that allows individual control of cone placement?
      Is there a function that allows individual control of cone absorption?
2. Should I create a cone mosaic which incorporates the mask into the
mosaic and then use the cMaskMosaic or should I apply the mask after
the scene has been applied to the cone mosaic matrix?   

%}

%%
ieInit;

%% Make a simple test scene

% 1 cpd harmonic
scene = sceneCreate('harmonic');
scene = sceneSet(scene,'fov',1);
sceneWindow(scene);

oi = oiCreate('wvf human');
oi = oiCompute(oi,scene);
oiWindow(oi);
%%  Let's try a rectangular mosaic, first

thisM = coneMosaicRect;
% class(thisM)

coneMask.row = linspace(-150,150,256)*1e-6;  % Units of meters
coneMask.col = linspace(-150,150,256)*1e-6;  %
% coneMask.img = rand(256,256);
coneMask.img = imageDeadLeaves(256,0.1);

imgLine = ones(256,256);
imgLine(:,25) = 0;
imgLine(:,97) = 0.2;
imgLine(:,185) = 0.5;
imgLine(:,215) = 0.0;
coneMask.img = imgLine;

ieNewGraphWin; imagesc(imgLine); colormap("gray"); axis image

thisM.coneMask = coneMask;

ieNewGraphWin; imagesc(coneMask.img); colormap("gray"); axis image

thisM.compute(oi);
thisM.window;

%%
chdir(fullfile(isetbioRootPath,'local'))
bv = imread('bloodvessels.png');
bv = sum(bv,3);
bv = imresize(bv,[256 256]);
ieNewGraphWin; imagesc(bv); colormap(gray); axis image
bv = ieScale(bv,0,1);
bv = 1- bv;
imwrite(bv,'coneMask_bloodvessels.png');




%%

thisM = cMosaic;  % Put in more parameters here.


%%