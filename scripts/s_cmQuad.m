%%
ieInit
chdir(fullfile(isetbioRootPath,'local'));

%% For later computations
scene = sceneCreate('rings rays');
fov = 1;
scene = sceneSet(scene,'fov',1);
oi = oiCreate; oi = oiCompute(oi,scene);
ieAddObject(oi);

cm = coneMosaic;
cm.setSizeToFOV(fov*0.8);
cm.emGenSequence(50);
cm.compute(oi);
cm.window('show','mean absorptions');

%% Calculate a move for the quadrature calculation

% Make the Gabor patches in quadrature phase
hparams = harmonicP; 
hparams.row = 32; hparams.col = 32;
hparams.GaborFlag = 0;
sQuad = imageHarmonic(hparams);
sQuad = sQuad - 1;
vcNewGraphWin; mesh(sQuad); colormap(gray)

hparams.ph = 0;
cQuad = imageHarmonic(hparams);
cQuad = cQuad - 1;
vcNewGraphWin; mesh(cQuad); colormap(gray)

%% Check that energy combination is invariant for lateral displacement
img = randn(32,32);
vcNewGraphWin; imagesc(img);

test1 = conv2(img,cQuad);
test2 = conv2(img,sQuad);
testA = test1.^2 + test2.^2;

img2 = circshift(img,1,2);
vcNewGraphWin; imagesc(img2);
test1 = conv2(img2,cQuad);
test2 = conv2(img2,sQuad);
testB = test1.^2 + test2.^2;

plot(testA(:),testB(:),'.')


%% Convolve the absorptions (with fixational eye movements) with sQuad/cQuad

[r,c,t] = size(cm.absorptions);
e1 = zeros(size(cm.absorptions));
e2 = zeros(size(cm.absorptions));
mEnergy = zeros(size(cm.absorptions));

for ii=1:t
    e1(:,:,ii) = conv2(cm.absorptions(:,:,ii),cQuad,'same');
    e2(:,:,ii) = conv2(cm.absorptions(:,:,ii),sQuad,'same');
end
ieMovie(e1);
% ieMovie(e2);

for ii=1:t
    mEnergy(:,:,ii) = e1(:,:,ii).^2 + e2(:,:,ii).^2;
end
ieMovie(mEnergy);


    



