%% t_coneEccentricitySize
%
%  Create cone mosaics at different eccentricities.  The goal is to
%  illustrate how the cone apertures change in size and how this impacts
%  properties such as total cone absorptions and the impact of eye
%  movements.
%
% See also
%
% 

%%
ieInit;

%% A simple harmonics stimulus

hparms = harmonicP('freq',4,'contrast',0.0);
scene = sceneCreate('harmonic',hparms);
% sceneWindow(scene);

%%
oi = oiCreate('human'); 
oi = oiCompute(oi,scene);
% oiWindow(oi);

%%  Build a human cone mosaic

cm = coneMosaicRect('center',[0 0]*1e-6,'eccentricityUnits','m');
cm.setSizeToFOV(1);

cm.emGenSequence(50);
cm.compute(oi);

% cm.plot('cone mosaic');
cm.window;

%%
cm = coneMosaicRect('center',[0 450]*1e-6,'eccentricityUnits','m');
cm.setSizeToFOV(1);

cm.emGenSequence(50);
cm.compute(oi);

% cm.plot('cone mosaic');
cm.window;

%%

cm = coneMosaicRect('center',[0 900]*1e-6,'eccentricityUnits','m');
cm.setSizeToFOV(1);

cm.emGenSequence(50);
cm.compute(oi);

% cm.plot('cone mosaic');
cm.window;