function mosaicPlot(ir, bp, sensor, irparams, cellType, ecc)
% Plots the cone, bipolar and RGC mosaic on top of one another
% 
%       mosaicPlot(ir, bp, sensor, irparams, cellType, ecc)
% 
% This is the mosaic plot function for an RGC model that has been generated
% by isetbio from an average of the parameters of a mosaic measured
% experimentally. The cones, bipolars and RGCs are plotted together in a
% way that allows the user to see their relative spatial scale.

% (c) 2016 JRG isetbio team
%% Plot the cone mosaic

figure;

spacing = 1e6*irGet(ir,'spacing');

cone_mosaic = sensor.human.coneType;

scale2 = size(cone_mosaic,2)/spacing;%1;%size(cone_mosaic,2)/40;
scale1 = scale2;
  
[xg yg] = meshgrid([1:size(cone_mosaic,1),1:size(cone_mosaic,2)]);
xg2 = xg(1:size(cone_mosaic,1),1:size(cone_mosaic,2)); yg2 = yg(1:size(cone_mosaic,1),1:size(cone_mosaic,2));
% scatter(yg2(:)./scale1,xg2(:)./scale2,40,4-cone_mosaic(:),'o','filled'); colormap jet; set(gca,'color',[0 0 0])

imagesc((1:size(cone_mosaic,2))./scale1, (1:size(cone_mosaic,1))./scale1, 4-cone_mosaic); colormap jet;

%% Plot the bipolar mosaic
bpRad = 1;%size(bp.sRFcenter,1);
circle_samples = 0:0.05:2*pi;
x_circle = bpRad*cos(circle_samples);
y_circle = bpRad*sin(circle_samples);
% 
szbp = size(bp.sRFcenter,1)-.5;
% figure;
hold on;
szmosaic = size(cone_mosaic);
ycind = 0; xcind = 0;
magjitter = 0.6;
for yc = 1:szbp*2:szmosaic(1)
    ycind = ycind+1;
    for xc = 1:szbp*2:szmosaic(2)
        xcind = xcind+1;
%         plot((szbp*y_circle+yc)./scale2,(szbp*x_circle+xc+(szbp)*mod(ycind+1,2))./scale1,'linewidth',1,'color','c');
        plot((magjitter*rand(1)+.5+szbp*x_circle+xc+1*(szbp)*mod(ycind+1,2))./scale1,(magjitter*rand(1)+.5+szbp*y_circle+yc)./scale2,'linewidth',1,'color','c');
    end
end


%% Plot the RGC mosaic

for i = 1:length(ir.mosaic{1}.cellLocation); 
    loc(i,:) = irparams.inputScale.*ir.mosaic{1}.cellLocation{i}; 
end;
% figure; scatter(loc(:,1),loc(:,2),100); axis equal

rfRad = 1*ir.mosaic{1}.rfDiaMagnitude/2;
circle_samples = 0:0.05:2*pi;
c1(1,:) = rfRad*cos(circle_samples);
c1(2,:) = rfRad*sin(circle_samples);

% figure; 
hold on;
magrand = 3.5;
for i = 1:length(ir.mosaic{1}.cellLocation)
%     plot((magrand*rand(1)+loc(i,2)+(c1(1,2:end)+(rand(1,length(c1(1,2:end)),1)-.5) )-1*rfRad-0*ir.mosaic{1}.cellLocation{1,1}(1))./scale1,(magrand*rand(1)+loc(i,1)+(c1(2,2:end)+(rand(1,length(c1(1,2:end)))-.5))-1*rfRad-0*ir.mosaic{1}.cellLocation{1,1}(2))./scale2,'m','linewidth',3); axis equal
    
    plot((magrand*rand(1)+loc(i,2)+c1(1,2:end)-1*rfRad-0*ir.mosaic{1}.cellLocation{1,1}(1))./scale1,(magrand*rand(1)+loc(i,1)+c1(2,2:end)-1*rfRad-0*ir.mosaic{1}.cellLocation{1,1}(2))./scale2,'m','linewidth',4); axis equal
%         plot(loc(i,2)+c1(1,2:end)+1*rfRad-1*ir.mosaic{1}.cellLocation{1,1}(1),loc(i,1)+c1(2,2:end)+rfRad-2*ir.mosaic{1}.cellLocation{1,1}(2),'m','linewidth',3); axis equal
%     plot(loc(i,1)+c1(1,2:end)+2*rfRad-1*ir.mosaic{1}.cellLocation{1,1}(1),loc(i,2)+c1(2,2:end)+rfRad-1*ir.mosaic{1}.cellLocation{1,1}(2),'m','linewidth',3); axis equal
%     plot(loc(i,2)+c1(2,2:end)-16.5,loc(i,1)+c1(1,2:end)-16.5,'m','linewidth',3); axis equal
    
end
axis equal
% axis([1 94 1 74])
title(sprintf('%s simulated mosaic at %1.1f\\circ Ecc',cellType(1:end-3),ecc/.3));
set(gca,'fontsize',16)
xlabel(sprintf('Distance (\\mum)')); 
ylabel(sprintf('Distance (\\mum)'));

% view(0,-90);