function color_var_disp_sc = rgbGaborColorOpponentNormalized(params)

% clear
% figure;

% s_iso = [0 0 1];
% lm_red = [0.71 -0.71 0];
% lms_orange = [0.14 -0.14 -0.98];
% lms_magenta = [0.14 -0.14 0.98];

s_iso = [0.05 0.10 0.99];
lm_red = [0.86 -0.50 -0.05];
lms_orange = [0.12 -0.19 -0.97];
lms_magenta = [0.12 -0.19 0.97];

sz = params.row;
im = zeros(sz,sz,3);

for color_val = params.color_val%1:4
    
switch color_val
    case 1
        color_select1 = s_iso;
    case 2
        color_select1 = lm_red;
        
    case 3
%         color_select1 = lms_magenta;
        color_select1 = s_iso + lm_red;
        color_select1 = color_select1./norm(color_select1);
    case 4
%         color_select1 = lms_orange;

        color_select1 = s_iso - lm_red;
        color_select1 = color_select1./norm(color_select1);
end

theta = 2*pi*(1:sz)./(sz/params.freq);
phi = params.ph; % 2*pi*(t-1)/10;

color_var = zeros(sz,sz,3);

for r = 1:sz
for i = 1:3
    
color_var(:,r,i) = color_select1(i)*cos(theta + phi);

end
end

% g = fspecial('gauss',64,.1);
gaussian_window = repmat(gausswin(sz,(sz*params.GaborFlag))*gausswin(sz,(sz*params.GaborFlag))',[1 1 3]);
color_var_disp = 0.5*(1+params.contrast*color_var.*gaussian_window);

scfac = sqrt(sum(reshape(color_var_disp.^2,sz*sz,3),2));
scfacrs = reshape(scfac,sz,sz);
scfacrsrm = repmat(scfacrs,1,1,3);
color_var_disp_sc = (sqrt(3)/2)*color_var_disp./scfacrsrm;

% figure;
% % subplot(2,2,color_val);
% imagesc(color_var_disp);
% 
% figure; hist(color_var_disp(:),50)

end


% scene = sceneFromFile(color_var_disp,'rgb')
% vcAddAndSelectObject(scene);
% sceneWindow
