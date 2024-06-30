% The Aguilar and Stiles rod saturation data.  
% 
% These data provide a quantitative evaluation of how much contrast one
% would need to see a rod signal in the periphery.
%
% Reference: Aguilar, M., and W. S. Stiles. 1954. “Saturation of the Rod
% Mechanism of the Retina at High Levels of Stimulation.” Optica Acta:
% International Journal of Optics 1 (1): 59–65.
%
% Stimuli:  See the paper for spectra.  9 deg target, 20 deg diameter
% field. Used test rays at the edge of the pupil for S-C effect advantage.
% Field went straight in.
%
% Summary:  Threshold is about 10 percent of the field intensity for most
% of the range.  Then at about 10 cd/m2 threshold contrast starts to
% increase, reaching about 100 percent at about 250 cd/m2. They could not
% set a threshold above this level.
%
% The characteristics of all the curves are similar : (a) an initial
% horizontal section (corresponding to the absolute threshold), (b) an
% ascending linear section covering about four log units of field intensity
% and with a gradient d log T/d log F of about 0.95, (c) a section of much
% increased gradient (about 2.7) which sets in fairly sharply at log F
% equal to about 2, (d) a final section of low gradient (about 0.5).
%
% They attribute (d) to the cones taking over.  
% 
% They also plot the Weber/Fechner fraction, as below.
%
% Their conclusion is that an uper limit on the saturation is about 3.6 log
% scotopic trolands, which is about 300 cd/m2.
%
%    10^(3.7-delta)
%
% Of note and to check:  A&S calculate that at 100 scotopic trolands each
% rod captures about 220 photons per second. (Page 64).  They were unaware
% that rods do not spike.
%
% Horiguchi et al. had a background level of about 2000 cd/m2.
% Others?
%
% NOTE
%
% BW used grabit with the PNG in this directory.  Figure 3 from their
% paper. 
%
% Subject A is correct w.r.t. y-axis
%
% Subjects B,C and D are displaced vertically by 0.5,1.0, and 1.5 log
% units.
%


%% If you can't find the data on your path, here they are.

% chdir(fullfile(isetbioRootPath,'data','datafiles','human','aguilarstiles'));

%%
load("asSubjectA.mat",'asSubjectA');
% I missed one point at negative infinity
asSubjectA = [-6 -4; asSubjectA];

load("asSubjectB.mat",'asSubjectB');
load("asSubjectC.mat",'asSubjectC');
load("asSubjectD.mat",'asSubjectD');

contrastA =  100* 10.^(asSubjectA(2:end,2) - asSubjectA(2:end,1));
contrastB =  100* 10.^(asSubjectB(2:end,2) - asSubjectB(2:end,1) - 0.5);
contrastC =  100* 10.^(asSubjectC(2:end,2) - asSubjectC(2:end,1) - 1.0);
contrastD =  100* 10.^(asSubjectD(2:end,2) - asSubjectD(2:end,1) - 1.5);

%% These are log-log plots, as in the original

ieNewGraphWin([],'tall');
hold on;
plot(asSubjectA(:,1),asSubjectA(:,2),'-o');
plot(asSubjectB(:,1),asSubjectB(:,2),'-o');
plot(asSubjectC(:,1),asSubjectC(:,2),'-o');
plot(asSubjectD(:,1),asSubjectD(:,2),'-o');
grid on;
xlabel('Log field intensity (scotopic trolands)');
ylabel('Log increment threshold (scotopic trolands)')
identityLine;

% Deviation from Weber's Law starts around 250 cd
xWeber = log10(250);
lWeber = line([ xWeber, xWeber],[-6 6],'Color','k'); 

% Saturation at about 2500
xSat = log10(4000);
lSat = line([ xSat, xSat],[-6 6],'Color','k'); 

%% Odd that subject A has such a low threshold compared to the others

ieNewGraphWin([],'tall');
hold on;
plot(asSubjectA(:,1),asSubjectA(:,2),'-o');
plot(asSubjectB(:,1),asSubjectB(:,2) - 0.5,'-o');
plot(asSubjectC(:,1),asSubjectC(:,2) - 1,'-o');
plot(asSubjectD(:,1),asSubjectD(:,2) - 1.5,'-o');
grid on;
xlabel('Log field intensity (scotopic trolands)');
ylabel('Log increment threshold (scotopic trolands)')

%% Plot vs cd/m2. 

%
% A-S say that 2000-5000 scotopic trolands corresponds to 120-300 cd/m2

delta = log10(5000/300);
asSubjectAcd = asSubjectA - delta;
asSubjectBcd = asSubjectB - delta;
asSubjectCcd = asSubjectC - delta;
asSubjectDcd = asSubjectD - delta;

ieNewGraphWin([],'tall');
hold on;
plot(asSubjectAcd(:,1),asSubjectAcd(:,2),'-o');
plot(asSubjectBcd(:,1) ,asSubjectBcd(:,2) - 0.5, '-o');
plot(asSubjectCcd(:,1),asSubjectCcd(:,2) - 1,'-o');
plot(asSubjectDcd(:,1),asSubjectDcd(:,2) - 1.5,'-o');
grid on;
xlabel('Log field intensity (cd/m^2)');
ylabel('Log increment threshold (cd/m^2)');


% Deviation from Weber's Law starts around 250 scotopic trolands
xWeber = log10(250);
lWeber = line([ xWeber - delta, xWeber - delta],[-6 6],'Color','k'); 

% Saturation at about 4000 scotop;ic trolands
xSat = log10(4000);
lSat = line([ xSat - delta, xSat - delta],[-6 6],'Color','k'); 

%% The contrast thresholds at different cd/m2 levels

ieNewGraphWin;
hold on;
plot(asSubjectAcd(2:end,1),contrastA)
plot(asSubjectBcd(2:end,1),contrastB)
plot(asSubjectCcd(2:end,1),contrastC)
plot(asSubjectDcd(2:end,1),contrastD)
grid on;

xlabel('Log field intensity (cd/m2)');
ylabel('Threshold contrast (percent)');
set(gca,'yscale','log')

% Subtract delta to make it cd/m2.

% Deviation from Weber's Law starts around 250 scotopic trolands
xWeber = log10(250);
lWeber = line([ xWeber - delta, xWeber - delta],[5 1000],'Color','k'); 

% Saturation at about 4000 scotop;ic trolands
xSat = log10(4000);
lSat = line([ xSat - delta, xSat - delta],[5 1000],'Color','k'); 

%% Evaluate the spline at new points (optional)

ieNewGraphWin([],'tall');
tiledlayout(2,2);
nexttile
xq = linspace(asSubjectAcd(2,1), asSubjectAcd(end,1), 100); % Creates num_points between min and max of x
s = spline(asSubjectAcd(2:end,1),asSubjectAcd(2:end,2),xq);
plot(asSubjectAcd(2:end,1),asSubjectAcd(2:end,2),'-o',xq,s,'--');
identityLine; 
xlabel('Log Field (cd/m2)'); ylabel('Log Threshold (cd/m2)'); grid on

nexttile
xq = linspace(asSubjectBcd(2,1), asSubjectBcd(end,1)-0.5, 100); % Creates num_points between min and max of x
s = spline(asSubjectBcd(2:end,1),asSubjectBcd(2:end,2)-0.5,xq);
plot(asSubjectBcd(2:end,1),asSubjectBcd(2:end,2)-0.5,'-o',xq,s,'--');
identityLine
xlabel('Log Field (cd/m2)'); ylabel('Log Threshold (cd/m2)'); grid on

nexttile
xq = linspace(asSubjectCcd(2,1), asSubjectCcd(end,1)-1, 100); % Creates num_points between min and max of x
s = spline(asSubjectCcd(2:end,1),asSubjectCcd(2:end,2)-1,xq);
plot(asSubjectCcd(2:end,1),asSubjectCcd(2:end,2)-1,'-o',xq,s,'--');
identityLine
xlabel('Log Field (cd/m2)'); ylabel('Log Threshold (cd/m2)'); grid on

nexttile
xq = linspace(asSubjectDcd(2,1), asSubjectDcd(end,1)- 1.5, 100); % Creates num_points between min and max of x
s = spline(asSubjectDcd(2:end,1),asSubjectDcd(2:end,2)-1.5,xq);
plot(asSubjectDcd(2:end,1),asSubjectDcd(2:end,2)-1.5,'-o',xq,s,'--');
identityLine
xlabel('Log Field (cd/m2)'); ylabel('Log Threshold (cd/m2)'); grid on

%%