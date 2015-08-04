function [cur] = riekeConeAdaptation(stimulus)
% Functional form of rieke cone adaptation model based on UW code.
%
% This function is temporary and for exploration ... See notes below
%
% Example:
%  stimulus = zeros(1,1000); stimulus(100:105) = 70000;
%  [cur] = riekeConeAdaptation(stimulus);
%  vcNewGraphWin; plot(cur); grid on
%
% See Also:
%
% HJ/BW ISETBIO Team, Copyright 2014


%% Programming Notes
%
% The initial m-file from UW (s_fredConeModel.m) had two sections that
% computed the current from the stimulus.  The formulae in the two sections
% differed at a key step - the computation of the current.  
%
% One group had computation as:
%
%   cur = -cgmp2cur * g.^h * 2 ./ (2 + cslow ./ cdark);
%
% The second had the computation as:
%
%   cur = -cgmp2cur * g.^h ./ (1 + cslow ./ cdark);
%
% There is a '2' that differs.  We don't know which is right.
% 
% Also, the second calculation in the original confounded the stimulus
% creation with the differential equations, making that form of the
% implementation highly specialized and not useful for general calculation.
%
% Finally, the powerpoint description of the differential equations was not
% precisely the same as either of the calculations in s_fredConeModel.  So,
% that needs to be sorted out also.
%
% The formulae on the slide look like this, and they are implemented in the
% riekeAdaptTemporal.  See the document uploaded by HJ for the derivation
% of the formula. There are a couple of differences.  In particular one
% difference is in how cur2ca is calculated (p.q). HJ says the slides
% resolve to this formula
%
%   q    = 2 * beta * cdark / (k * gdark^h);  % This is cur2ca
%
% but the s_fredConeModel has this formula
%
%   cur2ca = beta * cdark / (cgmp2cur * gdark^3);
%
%
% This function is intended to serve as a basis for discussion as we
% explore the calculations and resolve the differences.  In the end, we
% will use this to settle on the implementation in riekeAdaptTemporal.m.


%% Initialize parameters

sigma = 100;			% rhodopsin activity decay rate constant (1/sec)
phi = 50;				% phosphodiesterase activity decay rate constant (1/sec)
eta = 100;				% phosphodiesterase activation rate constant (1/sec)
r(1) = 0;				% initial condition for r
p(1) = eta/phi;			% initial condition for p
TimeStep = 0.001;		% time step 
gdark = 35;				% concentration of cGMP in darkness
cgmp2cur = 10e-3;		% constant relating cGMP to current
cdark = 0.5;			% dark calcium concentration
beta = 50;				% rate constant for calcium removal in 1/sec
betaSlow = 2;
hillcoef = 4;			% cooperativity
hillaffinity = 0.35;		% affinity
NumPts = length(stimulus);

clear g s c p cslow;

% This differs from our calculation by a factor of 2, see riekeInit.m
% q    = 2 * beta * cdark / (k * gdark^h);
cur2ca = beta * cdark / (cgmp2cur * gdark^3);				        % get q using steady state
smax = eta/phi * gdark * (1 + (cdark / hillaffinity)^hillcoef);		% get smax using steady state

% initial conditions
g(1) = gdark;
s(1) = gdark * eta/phi;		
c(1) = cdark;
p(1) = eta/phi;
cslow(1) = cdark;


%%  First loop from s_fredConeModel script

% Should zero out but allocate space for all the p,c, ... and so forth
clear g s c p cslow;

g(1) = gdark;
s(1) = gdark * eta/phi;		
c(1) = cdark;
p(1) = eta/phi;
cslow(1) = cdark;

h = 3;   % From the slide

% Solve difference equations
for pnt = 2:NumPts

    r(pnt) = r(pnt-1) + TimeStep * (-sigma * r(pnt-1) + stimulus(pnt));
    p(pnt) = p(pnt-1) + TimeStep * (r(pnt-1) + eta - phi * p(pnt-1));
	c(pnt) = c(pnt-1) + TimeStep * (cur2ca * cgmp2cur * g(pnt-1)^h - beta * c(pnt-1));
	cslow(pnt) = cslow(pnt-1) - TimeStep * (betaSlow * (cslow(pnt-1)-c(pnt-1)));
	s(pnt) = smax / (1 + (c(pnt) / hillaffinity)^hillcoef);
	g(pnt) = g(pnt-1) + TimeStep * (s(pnt-1) - p(pnt-1) * g(pnt-1));



end

% determine current change
cur = -cgmp2cur * g.^h * 2 ./ (2 + cslow ./ cdark);

% s_fredCOneModel has two different calculations for current.
% Here is the second calculation.  We don't know why they differ.
%
% cur = -cgmp2cur * g.^h ./ (1 + cslow ./ cdark);

return

%% This is the second loop ... not used in this function.

clear g s c p cslow;
g(1) = gdark;
s(1) = gdark * eta/phi;		
c(1) = cdark;
p(1) = eta/phi;
cslow(1) = cdark;

for pnt = 2:NumPts
	if (pnt <= NumPts/4)
        r(pnt) = 0; 
    end
    if ((pnt > NumPts/4) && (pnt < (NumPts/4+StmPts)))
        r(pnt) = r(pnt-1) + TimeStep * (-sigma * r(pnt-1) + StmAmp);
        
    end
    if (pnt >= (NumPts/4+StmPts))
        r(pnt) = r(pnt-1) + TimeStep * (-sigma * r(pnt-1));
    end
    
    p(pnt) = p(pnt-1) + TimeStep * (r(pnt-1) + eta - phi * p(pnt-1));
	c(pnt) = c(pnt-1) + TimeStep * (cur2ca * cgmp2cur * g(pnt-1)^3 - beta * c(pnt-1));
	cslow(pnt) = cslow(pnt-1) - TimeStep * (betaSlow * (cslow(pnt-1)-c(pnt-1)));
	s(pnt) = smax / (1 + (c(pnt) / hillaffinity)^hillcoef);
	g(pnt) = g(pnt-1) + TimeStep * (s(pnt-1) - p(pnt-1) * g(pnt-1));
end

% determine current change.
% Notice that this formula differs from the formula in the previous cell.
% We don't know why
cur2 = -cgmp2cur * g.^3 ./ (1 + cslow ./ cdark);



