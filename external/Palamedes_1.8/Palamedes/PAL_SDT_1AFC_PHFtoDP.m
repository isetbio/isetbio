%
% PAL_SDT_1AFC_PHFtoDP converts proportion hits and proportion false 
% alarms into d'(d-prime) and criterion C for a 1AFC 
% (one-alternative-forced-choice) task, e.g. a Yes/No or symmetric single-
% alternative task
%
% Syntax: [dP C lnB pC]=PAL_SDT_1AFC_PHFtoDP(pHF,{optional argument});
% 
% returns a scalar or N-length vector of d' ('dP'), criterion C ('C'), 
% criterion lnBeta ('lnB') and proportion correct ('pC'), for an Nx2 input 
% matrix of N proportion hits and proportion false alarms ('pHF') defined 
% in the range 0<p<1 (p=proportion).
%
% Optional argument: an N length vector 'R' of the ratios of noise (or 
% signal 1) to signal (or signal 2) SDs (standard deviations) 
% corresponding to each pair of proportion hits and false alarms, in the 
% range 0>R>inf.  If the option is not included R defaults to 1.
%
% Example:
%
% [dP C lnB pC]=PAL_SDT_1AFC_PHFtoDP([.5 .3; .7 .2; .9 .1],'ratioSD',...
%   [1 1 1])
%
% returns:
%
% dP =
%
%    0.5244
%    1.3660
%    2.5631
%
%
% C =
%
%     0.2622
%     0.1586
%          0
% 
% 
% lnB =
% 
%     0.1375
%     0.2167
%          0
% 
% 
% pC =
% 
%     0.6000
%     0.7500
%     0.9000
%
% The example input argument consists of a 3 x 2 matrix, with each
% row (demarcated by a semi-colon) consisting of a proportion of hits and 
% a corresponding proportion of false alarms, and a N=3 vector of R. 
% The columns in the output are the resulting N=3 vectors of dP, C, lnB 
% and pC.  Note that if R~=1, lnB is not calculated and returns a NaN
%
% Introduced: Palamedes version 1.0.0 (FK)
% Modified: Palamedes version 1.6.0, 1.6.3 (see History.m)


function [dP, C, lnB, pC] = PAL_SDT_1AFC_PHFtoDP(pHF,varargin)

R = ones(size(pHF(:,1))); % default

if ~isempty(varargin)
    valid = 0;
    if strncmpi(varargin{1}, 'ratioSDvalue',5)
        R = varargin{2};
        R=R';
        valid = 1;
    end
    if valid == 0
        warning('PALAMEDES:invalidOption','%s is not a valid option. Ignored.',varargin{1});
    end        
end   

zH = PAL_PtoZ(pHF(:,1));
zF = PAL_PtoZ(pHF(:,2));

k = sqrt(2./(1+R.^2));
dP = k.*(zH - R.*zF);
C = -R./(1+R).*(zH+zF);

lnB = -0.5.*(zH.^2-zF.^2); %default
for i = 1:length(lnB)
    if R(i,1) ~= 1
        lnB(i,1) = NaN;
    end
end

pC = (pHF(:,1)+(1.0-pHF(:,2)))./2;