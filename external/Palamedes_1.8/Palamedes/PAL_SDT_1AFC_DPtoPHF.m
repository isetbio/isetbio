%
% PAL_SDT_1AFC_DPtoPHF converts d'(d-prime) for a value of criterion C into 
% proportion hits and proportion false alarms for a 1AFC 
% (one-alternative-forced-choice) task, e.g. a Yes/No or symmetric single-
% alternative task
%
% Syntax: [pHF]=PAL_SDT_1AFC_DPtoPHFv2(dP,C,{optional argument});
% 
% returns a Nx2 matrix of N proportion hits and proportion false alarms 
% ('pHF') for a scalar or N-length vector of d' ('dP'), defined in the 
% range 0<d'<inf, criterion C ('C'), defined in the range -inf<C<inf. 
%
% Optional argument: an N-length vector 'R' of the ratios of noise (or 
% signal 1) to signal (or signal 2) SDs (standard deviations) 
% corresponding to each d-prime, in the range 0>R>inf.  If the option is 
% not included R defaults to 1.
%
% Example:
%
% [pHF]=PAL_SDT_1AFC_DPtoPHF([.0 1 10],[-2.5 0 4.0],'ratioSD',[1 1 1])
%
% pHF =
%
%    0.9938    0.9938
%    0.6915    0.3085
%    0.8413         0
%
% The example input arguments are three N=3 vectors of d', C and R.  
% The first column of the 3 x 2 matrix output gives the resulting 
% proportion of hits and the second column the corresponding proportion of 
% false alarms
%
% Introduced: Palamedes version 1.0.0 (FK)
% Modified: Palamedes version 1.6.0, 1.6.3 (see History.m)

function pHF = PAL_SDT_1AFC_DPtoPHF(dP,C,varargin)

R = ones(1,length(dP)); % default SD ratios

if ~isempty(varargin)
    valid = 0;
    if strncmpi(varargin{1}, 'ratioSDvalue',5)
        R = varargin{2};
        valid = 1;
    end
    if valid == 0
        warning('PALAMEDES:invalidOption','%s is not a valid option. Ignored.',varargin{1});
    end        
end   

k=sqrt(2./(1+R.^2));
zH=(dP./k).*(1-R./(1+R))-C;
zF=-(dP./k).*(1./(1+R))-C./R;

pHF(:,1)=PAL_ZtoP(zH);
pHF(:,2)=PAL_ZtoP(zF);