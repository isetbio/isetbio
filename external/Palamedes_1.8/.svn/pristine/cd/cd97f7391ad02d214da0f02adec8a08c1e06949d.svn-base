%
%PAL_MLDS_SimulateObserver  Generate simulated responses in MLDS setting
%
%syntax: NumGreater = PAL_MLDS_SimulateObserver(Stim, OutOfNum, ...
%   PsiValues, SDnoise)
%
%Input:
%   'Stim': Stimulus set (see e.g., PAL_MLDS_GenerateStimList)
%
%   'OutOfNum': Number of trials to be simulated for each row in 'Stim'
%
%   'PsiValues': PsiValues characterizing simulated observer.
%
%   'SDnoise': internal noise of simulated observer.
%
%Output:
%   'NumGreater': For each row in 'Stim' the number of 'greater than'
%       responses.
%
%Example:
%
%   Stim = PAL_MLDS_GenerateStimList(2, 4, 1, 2);
%   OutOfNum = ones(1,size(Stim,1));
%   PsiValues = [0 1/3 2/3 1];
%   SDnoise = .5;
%   NumGreater = PAL_MLDS_SimulateObserver(Stim, OutOfNum, PsiValues, ...
%       SDnoise)
%
%   might generate:
%
%   NumGreater =
%
%     0     1     0     1     0     0
%
%   Tip: use PAL_MLDS_GroupTrialsbyX to combine like trials.
%
%Introduced: Palamedes version 1.0.0 (NP)
%Modified: Palamedes version 1.6.3 (see History.m)

function NumGreater = PAL_MLDS_SimulateObserver(Stim, OutOfNum, PsiValues, SDnoise)

if size(Stim,2) == 4
    D = (PsiValues(Stim(:,2))-PsiValues(Stim(:,1)))-(PsiValues(Stim(:,4))-PsiValues(Stim(:,3)));
end
if size(Stim,2) == 3
    D = (PsiValues(Stim(:,2))-PsiValues(Stim(:,1)))-(PsiValues(Stim(:,3))-PsiValues(Stim(:,2)));
end
if size(Stim,2) == 2
    D = (PsiValues(Stim(:,1))-PsiValues(Stim(:,2)));
end

NumGreater = zeros(1,length(D));

for Level = 1:length(D)
    Z_D = D(Level)./SDnoise;
    pFirst = .5 + .5*(1-erfc(Z_D./sqrt(2)));
    Pos = rand(OutOfNum(Level),1);
    Pos(Pos < pFirst) = 1;
    Pos(Pos ~= 1) = 0;
    NumGreater(Level) = sum(Pos);
end