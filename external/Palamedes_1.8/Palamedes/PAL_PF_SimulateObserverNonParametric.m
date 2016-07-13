%
%PAL_PF_SimulateObserverNonParametric   Simulate observer characterized by
%   raw response proportions.
%
%syntax: NumPos = PAL_PF_SimulateObserverNonParametric(StimLevels, ...
%   NumPos, OutOfNum)
%
%Input:
%   'StimLevels': vector containing stimulus levels. Note that each entry
%       should be unique. Use PAL_PFML_GroupTrialsbyX to combine like 
%       stimulus levels.
%
%   'NumPos': Expected number of positive (e.g., 'correct') responses, for
%       each of the entries in 'StimLevels' (for most purposes, this would 
%       be the observed number of positive responses of a human observer).
%
%   'OutOfNum': Number of trials to be simulated for each of the entries in
%       'StimLevels'. 
%
%Output:
%   'NumPos': simulated number of positive responses for each entry in
%       'StimLevels'.
%
%Example: NumPos = PAL_PF_SimulateObserverNonParametric([1 3 10], ...
%   [25 60 125], [100 100 150]) might return:
%
%   NumPos =
%       
%       20  62  126
%
%Introduced: Palamedes version 1.0.0 (NP)

function NumPos = PAL_PF_SimulateObserverNonParametric(StimLevels, NumPos, OutOfNum)

I = NumPos ~= 0;

NumPos = NumPos(I);
OutOfNum = OutOfNum(I);
StimLevels = StimLevels(I);

pcorrect = NumPos./OutOfNum;
for Level = 1:length(StimLevels)
    Pos = rand(OutOfNum(Level),1);
    Pos(Pos < pcorrect(Level)) = 1;
    Pos(Pos ~= 1) = 0;
    NumPos(Level) = sum(Pos);
end

Temp = zeros(size(I));
Temp(I) = NumPos;
NumPos = Temp;