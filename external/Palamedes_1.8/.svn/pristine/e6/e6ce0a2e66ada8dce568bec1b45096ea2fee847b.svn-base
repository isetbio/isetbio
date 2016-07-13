%
%PAL_MLDS_GroupTrialsbyX    Combines like trials in MLDS setting
%
%syntax: [StimG NumGreaterG OutOfNumG] = PAL_MLDS_GroupTrialsbyX(Stim, ... 
%   NumGreater, OutOfNum)
%
%example:
%
%   Stim =  [1 2 3;
%            2 3 4;    
%            1 2 3]
%
%   NumGreater = [1 2 1];
%
%   OutOfNum = [1 3 3];
%
%   [StimG NumGreaterG OutOfNumG] = PAL_MLDS_GroupTrialsbyX(Stim, ...
%       NumGreater, OutOfNum) returns:
%StimG =
%
%     1     2     3
%     2     3     4
%
%NumGreaterG =
%
%     2     2
%
%OutOfNumG =
%
%     4     3
%
%Introduced: Palamedes version 1.0.0 (NP)

function [StimG, NumGreaterG, OutOfNumG] = PAL_MLDS_GroupTrialsbyX(Stim, NumGreater, OutOfNum)

if size(Stim,1) == 1
    Stim = Stim';
end

NumLevels = 1;

while size(Stim,1) ~= 0

    x = repmat(Stim(1,:), [size(Stim,1) 1]);
    xG(NumLevels,:) = Stim(1,:);
    match = sum(Stim == x,2) == size(Stim,2);
    OutOfNumG(NumLevels) = sum(OutOfNum(match));
    NumGreaterG(NumLevels) = sum(NumGreater(match));
    OutOfNum = OutOfNum(match == 0);
    NumGreater = NumGreater(match == 0);
    Stim = Stim(match == 0,:);
    NumLevels = NumLevels+1;
        
end

[StimG, I] = sortrows(xG);
NumGreaterG = NumGreaterG(I);
OutOfNumG = OutOfNumG(I);
if size(StimG,2) == 1
    StimG = StimG';
end