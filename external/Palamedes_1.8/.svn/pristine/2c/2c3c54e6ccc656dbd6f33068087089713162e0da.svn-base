%
%PAL_PFML_rangeTries    Set up a (number of conditions)x4 matrix of random
%   numbers to be used to create jitter on initial parameter search values
%   appropriate to the model that is to be fitted.
%
%Syntax: [multiplier] = PAL_PFML_rangeTries(MC, rangeTries)
%
%Internal function
%
%Introduced: Palamedes version 1.0.0 (NP)
%Modified: Palamedes version 1.1.0 (see History.m)

function [multiplier] = PAL_PFML_rangeTries(ArgsA, ArgsB, ArgsG, ArgsL, rangeTries)

if ~isstruct(ArgsA)
    if ischar(ArgsA)
        switch lower(ArgsA(1:4))
            case 'cons'
                ArgsA = ones(1,size(rangeTries,1));
            case 'unco'
                ArgsA = PAL_Contrasts(size(rangeTries,1),'iden');
            case 'fixe'
                ArgsA = [];
        end
    end
    if isempty(ArgsA)
        multiplier(1:size(rangeTries,1),1) = 0;
    else
        if size(ArgsA,1) >= 1 & size(ArgsA,1) < size(rangeTries,1);
            multiplier(1:size(rangeTries,1),1) = rand(1)-.5;
        else
            multiplier(1:size(rangeTries,1),1) = rand(size(rangeTries,1),1)-.5;
        end
    end
else
    multiplier(1:size(rangeTries,1),1) = 0;
end

if ~isstruct(ArgsB)
    if ischar(ArgsB)
        switch lower(ArgsB(1:4))
            case 'cons'
                ArgsB = ones(1,size(rangeTries,1));
            case 'unco'
                ArgsB = PAL_Contrasts(size(rangeTries,1),'iden');
            case 'fixe'
                ArgsB = [];
        end
    end
    if isempty(ArgsB)
        multiplier(1:size(rangeTries,1),2) = 0;
    else
        if size(ArgsB,1) >= 1 & size(ArgsB,1) < size(rangeTries,1);
            multiplier(1:size(rangeTries,1),2) = rand(1)-.5;
        else
            multiplier(1:size(rangeTries,1),2) = rand(size(rangeTries,1),1)-.5;
        end
    end
else
    multiplier(1:size(rangeTries,1),2) = 0;
end

if ~isstruct(ArgsG)
    if ischar(ArgsG)
        switch lower(ArgsG(1:4))
            case 'cons'
                ArgsG = ones(1,size(rangeTries,1));
            case 'unco'
                ArgsG = PAL_Contrasts(size(rangeTries,1),'iden');
            case 'fixe'
                ArgsG = [];
        end
    end
    if isempty(ArgsG)
        multiplier(1:size(rangeTries,1),3) = 0;
    else
        if size(ArgsG,1) >= 1 & size(ArgsG,1) < size(rangeTries,1);
            multiplier(1:size(rangeTries,1),3) = rand(1)-.5;
        else
            multiplier(1:size(rangeTries,1),3) = rand(size(rangeTries,1),1)-.5;
        end
    end
else
    multiplier(1:size(rangeTries,1),2) = 0;
end

if ~isstruct(ArgsL)
    if ischar(ArgsL)
        switch lower(ArgsL(1:4))
            case 'cons'
                ArgsL = ones(1,size(rangeTries,1));
            case 'unco'
                ArgsL = PAL_Contrasts(size(rangeTries,1),'iden');
            case 'fixe'
                ArgsL = [];
        end
    end
    if isempty(ArgsL)
        multiplier(1:size(rangeTries,1),4) = 0;
    else
        if size(ArgsL,1) >= 1 & size(ArgsL,1) < size(rangeTries,1);
            multiplier(1:size(rangeTries,1),4) = rand(1)-.5;
        else
            multiplier(1:size(rangeTries,1),4) = rand(size(rangeTries,1),1)-.5;
        end
    end
else
    multiplier(1:size(rangeTries,1),4) = 0;
end

