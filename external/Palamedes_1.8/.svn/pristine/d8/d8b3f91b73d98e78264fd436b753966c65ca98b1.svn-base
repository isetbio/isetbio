%
%PAL_PFML_TtoP      Transform theta parameters to PF parameters according 
%   to model specified in structure FM. More or less inverse function of
%   PAL_PFML_PtoT
%
%Syntax: params = PAL_PFML_TtoP(params, thetas, thetasID, FM)
%
%Internal function
%
%Introduced: Palamedes version 1.0.0 (NP)
%Modified: Palamedes version 1.1.0, 1.6.0 (see History.m)

function params = PAL_PFML_TtoP(params, thetas, thetasID, FM)

if ischar(FM.argsA)
    switch lower(FM.argsA(1:4))
        case 'cons'
            FM.argsA = ones(1,size(params,1));
        case 'unco'
            FM.argsA = PAL_Contrasts(size(params,1),'iden');
        case 'fixe'
            FM.argsA = [];
    end
end
if ischar(FM.argsB)
    switch lower(FM.argsB(1:4))
        case 'cons'
            FM.argsB = ones(1,size(params,1));
        case 'unco'
            FM.argsB = PAL_Contrasts(size(params,1),'iden');
        case 'fixe'
            FM.argsB = [];
    end
end
if ischar(FM.argsG)
    switch lower(FM.argsG(1:4))
        case 'cons'
            FM.argsG = ones(1,size(params,1));
        case 'unco'
            FM.argsG = PAL_Contrasts(size(params,1),'iden');
        case 'fixe'
            FM.argsG = [];
    end
end
if ischar(FM.argsL)
    switch lower(FM.argsL(1:4))
        case 'cons'
            FM.argsL = ones(1,size(params,1));
        case 'unco'
            FM.argsL = PAL_Contrasts(size(params,1),'iden');
        case 'fixe'
            FM.argsL = [];
    end
end

if isempty(FM.argsA)
    alphas = params(:,1)';
else
    if isstruct(FM.argsA)
        funcParams = zeros(1,length(FM.argsA.paramsValuesA));
        funcParams(FM.argsA.paramsFreeA == 1) = thetas(thetasID == 1);
        funcParams(FM.argsA.paramsFreeA == 0) = FM.argsA.paramsValuesA(FM.argsA.paramsFreeA == 0);
        alphas = FM.argsA.funcA(funcParams);
    else
        alphas = thetas(thetasID == 1)*FM.argsA;
        alphas(sum(FM.argsA~=0,1)==0) = params(sum(FM.argsA~=0,1)==0,1);               
    end
end
if isempty(FM.argsB)
    betas = params(:,2)';
else
    if isstruct(FM.argsB)
        funcParams = zeros(1,length(FM.argsB.paramsValuesB));
        funcParams(FM.argsB.paramsFreeB == 1) = thetas(thetasID == 2);
        funcParams(FM.argsB.paramsFreeB == 0) = FM.argsB.paramsValuesB(FM.argsB.paramsFreeB == 0);
        betas = FM.argsB.funcB(funcParams);
    else
        betas = thetas(thetasID == 2)*FM.argsB;    
        betas(sum(FM.argsB~=0,1)==0) = params(sum(FM.argsB~=0,1)==0,2);      
    end
end
if isempty(FM.argsG)
    gammas = params(:,3)';
else
    if isstruct(FM.argsG)
        funcParams = zeros(1,length(FM.argsG.paramsValuesG));
        funcParams(FM.argsG.paramsFreeG == 1) = thetas(thetasID == 3);
        funcParams(FM.argsG.paramsFreeG == 0) = FM.argsG.paramsValuesG(FM.argsG.paramsFreeG == 0);
        gammas = FM.argsG.funcG(funcParams);
    else
        gammas = thetas(thetasID == 3)*FM.argsG;  
        gammas(sum(FM.argsG~=0,1)==0) = params(sum(FM.argsG~=0,1)==0,3);        
    end
end
if isempty(FM.argsL)
    lambdas = params(:,4)';
else
    if isstruct(FM.argsL)
        funcParams = zeros(1,length(FM.argsL.paramsValuesL));
        funcParams(FM.argsL.paramsFreeL == 1) = thetas(thetasID == 4);
        funcParams(FM.argsL.paramsFreeL == 0) = FM.argsL.paramsValuesL(FM.argsL.paramsFreeL == 0);
        lambdas = FM.argsL.funcL(funcParams);
    else
        lambdas = thetas(thetasID == 4)*FM.argsL;    
        lambdas(sum(FM.argsL~=0,1)==0) = params(sum(FM.argsL~=0,1)==0,4);        
    end
end

params = [alphas' betas' gammas' lambdas'];