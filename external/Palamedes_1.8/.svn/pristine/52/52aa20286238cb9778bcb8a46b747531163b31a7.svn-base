%
%PAL_PFML_PtoT      Transform PF parameters according to model specified in
%   structure FM. ~inverse: PAL_PFML_TtoP
%
%Syntax: [thetas thetasID FM] = PAL_PFML_PtoT(params, FM)
%   'thetas' contains only free thetas, 'thetasID' identifies thetas as
%   specifying either alphas (1), betas (2), gammas (3), or lambdas (4).
%
%Internal function
%
%Introduced: Palamedes version 1.0.0 (NP)
%Modified: Palamedes version 1.1.0, 1.6.3 (see History.m)

function [thetas, thetasID, FM] = PAL_PFML_PtoT(params, FM)

thetas = [];
thetasID = [];

switch PAL_whatIs(FM.argsA)
    case 0
        thetas = [thetas];
        thetasID = [thetasID];
    case 1
        thetas = [thetas params(:,1)'/FM.argsA];
        thetasID = [thetasID ones(size(params(:,1)'/FM.argsA))];        
    case 2
        switch lower(FM.argsA(1:4))
            case {'unco'}
                FM.argsA = PAL_Contrasts(size(params,1),'iden');
                thetas = [thetas params(:,1)'/FM.argsA];
                thetasID = [thetasID ones(size(params(:,1)'/FM.argsA))];
            case {'cons'}
                FM.argsA = ones(1,size(params,1));
                thetas = [thetas params(:,1)'/FM.argsA];
                thetasID = [thetasID ones(size(params(:,1)'/FM.argsA))];                
            case {'fixe'}
                FM.argsA = [];
                thetas = [thetas];
                thetasID = [thetasID];
        end
    case 4
        thetas = [thetas FM.argsA.paramsValuesA(FM.argsA.paramsFreeA==1)];
        thetasID = [thetasID ones(1,sum(FM.argsA.paramsFreeA==1))];
end                        

switch PAL_whatIs(FM.argsB)
    case 0
        thetas = [thetas];
        thetasID = [thetasID];
    case 1
        thetas = [thetas params(:,2)'/FM.argsB];
        thetasID = [thetasID 2*ones(size(params(:,1)'/FM.argsB))];
    case 2
        switch lower(FM.argsB(1:4))
            case {'unco'}
                FM.argsB = PAL_Contrasts(size(params,1),'iden');
                thetas = [thetas params(:,2)'/FM.argsB];
                thetasID = [thetasID 2*ones(size(params(:,1)'/FM.argsB))];
            case {'cons'}
                FM.argsB = ones(1,size(params,1));
                thetas = [thetas params(:,2)'/FM.argsB];
                thetasID = [thetasID 2*ones(size(params(:,1)'/FM.argsB))];
            case {'fixe'}
                FM.argsB = [];
                thetas = [thetas];
                thetasID = [thetasID];
        end
    case 4
        thetas = [thetas FM.argsB.paramsValuesB(FM.argsB.paramsFreeB==1)];
        thetasID = [thetasID 2*ones(1,sum(FM.argsB.paramsFreeB==1))];
end                        

switch PAL_whatIs(FM.argsG)
    case 0
        thetas = [thetas];
        thetasID = [thetasID];
    case 1        
        thetas = [thetas params(:,3)'/FM.argsG];
        thetasID = [thetasID 3*ones(size(params(:,1)'/FM.argsG))];
    case 2
        switch lower(FM.argsG(1:4))
            case {'unco'}
                FM.argsG = PAL_Contrasts(size(params,1),'iden');
                thetas = [thetas params(:,3)'/FM.argsG];
                thetasID = [thetasID 3*ones(size(params(:,1)'/FM.argsG))];
            case {'cons'}
                FM.argsG = ones(1,size(params,1));
                thetas = [thetas params(:,3)'/FM.argsG];
                thetasID = [thetasID 3*ones(size(params(:,1)'/FM.argsG))];
            case {'fixe'}
                FM.argsG = [];
                thetas = [thetas];
                thetasID = [thetasID];
        end
    case 4
        thetas = [thetas FM.argsG.paramsValuesG(FM.argsG.paramsFreeG==1)];
        thetasID = [thetasID 3*ones(1,sum(FM.argsG.paramsFreeG==1))];
end                        

switch PAL_whatIs(FM.argsL)
    case 0
        thetas = [thetas];
        thetasID = [thetasID];
    case 1
        thetas = [thetas params(:,4)'/FM.argsL];
        thetasID = [thetasID 4*ones(size(params(:,1)'/FM.argsL))];
    case 2
        switch lower(FM.argsL(1:4))
            case {'unco'}
                FM.argsL = PAL_Contrasts(size(params,1),'iden');
                thetas = [thetas params(:,4)'/FM.argsL];
                thetasID = [thetasID 4*ones(size(params(:,1)'/FM.argsL))];
            case {'cons'}
                FM.argsL = ones(1,size(params,1));
                thetas = [thetas params(:,4)'/FM.argsL];
                thetasID = [thetasID 4*ones(size(params(:,1)'/FM.argsL))];
            case {'fixe'}
                FM.argsL = [];
                thetas = [thetas];
                thetasID = [thetasID];
        end
    case 4
        thetas = [thetas FM.argsL.paramsValuesL(FM.argsL.paramsFreeL==1)];
        thetasID = [thetasID 4*ones(1,sum(FM.argsL.paramsFreeL==1))];
end                        