%
% RGCmodels.BenardeteKaplan1997.figure7CenterSurroundFilterParams
%
function params = figure7CenterSurroundFilterParams()

    nL_tL_product = 48;

  
    params.centerIR.pVector(1) = 184.20;
    params.surroundIR.pVector(1)  = -117.9;

    params.centerIR.pVector(2) = 4.0/1000;
    params.surroundIR.pVector(2) = 4.03/1000;

    params.centerIR.pVector(3) = 0.70;
    params.surroundIR.pVector(3) = 0.41;

    params.centerIR.pVector(4) = 37.30/1000;
    params.surroundIR.pVector(4) = 42.49/1000;


    params.centerIR.pVector(6) = 32;
    params.surroundIR.pVector(6) = 93;

    params.centerIR.pVector(5) = nL_tL_product/params.centerIR.pVector(6)/1000;
    params.surroundIR.pVector(5) = nL_tL_product/params.surroundIR.pVector(6)/1000;

    params.centerIR.pVector(7) = 1;
    params.surroundIR.pVector(7) = 1;

    

end


