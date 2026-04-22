%
% RGCmodels.BenardeteKaplan1992.figure6CenterSurroundFilterParams
%
function params = figure6CenterSurroundFilterParams(whichCell)

    nL_tL_product = 48;

    switch (whichCell)
        case 'ON'
            params.centerIR.pVector(1) = 184.20;
            params.surroundIR.pVector(1)  = 125.33;

            params.centerIR.pVector(2) = 4.0/1000;
            params.surroundIR.pVector(2) = 4.0/1000;

            params.centerIR.pVector(3) = 0.69;
            params.surroundIR.pVector(3) = 0.56;

            params.centerIR.pVector(4) = 18.61/1000;
            params.surroundIR.pVector(4) = 33.28/1000;

            params.centerIR.pVector(6) = 38;
            params.surroundIR.pVector(6) = 124;

            %params.centerIR.pVector(5) = 1.23/1000;
            %params.surroundIR.pVector(5) = 0.42/1000;

            params.centerIR.pVector(5) = nL_tL_product/params.centerIR.pVector(6)/1000;
            params.surroundIR.pVector(5) = nL_tL_product/params.surroundIR.pVector(6)/1000;
 
            params.centerIR.pVector(7) = 1;
            params.surroundIR.pVector(7) = 1;

      case 'OFF'
            params.centerIR.pVector(1) = 114.12;
            params.surroundIR.pVector(1) = 74.57;

            params.centerIR.pVector(2) = 3.5/1000;
            params.surroundIR.pVector(2) = 3.5/1000;

            params.centerIR.pVector(3) = 0.82;
            params.surroundIR.pVector(3) = 0.72;

            params.centerIR.pVector(4) = 24.9/1000;
            params.surroundIR.pVector(4) = 49.81/1000;

            params.centerIR.pVector(6)= 25;
            params.surroundIR.pVector(6) = 83;

            %params.centerIR.pVector(5)= 2.12/1000;
            %params.surroundIR.pVector(5) = 0.76/1000;

            params.centerIR.pVector(5) = nL_tL_product/params.centerIR.pVector(6)/1000;
            params.surroundIR.pVector(5) = nL_tL_product/params.surroundIR.pVector(6)/1000;


            params.centerIR.pVector(7) = 1;
            params.surroundIR.pVector(7) = 1;

        otherwise
            error('No data for Banardete&Kaplan (1992a)''%s'' cell', whichCell);
    end

end


