%
% RGCmodels.PurpuraTranchinaKaplanShapley1990.table1FilterParams
%

function params = table1FilterParams(whichCellAndIllumination)

    switch (whichCellAndIllumination)
        case 'P8/25_@4600Trolands'
            params = [21  4   0.44   38.0/1000   3.0/1000   1.3/1000   10.0/1000   38/1000];

        case 'P8/25_@460Trolands'
            params = [21  4   0.40   52.0/1000   2.9/1000   1.3/1000   16.0/1000   43/1000];

        case 'P8/25_@150Trolands'
            params = [21  4   0.38   70.0/1000   3.7/1000   1.5/1000   26.0/1000   56/1000];

        case 'P8/25_@46Trolands'
            params = [21  4   0.40   42.0/1000   3.5/1000   3.6/1000   20.0/1000   58/1000];

        case 'P8/25_@15Trolands'
            params = [15  4   0.30   67.0/1000   2.8/1000   2.5/1000   38.0/1000   71/1000];


        case 'P24/9_@1900Trolands'
            params = [30  4   0.55   30.0/1000   2.0/1000   0.9/1000   9.6/1000    50/1000];

        case 'P24/9_@120Trolands'
            params = [21  4   0.44   45.0/1000   2.5/1000   1.3/1000   15.0/1000    63/1000];

        case 'P24/9_@3Trolands'
            params = [21  4   0.21  143/1000     2.7/10000   2.3/1000  52.0/1000    85/1000];

        case 'P24/9_@0.3Trolands'
            params = [30  4   0.14  562/1000    401/1000    2.1/1000   23/1000     133/1000];


        case 'P8/14_@4600Trolands'
            params = [24  0   0.22   90/1000      52/1000    1.8/1000   nan        40/1000];

        case 'P8/14_@460Trolands'
            params = [24  0   0.16   95/1000      53/1000    2.4/1000   nan        55/1000];

        case 'P26/10_@1900Trolands'
            params = [30  0   0.51   70/1000      36/1000    1.4/1000   nan        39/1000];

        case 'P26/10_@120Trolands'
            params = [21  0   0.29  107/1000      72/1000    2.5/1000   nan        49/1000];

        otherwise
            error('No data for PurpuraTranchinaKaplanShapley (1990)''%s'' cell', whichCellAndIllumination);

    end

    idx = [3 8 4 5 6 7 1 2];
    params = params(idx);

end


