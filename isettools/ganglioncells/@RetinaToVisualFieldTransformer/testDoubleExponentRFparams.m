function testDoubleExponentRFparams()

    r = -0.2:0.001:0.2;
    Kc = 1;
    Rc = 0.01;
    Rwide = 0.3;
    KsKcPeakRatio = 0.1;


    % From the 4 cells in Figure 6 of Packer & Dacey (2002)
    RnarrowToRwideRatios  = [152/515 170/718 115/902 221/1035];
    NWvolumeRatios        = [1.0     0.8     0.3     0.2];

    figure(1); clf
    for H1index = 1:4
        
        [Kwide, Knarrow, Rnarrow] = RetinaToVisualFieldTransformer.H1doubleExponentRFparams(...
            Kc, Rwide, KsKcPeakRatio, NWvolumeRatios(H1index), RnarrowToRwideRatios(H1index));

        centerRF = Kc * exp(-(r/Rc).^2);
        wideSurroundRF = Kwide * exp(-2.3*abs(r)/Rwide) ;
        narrowSurroundRF =  Knarrow * exp(-2.3*abs(r)/Rnarrow);

        subplot(1,4,H1index);
        plot(r, centerRF, 'r-', 'LineWidth', 1.5); hold on;
        plot(r, wideSurroundRF+narrowSurroundRF, 'b-', 'LineWidth', 1.5);
        plot(r, -narrowSurroundRF, 'k-', 'LineWidth', 1.5);
        plot(r, -wideSurroundRF, 'k--', 'LineWidth', 1.5);
        plot(r, centerRF-(wideSurroundRF+narrowSurroundRF), 'm-', 'LineWidth', 1.5);
    end
end
