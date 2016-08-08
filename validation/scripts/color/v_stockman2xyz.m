function varargout = v_stockman2xyz(varargin)
%
% Test that colorTransformMatrix does the right thing for Stockman-Sharpe to XYZ and back.
%
% 4/16/15  dhb  Added a few more comparisons.  Not really needed.  These
%               were part of trying to track down an issue that in the end was related to
%               sRGB subtleties.  See v_Colorimetry for a description of
%               those.
    varargout = UnitTest.runValidationRun(@ValidationFunction, nargout, varargin);
end

%% Actual validation code
function ValidationFunction(runTimeParams)

    %% Ihitialize
    ieInit;

    %% Get Stockman-Sharpe and XYZ functions
    wave = (400:5:700)';
    xyz = ieReadSpectra('XYZ',wave);
    stock = ieReadSpectra('stockman',wave);
    if (runTimeParams.generatePlots)
        vcNewGraphWin;
        plot(wave,stock)
        plot(wave,xyz)
    end
    UnitTest.validationData('wave', wave);
    UnitTest.validationData('xyz', xyz);
    UnitTest.validationData('stockman', stock);

    %% Make matrix that transforms XYZ to Stockman-Sharpe
    % Isetbio keeps sensitivities in columns and matrices
    % are applied to the right to transform.
    %
    % Note (sigh) that PTB keeps sensitivities in rows and
    % that matrices are appplied from the left, so if you
    % are going to use PTB routines you need to think about
    % this.
    T_xyz2ss = xyz \ stock;
    pred = xyz*T_xyz2ss;
    if (runTimeParams.generatePlots)
        plot(pred(:),stock(:), '.'); grid on; axis equal
        xlim([-1 2]); ylim([-1 2]);
        xlabel('Predicted Stockman-Sharpe Sensitivities');
        ylabel('Actual Stockman-Sharpe Sensitivities');
    end

    %% Make the one that goes the other way.
    T_ss2xyz = stock \ xyz;
    pred = stock*T_ss2xyz;
    if (runTimeParams.generatePlots)
        plot(pred(:),xyz(:),'.'); grid on; axis equal
        xlim([-1 2]); ylim([-1 2]);
        xlabel('Predicted XYZ Sensitivities');
        ylabel('Actual XYZ Sensitivities');
        title('ISETBIO actual XYZ and LMS->XYZ');
    end

    %% Get the same things out of colorTransformMatrix
    T1 = colorTransformMatrix('stockman 2 xyz');
    T2 = colorTransformMatrix('xyz 2 stockman');
    tolerance = 1e-3;
    quantityOfInterest = T1-T_ss2xyz;
    UnitTest.assertIsZero(quantityOfInterest,'Matrix T_ss2xyz returned by colorTransformMatrix',tolerance);
    
    quantityOfInterest = T2-T_xyz2ss;
    UnitTest.assertIsZero(quantityOfInterest,'Matrix T_xyz2ss returned by colorTransformMatrix',tolerance);
    UnitTest.validationData('T1', T1);
    UnitTest.validationData('T2', T2);

    %% Check for self inversion
    % This matrix should be close to but not exactly the identity matrix.  The
    % differnce comes in because Stockman-Sharpe and XYZ are not exact linear
    % transformatinos of each other.
    tolerance = 0.02;
    identityCheck = T1*T2;
    quantityOfInterest = identityCheck-eye(3);
    UnitTest.assertIsZero(quantityOfInterest,'Conversion self-inversion',tolerance);
    
    %% PTB can do these things too.
    %
    % Are they the same?   
    load T_xyz1931
    T_xyz = SplineCmf(S_xyz1931,T_xyz1931,wave);
    load T_cones_ss2
    T_cones = SplineCmf(S_cones_ss2,T_cones_ss2,wave);
    UnitTest.validationData('T_cones', T_cones);
    UnitTest.validationData('T_xyz', T_xyz);

    % Data agreement is good for XYZ.  For S-S it is good but not perfect.
    % To see this, zoom in on the peak of the cone plots.  
    %
    % Might be difference in underlying wavelength sampling and how they
    % were interpolated by someone sometime.  The PTB functions match the
    % CVRL data, as demonstrated below.
    %
    % The agreement would get a little better for the subsampled PTB functions if one
    % renormalized after interpolating.
    if (runTimeParams.generatePlots)
        figure; clf;
        subplot(1,2,1); hold on
        plot(wave,xyz,'k','LineWidth',4);
        subplot(1,2,2); hold on
        plot(wave,stock,'k','LineWidth',4);
        subplot(1,2,1);
        plot(wave,T_xyz','r','LineWidth',3);
        xlabel('Wavelength (nm');
        ylabel('Sensitivity');
        title('ISETBIO and PTB XYZ');
        subplot(1,2,2);
        plot(wave,T_cones','r','LineWidth',3);
        plot(SToWls(S_cones_ss2),T_cones_ss2','g','LineWidth',2);
        xlabel('Wavelength (nm');
        ylabel('Sensitivity');
        title('ISETBIO and PTB LMS');    
    end
    
    % PTB data matches CVRL
    [cvrlWls,cvrlLMS] = ReturnCVRLLMS;
    if (runTimeParams.generatePlots)      
        figure; clf; hold on
        plot(cvrlWls,cvrlLMS,'k','LineWidth',3);
        plot(SToWls(S_cones_ss2),T_cones_ss2','g');
        xlabel('Wavelength (nm');
        ylabel('Sensitivity');
        title('CVRL direct and PTB LMS');
    end
    UnitTest.validationData('cvrlWls',cvrlWls);
    UnitTest.validationData('T_cones_ss2',T_cones_ss2);
    
    % The PTB matrix is a little different numerically.
    % The plot below shows that applying each (isetbio/ptb)
    % matrix to its own cone sensitivities yields very similar
    % estimates of XYZ.
    T_ptbData = (T_cones'\T_xyz')';
    ptbPred = T_ptbData*T_cones;
    if (runTimeParams.generatePlots)
        figure; clf; hold on
        plot(wave,pred,'k','LineWidth',3)
        plot(wave,ptbPred,'r','LineWidth',2);
        xlabel('Wavelength (nm');
        ylabel('Sensitivity');
        title('ISETBIO and PTB LMS->XYZ');
    end
    UnitTest.validationData('T_cones_ss2',T_cones_ss2);
    UnitTest.validationData('T_ptbData',T_ptbData);
    UnitTest.validationData('ptbPred',ptbPred);

    %% End
end

function [cvrlWls,cvrlLMS] = ReturnCVRLLMS
% I downloaded the SS 2-deg fundamentals from CVRL, April 2015,
% and pasted them in here (linear energy units, 1 nm sampling).
% These are pasted here to two places, but there are more digits
% available on the CVRL site, and (as of Summer 2016) in the PTB
% data.
%
% This function just returns them.  
cvrlWls = [390
391
392
393
394
395
396
397
398
399
400
401
402
403
404
405
406
407
408
409
410
411
412
413
414
415
416
417
418
419
420
421
422
423
424
425
426
427
428
429
430
431
432
433
434
435
436
437
438
439
440
441
442
443
444
445
446
447
448
449
450
451
452
453
454
455
456
457
458
459
460
461
462
463
464
465
466
467
468
469
470
471
472
473
474
475
476
477
478
479
480
481
482
483
484
485
486
487
488
489
490
491
492
493
494
495
496
497
498
499
500
501
502
503
504
505
506
507
508
509
510
511
512
513
514
515
516
517
518
519
520
521
522
523
524
525
526
527
528
529
530
531
532
533
534
535
536
537
538
539
540
541
542
543
544
545
546
547
548
549
550
551
552
553
554
555
556
557
558
559
560
561
562
563
564
565
566
567
568
569
570
571
572
573
574
575
576
577
578
579
580
581
582
583
584
585
586
587
588
589
590
591
592
593
594
595
596
597
598
599
600
601
602
603
604
605
606
607
608
609
610
611
612
613
614
615
616
617
618
619
620
621
622
623
624
625
626
627
628
629
630
631
632
633
634
635
636
637
638
639
640
641
642
643
644
645
646
647
648
649
650
651
652
653
654
655
656
657
658
659
660
661
662
663
664
665
666
667
668
669
670
671
672
673
674
675
676
677
678
679
680
681
682
683
684
685
686
687
688
689
690
691
692
693
694
695
696
697
698
699
700
701
702
703
704
705
706
707
708
709
710
711
712
713
714
715
716
717
718
719
720
721
722
723
724
725
726
727
728
729
730
731
732
733
734
735
736
737
738
739
740
741
742
743
744
745
746
747
748
749
750
751
752
753
754
755
756
757
758
759
760
761
762
763
764
765
766
767
768
769
770
771
772
773
774
775
776
777
778
779
780
781
782
783
784
785
786
787
788
789
790
791
792
793
794
795
796
797
798
799
800
801
802
803
804
805
806
807
808
809
810
811
812
813
814
815
816
817
818
819
820
821
822
823
824
825
826
827
828
829
830];

cvrlLMS = [4.15E-04	3.68E-04	9.55E-03
5.03E-04	4.48E-04	1.15E-02
6.07E-04	5.44E-04	1.38E-02
7.32E-04	6.59E-04	1.66E-02
8.79E-04	7.96E-04	1.99E-02
1.05E-03	9.59E-04	2.38E-02
1.25E-03	1.15E-03	2.85E-02
1.49E-03	1.37E-03	3.40E-02
1.76E-03	1.63E-03	4.04E-02
2.06E-03	1.93E-03	4.79E-02
2.41E-03	2.27E-03	5.66E-02
2.80E-03	2.65E-03	6.67E-02
3.23E-03	3.08E-03	7.81E-02
3.71E-03	3.56E-03	9.12E-02
4.24E-03	4.10E-03	1.06E-01
4.83E-03	4.70E-03	1.22E-01
5.49E-03	5.37E-03	1.41E-01
6.22E-03	6.12E-03	1.61E-01
7.01E-03	6.94E-03	1.83E-01
7.85E-03	7.83E-03	2.07E-01
8.72E-03	8.79E-03	2.33E-01
9.62E-03	9.82E-03	2.60E-01
1.05E-02	1.09E-02	2.89E-01
1.15E-02	1.21E-02	3.19E-01
1.24E-02	1.33E-02	3.49E-01
1.34E-02	1.45E-02	3.81E-01
1.44E-02	1.59E-02	4.14E-01
1.54E-02	1.72E-02	4.47E-01
1.64E-02	1.87E-02	4.80E-01
1.75E-02	2.02E-02	5.13E-01
1.84E-02	2.17E-02	5.44E-01
1.94E-02	2.32E-02	5.72E-01
2.03E-02	2.47E-02	5.99E-01
2.12E-02	2.63E-02	6.25E-01
2.20E-02	2.79E-02	6.50E-01
2.29E-02	2.96E-02	6.74E-01
2.39E-02	3.13E-02	7.00E-01
2.49E-02	3.32E-02	7.26E-01
2.60E-02	3.52E-02	7.53E-01
2.71E-02	3.73E-02	7.78E-01
2.82E-02	3.95E-02	8.03E-01
2.93E-02	4.18E-02	8.25E-01
3.05E-02	4.42E-02	8.45E-01
3.17E-02	4.66E-02	8.65E-01
3.29E-02	4.92E-02	8.84E-01
3.41E-02	5.18E-02	9.04E-01
3.54E-02	5.45E-02	9.24E-01
3.66E-02	5.71E-02	9.44E-01
3.79E-02	5.97E-02	9.63E-01
3.91E-02	6.23E-02	9.79E-01
4.03E-02	6.48E-02	9.91E-01
4.13E-02	6.71E-02	9.98E-01
4.23E-02	6.94E-02	1.00E+00
4.32E-02	7.16E-02	9.99E-01
4.40E-02	7.37E-02	9.96E-01
4.49E-02	7.59E-02	9.92E-01
4.59E-02	7.81E-02	9.87E-01
4.68E-02	8.04E-02	9.83E-01
4.78E-02	8.26E-02	9.76E-01
4.89E-02	8.49E-02	9.68E-01
4.99E-02	8.71E-02	9.55E-01
5.09E-02	8.92E-02	9.39E-01
5.19E-02	9.13E-02	9.21E-01
5.29E-02	9.34E-02	9.01E-01
5.41E-02	9.57E-02	8.80E-01
5.53E-02	9.82E-02	8.60E-01
5.68E-02	1.01E-01	8.42E-01
5.84E-02	1.04E-01	8.26E-01
6.02E-02	1.08E-01	8.11E-01
6.23E-02	1.12E-01	7.98E-01
6.47E-02	1.16E-01	7.87E-01
6.74E-02	1.21E-01	7.77E-01
7.04E-02	1.27E-01	7.68E-01
7.36E-02	1.32E-01	7.60E-01
7.71E-02	1.38E-01	7.50E-01
8.07E-02	1.45E-01	7.38E-01
8.44E-02	1.51E-01	7.24E-01
8.81E-02	1.57E-01	7.08E-01
9.19E-02	1.63E-01	6.90E-01
9.57E-02	1.70E-01	6.69E-01
9.95E-02	1.76E-01	6.46E-01
1.03E-01	1.82E-01	6.22E-01
1.07E-01	1.88E-01	5.97E-01
1.11E-01	1.94E-01	5.70E-01
1.15E-01	2.00E-01	5.43E-01
1.19E-01	2.05E-01	5.16E-01
1.23E-01	2.11E-01	4.90E-01
1.27E-01	2.17E-01	4.63E-01
1.31E-01	2.23E-01	4.38E-01
1.36E-01	2.30E-01	4.14E-01
1.40E-01	2.36E-01	3.90E-01
1.45E-01	2.42E-01	3.68E-01
1.49E-01	2.49E-01	3.48E-01
1.54E-01	2.55E-01	3.28E-01
1.59E-01	2.62E-01	3.09E-01
1.64E-01	2.68E-01	2.90E-01
1.69E-01	2.75E-01	2.73E-01
1.74E-01	2.81E-01	2.56E-01
1.79E-01	2.88E-01	2.40E-01
1.85E-01	2.96E-01	2.25E-01
1.92E-01	3.04E-01	2.12E-01
1.99E-01	3.13E-01	2.00E-01
2.06E-01	3.23E-01	1.89E-01
2.14E-01	3.33E-01	1.79E-01
2.23E-01	3.45E-01	1.69E-01
2.33E-01	3.57E-01	1.61E-01
2.43E-01	3.70E-01	1.52E-01
2.54E-01	3.83E-01	1.45E-01
2.65E-01	3.97E-01	1.37E-01
2.77E-01	4.12E-01	1.30E-01
2.89E-01	4.28E-01	1.23E-01
3.02E-01	4.44E-01	1.16E-01
3.16E-01	4.61E-01	1.09E-01
3.30E-01	4.79E-01	1.02E-01
3.44E-01	4.97E-01	9.54E-02
3.60E-01	5.16E-01	8.89E-02
3.76E-01	5.35E-01	8.25E-02
3.92E-01	5.55E-01	7.65E-02
4.09E-01	5.75E-01	7.08E-02
4.26E-01	5.95E-01	6.55E-02
4.44E-01	6.16E-01	6.08E-02
4.62E-01	6.36E-01	5.66E-02
4.80E-01	6.57E-01	5.28E-02
4.99E-01	6.78E-01	4.92E-02
5.18E-01	6.98E-01	4.59E-02
5.36E-01	7.19E-01	4.28E-02
5.55E-01	7.40E-01	3.98E-02
5.74E-01	7.60E-01	3.69E-02
5.93E-01	7.80E-01	3.42E-02
6.11E-01	7.99E-01	3.16E-02
6.29E-01	8.17E-01	2.92E-02
6.45E-01	8.33E-01	2.70E-02
6.61E-01	8.48E-01	2.49E-02
6.76E-01	8.61E-01	2.29E-02
6.91E-01	8.74E-01	2.11E-02
7.05E-01	8.86E-01	1.94E-02
7.19E-01	8.97E-01	1.78E-02
7.32E-01	9.08E-01	1.64E-02
7.46E-01	9.18E-01	1.50E-02
7.58E-01	9.27E-01	1.37E-02
7.71E-01	9.36E-01	1.26E-02
7.82E-01	9.43E-01	1.15E-02
7.93E-01	9.50E-01	1.06E-02
8.04E-01	9.57E-01	9.68E-03
8.15E-01	9.63E-01	8.86E-03
8.26E-01	9.69E-01	8.09E-03
8.37E-01	9.75E-01	7.39E-03
8.48E-01	9.81E-01	6.74E-03
8.60E-01	9.86E-01	6.14E-03
8.71E-01	9.91E-01	5.59E-03
8.81E-01	9.95E-01	5.09E-03
8.90E-01	9.98E-01	4.63E-03
8.99E-01	1.00E+00	4.21E-03
9.07E-01	1.00E+00	3.83E-03
9.13E-01	9.99E-01	3.49E-03
9.19E-01	9.97E-01	3.17E-03
9.24E-01	9.94E-01	2.88E-03
9.28E-01	9.90E-01	2.62E-03
9.32E-01	9.86E-01	2.38E-03
9.36E-01	9.81E-01	2.16E-03
9.40E-01	9.77E-01	1.96E-03
9.45E-01	9.73E-01	1.78E-03
9.50E-01	9.70E-01	1.61E-03
9.56E-01	9.66E-01	1.46E-03
9.61E-01	9.62E-01	1.33E-03
9.66E-01	9.57E-01	1.20E-03
9.70E-01	9.50E-01	1.09E-03
9.73E-01	9.43E-01	9.90E-04
9.76E-01	9.35E-01	8.99E-04
9.79E-01	9.26E-01	8.16E-04
9.81E-01	9.18E-01	7.40E-04
9.84E-01	9.09E-01	6.72E-04
9.87E-01	9.01E-01	6.10E-04
9.90E-01	8.92E-01	5.53E-04
9.92E-01	8.83E-01	5.02E-04
9.94E-01	8.73E-01	4.56E-04
9.96E-01	8.63E-01	4.14E-04
9.98E-01	8.51E-01	3.76E-04
9.99E-01	8.39E-01	3.41E-04
1.00E+00	8.27E-01	3.10E-04
1.00E+00	8.14E-01	2.82E-04
1.00E+00	8.00E-01	2.56E-04
9.99E-01	7.86E-01	2.33E-04
9.98E-01	7.72E-01	2.12E-04
9.95E-01	7.56E-01	1.92E-04
9.92E-01	7.40E-01	1.75E-04
9.88E-01	7.23E-01	1.59E-04
9.83E-01	7.06E-01	1.45E-04
9.78E-01	6.88E-01	1.32E-04
9.74E-01	6.71E-01	1.20E-04
9.69E-01	6.53E-01	1.09E-04
9.66E-01	6.37E-01	9.97E-05
9.64E-01	6.20E-01	9.09E-05
9.61E-01	6.04E-01	8.29E-05
9.59E-01	5.88E-01	7.56E-05
9.56E-01	5.73E-01	6.90E-05
9.52E-01	5.57E-01	6.30E-05
9.47E-01	5.41E-01	5.75E-05
9.41E-01	5.25E-01	5.25E-05
9.35E-01	5.09E-01	4.80E-05
9.28E-01	4.93E-01	4.39E-05
9.20E-01	4.76E-01	4.02E-05
9.12E-01	4.60E-01	3.67E-05
9.04E-01	4.44E-01	3.36E-05
8.95E-01	4.27E-01	3.08E-05
8.86E-01	4.11E-01	2.82E-05
8.76E-01	3.95E-01	2.59E-05
8.66E-01	3.80E-01	2.37E-05
8.56E-01	3.64E-01	2.18E-05
8.45E-01	3.49E-01	2.00E-05
8.34E-01	3.34E-01	1.83E-05
8.23E-01	3.20E-01	1.69E-05
8.11E-01	3.06E-01	1.55E-05
8.00E-01	2.92E-01	1.42E-05
7.88E-01	2.78E-01	1.31E-05
7.75E-01	2.65E-01	1.21E-05
7.62E-01	2.52E-01	1.11E-05
7.48E-01	2.40E-01	1.02E-05
7.34E-01	2.28E-01	9.44E-06
7.20E-01	2.16E-01	8.71E-06
7.06E-01	2.05E-01	8.03E-06
6.91E-01	1.95E-01	7.42E-06
6.76E-01	1.84E-01	6.85E-06
6.61E-01	1.75E-01	6.33E-06
6.46E-01	1.65E-01	5.86E-06
6.31E-01	1.56E-01	5.42E-06
6.15E-01	1.48E-01	0.00E+00
6.00E-01	1.39E-01	0.00E+00
5.84E-01	1.31E-01	0.00E+00
5.69E-01	1.24E-01	0.00E+00
5.54E-01	1.17E-01	0.00E+00
5.39E-01	1.10E-01	0.00E+00
5.25E-01	1.03E-01	0.00E+00
5.10E-01	9.70E-02	0.00E+00
4.95E-01	9.11E-02	0.00E+00
4.80E-01	8.56E-02	0.00E+00
4.64E-01	8.03E-02	0.00E+00
4.48E-01	7.54E-02	0.00E+00
4.32E-01	7.07E-02	0.00E+00
4.16E-01	6.63E-02	0.00E+00
4.01E-01	6.21E-02	0.00E+00
3.85E-01	5.82E-02	0.00E+00
3.70E-01	5.44E-02	0.00E+00
3.56E-01	5.09E-02	0.00E+00
3.42E-01	4.76E-02	0.00E+00
3.28E-01	4.45E-02	0.00E+00
3.15E-01	4.16E-02	0.00E+00
3.02E-01	3.88E-02	0.00E+00
2.89E-01	3.62E-02	0.00E+00
2.77E-01	3.38E-02	0.00E+00
2.66E-01	3.14E-02	0.00E+00
2.55E-01	2.92E-02	0.00E+00
2.44E-01	2.71E-02	0.00E+00
2.34E-01	2.52E-02	0.00E+00
2.23E-01	2.34E-02	0.00E+00
2.13E-01	2.18E-02	0.00E+00
2.03E-01	2.03E-02	0.00E+00
1.93E-01	1.90E-02	0.00E+00
1.84E-01	1.77E-02	0.00E+00
1.74E-01	1.66E-02	0.00E+00
1.65E-01	1.54E-02	0.00E+00
1.56E-01	1.44E-02	0.00E+00
1.48E-01	1.34E-02	0.00E+00
1.40E-01	1.24E-02	0.00E+00
1.32E-01	1.15E-02	0.00E+00
1.25E-01	1.07E-02	0.00E+00
1.18E-01	9.93E-03	0.00E+00
1.11E-01	9.20E-03	0.00E+00
1.05E-01	8.52E-03	0.00E+00
9.87E-02	7.89E-03	0.00E+00
9.30E-02	7.30E-03	0.00E+00
8.76E-02	6.76E-03	0.00E+00
8.24E-02	6.25E-03	0.00E+00
7.75E-02	5.79E-03	0.00E+00
7.29E-02	5.36E-03	0.00E+00
6.85E-02	4.97E-03	0.00E+00
6.44E-02	4.61E-03	0.00E+00
6.04E-02	4.28E-03	0.00E+00
5.67E-02	3.98E-03	0.00E+00
5.32E-02	3.70E-03	0.00E+00
4.99E-02	3.44E-03	0.00E+00
4.67E-02	3.19E-03	0.00E+00
4.38E-02	2.97E-03	0.00E+00
4.10E-02	2.76E-03	0.00E+00
3.83E-02	2.56E-03	0.00E+00
3.58E-02	2.38E-03	0.00E+00
3.35E-02	2.21E-03	0.00E+00
3.13E-02	2.05E-03	0.00E+00
2.92E-02	1.90E-03	0.00E+00
2.72E-02	1.76E-03	0.00E+00
2.54E-02	1.64E-03	0.00E+00
2.37E-02	1.52E-03	0.00E+00
2.20E-02	1.41E-03	0.00E+00
2.05E-02	1.31E-03	0.00E+00
1.91E-02	1.21E-03	0.00E+00
1.77E-02	1.12E-03	0.00E+00
1.64E-02	1.04E-03	0.00E+00
1.53E-02	9.60E-04	0.00E+00
1.41E-02	8.88E-04	0.00E+00
1.31E-02	8.22E-04	0.00E+00
1.22E-02	7.61E-04	0.00E+00
1.13E-02	7.06E-04	0.00E+00
1.05E-02	6.55E-04	0.00E+00
9.78E-03	6.08E-04	0.00E+00
9.10E-03	5.65E-04	0.00E+00
8.47E-03	5.25E-04	0.00E+00
7.88E-03	4.89E-04	0.00E+00
7.33E-03	4.54E-04	0.00E+00
6.82E-03	4.22E-04	0.00E+00
6.34E-03	3.93E-04	0.00E+00
5.90E-03	3.65E-04	0.00E+00
5.48E-03	3.40E-04	0.00E+00
5.10E-03	3.16E-04	0.00E+00
4.74E-03	2.94E-04	0.00E+00
4.41E-03	2.73E-04	0.00E+00
4.09E-03	2.53E-04	0.00E+00
3.80E-03	2.35E-04	0.00E+00
3.52E-03	2.18E-04	0.00E+00
3.27E-03	2.03E-04	0.00E+00
3.03E-03	1.88E-04	0.00E+00
2.80E-03	1.74E-04	0.00E+00
2.60E-03	1.62E-04	0.00E+00
2.41E-03	1.50E-04	0.00E+00
2.23E-03	1.40E-04	0.00E+00
2.07E-03	1.30E-04	0.00E+00
1.92E-03	1.21E-04	0.00E+00
1.78E-03	1.12E-04	0.00E+00
1.66E-03	1.04E-04	0.00E+00
1.54E-03	9.71E-05	0.00E+00
1.43E-03	9.04E-05	0.00E+00
1.33E-03	8.42E-05	0.00E+00
1.23E-03	7.84E-05	0.00E+00
1.14E-03	7.29E-05	0.00E+00
1.06E-03	6.79E-05	0.00E+00
9.88E-04	6.33E-05	0.00E+00
9.18E-04	5.89E-05	0.00E+00
8.53E-04	5.49E-05	0.00E+00
7.94E-04	5.12E-05	0.00E+00
7.38E-04	4.78E-05	0.00E+00
6.87E-04	4.46E-05	0.00E+00
6.39E-04	4.16E-05	0.00E+00
5.95E-04	3.88E-05	0.00E+00
5.54E-04	3.62E-05	0.00E+00
5.15E-04	3.38E-05	0.00E+00
4.79E-04	3.15E-05	0.00E+00
4.46E-04	2.94E-05	0.00E+00
4.15E-04	2.75E-05	0.00E+00
3.86E-04	2.56E-05	0.00E+00
3.59E-04	2.39E-05	0.00E+00
3.34E-04	2.23E-05	0.00E+00
3.11E-04	2.09E-05	0.00E+00
2.90E-04	1.95E-05	0.00E+00
2.70E-04	1.83E-05	0.00E+00
2.52E-04	1.71E-05	0.00E+00
2.35E-04	1.61E-05	0.00E+00
2.19E-04	1.50E-05	0.00E+00
2.05E-04	1.41E-05	0.00E+00
1.91E-04	1.32E-05	0.00E+00
1.78E-04	1.23E-05	0.00E+00
1.66E-04	1.16E-05	0.00E+00
1.55E-04	1.08E-05	0.00E+00
1.44E-04	1.01E-05	0.00E+00
1.35E-04	9.49E-06	0.00E+00
1.26E-04	8.90E-06	0.00E+00
1.17E-04	8.34E-06	0.00E+00
1.10E-04	7.82E-06	0.00E+00
1.02E-04	7.34E-06	0.00E+00
9.56E-05	6.89E-06	0.00E+00
8.93E-05	6.46E-06	0.00E+00
8.35E-05	6.06E-06	0.00E+00
7.80E-05	5.69E-06	0.00E+00
7.29E-05	5.34E-06	0.00E+00
6.81E-05	5.01E-06	0.00E+00
6.36E-05	4.70E-06	0.00E+00
5.95E-05	4.41E-06	0.00E+00
5.56E-05	4.14E-06	0.00E+00
5.20E-05	3.89E-06	0.00E+00
4.87E-05	3.65E-06	0.00E+00
4.56E-05	3.43E-06	0.00E+00
4.27E-05	3.22E-06	0.00E+00
3.99E-05	3.03E-06	0.00E+00
3.74E-05	2.84E-06	0.00E+00
3.49E-05	2.67E-06	0.00E+00
3.27E-05	2.50E-06	0.00E+00
3.06E-05	2.35E-06	0.00E+00
2.86E-05	2.21E-06	0.00E+00
2.68E-05	2.08E-06	0.00E+00
2.51E-05	1.96E-06	0.00E+00
2.36E-05	1.84E-06	0.00E+00
2.21E-05	1.74E-06	0.00E+00
2.07E-05	1.63E-06	0.00E+00
1.94E-05	1.54E-06	0.00E+00
1.82E-05	1.45E-06	0.00E+00
1.71E-05	1.36E-06	0.00E+00
1.60E-05	1.29E-06	0.00E+00
1.50E-05	1.21E-06	0.00E+00
1.41E-05	1.14E-06	0.00E+00
1.32E-05	1.07E-06	0.00E+00
1.24E-05	1.01E-06	0.00E+00
1.17E-05	9.54E-07	0.00E+00
1.09E-05	8.99E-07	0.00E+00
1.03E-05	8.47E-07	0.00E+00
9.64E-06	7.99E-07	0.00E+00
9.05E-06	7.53E-07	0.00E+00
8.50E-06	7.10E-07	0.00E+00
7.98E-06	6.70E-07	0.00E+00
7.49E-06	6.32E-07	0.00E+00
7.04E-06	5.97E-07	0.00E+00
6.62E-06	5.63E-07	0.00E+00
6.22E-06	5.32E-07	0.00E+00
5.85E-06	5.03E-07	0.00E+00
5.50E-06	4.76E-07	0.00E+00
5.18E-06	4.50E-07	0.00E+00
4.87E-06	4.25E-07	0.00E+00
4.58E-06	4.02E-07	0.00E+00
4.31E-06	3.80E-07	0.00E+00
4.05E-06	3.59E-07	0.00E+00
3.81E-06	3.39E-07	0.00E+00
3.58E-06	3.21E-07	0.00E+00
3.37E-06	3.03E-07	0.00E+00
3.17E-06	2.86E-07	0.00E+00
2.98E-06	2.71E-07	0.00E+00
2.81E-06	2.56E-07	0.00E+00
2.64E-06	2.42E-07	0.00E+00
2.49E-06	2.29E-07	0.00E+00
2.34E-06	2.17E-07	0.00E+00
2.21E-06	2.05E-07	0.00E+00
2.08E-06	1.94E-07	0.00E+00
1.96E-06	1.84E-07	0.00E+00
1.85E-06	1.74E-07	0.00E+00
1.75E-06	1.65E-07	0.00E+00
1.65E-06	1.56E-07	0.00E+00
1.55E-06	1.48E-07	0.00E+00
1.46E-06	1.40E-07	0.00E+00
1.38E-06	1.33E-07	0.00E+00
1.30E-06	1.26E-07	0.00E+00
1.23E-06	1.19E-07	0.00E+00
1.16E-06	1.12E-07	0.00E+00
1.09E-06	1.06E-07	0.00E+00
1.03E-06	1.01E-07	0.00E+00
9.74E-07	9.53E-08	0.00E+00];
end

