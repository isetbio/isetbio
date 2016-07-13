function [testmovie, xval_mosaic_out] = loadDataRGCFigure2(experimentI, stimulusTestI, cellTypeI)
% Big switch statement to load data for t_rgcNaturalScenesFigure2. 
% 
% Inputs: 
%       the date of the experiment, 
%       the cell type (ON Parasol or OFF Parasol) and 
%       the stimulus (WN or NSEM).
% 
% Outputs:
%       the stimulus movie
%       recorded RGC spikes from the experiment
% 
% % from t_rgcNaturalScenesFigure2:
% switch experimentI
%     case 1; experimentID = '2013-08-19-6';
%     case 2; experimentID = '2012-08-09-3';
%     case 3; experimentID = '2013-10-10-0';
%     case 4; experimentID = '2012-09-27-3';
% end
% 
% switch cellTypeI
%     case 1; cellTypeI = 'On Parasol';
%     case 2; cellTypeI = 'Off Parasol';
% end
% 
% switch stimulusTestI
%     case 1; stimulusTest = 'WN';
%     case 2; stimulusTest = 'NSEM';
% end
% 
% 5/2016 JRG (c) isetbio team

% RDT initialization
rdt = RdtClient('isetbio');
rdt.crp('resources/data/rgc');

switch experimentI
    case 1 % 2013-08-19-6      
        
        switch stimulusTestI
            case 1 % WN
                
                % Binary white noise test movie
%                 data = rdt.readArtifact('testmovie_8pix_Identity_8pix', 'type', 'mat');                
%                 testmovie = data.testmovie;
                
                load(['/Users/james/Documents/MATLAB/'...
                    'akheitman/WN_mapPRJ/Stimuli/'...
                    'BW-8-1-0.48-11111_RNG_16807/testmovie_8pix_Identity_8pix.mat']);
                
                switch cellTypeI
                    case 1 % ON Parasol
                        load('isetbio misc/RDT Uploads/xval_mosaic_WN_ONParasol_2013_08_19_6.mat');
%                         data = rdt.readArtifact('xval_mosaic_WN_ONParasol_2013_08_19_6', 'type', 'mat');                
%                         xval_mosaic = data.xval_mosaic;
%                                                 
%                         data2 = rdt.readArtifact('goodind_2013_08_19_6_ONParasol', 'type', 'mat');
%                         goodind = data2.goodind;
                    case 2 % OFF Paraasol
%                         load('isetbio misc/RDT Uploads/xval_mosaic_WN_OFFParasol_2013_08_19_6.mat');
                        data = rdt.readArtifact('xval_mosaic_WN_OFFParasol_2013_08_19_6', 'type', 'mat');                
                        xval_mosaic = data.xval_mosaic;
                        
%                         load('/Users/james/Documents/MATLAB/isetbio misc/RDT uploads/goodind_2013_08_19_6_OFFParasol.mat')
                        data2 = rdt.readArtifact('goodind_2013_08_19_6_OFFParasol', 'type', 'mat');
                        goodind = data2.goodind;
                end
                
            case 2
                % NSEM test movie
                % These are small black and white van hatteren images with eye
                % movements superimposed.
%                 load(['/Users/james/Documents/MATLAB/'...
%                     'akheitman/NSEM_mapPRJ/Stimuli/'...
%                     'NSEM_eye-long-v2/testmovie_schemeA_8pix_Identity_8pix.mat']);
                
                data = rdt.readArtifact('testmovie_schemeA_8pix_Identity_8pix', 'type', 'mat');                
                testmovie = data.testmovie;
                
                switch cellTypeI
                    case 1 % ON Parasol
%                         load('isetbio misc/RDT Uploads/xval_mosaic_NSEM_ONParasol_2013_08_19_6.mat');
                        data = rdt.readArtifact('xval_mosaic_NSEM_ONParasol_2013_08_19_6', 'type', 'mat');                
                        xval_mosaic = data.xval_mosaic;
                        
                        data2 = rdt.readArtifact('goodind_2013_08_19_6_ONParasol', 'type', 'mat');
                        goodind = data2.goodind;

                    case 2 % OFF Paraasol
%                         load('isetbio misc/RDT Uploads/xval_mosaic_NSEM_OFFParasol_2013_08_19_6.mat');
                        data = rdt.readArtifact('xval_mosaic_NSEM_OFFParasol_2013_08_19_6', 'type', 'mat');                
                        xval_mosaic = data.xval_mosaic;
                        
                        data2 = rdt.readArtifact('goodind_2013_08_19_6_OFFParasol', 'type', 'mat');
                        goodind = data2.goodind;

                end
        end
% % % % % % % %         
    case 2 % 2012-08-09-3      
        
        switch stimulusTestI
            case 1 % WN
                
                % Binary white noise test movie
                data = rdt.readArtifact('testmovie_8pix_Identity_8pix', 'type', 'mat');                
                testmovie = data.testmovie;
                
                switch cellTypeI
                    case 1 % ON Parasol
                        load('isetbio misc/RDT Uploads/xval_mosaic_WN_ONParasol_2013_08_19_6.mat');
                    case 2 % OFF Paraasol
                        load('isetbio misc/RDT Uploads/xval_mosaic_WN_OFFParasol_2013_08_19_6.mat');
                end
                
            case 2
                % NSEM test movie
                % These are small black and white van hatteren images with eye
                % movements superimposed.
                load(['/Users/james/Documents/MATLAB/'...
                    'akheitman/NSEM_mapPRJ/Stimuli/'...
                    'NSEM_eye-long-v2/testmovie_schemeA_8pix_Identity_8pix.mat']);
                
                % data = rdt.readArtifact('testmovie_8pix_Identity_8pix', 'type', 'mat');                
                % testmovie = data.testmovie;
                
                switch cellTypeI
                    case 1 % ON Parasol
                        load('isetbio misc/RDT Uploads/xval_mosaic_NSEM_ONParasol_2013_08_19_6.mat');
                    case 2 % OFF Paraasol
                        load('isetbio misc/RDT Uploads/xval_mosaic_NSEM_OFFParasol_2013_08_19_6.mat');
                end
        end

% % % % % % % % % 
    case 3 % 2013-10-10-0     
        
        switch stimulusTestI
            case 1 % WN
                
                % Binary white noise test movie
                data = rdt.readArtifact('testmovie_8pix_Identity_8pix', 'type', 'mat');                
                testmovie = data.testmovie;
                
                switch cellTypeI
                    case 1 % ON Parasol
                        load('isetbio misc/RDT Uploads/xval_mosaic_WN_ONParasol_2013_08_19_6.mat');
                    case 2 % OFF Paraasol
                        load('isetbio misc/RDT Uploads/xval_mosaic_WN_OFFParasol_2013_08_19_6.mat');
                end
                
            case 2
                % NSEM test movie
                % These are small black and white van hatteren images with eye
                % movements superimposed.
                load(['/Users/james/Documents/MATLAB/'...
                    'akheitman/NSEM_mapPRJ/Stimuli/'...
                    'NSEM_eye-long-v2/testmovie_schemeA_8pix_Identity_8pix.mat']);
                
                % data = rdt.readArtifact('testmovie_8pix_Identity_8pix', 'type', 'mat');                
                % testmovie = data.testmovie;
                
                switch cellTypeI
                    case 1 % ON Parasol
                        load('isetbio misc/RDT Uploads/xval_mosaic_NSEM_ONParasol_2013_08_19_6.mat');
                    case 2 % OFF Paraasol
                        load('isetbio misc/RDT Uploads/xval_mosaic_NSEM_OFFParasol_2013_08_19_6.mat');
                end
        end

% % % % % % % % % 
    case 4 % 2012-09-27-3    
        
        switch stimulusTestI
            case 1 % WN
                
                % Binary white noise test movie
                data = rdt.readArtifact('testmovie_8pix_Identity_8pix', 'type', 'mat');                
                testmovie = data.testmovie;
                
                switch cellTypeI
                    case 1 % ON Parasol
                        load('isetbio misc/RDT Uploads/xval_mosaic_WN_ONParasol_2013_08_19_6.mat');
                    case 2 % OFF Paraasol
                        load('isetbio misc/RDT Uploads/xval_mosaic_WN_OFFParasol_2013_08_19_6.mat');
                end
                
            case 2
                % NSEM test movie
                % These are small black and white van hatteren images with eye
                % movements superimposed.
                load(['/Users/james/Documents/MATLAB/'...
                    'akheitman/NSEM_mapPRJ/Stimuli/'...
                    'NSEM_eye-long-v2/testmovie_schemeA_8pix_Identity_8pix.mat']);
                
                % data = rdt.readArtifact('testmovie_8pix_Identity_8pix', 'type', 'mat');                
                % testmovie = data.testmovie;
                
                switch cellTypeI
                    case 1 % ON Parasol
                        load('isetbio misc/RDT Uploads/xval_mosaic_NSEM_ONParasol_2013_08_19_6.mat');
                    case 2 % OFF Paraasol
                        load('isetbio misc/RDT Uploads/xval_mosaic_NSEM_OFFParasol_2013_08_19_6.mat');
                end
        end

end

goodind = 1:length(xval_mosaic);
for mosaicInd = 1:length(goodind)
    xval_mosaic_out{mosaicInd} = xval_mosaic{goodind(mosaicInd)};
end

testmovieshort = testmovie.matrix(:,:,1:1200); 


