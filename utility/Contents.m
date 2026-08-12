% UTILITY
%
% Top Level Files
%   isetbioDataPath                           - Return path to data files bundled with ISETBio
%   isetbioRootPath                           - Return the path to the root isetbio directory
%   temporalEquivEcc                          - (IN PROG.) Equivalent temporal eccentricity acc. to Kalmar/Chichilnisky
%   twoGammaResp                              - Create temporal response composed of the difference of 2 gamma functions
%   determineUnits                            - Function determines if unit other than default specified
%   localDropboxDir                           -
%   selectPlotSupport                         - Used with getMiddleMatrix to pull out the 'interesting' center of a plot
%
% FUNCTIONS (subdir functions)          Computational utility functions for isetbio
%   functions/coneEmpiricalDimFlash           - Return estimate of Schnapf et al. dim flash cone photocurrent response.
%   functions/hill                            - Compute hill function
%   functions/weberFechner                    - Calculate the Weber Fechner tvi function
%   functions/sumOfTwoDblExponentials         - Compute sum of two double exponential functions, with no mean offset.
%
% PSYCHOPHYSICS (subdir psychophysics)  Related to calculating psychophysical performance. The TSD
%   psychophysics/analyticpHitpFa             - [pHit,pFa] = analyticpHitpFa(signalMean,noiseMean,commonSd,rightCrit)
%   psychophysics/computeROCArea              - rocArea = computeROCArea(signalMean,noiseMean,commonSd,criteria)
%   psychophysics/computeDPrimeCritNorm       - [dprime,critNorm] = computeDPrimeCritNorm(pHit,pFa)
%   psychophysics/PoissonDecisionLogLikelihoood - Log-likelihood for Poisson response vectors
%   psychophysics/PoissonIdealObserverNAlternativeFC -
%   psychophysics/SupportVectorMachineObserverNAlternativeFC - SVM-based Monte-Carlo simulation of probability correct for N-alternative forced choice
%   psychophysics/TAFCFractionCorrectToDPrime - dPrime = TAFCFractionCorrectToDPrime(fractionCorrect)
%   psychophysics/analyticPoissonIdealObserver - ]  Get the ideal observer fraction correct in a TAFC task
%   psychophysics/dPrimeToTAFCFractionCorrect - fractionCorrect = dPrimeToTAFCPercentFraction(dPrime)
%
% APPLESILICONCPU (subdir AppleSiliconCPU)
%   AppleSiliconCPU/AppleSiliconParPoolManager - Create an AppleSiliconParPoolManager
%
% CONES (subdir cones)
%   cones/IsomerizationsFromLuminanceGeisler  - Estimate cone isomerizations from luminance level
%   cones/coneMeanIsomerizations              - Calculate spatial mean photon rate(R*/sec) for the 3 cone types in mosaic
%   cones/coneTypeLocations                   - Return row/col or indices of the three types of cones
%   cones/diameterForCircularApertureFromWidthForSquareAperture -
%   cones/sizeForSquareApertureFromDiameterForCircularAperture -
%
% MOSAIC (subdir mosaic)
%   mosaic/mosaicLoad                         - Load a stored cMosaic from the ISETBio SDR mosaic cache
%   mosaic/mosaicName                         - Return a standard mosaic name for data/cones file
%
% SDR (subdir sdr)
%   sdr/sdrLegacyFilenameMap                  - Legacy ISETBio filenames for the isetbio-mosaics SDR assets
%   sdr/sdrMosaicCacheRoot                    - Root of the local cache mirroring the isetbio-mosaics SDR deposit
%   sdr/sdrMosaicFetch                        - Return a verified local copy of one isetbio-mosaics asset
%   sdr/sdrMosaicRecords                      - Return the manifest records for one isetbio-mosaics collection
%   sdr/sdrMosaicVerify                       - True when a local file matches its isetbio-mosaics manifest record
%   sdr/sdrPrebakedCircuitFile                - Resolve a pre-baked ON-mRGC circuit filename to a local file
%
