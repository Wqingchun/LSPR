function LSPR (inputFilename,outputFilename,inputPath,outputPath,defaultPeriod,lower,upper)
%Analyze microarray data (evenly or unevenly sampled).
%
%   Usage:
%       LSPR (inputFilename,outputFilename,inputPath,outputPath,defaultPeriod,lower,upper)
%
%   INPUT: 
%       inputFilename          - input file name
%       outputFilename         - output file name
%       inputPath              - load input file from
%       outputPath             - save output file to
%       defaultPeriod          - use a default period (i.e. 24 for circadian microarray data) to do harmonic analysis when no periods could be detected in [lower,upper]
%       lower/upper            - endpoints of period range
%
%   OUTPUT: 
%       none
%
%   Examples:
%       LSPR('inputExample.txt','outputExample.txt','inputFolder/','outputFolder/',24,20,28)

%//////////////////////////////////////////////////////////////////////////
%   LSPR calls 6 functions when analyzing: 
%   1.  LSPR                 - main function to perform LSPR algorithm.
%   2.  loadFile             - load data from input file
%   3.  dataPreprocessing    - detrend or filter data
%   4.  spectrumEstimation   - apply Lomb-Scargle periodogram to detect periods
%   5.  harmonicRegression   - do harmonic analysis and compute AIC of harmonic models 
%   6.  saveFile             - save results to output file 
%
%  The other one auxiliary function in spectrumEstimation.m:
%  7.  lombFunction         - Lomb-Scargle periodogram function
%//////////////////////////////////////////////////////////////////////////
%   Copyright (C) 2010 Chen ZHANG and Rendong Yang.
%   $Revision Date: 2010/12/7 $
%//////////////////////////////////////////////////////////////////////////
%//  Authors:
%//        name            organization 					email
%//    --------------  ------------------------    ------------------------------
%//    Chen ZHANG         College of Science            zcreation@yahoo.cn
%//    Rendong Yang   College of Biological Sciences     cauyrd@gmail.com
%//
%//  Established Date:   2010/9/2
%//////////////////////////////////////////////////////////////////////////

% Local CONSTANTS
FINAL_VAR_NUM = 10;

% LOAD FILE
fprintf('Loading...\n');
expandInputPath = strcat(inputPath,inputFilename);
[probeNames,timepoints,microarrayData,totalGenesNum,validSign,stat] = loadFile(expandInputPath);
if stat == -1% Fail to load data
    exit;
end



% Step 1 : DATA PREPROCESSING
[detrendedData,filteredData] = dataPreprocessing(microarrayData);



% Step 2 : SPECTRUM ESTIMATION
mRows = size(microarrayData,1);
% Create two pages for each intermediate variable (followings) to save values:
%  - first page saves the intermediate variables of detrendedData after spectrum estimation
%  - second page saves the intermediate variables of filtereddData (been detrended and filter) after spectrum estimation 
peakSignLayers = zeros(mRows,1,2);% Mark if any peaks exists in [1/upper,1/lower]
periodsNumLayers = zeros(mRows,1,2);% Array of number of different oscillations 
periodsLayers = cell(mRows,1,2); % Array of detected periods

% Do spectrum estimation with detrendedData and filteredData, respectively
fprintf('Estimating periods with detrended data...\n');
[peakSignLayers(:,1,1),periodsLayers(:,1,1),periodsNumLayers(:,1,1)] = spectrumEstimation(detrendedData,lower,upper,defaultPeriod,timepoints);
fprintf('Estimating periods with detrended and flitered data...\n');
[peakSignLayers(:,1,2),periodsLayers(:,1,2),periodsNumLayers(:,1,2)] = spectrumEstimation(filteredData,lower,upper,defaultPeriod,timepoints);



% Step 3: HARMONIC REGRESSION
% Create two pages for each intermediate variable (followings) to save values:
%  - first page saves the intermediate variables of detrendedData after harmonic analysis
%  - second page saves the intermediate variables of filtereddData after harmonic analysis
rsquaresLayers = zeros(mRows,1,2);% Array of R squares
pvaluesLayers = zeros(mRows,1,2);% Array of pvalues
AICsLayers = zeros(mRows,1,2);%  Array of  AICs
amplitudesLayers = cell(mRows,1,2);% Array of amplitudes
phasesLayers = cell(mRows,1,2);% Array of phases

% Harmonic analysis with detrendedData and filteredData, respectively
fprintf('Do harmonic regression with detrended data...\n');
[AICsLayers(:,1,1),amplitudesLayers(:,1,1),phasesLayers(:,1,1),rsquaresLayers(:,1,1),pvaluesLayers(:,1,1)] = harmonicRegression(detrendedData,timepoints,periodsLayers(:,1,1),periodsNumLayers(:,1,1));
fprintf('Do harmonic regression with detrended and filtered data...\n');
[AICsLayers(:,1,2),amplitudesLayers(:,1,2),phasesLayers(:,1,2),rsquaresLayers(:,1,2),pvaluesLayers(:,1,2)] = harmonicRegression(detrendedData,timepoints,periodsLayers(:,1,2),periodsNumLayers(:,1,2));

% Choose a better one from two harmonic models according to Akaike information criterion
% (AIC) and compute corresponding qvalues and FDR-BHs
Results = cell(mRows,FINAL_VAR_NUM);% Save output variables
%--------------------------------------------------------------------------
%      1                   2                   3                   4
% filter type           method      number of oscillations      period
%      5                   6                   7                   8
%  amplitude            phase              R-square             pvalue
%      9                  10
%    qvalue             FDR-BH
%--------------------------------------------------------------------------
for iRow=1:mRows
    [~,page] = min(AICsLayers(iRow,1,:));% Find the better model
    if(page == 1)
        Results{iRow,1} = -1; % Page of detrended data
    else
        Results{iRow,1} = 1;% Page of detrended and filtered data
    end
    if peakSignLayers(iRow,1,page) == 1
        Results{iRow,2}='LSPR';
    else
        Results{iRow,2}='Default';
    end
    Results{iRow,3} = periodsNumLayers(iRow,1,page);  % Number of oscillations
    Results{iRow,4} = periodsLayers{iRow,1,page};  % Detected period
    Results{iRow,5} = amplitudesLayers{iRow,1,page}; % Amplitude
    Results{iRow,6} = phasesLayers{iRow,1,page}; % Phase
    Results{iRow,7} = rsquaresLayers(iRow,1,page); % R-square
    Results{iRow,8} = pvaluesLayers(iRow,1,page);  % Pvalue
end

% Calculate FDR (FDR-BHs, qvalues)
pvalues = cell2mat(Results(:,8));
qvalues=[];
% In some special cases, qvalues cannot be calculated (represented as NaNs).
% For example, when microarray data contains constant or monotonic time-series,
% so we suggest to remove them before analyzing
try
    [~,qvalues] = mafdr(pvalues);
catch
    qvalues(1:length(pvalues),1) = NaN;% Fail to generate qvalues
end
[FdrBHs] = mafdr(pvalues, 'BHFDR', true);
for iRow=1:mRows
    Results{iRow,9}=qvalues(iRow); % qvalues
    Results{iRow,10}=FdrBHs(iRow); % FDR-BH
end



% SAVE RESULT
fprintf('Saving...\n');
expandOutputPath = strcat(outputPath,outputFilename);
saveFile(expandOutputPath,Results,totalGenesNum,probeNames,validSign);