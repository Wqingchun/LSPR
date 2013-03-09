function saveFile(expandOutputPath,Results,totalGenesNum,probeNames,validSign)
%Save data to output file
%
%   Usage:
%      saveFile(expandOutputPath,Results,totalGenesNum,probeNames,validSign)
%
%   INPUT:
%       expandOutputPath    - expand path to load the output file
%       Results             - a matrix saves all computing variables, see more in LSPR.m
%       totalGenesNum       - number of genes contained in the input file
%       probeNames          - a (N+1)x1 cell vector contains a tag and N probe names, N is the number of genes that can be analyzed
%       validSign           - a Nx1 vector to mark if one gene expression profile can be analyzed 
%
%   OUTPUT:
%       none

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

outputFid = fopen(expandOutputPath,'w');

% Generate the first line
firstLine = cell(1,10);
firstLine{1,1} = 'filter type';
% '-1' - non-filtered
% ' 1' - filtered
firstLine{1,2} = 'method';
% 'LSPR'    - harmonic analysis with detected periods
% 'default' - harmonic analysis with default period in [lower,upper]
firstLine{1,3} = 'number of oscillations';% Number of oscillations
firstLine{1,4} = 'period'; % Detected periods
firstLine{1,5} = 'amplitude'; % Amplitude for single cosine model
firstLine{1,6} = 'phase'; % Phase for single cosine model
firstLine{1,7} = 'R square';% R square of regression curve
firstLine{1,8} = 'pvalue';% p value after harmonic regression
firstLine{1,9} = 'qvalue';% q value
firstLine{1,10}= 'FDR-BH';% Fdr-BH

% Print probe names
fprintf(outputFid,'%s\t',cell2mat(probeNames{1}));
for m=1:size(firstLine,2)
    fprintf(outputFid,'%s\t',firstLine{1,m});
end
fprintf(outputFid,'\n');

% Print other variables
rowCounter = 1;
for l=1:totalGenesNum
    fprintf(outputFid,'%s\t',cell2mat(probeNames{l+1}));
    if validSign(l) == -1 % No output variables
        fprintf(outputFid,'NaN\t');
        fprintf(outputFid,'NaN\t');
        fprintf(outputFid,'NaN\t');
        fprintf(outputFid,'NaN\t');
        fprintf(outputFid,'NaN\t');
        fprintf(outputFid,'NaN\t');
        fprintf(outputFid,'NaN\t');
        fprintf(outputFid,'NaN\t');
        fprintf(outputFid,'NaN\t');
        fprintf(outputFid,'NaN\t');
        fprintf(outputFid,'\n');
        continue;
    end  
    % Print filter type, method and number of oscillations
    fprintf(outputFid,'%g\t',Results{rowCounter,1});
    fprintf(outputFid,'%s\t',Results{rowCounter,2});
    fprintf(outputFid,'%g\t',Results{rowCounter,3});
    % Print periods
    if ~isempty(Results{rowCounter,4});
        fprintf(outputFid,'%g',Results{rowCounter,4}(1));
        i = Results{rowCounter,3};
        while i  > 1 % More than one period
            fprintf(outputFid,',%g',Results{rowCounter,4}(i));
            i = i-1;
        end
    end
    fprintf(outputFid,'\t');
    % Print amplitudes
    if ~isempty(Results{rowCounter,5});
        fprintf(outputFid,'%g',Results{rowCounter,5}(1));
        i = Results{rowCounter,3};
        while i  > 1 % More than one amplitude
            fprintf(outputFid,',%g',Results{rowCounter,5}(i));
            i = i-1;
        end
    end
    fprintf(outputFid,'\t'); 
    % Print phases
    if ~isempty(Results{rowCounter,6});
        fprintf(outputFid,'%g',Results{rowCounter,6}(1));
        i = Results{rowCounter,3};
        while i  > 1 % More than one phase
            fprintf(outputFid,',%g',Results{rowCounter,6}(i));
            i = i-1;
        end
    end
    fprintf(outputFid,'\t');
    % Print R square, pvalue, qvalue and FDR-BH
    fprintf(outputFid,'%g\t',Results{rowCounter,7});
    fprintf(outputFid,'%g\t',Results{rowCounter,8});
    fprintf(outputFid,'%g\t',Results{rowCounter,9});
    fprintf(outputFid,'%g\t',Results{rowCounter,10});
    fprintf(outputFid,'\n');
    rowCounter = rowCounter +1;
end
fclose(outputFid);
end