function  [probeNames,timepoints,microarrayData,totalGenesNum,validSign,stat] = loadFile(expandInputPath)
%Load data from input file
%
%   Usage:
%       [probeNames,timepoints,microarrayData,totalGenesNum,validSign,stat] = loadFile(expandInputPath)
%
%   INPUT£º
%       expandInputPath     - expand path to load the input file
%
%   OUTPUT:
%       probeNames          - a (N+1)x1 cell vector contains a tag and N probe names, N is the number of genes that can be analyzed
%       timepoints          - a vector contains all sample time points
%       microarrayData      - a NxM matrix representing N genes (probes) with M expression values
%       totalGenesNum       - number of genes contained in the input file
%       validSign           - a Nx1 vector to mark if one gene expression profile can be analyzed 
%       stat                - mark if load data from input file successfully

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

% Local variables
totalGenesNum = 0;
probeNames = cell(0,0);
nMicroarrayData = 0;
validSign = [];% Mark if one profile can be analyzed
microarrayData = [];
TAB = 9;% Horizontal tab (ASCII 9). Only tab delimited text file is allowed
glue = strcat('NaN',TAB);
inputFid = fopen(expandInputPath);% Get file ID

% Get first line
firstLineStr = fgetl(inputFid);
firstLines = strread(firstLineStr,'%s');
nTimePoints = length(firstLines)-1;
probeNames{1} = firstLines(1);

% Extract time points from the first line
for iTimePoint=1:nTimePoints
    timepoints(iTimePoint) = str2num(firstLines{iTimePoint+1});
end

% Check the validation of each profile
while ~feof(inputFid)
    totalGenesNum = totalGenesNum + 1;
    dataLines = fgetl(inputFid);
    strings = '';
    
    % Load data from each line
    for idataLine = 1: length(dataLines)
        if(dataLines(idataLine) ~= TAB)% Not a horizontal tab
            strings = [strings dataLines(idataLine)];
        else
            % If a dataline contains some strings formatted like TAB TAB, 
            % then a missing value will exists between TAB and TAB, 
            % so we insert a NaN in place of this missing value 
            if dataLines(idataLine)-dataLines(idataLine-1) == 0
                strings = [strings glue];% Replace with NaN
            else
                strings = [strings TAB];
            end
        end
    end
    
    tempRow = strread(strings,'%s');
    
    % A missing value exists at the end of the line
    if length(tempRow) ~= (nTimePoints + 1 )
        tempRow{nTimePoints + 1} = 'NaN';
    end
    
    % Count the number of missing values in each row
    nNaNs = 0;
    probeNames{totalGenesNum+1} = tempRow(1);
    for iTempRow = 2:length(tempRow)
        if isnan(str2double(tempRow{iTempRow}))
            nNaNs = nNaNs + 1;
        end
    end
    
    validSign(totalGenesNum)= 1;% Label invalid rows with -1
    
    % Ignore those time-series (rows) whose values are missing more than
    % 50% of the size of sample time points
    if (nTimePoints-nNaNs) < fix(0.5 * nTimePoints)
        validSign(totalGenesNum) = -1;% 
        continue;
    end
    
    % Save to microarrayData
    if validSign(totalGenesNum) == 1
        nMicroarrayData = nMicroarrayData + 1;
        for iTimepoint = 1:nTimePoints
            microarrayData(nMicroarrayData,iTimepoint) = str2double(tempRow{iTimepoint+1});
        end
    end
end

stat = 1; % Succeed in loading data from the input file
% Error: fail to load data
if isempty(microarrayData)
    stat = -1;
end

fclose(inputFid);
end