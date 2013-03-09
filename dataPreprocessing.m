function [detrendedData,filteredData] = dataPreprocessing(inputData)
%Detrend or filter data
%
%   Usage:
%       [detrendedData,filteredData] = dataPreprocessing(inputData)
%
%   INPUT£º
%       inputData           - a NxM matrix representing N genes (probes) with M expression values over time 
%
%   OUTPUT:
%       detrendedData       - a NxM matrix representing N genes (probes) with M detrended values 
%       filteredData        - a NxM matrix representing N genes (probes) with M detrended and filtered values 

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

mRows = size(inputData,1);

% Detrending
fprintf('Detrending...\n');
detrendedData = zeros(mRows, size(inputData,2));
% Remove linear trends of inputData
for iRow=1:mRows
    getRow = inputData(iRow,:);% Get valid data from each row
    detrendedData(iRow,:) = getRow;
    index = find(~isnan(getRow));
    tempRow = getRow(index);    
    % Find constant and monotonic time-series, then adjust their detrended
    % values, if not, analyzing them may lead to some unexplainable results
    if find(detrend(tempRow)>1e-6)
        detrendedData(iRow,index) = detrend(tempRow);
    else
        detrendedData(iRow,index) = 0;% Constant and monotonic time-series
    end
end



% Filtering
fprintf('Filtering...\n');
filteredData = zeros(mRows,size(inputData,2));
% Use Savitzky-Golay smoothing filter
for iRow=1:mRows
    getRow= detrendedData(iRow,:);
    filteredData(iRow,:) =  getRow;
    index = find(~isnan(getRow));
    tempRow = getRow(index);
    % Use try-catch structure to choose arguments of S-G filter. In some cases, 
    % some time-series cannot be filtered by sgolayfilt with 4th order, so 
    % it is necessary to replace with a smaller polynomial order. 
    try
        filteredData(iRow,index) =sgolayfilt(tempRow,4,11);
    catch
        filteredData(iRow,index) =sgolayfilt(tempRow,2,5);
    end
end
end