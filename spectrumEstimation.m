function [peakSign,periods,periodsNum] = spectrumEstimation(inputData,lower,upper,defaultPeriod,timepoints)
%Detect periods of each time-series
%
%   Usage:
%       [peakSign,periods,periodsNum] = spectrumEstimation(inputData,inputData,upper,defaultPeriod,timepoints)
%
%   INPUT:
%       inputData      - a NxM matrix representing N genes (probes) with M expression values over time 
%       lower/upper    - endpoints of period range, LSPR try to detect periods of each genes in [lower,upper] 
%       defaultPeriod  - a given default period of microarry time-series(i.e.24 for circadian microarray data)
%       timepoints     - sample time points
%
%   OUTPUT:
%       peakSign       - a Nx1 cell vector contains probe names, N is the number of genes that can be analyzed
%       periods        - a Nx1 cell vector, N stands for genes (probes) numbers, each cell contains the periods detected by L-S
%       periodsNum     - a Nx1 vector contains the number of oscillations of each genes

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

fprintf('Lomb-Scargle periodogram analyzing ...\n');

nInputData = size(inputData,1);
periods = cell(nInputData,1);
periodsNum = [];
minFreq = 1/upper;
maxFreq = 1/lower;
peakSign = zeros(nInputData,1);

% Generate test frequencies
difference = zeros(1,length(timepoints)-1);
for m = 1 :(length(timepoints)-1)
    difference(m) = timepoints(m+1)-timepoints(m);
end
% The maximum of timepoints should be larger than possible periods of the inputData
testPeriods = min(difference):0.05:timepoints(end);% SampleRate: every 3 mins
testFreqs = sort((1./testPeriods),'descend');

for iInputData = 1:nInputData
    % Display the percentage of completion
    if(mod(iInputData,10) == 0)
        fprintf('     %3.2f percent complete...\n', 100*iInputData/nInputData);
    end
    
    % Remove invalid data
    getRow = inputData(iInputData,:);
    index = find(~isnan(getRow));
    tempRow = getRow(index);
    tempTimePoints = timepoints(index);
    
    % Lomb-Scargle periodogram analysis
    [powers,~] = lombFunction(tempTimePoints',tempRow',testFreqs);
    if findpeaks(powers)
        [~,locs] = findpeaks(powers);
        frequencies = testFreqs(locs);
        
        % Find peaks in [1/upper,1/lower]
        peakIndex = find((frequencies <= maxFreq)&(frequencies >= minFreq));
        if peakIndex
            periods{iInputData,1} = 1./frequencies(peakIndex);
            periodsNum(iInputData,1) = length(peakIndex);
            peakSign(iInputData,1) = 1;
        else
            % No peak shows in [1/upper,1/lower]
            periods{iInputData,1}= defaultPeriod;
            periodsNum(iInputData,1)= 1;
            peakSign(iInputData,1) = -1;
        end
    else
        % No peak shows in testFreqs
        periods{iInputData,1} = defaultPeriod;
        periodsNum(iInputData,1) = 1;
        peakSign(iInputData,1) = -1;
    end
end
end
%% Subfunction:
%%
function [Pn, Prob] = lombFunction(t, y, f)
%   
%   [Pn, Prob] = lomb(t, y, freq)
%		Uses Lomb's method to compute normalized
%		periodogram values "Pn" as a function of
%		supplied vector of frequencies "freq" for
%		input vectors "t" (time) and "y" (observations).
%		Also returned is probability "Prob" of same
%		length as Pn (and freq) that the null hypothesis
%		is valid.
%		x and y must be the same length.

% This file is an adaptation of lomb.m obtained from a forgotten source.
% We found an alternative source: http://w3eos.whoi.edu/12.747/notes/lect07/lomb.m

%		check inputs
if length(t) ~= length(y); error('t and y not same length');
	exit; end;
%	subtract mean, compute variance, initialize Pn
z = y - mean(y);
var = std(y);
N=length(f);
Pn=zeros(size(f));
%	now do main loop for all frequencies
for i=1:length(f)
    w=2*pi*f(i);
    if w > 0 
       twt = 2*w*t;
       tau = atan2(sum(sin(twt)),sum(cos(twt)))/2/w;
       wtmt = w*(t - tau);
       Pn(i) = (sum(z.*cos(wtmt)).^2)/sum(cos(wtmt).^2) + ...
		(sum(z.*sin(wtmt)).^2)/sum(sin(wtmt).^2);
     else
	Pn(i) = (sum(z.*t).^2)/sum(t.^2);
     end
end
%	and normalize by variance, compute probs
Pn=Pn/2/var.^2;
Prob = 1-(1-exp(-Pn)).^N;
for i=1:length(Pn)		% accomodate possible roundoff error
    if Prob(i) < .001
	Prob(i) = N*exp(-Pn(i));
    end
end
end