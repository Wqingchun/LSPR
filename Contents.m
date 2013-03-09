%% Introduction to LSPR Package
% LSPR is a Matlab package to detect periodic expression profiles in
% microarray time-series.
%% Package content
% There are 6 functions in the LSPR algorithm. 
%
%  1.  LSPR               - main function to perform LSPR algorithm.
%  2.  loadFile           - load microarray data from input text file
%  3.  dataPreprocessing  - preprocess input data with detrend or filter procedure
%  4.  spectrumEstimation - apply Lomb-Scargle periodogram to detect periods
%  5.  harmonicRegression - do harmonic analysis and compute AIC value of harmonic models 
%  6.  saveFile           - save computing results to output text file 
%
% The other one auxiliary function in spectrumEstimation.m:
%  7.  lombFunction       - lomb-Scargle periodogram function
%
% Each function is provided with a help and examples of its usage and can be
% displayed by typing |help function name| in the command window. E.g.|help LSPR|.

% Copyright (c) 2010, Chen ZHANG and Rendong Yang
% <http://bioinformatics.cau.edu.cn/LSPR>
% For comments, bugs and suggestions, please contact:
% - Chen ZHANG
% - zcreation@yahoo.cn