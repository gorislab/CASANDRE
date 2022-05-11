%% plotAdlerMa2018

% Code to  visualize CASANDRE fits to Adler & Ma (2018)
% Expt 1 data.

% Start with clean slate
clearvars;
clc;
close all;

% Set Paths
thisPath    = fullfile(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(thisPath,'..')));

% Choose example data-set  (two tasks, many observers)
plotExp  = 1;  % [1]
plotTask = 1;  % [1, 2]
plotObs  = 6;  % [1,...]

% Specify some plotting variables
nSamples    = 100;      % Resolution of predicted psychometric and confidence functions
nPoints     = [9, 11];   % Number of points in the observed psychometric and confidence functions [task 1, task 2]

% Calculation precision
sampleRate    = 100;        % Higher values produce slower, more precise estimates.
delta         = 5;          % Number of standard deviations below and above mean, used to compute confidence variable distributions
calcPrecision = [sampleRate, delta];
asymFlag      = 1;


set(figure(1), 'OuterPosition', [100 100 2000 1000])
set(figure(2), 'OuterPosition', [100 100 2000 1000])

iE = plotExp;
iT = plotTask;
iO = plotObs;

% Specify load name associated with each experiment
if ~asymFlag
    loadName = strcat('Adler_2018_Expt', int2str(iE), '_sorted.mat');
elseif asymFlag
    loadName = strcat('Adler_2018_Expt', int2str(iE), '_sorted_ASYM_CC.mat');
end

% Load data
load(loadName);

% Step 1: unpack data
nObs(iE)  = numel(trials);
nTasks    = numel(trials{1});
nRel      = numel(trials{1}{2});
nConfCrit = 3;

for iR = 1:nRel
    
    % Set plot color
    col = [1-iR/nRel 0 iR/nRel];
    
    % Set experiment parameters
    stimCat{iR}   = trials{iO}{iT}{iR}(:,1);
    stimValue{iR} = trials{iO}{iT}{iR}(:,2);   % The different stimulus conditions in units of stimulus magnitude (e.g., orientation in degrees)
    choice        = trials{iO}{iT}{iR}(:,3);
    
    % Compute observed psychometric and confidence functions (this requires grouping of stimuli)
    edges = prctile(stimValue{iR}, linspace(0, 100, nPoints(iT) + 1));
    [nTrials, edges, binInd] = histcounts(stimValue{iR}, edges);
    for iP = 1:nPoints(iT)
        oriPF(iP) = mean(stimValue{iR}(binInd == iP));
        obsPF(iP) = mean((choice(binInd == iP)) > 0);
        obsCF(iP) = mean(abs(choice(binInd == iP)));
        
        % Psychometric for low and high confidence trials
        indHighConf    = abs(choice) >= 3;
        obsPF_HC(iP)   = mean((choice(binInd == iP & indHighConf)) > 0);
        obsPF_LC(iP)   = mean((choice(binInd == iP & ~indHighConf)) > 0);
        nTrials_HC(iP) = sum(binInd == iP & indHighConf);
        nTrials_LC(iP) = sum(binInd == iP & ~indHighConf);
    end
    
    % First specify function arguments
    stimPlot = linspace(min(stimValue{iR}), max(stimValue{iR}), nSamples);
    
    if plotTask == 1
        % Required order for getLlhChoice: [guess rate, stim sens, stim crit, meta uncertainty, conf criteria]
        params   = bestParamEst{iO}{iT}([1, 2+iR, 2+nRel+iR, 2, 3+(2*nRel):end]);
        paramsNM = [params(1:3), 0.01, params(5:end)];
        paramsNG = [0, params(2:end)];

        % Get model predictions
        fitLlh   = getLlhChoice(stimPlot, params, calcPrecision, asymFlag);
        fitLlhNG = getLlhChoice(stimPlot, paramsNG, calcPrecision, asymFlag);
        
        % The PF and CF predicted on the basis of the likelihood functions
        predPF = sum(fitLlh(size(fitLlh, 1)/2+1:end, :));
        predCF = fitLlh' * [4 3 2 1 1 2 3 4]'; % multiply choice-likelihood by confidence ratings
        
        % Predicted CF without meta-uncertainty
        predCFnoMN  = getLlhChoice(stimPlot, paramsNM, calcPrecision, asymFlag)' * [4 3 2 1 1 2 3 4]';
        
    elseif plotTask == 2
        % Required order for getLlhChoiceTaskB: [guess rate, stim sens, stim crit low, stim crit high, meta uncertainty, conf criteria]
        params   = bestParamEst{iO}{iT}([1, 2+iR, 2+nRel+iR, 2+(2*nRel)+iR, 2, 3+(3*nRel):end]);
        paramsNM = [params(1:4), 0.01, params(6:end)];
        paramsNG = [0, params(2:end)];
        
        % Get model predictions
        fitLlh   = getLlhChoiceTaskB(stimPlot, params, calcPrecision, asymFlag);
        fitLlhNG = getLlhChoiceTaskB(stimPlot, paramsNG, calcPrecision, asymFlag);
        
        % The PF and CF predicted on the basis of the likelihood functions
        predPF = sum(fitLlh(size(fitLlh, 1)/2+1:end, :));
        predCF = fitLlh' * [4 3 2 1 1 2 3 4]';
        
        % Predicted CF without meta-uncertainty
        predCFnoMN  = getLlhChoiceTaskB(stimPlot, paramsNM, calcPrecision, asymFlag)' * [4 3 2 1 1 2 3 4]';
    end
    
    % Predicted PF without guesses, split out for low and high confidence trials
    predPF_HC = sum(fitLlhNG(7:8, :))./sum(fitLlhNG([1:2, 7:8], :));
    predPF_LC = sum(fitLlhNG(5:6, :))./sum(fitLlhNG([3:4, 5:6], :));
    
    % Plot PF
    figure(1)
    subplot(4,6,iR)
    plot(stimPlot, predPF, '-', 'linewidth', 2, 'color', col)
    hold on, box off, axis square
    axis([-20 20 0 1])
    xlabel('Stimulus value')
    ylabel('Proportion category 1')
    
    for iP = 1:nPoints(iT)
        plot(oriPF(iP), obsPF(iP), 'ko', 'markerfacecolor', col, 'markersize', round(nTrials(iP)/5)+1)
    end
    
    if iR == 6
        legend('Full model', 'location', 'NorthWest')
    end
    
    % Plot CF
    figure(1)
    subplot(4,6, nRel+iR)
    plot(stimPlot, predCF, '-', 'linewidth', 2, 'color', col)
    hold on, box off, axis square
    plot(stimPlot, predCFnoMN, 'k--', 'linewidth', 2)
    axis([-20 20 1 4])
    xlabel('Stimulus value')
    ylabel('Mean confidence level')
    
    for iP = 1:nPoints(iT)
        plot(oriPF(iP), obsCF(iP), 'ko', 'markerfacecolor', col, 'markersize', round(nTrials(iP)/5)+1)
    end
    
    if iR == 6
        legend('Full model', 'meta uncertainty = 0', 'location', 'NorthWest')
    end
    
    % Plot PF split out for high and low confidence trials
    figure(1)
    subplot(4,6,2*nRel+iR)
    plot(stimPlot, predPF_HC, '-', 'linewidth', 2, 'color', [0 1 0])
    hold on, box off, axis square
    plot(stimPlot, predPF_LC, '-', 'linewidth', 2, 'color', [1 0 0])
    axis([-20 20 0 1])
    xlabel('Stimulus value')
    ylabel('Proportion category 1')
    
    for iP = 1:nPoints(iT)
        plot(oriPF(iP), obsPF_HC(iP), 'ko', 'markerfacecolor', [0 1 0], 'markersize', round(nTrials_HC(iP)/5)+1)
        plot(oriPF(iP), obsPF_LC(iP), 'ko', 'markerfacecolor', [1 0 0], 'markersize', round(nTrials_LC(iP)/5)+1)
    end
    
    if iR == 6
        legend('High confidence', 'Low confidence', 'location', 'NorthWest')
    end

    % Plot confidence vs performance
    figure(2)
    subplot(1,3,1)
    plot(predPF, predCF, '-', 'linewidth', 2, 'color', col)
    hold on, box off, axis square
    for iP = 1:nPoints(iT)
        plot(obsPF(iP), obsCF(iP), 'ko', 'markerfacecolor', col, 'markersize', round(nTrials(iP)/5)+1)
    end
    axis([0 1 1 4])
    xlabel('Proportion category 1')
    ylabel('Mean confidence level')
    
end

