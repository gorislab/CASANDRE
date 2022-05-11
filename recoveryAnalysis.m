%% recoveryAnalysis

% Script to generate and fit data with CASANDRE process model. In the paper
% associated with CASADRE we explored the recoverability of model parameters
% with changing the unique number of stimulus values (parameter: stimValue)
% and repetitions per stimulus value (parameter: stimReps).

% Parameters used to generate figure 4A & C: 
% guessRate   = 0;
% stimSens    = 1;
% stimCrit    = 0
% uncMeta     = [0.2 04 0.8 1.6 3.2];
% confCrit    = 0.75;
% asymFlag    = 0;

close all;
clearvars;
clc;

% Set experiment parameters
stimValue = linspace(-3, 3, 11);   % The different stimulus conditions in units of stimulus magnitude (e.g., orientation in degrees)
stimReps  = 200;                   % The number of repeats per stimulus

% Set model parameters
guessRate   = 0.000;                % The fraction of guesses
stimSens    = .5;                   % Stimulus sensitvity parameter, higher values produce a steeper psychometric function, strictly positive
stimCrit    = 0;                    % The sensory decision criterion in units of stimulus magnitude (e.g., orientation in degrees)
uncMeta     = .5;                   % Meta-uncertainty: the second stage noise parameter, only affects confidence judgments, strictly positive
confCrit    = [.75 1];              % The confidence criteria, unitless (can include more than 1)
asymFlag    = 0;                    % If set to 1, it allows for asymmetrical confidence criteria and confCrit needs two times as many elements    
modelParams = [guessRate, stimSens, stimCrit, uncMeta, confCrit];
modelParamsLabel = [{'guessRate', 'stimSens', 'stimCrit', 'uncMeta'} repmat({'confCrit'},1,numel(confCrit))];

% Set calulation precision
calcPrecision = 100;                % Higher values produce slower, more precise estimates. Precision saturates after ~25

% Get model predictions
[choiceLlh] = getLlhChoice(stimValue, modelParams, calcPrecision, asymFlag);

% Simulate choice data
randNumbers = rand(stimReps, numel(stimValue));
criteria    = cumsum(choiceLlh);

for iX = 1:size(criteria, 1)
    if iX == 1
        n{iX} = sum(randNumbers <= criteria(1,:));
    elseif (iX > 1 && iX < size(criteria, 1))
        n{iX} = sum((randNumbers > criteria(iX-1,:)) & (randNumbers <= criteria(iX,:)));
    elseif iX == size(criteria, 1)
        n{iX} = sum(randNumbers > criteria(end-1,:));
    end
end
nChoice  = cell2mat(n');

% Fit simulated data
options  = optimset('Display', 'off', 'Maxiter', 10^5, 'MaxFuneval', 10^5);
obFun    = @(paramVec) giveNLL(paramVec, stimValue, nChoice, calcPrecision, asymFlag);
startVec = [.01 1 -0.1 0.5 sort(2*rand(1,numel(confCrit)))];

% Search bounds:
LB          = zeros(numel(startVec),1);
UB          = zeros(numel(startVec),1);

LB(1,1)     = 0;                        UB(1,1)        = 0.1;                   % Guess rate
LB(2,1)     = 0;                        UB(2,1)        = 10;                    % Stimulus sensitivity
LB(3,1)     = -3;                       UB(3,1)        = 3;                     % Stimulus criterion
LB(4,1)     = 0.01;                     UB(4,1)        = 5;                     % Meta uncertainty 
LB(5:end,1) = 0;                        UB(5:end,1)    = 5;                     % Confidence criteria
paramEst    = fmincon(obFun, startVec, [], [], [], [], LB, UB, [], options);

% Computations for plotting
stimPlot = linspace(stimValue(1), stimValue(end), 100);
genLlh   = getLlhChoice(stimPlot, modelParams,calcPrecision, asymFlag);         % Likelihood of choices for ground truth parameters
fitLlh   = getLlhChoice(stimPlot, paramEst,calcPrecision, asymFlag);            % Likelihood of choices for fit parameters

genPF = sum(genLlh(size(genLlh, 1)/2+1:end, :));                                % Psychometric function predicted by generating parameters
fitPF = sum(fitLlh(size(genLlh, 1)/2+1:end, :));                                % Psychometric function predicted by model fit
obsPF = sum(nChoice(size(genLlh, 1)/2+1:end, :))/stimReps;                      % Observed psychometric function

%% Plot results
set(figure(1), 'OuterPosition', [100 100 500 1000])
subplot(3,1,1)
plot([-.5 ceil(max([modelParams, paramEst]))], [-.5 ceil(max([modelParams, paramEst]))], 'k--')
hold on, box off
for iP = 1:numel(modelParams)
       paramPlot(iP) = plot(modelParams(iP), paramEst(iP), 'o', 'markerfacecolor', [.8 .8 .8], 'markersize', 10, 'linewidth', 1);
       paramPlotLegend{iP} = modelParamsLabel{iP};
end
legend(paramPlot, paramPlotLegend,'location', 'southeast')
xlabel('Ground truth')
ylabel('Parameter estimate')

subplot(3,1,2)
plot(stimPlot, genPF, 'r-', 'linewidth', 2)
hold on, box off
plot(stimPlot, fitPF, 'k--', 'linewidth', 2)
plot(stimValue, obsPF, 'ko', 'linewidth', 1, 'markerfacecolor', [1 0 0], 'markersize', 12)
axis([-3 3 0 1])
xlabel('Stimulus value')
ylabel('Proportion clockwise')
legend('ground truth', 'model fit', 'observations', 'Location', 'NorthWest')

subplot(3,1,3)
plot(stimPlot, genLlh', 'r-', 'linewidth', 2)
hold on, box off
plot(stimPlot, fitLlh, 'k--', 'linewidth', 2)
xlabel('Stimulus value')
ylabel('Probability')
keyboard


function [NLL] = giveNLL(paramVec, stimValue, nChoice, calcPrecision, asymFlag)
choiceLlh = getLlhChoice(stimValue, paramVec,calcPrecision, asymFlag);
NLL       = -sum(sum(nChoice.*log(choiceLlh)));
end

