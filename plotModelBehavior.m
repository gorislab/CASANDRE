%% plotModelBehavior

% This script illustrates the predictions of a two-stage model of
% perceptual decision-making in a (2x2)AFC-task. In the first stage, a
% sensory stimulus gives rise to a 1-D normally distributed internal
% experience, which is compared with a fixed decision criterion. In the
% second stage, this internal sensation is (1) subtracted from the decision
% criterion, and (2) divided by the 'estimated sensory uncertainty' (i.e.
% the cross-repeats standard deviation of the internal experience). This
% term is on average accurate, but fluctuates from trial to trial according
% to a lognormal distribution. The ratio represents a confidence estimate
% and is in turn compared with a fixed confidence criterion, giving rise to
% 4 response alternatives in total ('category 1 - high confidence',
% 'category 1 - low confidence', 'category 2 - low confidence', and
% 'category 2 - high confidence').

% Parameters used to generate figure 2; order [ stimCrit, stimSens, confCrit, uncMeta, guessRate]
% 2A [0 1.75 1.75 .5 0];  [1 1.75   1.75 .5 0]
% 2B [0 1.75 1.75 .5 0];  [0 0.5833 1.75 .5 0]
% 2C [0 1.75 1.75 .5 0];  [0 0.75   1.75 .5 0]
% 2D [0 1.75 1.75 .5 0];  [0 1.75   1.75 2 0]
% 2E [0 1    1.75 .25 0]; [0 1      1.75 4 0]
% 2F [0 1    0.75 .25 0]; [0 1      0.75 4 0]; [0 1 0.75 Inf 0]
% 2G [0 1    0.5 .25 0];  [0 1      1.75 0.25 0]

%% Clear variables
clearvars; close all

%% Set Paths
thisPath    = fullfile(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(thisPath,'..')));

%% Stimulus parameters
stimVal = linspace(-3, 3, 49);   % The different stimulus conditions in units of stimulus magnitude (e.g., orientation in degrees)

%% Model parameters
% Don't touch
noiseSens   = 1;        % The sensory noise; if set to 1, then sensory variable and confidence variable distributions can be compared directly

% These can be changed
stimSens    = .5;       % Stimulus sensitvity parameter, higher values produce a steeper psychometric function, strictly positive
uncMeta     = .01;      % Meta-uncertainty: the second stage noise parameter, only affects confidence judgments, strictly positive
stimCrit    = 0;       % The sensory decision criterion in units of stimulus magnitude (e.g., orientation in degrees)
confCrit    = [.5];    % The confidence criterion, unitless [ use two values if asymFlag = 1]
asymFlag    = 0;       % If asymFlag = 0 - will use symmetrical confidence criteria; if asymFlag = 1 - will use asymmetrical confidence criteria and requires two value    
guessRate   = 0.00;    % The fraction of guesses, distributed uniformly over all response alternatives

%% Calculation precision
sampleRate = 100;       % Higher values produce slower, more precise estimates.
delta      = 5;         % Number of standard deviations below and above mean, used to compute confidence variable distributions

%% Plotting settings
stimPlotInd = 33;%[17, 25, 33]; % Indices for plotting

%% Specify function arguments
modelParams               = [guessRate, stimSens, stimCrit, uncMeta, confCrit];
modelParamsNM             = [guessRate, stimSens, stimCrit, 0.01, confCrit];        % Model parameters without meta-uncertainty (NM)
modelParamsNG             = [0, stimSens, stimCrit, uncMeta, confCrit];             % Model parameters with guess-rate fixed at zero (NG)
calcPrecision             = [sampleRate, delta];
probC1eval1               = double(stimVal < 0);                                    % The probability of getting rewarded for a 'C1' (category 1) response  
probC1eval1(stimVal == 0) = 0.5;

%% Get model predictions
[sensVarX, sensVarPdf]     = sensVarpdf(stimVal, modelParams, calcPrecision);                               % Get primary decision (sensory) variable
[confVarX, confVarPdf]     = confVarpdf(stimVal, modelParams, calcPrecision);                     % Get confidence variable
[choiceLlh]                = getLlhChoice(stimVal, modelParams, calcPrecision, asymFlag);                   % Get choice likelihood for given set of model parameters
[choiceLlhNG]              = getLlhChoice(stimVal, modelParamsNG, calcPrecision, asymFlag);                 % Get choice likelihood for given set of model parameters, meta-uncertainty = 0

%% Rename llh for semantic clarity 
llhHiC1 = choiceLlh(1,:);
llhLoC1 = choiceLlh(2,:);
llhLoC2 = choiceLlh(3,:);
llhHiC2 = choiceLlh(4,:);

%% The PF and CF predicted on the basis of the likelihood functions
predPF = llhLoC2 + llhHiC2;
predCF = llhHiC1 + llhHiC2;

%% PF split out for high and low confidence trials (effect of guesses removed)
predPF_HC = choiceLlhNG(4,:)./(choiceLlhNG(1,:) + choiceLlhNG(4,:));
predPF_LC = choiceLlhNG(3,:)./(choiceLlhNG(2,:) + choiceLlhNG(3,:));

%% Predicted CF, split out for correct and error trials (effect of guesses removed)
muConfC1        = (choiceLlhNG(1:2,:)./repmat(sum(choiceLlhNG(1:2,:)), [2 1]))' * [1 0]';   % The expected confidence given that C1 is the chosen response
muConfC2        = (choiceLlhNG(3:4,:)./repmat(sum(choiceLlhNG(3:4,:)), [2 1]))' * [0 1]';   % The expected confidence given that C2 is the chosen response
predCF_correct  = muConfC1.*probC1eval1' + muConfC2.*(1 - probC1eval1');
predCF_error    = muConfC1.*(1 - probC1eval1') + muConfC2.*probC1eval1';



%% Plot some output
%% Figure 1
set(figure(1), 'OuterPosition', [100 100 500 1000])
sensCrit    = stimCrit*stimSens;
if ~asymFlag
    confCritPos = sensCrit+confCrit;
    confCritNeg = sensCrit-confCrit;
elseif asymFlag
    confCritPos = sensCrit+confCrit(2);
    confCritNeg = sensCrit-confCrit(1);
end

% Panel 1 - sensory variable 
subplot(4,1,1)
plot(sensVarX, sensVarPdf(:,stimPlotInd), 'm-', 'linewidth', 2)
hold on, box off
plot([sensCrit sensCrit], [0 1], 'k--')
xlabel('Sensory Variable')
ylabel('Probability density')
axis([-5 5 0 1])
    
% Panel 2 - confidence variable 
subplot(4,1,2)
plot(sensVarX, sensVarPdf(:,stimPlotInd), 'm--', 'linewidth', 2)
hold on, box off
plot(confVarX, confVarPdf(:,stimPlotInd), 'b-', 'linewidth', 2)
plot([sensCrit sensCrit], [0 1], 'k--')
plot([confCritPos confCritPos], [0 1], 'r--')
plot([confCritNeg confCritNeg], [0 1], 'r--')
xlabel('Sensory/Confidence Variable')
ylabel('Probability density')
axis([-5 5 0 1])

% Panel 3 - confidence as a function of choice consistancy
subplot(4,1,3)
plot(predPF,predCF , 'k-', 'linewidth', 2)
hold on, box off
axis([0 1 0 1])
xlabel('Proportion category 1')
ylabel('Mean confidence level')

% Panel 4 - the response likelihood functions
subplot(4,1,4)
plot(stimVal, llhHiC2, 'c-', 'linewidth', 2)
hold on, box off
plot(stimVal, llhLoC2, 'c--', 'linewidth', 2)
plot(stimVal, llhLoC1, 'k--', 'linewidth', 2)
plot(stimVal, llhHiC1, 'k-', 'linewidth', 2)
axis([-3 3 0 1])
xlabel('Stimulus value')
ylabel('Probability')


%% Figure 2
set(figure(2), 'OuterPosition', [700 100 500 1000])

% Panel 1 - Psychometric function
subplot(4,1,1)
plot(stimVal, predPF, 'm-', 'linewidth', 2)
hold on, box off
plot(stimVal(stimPlotInd), predPF(stimPlotInd), 'ko', 'markerfacecolor', 'm', 'markersize', 10)
plot([stimCrit stimCrit], [0 1], 'k--')
axis([-3 3 0 1])
xlabel('Stimulus value')
ylabel('Proportion category 2')

% Panel 2 - Confidence function
subplot(4,1,2)
plot(stimVal, predCF, 'b-', 'linewidth', 2)
hold on, box off
plot(stimVal(stimPlotInd), predCF(stimPlotInd), 'ko', 'markerfacecolor', 'b', 'markersize', 10)
plot([stimCrit stimCrit], [0 1], 'k--')
plot([confCritPos confCritPos], [0 1], 'r--')
plot([confCritNeg confCritNeg], [0 1], 'r--')
axis([-3 3 0 1])
xlabel('Stimulus value')
ylabel('Proportion high confidence')

% Panel 3 - the psychometric function for high and low confidence trials
subplot(4,1,3)
plot(stimVal, predPF_HC, 'g--', 'linewidth', 2)
hold on, box off
plot(stimVal, predPF_LC, 'r--', 'linewidth', 2)
plot([stimCrit stimCrit], [0 1], 'k--')
axis([-3 3 0 1])
legend('High confidence', 'Low confidence', 'Location', 'NorthWest')
xlabel('Stimulus value')
ylabel('Proportion category 2')

% Panel 4 - the confidence function for correct and error trials
subplot(4,1,4)
plot(stimVal, predCF_correct, 'g-', 'linewidth', 2)
hold on, box off
plot(stimVal, predCF_error, 'r-', 'linewidth', 2)
plot([stimCrit stimCrit], [0 1], 'k--')
plot([confCritPos confCritPos], [0 1], 'r--')
plot([confCritNeg confCritNeg], [0 1], 'r--')
axis([-3 3 0 1])
legend('Correct', 'Error', 'Location', 'NorthWest')
xlabel('Stimulus value')
ylabel('Proportion high confidence')

%%
function [sensVarX, sensVarPdf] = sensVarpdf(stimValue, modelParams, calcPrecision)
%% sensVarpdf
% calculates the probability density function of the primary ('sens' or sensory) decsion varaible

% Decode function arguments
stimVal     = stimValue;                   % The different stimulus conditions in units of stimulus magnitude (e.g., orientation in degrees)
noiseSens   = 1;                           % If the sensory noise is set to 1, then distributions of decision variable and confidence variable can be compared directly
stimSens    = modelParams(2);              % Stimulus sensitvity parameter, higher values produce a steeper psychometric function, strictly positive
stimCrit    = modelParams(3);              % The sensory decision criterion in units of stimulus magnitude (e.g., orientation in degrees)
sampleRate  = 2 + 2*calcPrecision(1);      % To match confVarpdf
delta       = 3;

% Step 0 - rescale sensory representation by sensitivity parameter
sensMean  = stimVal*stimSens;
sensCrit  = stimCrit*stimSens;

% Step 1 - sample the sensory domain representation
sensVal = linspace(sensMean(1) - delta*noiseSens, sensMean(end) + delta*noiseSens, sampleRate);

% Step 2 - Probability density of sensory variable
sensVarX   = sensVal - sensCrit;
sensVarPdf = normpdf(repmat(sensVarX', [1 numel(stimVal)]), sensMean, noiseSens);
end

function [confVarX, confVarPdf] = confVarpdf(stimValue, modelParams, calcPrecision)
%% confVarpdf
% calculates the probability density function of the secondary (confidence)  varaible

% Decode function arguments
stimVal    = stimValue;                   % The different stimulus conditions in units of stimulus magnitude (e.g., orientation in degrees)
stimSens   = modelParams(2);              % Stimulus sensitvity parameter, higher values produce a steeper psychometric function, strictly positive
stimCrit   = modelParams(3);              % The sensory decision criterion in units of stimulus magnitude (e.g., orientation in degrees)
sampleRate = calcPrecision(1);
delta      = calcPrecision(2);

% Step 0 - rescale sensory representation by sensitivity parameter
sensCrit  = stimCrit*stimSens;

% Step 1 - reset confidence criteria for specific purpose of this function
sampleVec   = repmat(delta/(sampleRate), [1 sampleRate+1]);
modelParams = [0, modelParams(2:4), sampleVec];

% Step 2 - compute the vector over which the density function is evaluated
x        = sort([-cumsum(sampleVec), 0, cumsum(sampleVec)]);
confVarX = sensCrit + x(1:end-1) + diff(x(1:2)/2);

% Step 3 - construct the cumulative density function
confVarCdf  = cumsum(getLlhChoice(stimVal, modelParams, calcPrecision, 0));

% Step 4 - convert this into the probability density function
confVarPdf  = (numel(x)-1)/(x(end)-x(1)) * diff(confVarCdf(1:end-1,:));
end


