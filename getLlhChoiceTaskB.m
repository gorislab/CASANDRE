function [choiceLlh] = getLlhChoiceTaskB(stimValue, modelParams, calcPrecision, asymFlag)
%% getLlhChoiceTaskB
% Takes as input: 
% stimulusValue - different stimulus conditions in units of stimulus magnitude
% modelParams   - CASANDRE paramters in order [guessRate, stimSens, stimCritLow, stimCritHigh, uncMeta, confCrit]
% calcPrecision - calcPrecision(1) is the sample rate; 25 + recommended
% asymFlag      - fit CASANDRE with or without asymmetrical confidence criteria 

% Output:
% choiceLlh     - likelihood of each choice (2 * N confidence levels x N stimValues) 

% Decode function arguments
stimVal      = stimValue;                   % The different stimulus conditions in units of stimulus magnitude (e.g., orientation in degrees)
noiseSens    = 1;                           % If the sensory noise is set to 1, then distributions of decision variable and confidence variable can be compared directly
guessRate    = modelParams(1);              % The fraction of guesses
stimSens     = modelParams(2);              % Stimulus sensitvity parameter, higher values produce a steeper psychometric function, strictly positive
stimCritLow  = modelParams(3);              % The first sensory decision criterion in units of stimulus magnitude (e.g., orientation in degrees)
stimCritHigh = modelParams(4);              % The second sensory decision criterion in units of stimulus magnitude (e.g., orientation in degrees)
uncMeta      = modelParams(5);              % Meta-uncertainty: the second stage noise parameter, only affects confidence judgments, strictly positive
confCrit     = cumsum(modelParams(6:end));  % The confidence criteria, unitless

% Set calculation precision
sampleRate = calcPrecision(1);             % Higher values produce slower, more precise estimates. Precision saturates after ~25

%% Compute model prediction
if asymFlag == 0
    nRespOpt = 2*(numel(confCrit)+1);
elseif asymFlag == 1
    nRespOpt = 2*(numel(confCrit)/2+1);
end

% Step 0 - rescale sensory representation by sensitivity parameter
sensMean     = stimVal*stimSens;
sensCritLow  = stimCritLow*stimSens;
sensCritHigh = stimCritHigh*stimSens;
sensCritDist = sensCritHigh - sensCritLow;       % The distance between both perceptual criteria

for iC = 1:numel(stimVal)
    
    %% Compute llh of each response alternative
    % Step 1 - sample decision variable denominator in steps of constant cumulative density
    muLogN    = log((noiseSens.^2)./sqrt(uncMeta.^2 + noiseSens.^2));
    sigmaLogN = sqrt(log((uncMeta.^2)./(noiseSens.^2) + 1));
    dv_Den_x  = logninv(linspace(.5/sampleRate, 1-(.5/sampleRate), sampleRate), muLogN, sigmaLogN);
    
    % Step 2 - compute choice distribution under each scaled sensory distribution
    % Crucial property: linear transformation of normal variable is itself normal variable
    % Trick: we take inverse of denominator to work with products instead of ratios
    muCritLow  = (1./dv_Den_x').*(sensMean(iC) - sensCritLow);
    muCritHigh = (1./dv_Den_x').*(sensMean(iC) - sensCritHigh);
    sigma      = (1./dv_Den_x').*noiseSens;
    
    % Include mid point between both perceptual criteria
    % The mid-point determines which perceptual criterion informs distance
    % estimate for confidence computation
    xCritLow  = sort([-confCrit, 0, confCrit, sensCritDist/2]);
    xCritHigh = sort([-confCrit, 0, confCrit, -sensCritDist/2]);
    try if asymFlag == 1; confCrit = modelParams(6:end); xCritLow  = sort([-cumsum(confCrit(numel(confCrit)/2 + 1:end)), 0, cumsum(confCrit(1:numel(confCrit)/2)), sensCritDist/2]); end; catch; end
    try if asymFlag == 1; confCrit = modelParams(6:end); xCritHigh = sort([-cumsum(confCrit(1:numel(confCrit)/2)), 0, cumsum(confCrit(numel(confCrit)/2 + 1:end)), -sensCritDist/2]); end; catch; end
    xCritLow  = unique(xCritLow(xCritLow <= sensCritDist/2));
    xCritHigh = unique(xCritHigh(xCritHigh >= -sensCritDist/2));
    PCritLow  = normcdf(repmat(xCritLow, [sampleRate 1]), repmat(muCritLow, [1 numel(xCritLow)]), repmat(sigma, [1 numel(xCritLow)]));
    PCritHigh = normcdf(repmat(xCritHigh, [sampleRate 1]), repmat(muCritHigh, [1 numel(xCritHigh)]), repmat(sigma, [1 numel(xCritHigh)]));
    
    
    % Step 3 - average across all scaled sensory distributions to get likelihood functions
    ratio_dist_pCritLow  = mean(PCritLow);
    ratio_dist_pCritHigh = mean(PCritHigh);
    
    for iX = 1:numel(xCritLow)
        if iX == 1
            llhC{iX}(iC) = (ratio_dist_pCritLow(1) + (1 - ratio_dist_pCritHigh(numel(xCritHigh))));
        else
            llhC{iX}(iC) = (ratio_dist_pCritLow(iX) - ratio_dist_pCritLow(iX-1) + ...
                ratio_dist_pCritHigh(numel(xCritHigh)+2-iX) - ratio_dist_pCritHigh(numel(xCritHigh)+1-iX));
        end
    end
end
choiceLlh = cell2mat(llhC');

if size(choiceLlh,1) < nRespOpt
    choiceLlh(size(choiceLlh,1)+1:nRespOpt,:) = 0;
end

% Approximation is not perfect due to use of two different perceptual criteria, this trick ensures that llh functions sum to one
scalar = repmat((1 - sum(choiceLlh(1:nRespOpt/2,:))), [nRespOpt/2 1])./max(.00001, repmat(sum(choiceLlh(nRespOpt/2 + 1:end,:)), [nRespOpt/2 1]));
choiceLlh(nRespOpt/2 + 1:end,:) = scalar .* choiceLlh(nRespOpt/2 + 1:end,:);

% Effect of guessing
choiceLlh = guessRate/nRespOpt + (1 - guessRate) * choiceLlh;

% Inconvenient step, needed because of inconsistency in which category was
% labeled "stimulus 1" across scripts.
choiceLlh = flipud(choiceLlh);

