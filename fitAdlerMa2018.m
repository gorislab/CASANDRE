function [] = fitAdlerMa2018()

% Start with clean slate
clearvars;
clc;

% Set Paths
thisPath    = fullfile(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(thisPath,'..')));

% Specify data-set
expt = [1];

% Set parameters
nRuns         = 3;         % Number of restarts for model fit
sampleRate    = 100;       % Higher values produce slower, more precise estimates.
delta         = 5;         % Number of standard deviations below and above mean, used to compute confidence variable distributions
calcPrecision = [sampleRate, delta];
asymFlag      = 1;         % if 1, fit CASANDRE with asymmetrical confidence criteria; if 0, fit CASANDRE with symmetrical confidence criteria 

% Set options
options               = optimoptions('fmincon');
options.MaxIterations = 25;
options.Display       = 'off';

% Load data
loadName = strcat('Adler_2018_Expt', int2str(expt), '_sorted.mat');
load(loadName, 'trials');

% Step 1: unpack data
nObs      = numel(trials);
nRel      = numel(trials{1}{2});
nConfCrit = 3;              % Number of confidence criteria if asymFlag = 0;

for iL = 1:nRuns
    for iO = 1:nObs
        for iT = 1 % 1:2 % task index
            taskID = iT;
            
            fprintf('fitting expt %d, observer %d, task %d, loop %d... \n', expt, iO, iT, iL)
            
            if taskID == 1
                nParams = 2 + nRel + nRel + nConfCrit; % [Guess rate, meta-noise], [stimulus sensitivity], [stimulus criterion], [confidence criteria]
            elseif taskID == 2
                nParams = 2 + nRel + 2*nRel + nConfCrit; % [Guess rate, meta-noise], [stimulus sensitivity], [stimulus criteria], [confidence criteria]
            end
            
            if asymFlag == 1
               nParams = nParams + nConfCrit; %account for additional parameters if using asymmetrical confidence critiera 
            end

            for iR = 1:nRel                
                % Set experiment parameters
                stimValue{iR} = trials{iO}{iT}{iR}(:,2);   % The different stimulus conditions in units of stimulus magnitude (e.g., orientation in degrees)
                choice        = trials{iO}{iT}{iR}(:,3);   % Subject choice (sign = perceptual choice; value = confidence rating)
                
                % Recode responses
                respOptions   = [-(nConfCrit+1):-1, 1:(nConfCrit+1)];
                choiceVec{iR} = zeros(2*(nConfCrit+1), numel(stimValue{iR}));
                
                for iC = 1:numel(respOptions)
                    selTrials = (choice == respOptions(iC));
                    choiceVec{iR}(iC, selTrials) = 1;
                end
            end
            
            % Fit data-set (i.e., 1 obs, 1 task, all levels of stimulus uncertainty)
            % Set parameter bounds
            LB  = zeros(nParams,1);
            UB  = zeros(nParams,1);
            
            LB(1,1)                 = 0;           UB(1,1)                 = 0.1;       % Guess rate
            LB(2,1)                 = 0.1;         UB(2,1)                 = 5;         % Meta uncertainty
            LB(3:2+nRel,1)          = 0.005;       UB(3:2+nRel,1)          = 1;         % Stimulus sensitivity
            
            if taskID == 1
                LB(3+nRel:2+(2*nRel),1) = -20;         UB(3+nRel:2+(2*nRel),1) = 20;       % Stimulus criterion
                LB(3+(2*nRel):end,1)    = 0;           UB(3+(2*nRel):end, 1)   = 10;       % Confidence criteria
            elseif taskID == 2
                LB(3+nRel:2+(2*nRel),1)     = -20;     UB(3+nRel:2+(2*nRel),1)     = 0;     % Stimulus criterion 1, negative
                LB(3+(2*nRel):2+(3*nRel),1) = 0;       UB(3+(2*nRel):2+(3*nRel),1) = 20;    % Stimulus criterion 2, positive
                LB(3+(3*nRel):end,1)        = 0;       UB(3+(3*nRel):end, 1)       = 10;    % Confidence criteria
            end
            
            % Define objective function
            obFun = @(paramVec) giveNLL(paramVec, stimValue, choiceVec, taskID, calcPrecision, asymFlag);
            
            % Initialize parameters
            try load(loadName, 'bestParamEst');
                    startVec = bestParamEst{iO}{iT};
            catch
                if taskID == 1
                    startVec = [.01 0.5 0.6 0.4 0.25 0.1 0.05 0.02 0 0 0 0 0 0 rand(1, nConfCrit)];
                elseif taskID == 2
                    startVec = [.01 0.5 0.6 0.4 0.25 0.1 0.05 0.02 -5 -6 -7 -8 -9 -10 5 6 7 8 9 10 rand(1, nConfCrit)];
                end
                if asymFlag == 1
                    startVec = [startVec, rand(1, nConfCrit)];
                end
            end
                           
            % Introduce some noise
             startVec = startVec .* .9+.2*rand(size(startVec));

            % Make sure start is within bounds
            startVec = max(min(startVec, UB' - 0.001), LB' + 0.001);
            
            % Fit model
            [paramEst{iO}{iT}, NLL{iO}{iT}] = fmincon(obFun, startVec, [], [], [], [], LB, UB, [], options);
            
            % Fit model without meta-uncertainty
            LB(2,1) = 0.001;    UB(2,1) = 0.002;   startVec(2) = 0.0015;
            [paramEst_NM{iO}{iT}, NLL_NM{iO}{iT}] = fmincon(obFun, startVec, [], [], [], [], LB, UB, [], options);
            
            % Compare current fits with best fits
            try load(loadName, 'bestParamEst', 'bestNLL');
                if bestNLL{iO}{iT} <= NLL{iO}{iT}
                else
                    bestNLL{iO}{iT}      = NLL{iO}{iT};
                    bestParamEst{iO}{iT} = paramEst{iO}{iT};
                    fprintf('New best fit found for full model! \n')
                end
            catch
                bestNLL{iO}{iT}      = NLL{iO}{iT};
                bestParamEst{iO}{iT} = paramEst{iO}{iT};
            end
            
            try load(loadName, 'bestParamEst_NM', 'bestNLL_NM');
                if bestNLL_NM{iO}{iT} <= NLL_NM{iO}{iT}
                else
                    bestNLL_NM{iO}{iT}      = NLL_NM{iO}{iT};
                    bestParamEst_NM{iO}{iT} = paramEst_NM{iO}{iT};
                    fprintf('New best fit found for meta-uncertainty-less model! \n')
                end
            catch
                bestNLL_NM{iO}{iT}      = NLL_NM{iO}{iT};
                bestParamEst_NM{iO}{iT} = paramEst_NM{iO}{iT};
            end
            
            % Save fits
            save(loadName, 'trials', 'bestParamEst', 'bestNLL', 'bestParamEst_NM', 'bestNLL_NM');
        end
    end
end
end


function [NLL] = giveNLL(paramVec, stimValue, choiceVec, taskID, calcPrecision, asymFlag)

% Unpack data
nRel = numel(stimValue);  % n levels of stimulus uncertainty

for iR = 1:nRel
    % First select appropriate subset of parameters, then get NLL for each level of stimulus uncertainty
    if taskID == 1
        % Required order for getLlhChoice: [guess rate, stim sens, stim crit, meta-uncertainty, conf criteria]
        params    = paramVec([1, 2+iR, 2+nRel+iR, 2, 3+(2*nRel):end]);
        choiceLlh = getLlhChoice(stimValue{iR}, params, calcPrecision, asymFlag);
    elseif taskID == 2
        % Required order for getLlhChoiceTaskB: [guess rate, stim sens, stim crit low, stim crit high, meta-uncertainty, conf criteria]
        params    = paramVec([1, 2+iR, 2+nRel+iR, 2+(2*nRel)+iR, 2, 3+(3*nRel):end]);
        choiceLlh = getLlhChoiceTaskB(stimValue{iR}, params, calcPrecision, asymFlag);
    end    
    NLLvec(iR) = -sum(sum(choiceVec{iR}.*log(choiceLlh)));
end

% Compute NLL across entire data-set
NLL = sum(NLLvec);
end

