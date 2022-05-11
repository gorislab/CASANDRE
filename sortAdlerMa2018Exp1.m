%% sortAdlerMa2018Expt1

% This script sorts the data from their Experiment 1. This sorting script
% provides a general template for organizing data into the format used for
% fitting CASANDRE.

% Citation: WT Adler, WJ Ma. (2018). Comparing Bayesian and non-Bayesian 
% accounts of human confidence reports. PLOS Computational Biology. 14(11):
% e1006572. https://doi.org/10.1371/journal. pcbi.1006572

% Stimulus and task:
% During each session, each subject completed two orientation 
% categorization tasks, Tasks A and B. The stimuli were drifting Gabors for
% some subjects and ellipses for other subjects. On each trial, a category 
% C was selected randomly (both categories were equally probable), and a 
% stimulus s was drawn from the corresponding stimulus distribution and 
% displayed. The subject categorized the stimulus and simultaneously 
% reported their confidence on a 4-point scale, with a single button press. 
% The categories were defined by normal distributions on orientation, which
% differed by task. In Task A, the distributions had different means (±mu_C)
% and the same standard deviation (sigma_C); leftward-tilting stimuli were more 
% likely to be from category 1. In Task B, the distributions had the same 
% mean (0°) and different standard deviations (sigma_1, sigma_2); stimuli around 
% the horizontal were more likely to be from category 1.

% Please see these links for more information about these data: 
% https://github.com/wtadler/confidence/
% https://osf.io/s46pr/


% Start with clean slate
clearvars

% Set Paths
thisPath    = fullfile(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(thisPath,'..')));

% Load data
cd(thisPath)
load('Adler_2018_Expt1-3.mat');

%% Preprocess data
% Step 0: native structure
% Column 1 = subject index (1 to 11)
% Column 2 = stimulus category (1 or 2)
% Column 3 = response category (1 or 2)
% Column 4 = confidence level (1 to 4)
% Column 5 = RT
% Column 6 = stimulus uncertainty (higher is more difficult)
% Column 7 = stimulus orientation in visual degrees
% Column 8 = task (A or B)

dataRaw = table2cell(Expt1);
dataMat = cell2mat(dataRaw(:,1:7));

% Step 1: convert to one 2-D matrix per observer (trial x [stimulus strength, 2-D choice]) 
obsInd  = dataMat(:,1);
obsList = unique(obsInd);
relList = unique(dataMat(:,6)); % 6 levels of stimulus reliability

for iO = 1:numel(obsList)   % n observers
    for iT = [1 2]          % 2 tasks    
        if iT == 1
            taskInd = cell2mat(dataRaw(:,8)) == 'A';
        else
            taskInd = cell2mat(dataRaw(:,8)) == 'B';
        end
            
        for iR = 1:numel(relList)
            trialList = find(obsInd == obsList(iO) & dataMat(:,6) == relList(iR) & taskInd == 1);
            respInd   = round(dataMat(trialList, 3) - 1.5);

            trials{iO}{iT}{iR}(:,1) = dataMat(trialList, 2);
            trials{iO}{iT}{iR}(:,2) = dataMat(trialList, 7);
            trials{iO}{iT}{iR}(:,3) = dataMat(trialList, 4).*respInd;
        end
    end
end


% Save sorted data
saveName = 'Adler_2018_Expt1_sorted';
save(saveName, 'trials');



