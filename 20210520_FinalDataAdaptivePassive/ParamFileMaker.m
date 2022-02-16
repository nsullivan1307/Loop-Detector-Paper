% 2021-05-22
% This file makes csv data files of simulations for both passive and adaptive 
% strategies for different values of eta, the loop efficiency.

paramFolder = 'Params/';
dataFolder = 'Data/';
figureFolder = 'Figures/';

paramStatus = mkdir(paramFolder);
dataStatus = mkdir(dataFolder);
figureStatus = mkdir(figureFolder);

%Nmax = 50; %maximum number of photons that might be sent into the system (assumption)

%N0Array = [1, 2, 3, 5, 10, 15, 20, 25, 30, 35, 40, 45, 49, 50]; %Actual starting photon number

Nmax = 100; %maximum number of photons that might be sent into the system (assumption)

N0Array = [1, 2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 95, 100]; %Actual starting photon number

Rounds = 400; %How many round-trips to simulate

% Outcoupling Rate
epsArray = 10.^(-2:0.1:-1);

% Loop Efficiency
etaArray = 1 - 10.^(-3:1:-1);
%etaArray = [0.99];

% Dark Counts
nuArray = [1e-6]; 

% Detector Efficiency
gammaArray = [0.90];

Trials = 1000; % How many full trials to simulate

csvwrite(strcat(paramFolder, 'Params.csv'), [Nmax, Trials, Rounds]);
csvwrite(strcat(paramFolder, 'epsArray.csv'), epsArray);
csvwrite(strcat(paramFolder, 'N0Array.csv'), N0Array);
csvwrite(strcat(paramFolder, 'etaArray.csv'), etaArray);
csvwrite(strcat(paramFolder, 'nuArray.csv'), nuArray);
csvwrite(strcat(paramFolder, 'gammaArray.csv'), gammaArray);

