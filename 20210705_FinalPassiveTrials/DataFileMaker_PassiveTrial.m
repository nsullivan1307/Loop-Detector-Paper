% 2021-05-20
% This file makes csv data files of simulations for both passive and adaptive 
% strategies for different values of eta, the loop efficiency.

paramFolder = 'Params/';
dataFolder = 'Data/';

Params = readmatrix(strcat(paramFolder, 'Params.csv'));
N0Array = readmatrix(strcat(paramFolder, 'N0Array.csv'));
epsArray = readmatrix(strcat(paramFolder, 'epsArray.csv'));
etaArray = readmatrix(strcat(paramFolder, 'etaArray.csv'));
nuArray = readmatrix(strcat(paramFolder, 'nuArray.csv'));
gammaArray = readmatrix(strcat(paramFolder, 'gammaArray.csv'));

Nmax = Params(1);
Trials = Params(2);
Rounds = Params(3);

LN0 = length(N0Array);
Leps = length(epsArray);
Leta = length(etaArray);
Lnu = length(nuArray);
Lgamma = length(gammaArray);

CheckEpsilon = 1:Rounds;
EpsTries = 10.^(-3:0.1:0);

for etaIndex = 1:length(etaArray)
    eta = etaArray(etaIndex);
    eta %#ok<*NOPTS>
    % This matrix characterizes how ProbMat will evolve due to the passive loss of the loop
    PassiveLossMatrix = lossmatrix_nologbinom ( Nmax,1-eta );
    
    for nuIndex = 1:length(nuArray)
        nu = nuArray(nuIndex);
        nu
        
        for gammaIndex = 1:length(gammaArray)
            gamma = gammaArray(gammaIndex);
            gamma
            
            % The probability of a certain combination of photons before and after a
            % given detector outcoupling step leading to a detection event
            ClickProbMatrix = zeros(Nmax+1);
            
            for j = 1:Nmax+1
                for k = 1:Nmax+1
                    nloss = k-j;
                    if nloss >= 0
                        PClick = nu + (1-nu) * (1 - (1-gamma)^nloss);
                        ClickProbMatrix(j,k) = PClick;
                    end
                end
            end
            
            
            for N0Index = 1:length(N0Array)
                N0 = N0Array(N0Index);
                N0
                
                % Beginning of Passive Approach
                for epsIndex = 1:length(epsArray)
                    eps = epsArray(epsIndex);
                    eps
                    
                    passivetrialsuffix = sprintf('PassiveTrialN%deps%deta%dnu%dgamma%d.csv',N0Index,epsIndex,etaIndex,nuIndex,gammaIndex);
                    fnamestats = strcat(dataFolder, 'Stats', passivetrialsuffix);
                    fnameinit = strcat(dataFolder, 'Init', passivetrialsuffix);
                    fnameloop = strcat(dataFolder, 'Loop', passivetrialsuffix);
                    
                    % Likelihood matrix reconstruction based on the measurement record in
                    % ClickRecord.
                    
                    % These are the transition matrices that characterize how the matrix
                    % ProbMat will evolve due to loss/outcoupling from detection (this is
                    % independent of the click record)
                    DetectionLossMatrix = lossmatrix_nologbinom ( Nmax,eps );
                    
                    %Transition matrices conditional on detection (or absence thereof) using
                    %the fact that P(A and B) = P(A) * P(B assuming A)
                    DetectionLossMatrixClick = DetectionLossMatrix .* ClickProbMatrix;
                    DetectionLossMatrixNoClick = DetectionLossMatrix .* (1 - ClickProbMatrix);
                    
                    DetectionLossMatrixClickWithPassiveLoss = DetectionLossMatrixClick*PassiveLossMatrix;
                    DetectionLossMatrixNoClickWithPassiveLoss = DetectionLossMatrixNoClick*PassiveLossMatrix;
                    
                    ProbEstimatesVsTrialPassive = zeros(Trials, Nmax+1);
                    
                    
                    % STATS:
                    % 1 : Click Record
                    % 2 : Current Photons in Loop
                    % 3 : Estimate of Photons in Loop (Mean)
                    % 4 : Estimate of Photons in Loop (MLE)
                    % 5 : Estimate of Photons in Loop (Density Matrix)
                    % 6 : Estimate of Initial Photons (Mean)
                    % 7 : Estimate of Initial Photons (MLE)
                    
                    ClickRecord = zeros(Rounds,1);  %This is the click record of the detector
                    Loop = zeros(Rounds,1);
                    LoopEstMean = zeros(Rounds,1);
                    LoopEstMLE = zeros(Rounds,1);
                    LoopEstDensity = zeros(Rounds,1);
                    InitEstMean = zeros(Rounds,1);
                    InitEstMLE = zeros(Rounds,1);
                    
                    % START OF SIMULATION

                    NumMarbles = zeros(2*Rounds,1); %In each round, there are two steps (loop loss and detection). This is how many photons are left at each sub-step
                    NumRemoved = zeros(2*Rounds,1); %In each round, there are two steps (loop loss and detection). This is how many photons are removed at each sub-step

                    %columns are the original photon number, rows are the current photon
                    %number. This will get update to match what the player believes to be the
                    %probabilities of different possible current photon numbers, given the
                    %experimental record, and assuming each possible initial photon number.
                    ProbInitial = ones(1,Nmax+1)/(Nmax+1);
                    ProbInitial = ProbInitial/sum(ProbInitial);

                    ProbMat = diag(ProbInitial); 

                    CurrentNum = N0; %current number of photons in the loop
                    
                    AvailableInformation = zeros(Rounds,1);
                    KnownInformation = zeros(Rounds,1);

                    % Bayesian estimate of the probability distribution at each step
                    InitEstimatesVsRoundPassive = zeros(Rounds,Nmax+1);
                    LoopEstimatesVsRoundPassive = zeros(Rounds,Nmax+1);

                    for j = 1:Rounds    
                        Loop(j) = CurrentNum;
                        LoopEstMean(j) = sum((0:Nmax).*sum(ProbMat,2)');
                        [~,index] = max(sum(ProbMat,2));
                        LoopEstMLE(j) = index - 1;
                        LoopEstDensity(j) = sum((0:Nmax).*ProbMat(:,N0+1)')/sum(ProbMat(:,N0+1));
                        InitEstMean(j) = sum((0:Nmax).*sum(ProbMat,1));
                        [~,index] = max(sum(ProbMat,1));
                        InitEstMLE(j) = index - 1;
                        
                        PN = sum(ProbMat,2);
                        Divergence = ProbMat.*(log(ProbMat) - log(PN) - log(ProbInitial));
                        Divergence(isnan(Divergence)) = 0;
                        AvailableInformation(j) = sum(sum(Divergence)); 
                        PN0 = sum(ProbMat,1);
                        Divergence = ProbMat.*(log(PN0) - log(ProbInitial));
                        Divergence(isnan(Divergence)) = 0;
                        KnownInformation(j) = sum(sum(Divergence));

                        NumRemoved(2*j-1) = sum(rand(CurrentNum,1)<1 - eta); %Loss due to loop efficiency
                        NumMarbles(2*j-1) = CurrentNum - NumRemoved(2*j-1);
                        CurrentNum = NumMarbles(2*j-1);
                        NumRemoved(2*j) = sum(rand(CurrentNum,1)<eps); %Loss due to detection
                        NumMarbles(2*j) = CurrentNum - NumRemoved(2*j);
                        CurrentNum = NumMarbles(2*j);
                        PClick = nu + (1-nu) * (1 - (1-gamma)^NumRemoved(2*j)); %Probability of the detector clicking with NumRemoved(2*j) photons hitting it
                        ClickRecord(j) = (rand<PClick); %Generate a click with probability PClick

                        ThisClick = ClickRecord(j);

                        %Perform Bayesian update
                        if ThisClick
                            ProbMat = DetectionLossMatrixClickWithPassiveLoss*ProbMat;
                        else
                            ProbMat = DetectionLossMatrixNoClickWithPassiveLoss*ProbMat;
                        end

                        ProbMat = ProbMat/sum(sum(ProbMat));

                        InitEstimatesVsRoundPassive(j,:) = sum(ProbMat,1);
                        LoopEstimatesVsRoundPassive(j,:) = sum(ProbMat,2)';
                    end
                        
                    % END OF SIMULATION
                    
                    TrialStats = [
                        ClickRecord'
                        Loop'
                        LoopEstMean'
                        LoopEstMLE'
                        LoopEstDensity'
                        InitEstMean'
                        InitEstMLE'
                        AvailableInformation'
                        KnownInformation']';
                        
                    writematrix(TrialStats, fnamestats);
                    writematrix(InitEstimatesVsRoundPassive, fnameinit);
                    writematrix(LoopEstimatesVsRoundPassive, fnameloop);
                end    
            end    
        end    
    end    
end




