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
EpsTries = 10.^(-3:0.1:-0.1);

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
            
            % This matrix characterizes how ProbMat will evolve due to the passive loss of the loop
            PassiveLossMatrix = lossmatrix_nologbinom ( Nmax,1-eta );
            
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
            
            DetectionLossMatrixClickWithPassiveLossArray = zeros([Nmax+1,Nmax+1,length(EpsTries)]);
            DetectionLossMatrixNoClickWithPassiveLossArray = zeros([Nmax+1,Nmax+1,length(EpsTries)]);
            
            for EpsIndex = 1:length(EpsTries)
                epsTemp = EpsTries(EpsIndex);
                DetectionLossMatrixTemp = lossmatrix_nologbinom ( Nmax,epsTemp );
                
                %Transition matrices conditional on detection (or absence thereof) using
                %the fact that P(A and B) = P(A) * P(B assuming A)
                DetectionLossMatrixClickTemp = DetectionLossMatrixTemp .* ClickProbMatrix;
                DetectionLossMatrixNoClickTemp = DetectionLossMatrixTemp .* (1 - ClickProbMatrix);
                
                DetectionLossMatrixClickWithPassiveLossArray(:,:,EpsIndex) = DetectionLossMatrixClickTemp*PassiveLossMatrix;
                DetectionLossMatrixNoClickWithPassiveLossArray(:,:,EpsIndex) = DetectionLossMatrixNoClickTemp*PassiveLossMatrix;
            end
            
            DetectionClickArrayPermute = permute(DetectionLossMatrixClickWithPassiveLossArray,[3 1 2]);
            DetectionNoClickArrayPermute = permute(DetectionLossMatrixNoClickWithPassiveLossArray,[3 1 2]);
            
            for N0Index = 1:length(N0Array)
                N0 = N0Array(N0Index);
                N0
                
                % Beginning of Adaptive Approach
                adaptivesuffix = sprintf('AdaptiveN%deta%dnu%dgamma%d.csv',N0Index,etaIndex,nuIndex,gammaIndex);
                fnameprob = strcat(dataFolder, 'Prob', adaptivesuffix);
                fnameclick = strcat(dataFolder, 'Click', adaptivesuffix);
                fnameeps = strcat(dataFolder, 'Eps', adaptivesuffix);
                fnamerounds = strcat(dataFolder, 'Rounds', adaptivesuffix);
                fnamecurrent = strcat(dataFolder, 'CurrentPhotons', adaptivesuffix);
                fnamecurrentest = strcat(dataFolder, 'CurrentPhotonsEst', adaptivesuffix);
                fnamecurrentdensity = strcat(dataFolder, 'CurrentPhotonsDensity', adaptivesuffix);
                
                ProbEstimatesVsTrialAdpative = zeros(Trials, Nmax+1);
                
                ClickRecord = zeros(Trials, Rounds);  %This is the click record of the detector
                EpsRecord = zeros(Trials, Rounds);
                SimulatedRounds = Rounds*ones(Trials,1);
                CurrentPhotons = zeros(Trials, Rounds);
                CurrentPhotonsEst = zeros(Trials, Rounds);
                CurrentPhotonsDensity = zeros(Trials, Rounds);
                
                for n = 1:Trials
                    
                    NumMarbles = zeros(2*Rounds,1); %In each round, there are two steps (loop loss and detection). This is how many photons are left at each sub-step
                    NumRemoved = zeros(2*Rounds,1); %In each round, there are two steps (loop loss and detection). This is how many photons are removed at each sub-step
                    
                    %columns are the original photon number, rows are the current photon
                    %number. This will get update to match what the player believes to be the
                    %probabilities of different possible current photon numbers, given the
                    %experimental record, and assuming each possible initial photon number.
                    ProbInitial = ones(1,Nmax+1)/(Nmax+1);
                    ProbInitial = ProbInitial/sum(ProbInitial);
                    
                    ProbMat = diag(ProbInitial); 
                    
                    AvailableInformation = zeros(Rounds,1);
                    KnownInformation = zeros(Rounds,1);
                    
                    IAvailable = zeros(length(EpsTries),length(CheckEpsilon));
                    IGain = zeros(length(EpsTries),length(CheckEpsilon));
                    InformationRatio = zeros(length(EpsTries),length(CheckEpsilon));
                    OptimalEps = zeros(length(CheckEpsilon),1);
                    CheckIndex = 0;
                    
                    EndOfSim = false;
                    
                    CurrentNum = N0; %current number of photons in the loop
                    
                    % Bayesian estimate of the probability distribution at each step
                    ProbEstimatesVsRoundAdaptive = zeros(Rounds,Nmax+1);
                    
                    for j = 1:Rounds
                        if EndOfSim == false
                            PN = sum(ProbMat,2);
                            Divergence = ProbMat.*(log(ProbMat) - log(PN) - log(ProbInitial));
                            Divergence(isnan(Divergence)) = 0;
                            AvailableInformation(j) = sum(sum(Divergence)); 
                            PN0 = sum(ProbMat,1);
                            Divergence = ProbMat.*(log(PN0) - log(ProbInitial));
                            Divergence(isnan(Divergence)) = 0;
                            KnownInformation(j) = sum(sum(Divergence));
                            
                            PreviousAvailableInformation = AvailableInformation(j);
                            PreviousKnownInformation = KnownInformation(j);
                            
                            
                            PreviousEstimate = PN0;
                            
                            
                            if ismember(j,CheckEpsilon)
                                CheckIndex = CheckIndex + 1;
                                IAvailable = zeros(length(EpsTries),1);
                                IGain = zeros(length(EpsTries),1);
                                for EpsIndex = 1:length(EpsTries)
                                    epsTemp = EpsTries(EpsIndex);
                                    
                                    ProbMatClick = DetectionLossMatrixClickWithPassiveLossArray(:,:,EpsIndex)*ProbMat;
                                    ProbMatNoClick = DetectionLossMatrixNoClickWithPassiveLossArray(:,:,EpsIndex)*ProbMat;
                                    
                                    ProbClick = sum(sum(ProbMatClick));
                                    ProbNoClick = sum(sum(ProbMatNoClick));
                                    
                                    ProbMatClick = ProbMatClick/ProbClick;
                                    ProbMatNoClick = ProbMatNoClick/ProbNoClick;
                                    
                                    PNClick = sum(ProbMatClick,2);
                                    DivergenceClick = ProbMatClick.*(log(ProbMatClick) - log(PNClick) - log(ProbInitial));
                                    DivergenceClick(isnan(DivergenceClick)) = 0;
                                    IClick = sum(sum(DivergenceClick));
                                    PN0Click = sum(ProbMatClick,1);
                                    DivergenceClick = ProbMatClick.*(log(PN0Click) - log(ProbInitial));
                                    DivergenceClick(isnan(DivergenceClick)) = 0;
                                    IGainClick = sum(sum(DivergenceClick));
                                    
                                    PNNoClick = sum(ProbMatNoClick,2);
                                    DivergenceNoClick = ProbMatNoClick.*(log(ProbMatNoClick) - log(PNNoClick) - log(ProbInitial));
                                    DivergenceNoClick(isnan(DivergenceNoClick)) = 0;
                                    INoClick = sum(sum(DivergenceNoClick)); 
                                    PN0NoClick = sum(ProbMatNoClick,1);
                                    DivergenceNoClick = ProbMatNoClick.*(log(PN0NoClick) - log(ProbInitial));
                                    DivergenceNoClick(isnan(DivergenceNoClick)) = 0;
                                    IGainNoClick = sum(sum(DivergenceNoClick));
                                    
                                    IAvailable(EpsIndex) = IClick*ProbClick + INoClick*ProbNoClick;
                                    IGain(EpsIndex) = IGainClick*ProbClick + IGainNoClick*ProbNoClick;
                                    
                                end
                                
                                InformationRatio = (IGain-PreviousKnownInformation)./(IAvailable-PreviousAvailableInformation);
                                FinerEpsilon = 10.^(-3:0.01:0);
                                FinerInformationRatio = spline(EpsTries,InformationRatio,FinerEpsilon);
                                [V,I] = min(FinerInformationRatio);
                                OptimalEps = FinerEpsilon(I);
                                
                                if true
                                    eps = OptimalEps;
                                    
                                    if (eps >= EpsTries(length(EpsTries)-1)-0.01)
                                        EndOfSim = true;
                                        SimulatedRounds(n) = CheckIndex;
                                        eps = EpsTries(length(EpsTries)-1)-0.01;
                                    end
                                    
                                    DetectionLossMatrixClickWithPassiveLoss = interp1(EpsTries, DetectionClickArrayPermute, eps);
                                    DetectionLossMatrixNoClickWithPassiveLoss = interp1(EpsTries, DetectionNoClickArrayPermute, eps);
                                end
                                
                            end
                            
                            EpsRecord(n,j) = eps;
                            
                            CurrentPhotons(n,j) = CurrentNum;
                            CurrentPhotonsEst(n,j) = sum((0:Nmax).*sum(ProbMat,2)');
                            CurrentPhotonsDensity(n,j) = sum((0:Nmax).*ProbMat(:,N0+1)')/sum(ProbMat(:,N0+1));
                            
                            NumRemoved(2*j-1) = sum(rand(CurrentNum,1)<1 - eta); %Loss due to loop efficiency
                            NumMarbles(2*j-1) = CurrentNum - NumRemoved(2*j-1);
                            CurrentNum = NumMarbles(2*j-1);
                            NumRemoved(2*j) = sum(rand(CurrentNum,1)<eps); %Loss due to detection
                            NumMarbles(2*j) = CurrentNum - NumRemoved(2*j);
                            CurrentNum = NumMarbles(2*j);
                            PClick = nu + (1-nu) * (1 - (1-gamma)^NumRemoved(2*j)); %Probability of the detector clicking with NumRemoved(2*j) photons hitting it
                            ClickRecord(n,j) = (rand<PClick); %Generate a click with probability PClick
                            
                            %ProbMat = PassiveLossMatrix*ProbMat; %Account for the loss due to the loop. 
                            ThisClick = ClickRecord(n,j);
                            
                            %Perform Bayesian update
                            if ThisClick
                                ProbMat = reshape(DetectionLossMatrixClickWithPassiveLoss, Nmax+1, Nmax+1)*ProbMat;
                            else
                                ProbMat = reshape(DetectionLossMatrixNoClickWithPassiveLoss, Nmax+1, Nmax+1)*ProbMat;
                            end
                            
                            ProbMat = ProbMat/sum(sum(ProbMat));
                            
                            ProbEstimatesVsRoundAdaptive(j,:) = sum(ProbMat,1);
                        end
                    end
                    
                    ProbEstimatesVsTrialAdpative(n,:) = ProbEstimatesVsRoundAdaptive(SimulatedRounds(n),:);
                    
                end
                
                csvwrite(fnameprob, ProbEstimatesVsTrialAdpative);
                csvwrite(fnameclick, ClickRecord);
                csvwrite(fnameeps, EpsRecord);
                csvwrite(fnamerounds, SimulatedRounds);
                csvwrite(fnamecurrent, CurrentPhotons);
                csvwrite(fnamecurrentest, CurrentPhotonsEst);
                csvwrite(fnamecurrentdensity, CurrentPhotonsDensity);
                
            end
        end    
    end    
end




