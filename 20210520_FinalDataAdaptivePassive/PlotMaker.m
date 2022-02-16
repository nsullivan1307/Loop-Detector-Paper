% 2021-05-20
% This file makes plots comparing passive and adaptive strategies for different
% values of eta, the loop efficiency.

paramFolder = 'Params/';
dataFolder = 'Data/';
figureFolder = 'Figures/';

Params = csvread(strcat(paramFolder, 'Params.csv'));
N0Array = csvread(strcat(paramFolder, 'N0Array.csv'));
epsArray = csvread(strcat(paramFolder, 'epsArray.csv'));
etaArray = csvread(strcat(paramFolder, 'etaArray.csv'));
nuArray = csvread(strcat(paramFolder, 'nuArray.csv'));
gammaArray = csvread(strcat(paramFolder, 'gammaArray.csv'));

Nmax = Params(1);
Trials = Params(2);
Rounds = Params(3);

LN0 = length(N0Array);
Leps = length(epsArray);
Leta = length(etaArray);
Lnu = length(nuArray);
Lgamma = length(gammaArray);

%%
% Loading stats data into arrays

% 1: MS Error
% 2: Mean of Estimate Variance
% 3: Estimate Mean
% 4: Standard Deviation of Estimate Mean
% 5: Maximum Likelihood Estimate
% 6: MLE MS Error
% 7: Number of Rounds (before estimated photons < 0.5)

VarianceStatsAdaptive = zeros([LN0, Leta, Lnu, Lgamma, 6]);
VarianceStatsPassive = zeros([LN0, Leps, Leta, Lnu, Lgamma, 6]);

for etaIndex = 1:length(etaArray)
    eta = etaArray(etaIndex);
    eta %#ok<*NOPTS>
    for nuIndex = 1:length(nuArray)
        nu = nuArray(nuIndex);
        nu
        for gammaIndex = 1:length(gammaArray)
            gamma = gammaArray(gammaIndex);
            gamma
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
                fnamecurrentmle = strcat(dataFolder, 'CurrentPhotonsMLE', adaptivesuffix);
                fnamecurrentdensity = strcat(dataFolder, 'CurrentPhotonsDensity', adaptivesuffix);
                
                ProbEstimatesVsTrial = csvread(fnameprob);
                EpsVsTrial = csvread(fnameeps);
                RoundNumberVsTrial = csvread(fnamerounds);
                ClickRecord = csvread(fnameclick);
                CurrentPhotons = csvread(fnamecurrent);
                CurrentPhotonsEst = csvread(fnamecurrentest);
                CurrentPhotonsMLE = csvread(fnamecurrentmle);
                CurrentPhotonsDensity = csvread(fnamecurrentdensity);
                
                StatEstimatesVsTrial = zeros(Trials, 5);
                
                Count = ones(Trials,3);
                for j = 1:Trials
                    StatEstimatesVsTrial(j,1) = sum((0:1:Nmax).*ProbEstimatesVsTrial(j,:)); % Mean estimate
                    StatEstimatesVsTrial(j,2) = sum((0:1:Nmax).^2.*ProbEstimatesVsTrial(j,:)) - StatEstimatesVsTrial(j,1)^2; % Variance of estimate
                    StatEstimatesVsTrial(j,3) = sum(((0:1:Nmax) < N0).*ProbEstimatesVsTrial(j,:)) + 0.5*ProbEstimatesVsTrial(j,N0+1);
                    [~,index] = max(ProbEstimatesVsTrial(j,:));
                    StatEstimatesVsTrial(j,4) = index + 1;
                    StatEstimatesVsTrial(j,5) = sum(CurrentPhotonsEst(j,:) > 0.5);
                end
                
                StatEstimatesVsTrial = StatEstimatesVsTrial(~isnan(StatEstimatesVsTrial(:,3)),:);
                Total = length(StatEstimatesVsTrial(:,1));
                Mean = sum(StatEstimatesVsTrial(:,1))/Total;
                MSErrorMean = sum((StatEstimatesVsTrial(:,1) - N0).^2)/Total;
                VarMean = sum((StatEstimatesVsTrial(:,1) - Mean).^2)/Total;
                SDMean = sqrt(VarMean);
                MeanVar = sum(StatEstimatesVsTrial(:,2))/Total;
                MLE = StatEstimatesVsTrial(j,4);
                MSErrorMLE = sum((StatEstimatesVsTrial(:,4) - N0).^2)/Total;
                UsefulRounds = sum(StatEstimatesVsTrial(:,5))/Total;
                
                VarianceStatsAdaptive(N0Index, etaIndex, nuIndex, gammaIndex, 1) = MSErrorMean;
                VarianceStatsAdaptive(N0Index, etaIndex, nuIndex, gammaIndex, 2) = MeanVar;
                VarianceStatsAdaptive(N0Index, etaIndex, nuIndex, gammaIndex, 3) = Mean;
                VarianceStatsAdaptive(N0Index, etaIndex, nuIndex, gammaIndex, 4) = SDMean;
                VarianceStatsAdaptive(N0Index, etaIndex, nuIndex, gammaIndex, 5) = MLE;
                VarianceStatsAdaptive(N0Index, etaIndex, nuIndex, gammaIndex, 6) = MSErrorMLE;
                VarianceStatsAdaptive(N0Index, etaIndex, nuIndex, gammaIndex, 7) = UsefulRounds;
                
                
                % Beginning of Passive Approach
                for epsIndex = 1:length(epsArray)
                    eps = epsArray(epsIndex);
                    eps
                    
                    passivesuffix = sprintf('PassiveN%deps%deta%dnu%dgamma%d.csv',N0Index,epsIndex,etaIndex,nuIndex,gammaIndex);
                    fnameprob = strcat(dataFolder, 'Prob', passivesuffix);
                    fnameclick = strcat(dataFolder, 'Click', passivesuffix);
                    fnamecurrent = strcat(dataFolder, 'CurrentPhotons', passivesuffix);
                    fnamecurrentest = strcat(dataFolder, 'CurrentPhotonsEst', passivesuffix);
                    fnamecurrentmle = strcat(dataFolder, 'CurrentPhotonsMLE', passivesuffix);
                    fnamecurrentdensity = strcat(dataFolder, 'CurrentPhotonsDensity', passivesuffix);
                    
                    CurrentPhotonsEst = csvread(fnamecurrentest);
                    
                    ProbEstimatesVsTrial = csvread(fnameprob);
                    
                    StatEstimatesVsTrial = zeros(Trials, 5);
                    
                    Count = ones(Trials,3);
                    for j = 1:Trials
                        StatEstimatesVsTrial(j,1) = sum((0:1:Nmax).*ProbEstimatesVsTrial(j,:)); % Mean estimate
                        StatEstimatesVsTrial(j,2) = sum((0:1:Nmax).^2.*ProbEstimatesVsTrial(j,:)) - StatEstimatesVsTrial(j,1)^2; % Variance of estimate
                        StatEstimatesVsTrial(j,3) = sum(((0:1:Nmax) < N0).*ProbEstimatesVsTrial(j,:)) + 0.5*ProbEstimatesVsTrial(j,N0+1);
                        [~,index] = max(ProbEstimatesVsTrial(j,:));
                        StatEstimatesVsTrial(j,4) = index + 1;
                        StatEstimatesVsTrial(j,5) = sum(CurrentPhotonsEst(j,:) > 0.5);
                    end
                    
                    
                    StatEstimatesVsTrial = StatEstimatesVsTrial(~isnan(StatEstimatesVsTrial(:,3)),:);
                    Total = length(StatEstimatesVsTrial(:,1));
                    Mean = sum(StatEstimatesVsTrial(:,1))/Total;
                    MSErrorMean = sum((StatEstimatesVsTrial(:,1) - N0).^2)/Total;
                    VarMean = sum((StatEstimatesVsTrial(:,1) - Mean).^2)/Total;
                    SDMean = sqrt(VarMean);
                    MeanVar = sum(StatEstimatesVsTrial(:,2))/Total;
                    MLE = StatEstimatesVsTrial(j,4);
                    MSErrorMLE = sum((StatEstimatesVsTrial(:,4) - N0).^2)/Total;
                    UsefulRounds = sum(StatEstimatesVsTrial(:,5))/Total;
                    
                    
                    
                    VarianceStatsPassive(N0Index, epsIndex,etaIndex, nuIndex, gammaIndex, 1) = MSErrorMean;
                    VarianceStatsPassive(N0Index, epsIndex,etaIndex, nuIndex, gammaIndex, 2) = MeanVar;
                    VarianceStatsPassive(N0Index, epsIndex,etaIndex, nuIndex, gammaIndex, 3) = Mean;
                    VarianceStatsPassive(N0Index, epsIndex,etaIndex, nuIndex, gammaIndex, 4) = SDMean;
                    VarianceStatsPassive(N0Index, epsIndex,etaIndex, nuIndex, gammaIndex, 5) = MLE;
                    VarianceStatsPassive(N0Index, epsIndex,etaIndex, nuIndex, gammaIndex, 6) = MSErrorMLE;
                    VarianceStatsPassive(N0Index, epsIndex,etaIndex, nuIndex, gammaIndex, 7) = UsefulRounds;
                    
                end 
            end
        end
    end
end

%%
% Individual trial stats

CURRENT_PHOTONS_EPSILON_VS_ROUND = 1;
CURRENT_PHOTONS_VS_EPSILON = 2;
CURRENT_PHOTONS_VS_ROUND = 3;

FigureOption = CURRENT_PHOTONS_EPSILON_VS_ROUND;

for etaIndex = 1:length(etaArray)
    eta = etaArray(etaIndex);
    eta %#ok<*NOPTS>
    for nuIndex = 1:length(nuArray)
        nu = nuArray(nuIndex);
        nu
        for gammaIndex = 1:length(gammaArray)
            gamma = gammaArray(gammaIndex);
            gamma
            for N0Index = 1:length(N0Array)
                N0 = N0Array(N0Index);
                N0
                if FigureOption == CURRENT_PHOTONS_EPSILON_VS_ROUND
                  FigW = 8;
                  FigH = 9;
                  PosVector = [0,0,500,400];
                  legfont = 8;

                  figure
                  hold on
                  set(gcf,'PaperSize',[FigW FigH],'PaperPosition',[0 0 FigW FigH]);
                  set(gca,'FontSize',12,'Linewidth',1)
                  grid on
                  grid minor
                  set(gca,'YDir','normal')
                  %semilogx(epsArray, VarianceStats(:,1),'-','Linewidth',2,'Color',[1 0 0])
                  %semilogx(epsArray, VarianceStats(:,2),'-','Linewidth',2,'Color',[0 1 0])
                  %semilogx(epsArray, VarianceStats(:,3),'-','Linewidth',2,'Color',[0 0 1])
                  %plot(1:Rounds, sum(CurrentPhotons,1)/(N0*Trials),'-','Linewidth',2,'Color',[1 0 0])
                  %plot(1:Rounds, sum(CurrentPhotonsEst,1)/(N0*Trials),'-','Linewidth',2,'Color',[0 1 0])
                  %plot(1:Rounds, sum(ClickRecord,1)/Trials,'-','Linewidth',2,'Color',[0 0 1])
                  %plot(1:Rounds, sum(EpsVsTrial,1)/Trials, '--','Linewidth',2,'Color',[0.8 0.3 1])
                  plot(1:Rounds, CurrentPhotons(1,:)/N0,'-','Linewidth',2,'Color',[1 0 0])
                  plot(1:Rounds, CurrentPhotonsEst(1,:)/N0,'-','Linewidth',2,'Color',[0 1 0])
                  plot(1:Rounds, CurrentPhotonsDensity(1,:)/N0,'-','Linewidth',2,'Color',[0 1 1])
                  plot(1:Rounds, ClickRecord(1,:),'.','Linewidth',2,'Color',[0 0 1])
                  plot(1:Rounds, EpsVsTrial(1,:).*CurrentPhotonsEst(1,:),'-','Linewidth',2,'Color',[1 0 1])
                  plot(1:Rounds, EpsVsTrial(1,:), '--','Linewidth',2,'Color',[0.8 0.3 1])
                  leg = {'Current Photons','Estimated Current Photons', 'Current Photons from Density Matrix', 'Average Clicks', 'Epsilon *Current Photons', 'Average Epsilon'};
                  legend(leg,'FontSize',legfont,'Location','SouthOutside');
                  xlabel('Scale')
                  ylabel('Round')
                  title(strcat('Current Photons, Epsilon vs Round N=', num2str(N0)))
                  fname = strcat(figureFolder, 'CurrentPhotonsEpsilonvsRoundN', num2str(N0), 'Eta', num2str(etaIndex), 'Nu', num2str(nuIndex), 'Gamma', num2str(gammaIndex));
                  %fname = 'LogLogPlotRMSErrorVariancevsEpsilonEta0.99Nu1e-6Gamma0.9N3';
                  print3(gcf,fname)
                  close
                end
                
                if FigureOption == CURRENT_PHOTONS_VS_EPSILON
                  FigW = 8;
                  FigH = 9;
                  PosVector = [0,0,500,400];
                  legfont = 8;

                  figure
                  hold on
                  set(gcf,'PaperSize',[FigW FigH],'PaperPosition',[0 0 FigW FigH]);
                  set(gca,'FontSize',12,'Linewidth',1)
                  grid on
                  grid minor
                  set(gca,'YDir','normal')
                  %semilogx(epsArray, VarianceStats(:,1),'-','Linewidth',2,'Color',[1 0 0])
                  %semilogx(epsArray, VarianceStats(:,2),'-','Linewidth',2,'Color',[0 1 0])
                  %semilogx(epsArray, VarianceStats(:,3),'-','Linewidth',2,'Color',[0 0 1])
                  %plot(1:Rounds, sum(CurrentPhotons,1)/(N0*Trials),'-','Linewidth',2,'Color',[1 0 0])
                  %plot(1:Rounds, sum(CurrentPhotonsEst,1)/(N0*Trials),'-','Linewidth',2,'Color',[0 1 0])
                  %plot(1:Rounds, sum(ClickRecord,1)/Trials,'-','Linewidth',2,'Color',[0 0 1])
                  %plot(1:Rounds, sum(EpsVsTrial,1)/Trials, '--','Linewidth',2,'Color',[0.8 0.3 1])
                  for j = 1:Trials
                      plot(EpsVsTrial(j,:), CurrentPhotons(j,:),'-','Linewidth',1,'Color',[0 0.7 0.3])
                      plot(EpsVsTrial(j,:), CurrentPhotonsEst(j,:),'-','Linewidth',1,'Color',[0 0.3 0.7])
                      plot(EpsVsTrial(j,:), CurrentPhotonsDensity(j,:),'-','Linewidth',1,'Color',[0.3 0.7 0])
                  end
                  leg = {'Current Photons','Estimated Current Photons', 'Current Photons from Density Matrix'};
                  legend(leg,'FontSize',legfont,'Location','SouthOutside');
                  xlabel('Epsilon')
                  ylabel('Photons')
                  title(strcat('Current Photons Vs Epsilon N=', num2str(N0)))
                  fname = strcat(figureFolder, 'CurrentPhotonsVsEpsilonN', num2str(N0), 'Eta', num2str(etaIndex), 'Nu', num2str(nuIndex), 'Gamma', num2str(gammaIndex));
                  %fname = 'LogLogPlotRMSErrorVariancevsEpsilonEta0.99Nu1e-6Gamma0.9N3';
                  print3(gcf,fname)
                  close
                end 
                
                if FigureOption == CURRENT_PHOTONS_VS_ROUND
                  FigW = 8;
                  FigH = 9;
                  PosVector = [0,0,500,400];
                  legfont = 8;

                  figure
                  hold on
                  set(gcf,'PaperSize',[FigW FigH],'PaperPosition',[0 0 FigW FigH]);
                  set(gca,'FontSize',12,'Linewidth',1)
                  grid on
                  grid minor
                  set(gca,'YDir','normal')
                  %semilogx(epsArray, VarianceStats(:,1),'-','Linewidth',2,'Color',[1 0 0])
                  %semilogx(epsArray, VarianceStats(:,2),'-','Linewidth',2,'Color',[0 1 0])
                  %semilogx(epsArray, VarianceStats(:,3),'-','Linewidth',2,'Color',[0 0 1])
                  %plot(1:Rounds, sum(CurrentPhotons,1)/(N0*Trials),'-','Linewidth',2,'Color',[1 0 0])
                  %plot(1:Rounds, sum(CurrentPhotonsEst,1)/(N0*Trials),'-','Linewidth',2,'Color',[0 1 0])
                  %plot(1:Rounds, sum(ClickRecord,1)/Trials,'-','Linewidth',2,'Color',[0 0 1])
                  %plot(1:Rounds, sum(EpsVsTrial,1)/Trials, '--','Linewidth',2,'Color',[0.8 0.3 1])
                  for j = 1:Trials
                      plot(1:Rounds, CurrentPhotons(j,:),'-','Linewidth',1,'Color',[0 0.7 0.3])
                      plot(1:Rounds, CurrentPhotonsEst(j,:),'-','Linewidth',1,'Color',[0 0.3 0.7])
                      plot(1:Rounds, CurrentPhotonsDensity(j,:),'-','Linewidth',1,'Color',[0.3 0.7 0])
                  end
                  leg = {'Current Photons','Estimated Current Photons', 'Current Photons from Density Matrix'};
                  legend(leg,'FontSize',legfont,'Location','SouthOutside');
                  xlabel('Rounds')
                  ylabel('Photons')
                  title(strcat('Current Photons Vs Round N=', num2str(N0)))
                  fname = strcat(figureFolder, 'CurrentPhotonsVsRoundN', num2str(N0), 'Eta', num2str(etaIndex), 'Nu', num2str(nuIndex), 'Gamma', num2str(gammaIndex));
                  %fname = 'LogLogPlotRMSErrorVariancevsEpsilonEta0.99Nu1e-6Gamma0.9N3';
                  print3(gcf,fname)
                  close
                end  
            end
        end    
    end    
end

%%
NPick = 20;

for etaIndex = 1:length(etaArray)
    eta = etaArray(etaIndex);
    eta %#ok<*NOPTS>
    for nuIndex = 1:length(nuArray)
        nu = nuArray(nuIndex);
        nu
        for gammaIndex = 1:length(gammaArray)
            gamma = gammaArray(gammaIndex);
            gamma
            for N0Index = 1:length(N0Array)
                N0 = N0Array(N0Index);
                N0
                if N0 == NPick
                  FigW = 8;
                  FigH = 9;
                  PosVector = [0,0,500,400];
                  legfont = 8;

                  figure
                  hold on
                  set(gcf,'PaperSize',[FigW FigH],'PaperPosition',[0 0 FigW FigH]);
                  set(gca,'FontSize',12,'Linewidth',1)
                  grid on
                  grid minor
                  set(gca,'YDir','normal')
                  %semilogx(epsArray, VarianceStats(:,1),'-','Linewidth',2,'Color',[1 0 0])
                  %semilogx(epsArray, VarianceStats(:,2),'-','Linewidth',2,'Color',[0 1 0])
                  %semilogx(epsArray, VarianceStats(:,3),'-','Linewidth',2,'Color',[0 0 1])
                  %plot(1:Rounds, sum(CurrentPhotons,1)/(N0*Trials),'-','Linewidth',2,'Color',[1 0 0])
                  %plot(1:Rounds, sum(CurrentPhotonsEst,1)/(N0*Trials),'-','Linewidth',2,'Color',[0 1 0])
                  %plot(1:Rounds, sum(ClickRecord,1)/Trials,'-','Linewidth',2,'Color',[0 0 1])
                  %plot(1:Rounds, sum(EpsVsTrial,1)/Trials, '--','Linewidth',2,'Color',[0.8 0.3 1])
                  
                  %plot(1:Rounds, CurrentPhotons(1,:)/N0,'-','Linewidth',2,'Color',[1 0 0])
                  %plot(1:Rounds, CurrentPhotonsEst(1,:)/N0,'-','Linewidth',2,'Color',[0 1 0])
                  %plot(1:Rounds, ClickRecord(1,:),'.','Linewidth',2,'Color',[0 0 1])
                  plot(1:Rounds, 1 - (1-nu)*(1-EpsVsTrial(1,:)*gamma).^CurrentPhotonsEst(1,:),'-','Linewidth',2,'Color',[0 0.5 0.5])
                  plot(1:Rounds, EpsVsTrial(1,:), '-','Linewidth',2,'Color',[0.8 0.3 1])
                  %leg = {'Current Photons','Estimated Current Photons', 'Average Clicks', 'Click Probability', 'Average Epsilon'};
                  leg = {'Click Probability', 'Epsilon'};
                  legend(leg,'FontSize',legfont,'Location','SouthOutside');
                  xlabel('Round')
                  ylabel('Click Probability')
                  title(strcat('Current Photons, Epsilon vs Round N=', num2str(N0), ' Eta=', num2str(eta)))
                  fname = strcat(figureFolder, 'ClickProbvsRoundN', num2str(N0), sprintf('eta%dnu%dgamma%d.csv',etaIndex,nuIndex,gammaIndex));
                  %fname = 'LogLogPlotRMSErrorVariancevsEpsilonEta0.99Nu1e-6Gamma0.9N3';
                  print3(gcf,fname)
                  close
                end
            end
        end
    end
end

%%
% Stats vs. N

MS_ERROR_MEAN_ALL_PASSIVE = 1;
MS_ERROR_MEAN_OPTIMAL_PASSIVE = 2;
MS_ERROR_MLE_ALL_PASSIVE = 3;
MS_ERROR_MLE_OPTIMAL_PASSIVE = 4;
MS_ERROR_MEAN_ALL_PASSIVE_LL = 5;
MS_ERROR_MEAN_OPTIMAL_PASSIVE_LL = 6;
USEFUL_ROUNDS = 7

FigureOption = USEFUL_ROUNDS;

epsDisplay = [1 4 8 11];
epsArrayDisplay = epsArray(epsDisplay);
N0Display = [1 2 3 4 5 7 10 12 14];

for etaIndex = 1:length(etaArray)
    eta = etaArray(etaIndex);
    eta %#ok<*NOPTS>
    for nuIndex = 1:length(nuArray)
        nu = nuArray(nuIndex);
        nu
        for gammaIndex = 1:length(gammaArray)
            gamma = gammaArray(gammaIndex);
            gamma
            suffix = sprintf('Nm%deta%dnu%dgamma%d',Nmax,etaIndex,nuIndex,gammaIndex)
            
            if FigureOption == MS_ERROR_MEAN_ALL_PASSIVE
                FigW = 8;
                FigH = 9;
                PosVector = [0,0,500,400];
                legfont = 8;

                figure
                hold on
                set(gcf,'PaperSize',[FigW FigH],'PaperPosition',[0 0 FigW FigH]);
                set(gca,'FontSize',12,'Linewidth',1)
                grid on
                grid minor
                set(gca,'YDir','normal')
                for epsIndex = epsDisplay
                    plot(N0Array(N0Display), VarianceStatsPassive(N0Display,epsIndex, etaIndex, nuIndex, gammaIndex, 1),'-','Linewidth',1,'Color',[0 1-epsIndex/length(epsArray) epsIndex/length(epsArray)])
                end
                plot(N0Array(N0Display), VarianceStatsAdaptive(N0Display,etaIndex, nuIndex, gammaIndex, 1),'-','Linewidth',2,'Color',[1 0 0])
                plot(N0Array(N0Display), N0Array(N0Display),'-','Linewidth',2,'Color',[0 0 0])
                
                leg = num2cell(["eps="+string(epsArrayDisplay) 'Adaptive' 'Shot Noise']);
                legend(leg,'FontSize',legfont,'Location','Northwest');
                
                xlabel('N0')
                ylabel('MS Error')
                title(strcat('MS Error of Mean vs. N eta=', num2str(eta)))
                fname = strcat(figureFolder, 'MSErrorMeanvsN0', suffix);
                print3(gcf,fname)
                close
            end  
            if FigureOption == MS_ERROR_MEAN_OPTIMAL_PASSIVE
                FigW = 8;
                FigH = 9;
                PosVector = [0,0,500,400];
                legfont = 8;

                figure
                hold on
                set(gcf,'PaperSize',[FigW FigH],'PaperPosition',[0 0 FigW FigH]);
                set(gca,'FontSize',12,'Linewidth',1)
                grid on
                grid minor
                set(gca,'YDir','normal')
                OptMSError = sum((VarianceStatsPassive(:,1, etaIndex, nuIndex, gammaIndex, 1)./N0Array).^2);
                OptEpsIndex = 1;
                for epsIndex = epsDisplay
                    NewMSError = sum((VarianceStatsPassive(:,epsIndex, etaIndex, nuIndex, gammaIndex, 1)./N0Array).^2);
                    if NewMSError < OptMSError
                      OptMSError = NewMSError;
                      OptEpsIndex = epsIndex;
                    end
                end
                plot(N0Array(N0Display), VarianceStatsPassive(N0Display,OptEpsIndex, etaIndex, nuIndex, gammaIndex, 1),'-','Linewidth',1,'Color',[0 1-OptEpsIndex/length(epsArray) OptEpsIndex/length(epsArray)])
                
                plot(N0Array(N0Display), VarianceStatsAdaptive(N0Display,etaIndex, nuIndex, gammaIndex, 1),'-','Linewidth',2,'Color',[1 0 0])
                plot(N0Array(N0Display), N0Array(N0Display),'-','Linewidth',2,'Color',[0 0 0])
                leg = {strcat('Passive MS Error eps=', num2str(epsArray(OptEpsIndex))),'Adaptive MS Error','Shot Noise'};
                legend(leg,'FontSize',legfont,'Location','SouthOutside');
                xlabel('N0')
                ylabel('MS Error')
                title(strcat('MS Error of Mean vs. N eta=', num2str(eta)))
                fname = strcat(figureFolder, 'OptPassiveMSErrorMeanvsN0', suffix);
                
                print3(gcf,fname)
                close
            end   
            if FigureOption == MS_ERROR_MLE_ALL_PASSIVE
                FigW = 8;
                FigH = 9;
                PosVector = [0,0,500,400];
                legfont = 8;

                figure
                hold on
                set(gcf,'PaperSize',[FigW FigH],'PaperPosition',[0 0 FigW FigH]);
                set(gca,'FontSize',12,'Linewidth',1)
                grid on
                grid minor
                set(gca,'YDir','normal')
                for epsIndex = epsDisplay
                    plot(N0Array(N0Display), VarianceStatsPassive(N0Display,epsIndex, etaIndex, nuIndex, gammaIndex, 6),'-','Linewidth',1,'Color',[0 1-epsIndex/length(epsArray) epsIndex/length(epsArray)])
                end
                plot(N0Array(N0Display), VarianceStatsAdaptive(N0Display,etaIndex, nuIndex, gammaIndex, 6),'-','Linewidth',2,'Color',[1 0 0])
                plot(N0Array(N0Display), N0Array(N0Display),'-','Linewidth',2,'Color',[0 0 0])
                leg = num2cell(["eps="+string(epsArrayDisplay) 'Adaptive' 'Shot Noise']);
                legend(leg,'FontSize',legfont,'Location','Northwest');
                xlabel('N0')
                ylabel('MS Error')
                title(strcat('MS Error of MLE vs. N eta=', num2str(eta)))
                fname = strcat(figureFolder, 'MSErrorMLEvsN0', suffix);
                %fname = 'LogLogPlotRMSErrorVariancevsEpsilonEta0.99Nu1e-6Gamma0.9N3';
                print3(gcf,fname)
                close
            end  
            if FigureOption == MS_ERROR_MLE_OPTIMAL_PASSIVE
                FigW = 8;
                FigH = 9;
                PosVector = [0,0,500,400];
                legfont = 8;

                figure
                hold on
                set(gcf,'PaperSize',[FigW FigH],'PaperPosition',[0 0 FigW FigH]);
                set(gca,'FontSize',12,'Linewidth',1)
                grid on
                grid minor
                set(gca,'YDir','normal')
                OptMSError = sum((VarianceStatsPassive(:,1, etaIndex, nuIndex, gammaIndex, 6)./N0Array).^2);
                OptEpsIndex = 1;
                for epsIndex = epsDisplay
                    NewMSError = sum((VarianceStatsPassive(:,epsIndex, etaIndex, nuIndex, gammaIndex, 6)./N0Array).^2);
                    if NewMSError < OptMSError
                      OptMSError = NewMSError;
                      OptEpsIndex = epsIndex;
                    end
                end
                plot(N0Array(N0Display), VarianceStatsPassive(N0Display,OptEpsIndex, etaIndex, nuIndex, gammaIndex, 6),'-','Linewidth',1,'Color',[0 1-OptEpsIndex/length(epsArray) OptEpsIndex/length(epsArray)])
                
                plot(N0Array(N0Display), VarianceStatsAdaptive(N0Display,etaIndex, nuIndex, gammaIndex, 1),'-','Linewidth',2,'Color',[1 0 0])
                plot(N0Array(N0Display), N0Array(N0Display),'-','Linewidth',2,'Color',[0 0 0])
                leg = {strcat('Passive MS Error eps=', num2str(epsArray(OptEpsIndex))),'Adaptive MS Error','Shot Noise'};
                legend(leg,'FontSize',legfont,'Location','SouthOutside');
                xlabel('N0')
                ylabel('MS Error')
                title(strcat('MS Error of MLE vs. N eta=', num2str(eta)))
                fname = strcat(figureFolder, 'OptPassiveMSErrorMLEvsN0', suffix);
                
                print3(gcf,fname)
                close
            end   
            if FigureOption == MS_ERROR_MEAN_ALL_PASSIVE_LL
                FigW = 4;
                FigH = 4;
                PosVector = [0,0,500,400];
                legfont = 8;

                figure
                hold on
                set(gcf,'PaperSize',[FigW FigH],'PaperPosition',[0 0 FigW FigH]);
                set(gca,'FontSize',12,'Linewidth',1)
                grid on
                grid minor
                set(gca,'YDir','normal')
                set(gca,'YScale','log')
                set(gca,'XScale','log')
                for epsIndex = epsDisplay
                    plot(N0Array(N0Display), VarianceStatsPassive(N0Display,epsIndex, etaIndex, nuIndex, gammaIndex, 1),'-s','Linewidth',1,'Color',[0 1-epsIndex/length(epsArray) epsIndex/length(epsArray)])
                end
                plot(N0Array(N0Display), VarianceStatsAdaptive(N0Display,etaIndex, nuIndex, gammaIndex, 1),'-s','Linewidth',2,'Color',[1 0 0])
                plot(N0Array(N0Display), N0Array(N0Display),'-','Linewidth',2,'Color',[0 0 0])
                leg = num2cell(["eps="+string(epsArrayDisplay) 'Adaptive' 'Shot Noise']);
                legend(leg,'FontSize',legfont,'Location','Northwest');
                xlabel('N0')
                ylabel('MS Error')
                %title(strcat('MS Error of Mean vs. N eta=', num2str(eta)))
                fname = strcat(figureFolder, 'LogLogMSErrorMeanvsN0', suffix);
                %fname = 'LogLogPlotRMSErrorVariancevsEpsilonEta0.99Nu1e-6Gamma0.9N3';
                print3(gcf,fname)
                close
            end  
            if FigureOption == MS_ERROR_MEAN_OPTIMAL_PASSIVE_LL
                FigW = 8;
                FigH = 9;
                PosVector = [0,0,500,400];
                legfont = 8;

                figure
                hold on
                set(gcf,'PaperSize',[FigW FigH],'PaperPosition',[0 0 FigW FigH]);
                set(gca,'FontSize',12,'Linewidth',1)
                grid on
                grid minor
                set(gca,'YDir','normal')
                OptMSError = sum((VarianceStatsPassive(:,1, etaIndex, nuIndex, gammaIndex, 1)./N0Array).^2);
                OptEpsIndex = 1;
                for epsIndex = epsDisplay
                    NewMSError = sum((VarianceStatsPassive(:,epsIndex, etaIndex, nuIndex, gammaIndex, 1)./N0Array).^2);
                    if NewMSError < OptMSError
                      OptMSError = NewMSError;
                      OptEpsIndex = epsIndex;
                    end
                end
                plot(N0Array(N0Display), VarianceStatsPassive(N0Display,OptEpsIndex, etaIndex, nuIndex, gammaIndex, 1),'-','Linewidth',1,'Color',[0 1-OptEpsIndex/length(epsArray) OptEpsIndex/length(epsArray)])
                
                plot(N0Array(N0Display), VarianceStatsAdaptive(N0Display,etaIndex, nuIndex, gammaIndex, 1),'-','Linewidth',2,'Color',[1 0 0])
                plot(N0Array(N0Display), N0Array(N0Display),'-','Linewidth',2,'Color',[0 0 0])
                leg = {strcat('Passive MS Error eps=', num2str(epsArrayDisplay(OptEpsIndex))),'Adaptive MS Error','Shot Noise'};
                legend(leg,'FontSize',legfont,'Location','SouthOutside');
                xlabel('N0')
                ylabel('MS Error')
                title(strcat('MS Error of Mean vs. N eta=', num2str(eta)))
                fname = strcat(figureFolder, 'LogLogOptPassiveMSErrorMeanvsN0', suffix);
                
                print3(gcf,fname)
                close
            end 
            if FigureOption == USEFUL_ROUNDS
                FigW = 4;
                FigH = 4;
                PosVector = [0,0,500,400];
                legfont = 8;

                figure
                hold on
                set(gcf,'PaperSize',[FigW FigH],'PaperPosition',[0 0 FigW FigH]);
                set(gca,'FontSize',12,'Linewidth',1)
                grid on
                grid minor
                set(gca,'YDir','normal')
                set(gca,'YScale','log')
                set(gca,'XScale','log')
                for epsIndex = epsDisplay
                    plot(N0Array(N0Display), VarianceStatsPassive(N0Display,epsIndex, etaIndex, nuIndex, gammaIndex, 7),'-s','Linewidth',1,'Color',[0 1-epsIndex/length(epsArray) epsIndex/length(epsArray)])
                end
                plot(N0Array(N0Display), VarianceStatsAdaptive(N0Display,etaIndex, nuIndex, gammaIndex, 7),'-s','Linewidth',2,'Color',[1 0 0])
                leg = num2cell(["eps="+string(epsArrayDisplay) 'Adaptive']);
                legend(leg,'FontSize',legfont,'Location','Northwest');
                xlabel('N0')
                ylabel('Rounds')
                ylim([0 Rounds])
                %title(strcat('MS Error of Mean vs. N eta=', num2str(eta)))
                fname = strcat(figureFolder, 'UsefulRoundsvsN0', suffix);
                print3(gcf,fname)
                close
            end  
        end
    end
end

%%
% Prob Distribution for values of N0

ADAPTIVE = 1;

FigureOption = ADAPTIVE;

N0Display = [4 8 12];
N0ArrayDisplay = N0Array(N0Display);

for etaIndex = 1:length(etaArray)
    eta = etaArray(etaIndex);
    eta %#ok<*NOPTS>
    for nuIndex = 1:length(nuArray)
        nu = nuArray(nuIndex);
        nu
        for gammaIndex = 1:length(gammaArray)
            gamma = gammaArray(gammaIndex);
            gamma
            
            suffix = sprintf('Nm%deta%dnu%dgamma%d',Nmax,etaIndex,nuIndex,gammaIndex);
            ProbEstimatesVsN0 = zeros(Nmax+1, length(N0Array));
            
            for N0Index = 1:length(N0Array)
                N0 = N0Array(N0Index);
                N0
                % Beginning of Adaptive Approach
                adaptivesuffix = sprintf('AdaptiveN%deta%dnu%dgamma%d.csv',N0Index,etaIndex,nuIndex,gammaIndex);
                fnameprob = strcat(dataFolder, 'Prob', adaptivesuffix);
                ProbEstimatesVsTrial = csvread(fnameprob);
                ProbEstimatesVsN0(:,N0Index) = ProbEstimatesVsTrial(1,:);
            end

            if FigureOption == ADAPTIVE
                FigW = 4;
                FigH = 4;
                PosVector = [0,0,500,400];
                legfont = 8;

                figure
                hold on
                set(gcf,'PaperSize',[FigW FigH],'PaperPosition',[0 0 FigW FigH]);
                set(gca,'FontSize',12,'Linewidth',1)
                grid on
                grid minor
                set(gca,'YDir','normal')
                for N0Index = N0Display
                    bar(0:Nmax, ProbEstimatesVsN0(:,N0Index))
                end

                leg = num2cell(["N0="+string(N0ArrayDisplay)]);
                legend(leg,'FontSize',legfont,'Location','Northwest');

                xlabel('N0')
                ylabel('Probability')
                %title(strcat('Bayesian Estimates for various N0, eta=', num2str(eta)))
                fname = strcat(figureFolder, 'BarPlotBayesianEst', suffix);
                print3(gcf,fname)
                close
            end  
        end
    end
end

%%
% Stats vs. N and eta
% Want to make graph 
% INCOMPLETE

for nuIndex = 1:length(nuArray)
    nu = nuArray(nuIndex);
    nu
    for gammaIndex = 1:length(gammaArray)
        gamma = gammaArray(gammaIndex);
        gamma
        if true
            FigW = 8;
            FigH = 9;
            PosVector = [0,0,500,400];
            legfont = 8;

            figure
            hold on
            set(gcf,'PaperSize',[FigW FigH],'PaperPosition',[0 0 FigW FigH]);
            set(gca,'FontSize',12,'Linewidth',1)
            grid on
            grid minor
            set(gca,'YDir','normal')
            for epsIndex = 1:length(epsArray)
                plot(N0Array, VarianceStatsPassive(:,epsIndex, etaIndex, nuIndex, gammaIndex, 1),'-','Linewidth',1,'Color',[0 1-epsIndex/length(epsArray) epsIndex/length(epsArray)])
            end
            plot(N0Array, VarianceStatsAdaptive(:,etaIndex, nuIndex, gammaIndex, 1),'-','Linewidth',2,'Color',[1 0 0])
            plot(N0Array, N0Array,'-','Linewidth',2,'Color',[0 0 0])
            %leg = {'Current Photons','Estimated Current Photons'};
            %legend(leg,'FontSize',legfont,'Location','SouthOutside');
            xlabel('N0')
            ylabel('Photons')
            title(strcat('MS Error vs. N eta=', num2str(eta)))
            fname = strcat(figureFolder, 'MSErrorvsN0', sprintf('eta%dnu%dgamma%d.csv',etaIndex,nuIndex,gammaIndex));
            %fname = 'LogLogPlotRMSErrorVariancevsEpsilonEta0.99Nu1e-6Gamma0.9N3';
            print3(gcf,fname)
            close
        end  
        if true
            FigW = 8;
            FigH = 9;
            PosVector = [0,0,500,400];
            legfont = 8;

            figure
            hold on
            set(gcf,'PaperSize',[FigW FigH],'PaperPosition',[0 0 FigW FigH]);
            set(gca,'FontSize',12,'Linewidth',1)
            grid on
            grid minor
            set(gca,'YDir','normal')
            OptMSError = sum((VarianceStatsPassive(:,1, etaIndex, nuIndex, gammaIndex, 1)./N0Array).^2);
            OptEpsIndex = 1;
            for epsIndex = 1:length(epsArray)
                NewMSError = sum((VarianceStatsPassive(:,epsIndex, etaIndex, nuIndex, gammaIndex, 1)./N0Array).^2);
                if NewMSError < OptMSError
                  OptMSError = NewMSError;
                  OptEpsIndex = epsIndex;
                end
            end
            plot(N0Array, VarianceStatsPassive(:,OptEpsIndex, etaIndex, nuIndex, gammaIndex, 1),'-','Linewidth',1,'Color',[0 1-OptEpsIndex/length(epsArray) OptEpsIndex/length(epsArray)])

            plot(N0Array, VarianceStatsAdaptive(:,etaIndex, nuIndex, gammaIndex, 1),'-','Linewidth',2,'Color',[1 0 0])
            plot(N0Array, N0Array,'-','Linewidth',2,'Color',[0 0 0])
            leg = {strcat('Passive MS Error eps=', num2str(epsArray(OptEpsIndex))),'Adaptive MS Error','Shot Noise'};
            legend(leg,'FontSize',legfont,'Location','SouthOutside');
            xlabel('N0')
            ylabel('Photons')
            title(strcat('MS Error vs. N eta=', num2str(eta)))
            fname = strcat(figureFolder, 'OptPassiveMSErrorvsN0', sprintf('eta%dnu%dgamma%d.csv',etaIndex,nuIndex,gammaIndex));

            print3(gcf,fname)
            close
        end   
    end
end