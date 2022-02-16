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

VarianceStatsAdaptive = zeros([LN0, Leta, Lnu, Lgamma, 5]);
VarianceStatsPassive = zeros([LN0, Leps, Leta, Lnu, Lgamma, 5]);

set(0,'DefaultFigureWindowStyle','docked')
figure

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
                    
                adaptivetrialsuffix = sprintf('AdaptiveTrialN%deta%dnu%dgamma%d.csv',N0Index,etaIndex,nuIndex,gammaIndex);
                fnamestats = strcat(dataFolder, 'Stats', adaptivetrialsuffix);
                fnameinit = strcat(dataFolder, 'Init', adaptivetrialsuffix);
                fnameloop = strcat(dataFolder, 'Loop', adaptivetrialsuffix);

                TrialStats = readmatrix(fnamestats);
                InitEstimatesVsRoundPassive = readmatrix(fnameinit);
                LoopEstimatesVsRoundPassive = readmatrix(fnameloop);

                % MAKING THE FIGURE
                FigW = 8;
                FigH = 9;
                PosVector = [0,0,500,400];
                legfont = 8;

                hold on
                set(gcf,'PaperSize',[FigW FigH],'PaperPosition',[0 0 FigW FigH]);
                set(gca,'FontSize',12,'Linewidth',1)
                %grid on
                %grid minor
                ylim([1 Nmax])
                xlim([1 Rounds])
                set(gca,'YDir','normal')
                imagesc(1:1:Rounds,0:1:Nmax,InitEstimatesVsRoundPassive')
                grid off
                colormap parula
                plot(1:Rounds, N0*ones(Rounds,1),'-','Linewidth',3,'Color',[0 0 0])
                plot(1:Rounds, TrialStats(:,7),'-','Linewidth',2,'Color',[0 0.5 1])
                plot(1:Rounds, TrialStats(:,8),'-','Linewidth',2,'Color',[0 1 0.5])
                leg = {
                    'Init Photons'
                    'Init Photons Mean Estimate'
                    'Init Photons MLE'};
                legend(leg,'FontSize',legfont,'Location','SouthOutside');
                xlabel('Rounds')
                ylabel('Photons')
                title(strcat('Adaptive Trial Init N=', int2str(N0)))
                fname = strcat(figureFolder, 'AdaptiveTrialInit', sprintf('N%deta%dnu%dgamma%d',N0,etaIndex,nuIndex,gammaIndex));

                print3(gcf,fname)
                hold off

                FigW = 8;
                FigH = 9;
                PosVector = [0,0,500,400];
                legfont = 8;

                %set(0,'DefaultFigureWindowStyle','docked')
                hold on
                set(gcf,'PaperSize',[FigW FigH],'PaperPosition',[0 0 FigW FigH]);
                set(gca,'FontSize',12,'Linewidth',1)
                %grid on
                %grid minor
                ylim([1 Nmax])
                xlim([1 Rounds])
                set(gca,'YDir','normal')
                imagesc(1:1:Rounds,0:1:Nmax,LoopEstimatesVsRoundPassive')
                grid off
                colormap parula
                plot(1:Roundks, TrialStats(:,3),'-','Linewidth',3,'Color',[0 0 0])
                plot(1:Rounds, TrialStats(:,4),'-','Linewidth',2,'Color',[0 0.5 1])
                plot(1:Rounds, TrialStats(:,5),'-','Linewidth',2,'Color',[0 1 0.5])
                leg = {
                    'Loop Photons'
                    'Loop Photons Mean Estimate'
                    'Loop Photons MLE'};
                legend(leg,'FontSize',legfont,'Location','SouthOutside');
                xlabel('Rounds')
                ylabel('Photons')
                title(strcat('Adaptive Trial Loop N=', int2str(N0)))
                fname = strcat(figureFolder, 'AdaptiveTrialLoop', sprintf('N%deta%dnu%dgamma%d',N0,etaIndex,nuIndex,gammaIndex));

                print3(gcf,fname)
                close
                
                FigW = 8;
                FigH = 9;
                PosVector = [0,0,500,400];
                legfont = 8;

                %set(0,'DefaultFigureWindowStyle','docked')
                hold on
                set(gcf,'PaperSize',[FigW FigH],'PaperPosition',[0 0 FigW FigH]);
                set(gca,'FontSize',12,'Linewidth',1)
                grid on
                grid minor
                %ylim([1 Nmax])
                xlim([1 Rounds])
                set(gca,'YDir','normal')
                plot(1:Rounds, TrialStats(:,9),'-','Linewidth',2,'Color',[0 0 0])
                plot(1:Rounds, TrialStats(:,10),'-','Linewidth',2,'Color',[0 0.5 1])
                leg = {
                    'Available Information'
                    'Known Information'};
                legend(leg,'FontSize',legfont,'Location','SouthOutside');
                xlabel('Rounds')
                ylabel('Information')
                title(strcat('Available and Known Information N=', int2str(N0)))
                fname = strcat(figureFolder, 'AdaptiveTrialInformation', sprintf('N%deta%dnu%dgamma%d',N0,etaIndex,nuIndex,gammaIndex));

                print3(gcf,fname)
                close
            end    
        end    
    end    
end