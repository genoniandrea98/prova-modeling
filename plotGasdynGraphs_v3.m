%% main function
function varargout = plotGasdynGraph(varargin)

% =========================================================================
% initialization

close all
clear
clc
fclose('all');
format long;
varargout={};
color = ['k', 'b', 'r', 'g', 'c', 'm', 'y', 'k', 'b', 'r', 'g', 'c', 'm', 'y'];

global nRegimes

% =========================================================================

% confCase =  { 'base', 'confA', 'confB', 'confC', 'confC1', 'confD', 'confE', 'confF', 'confG1', 'confH1'}
confCase =  { 'REFERENCE CI','1 CYL LPG OPTIMIZED','1 CYL LPG BASE'}
% confCase =  { 'confA', 'confB'}
% confCase =  { 'base', 'confA', 'confC', 'confC1'}
% confCase =  { 'REFERENCE_CI','1 CYL BOOST','TEST 30 5 VALVE 1'}
% confCase =  { 'confA', 'confF'}
% confCase =  { 'confG', 'confG1'}
% confCase =  { 'confF', 'confG1', 'confH1'}

% plot point:
probe = [1, 2];  %asp, sca
plotRegimes = [1000, 3000, 4000, 6000];
% intakeTiming  = [341, 570];
% exhaustTiming = [151, 380];
intakeTiming  = [360, 540];
exhaustTiming = [180, 360];

% =========================================================================

h = waitbar(0,'Reading files: please wait...');
for iCase = 1:length(confCase)
    
    gasdynCase = confCase(iCase);
    
    [...
        regime, ...
        volEfficiency, ...
        BMEP, ...
        BSFC, ...        
        indEfficiency, ...
        brakeTorque, ...
        brakePower, ...
        pumpingPower ...
    ] = ...
    readGasdynRegimesFile(gasdynCase{1});
    confReg(iCase, :) = regime(:);
    confVolEfficiency(iCase, :) = volEfficiency(:);
    confBMEP(iCase, :) = BMEP(:);
    confBSFC(iCase, :) = BSFC(:);    
    confIndEfficiency(iCase, :) = indEfficiency(:);
    confBrakeTorque(iCase, :) = brakeTorque(:);
    confBrakePower(iCase, :) = brakePower(:);
    confPumpingPower(iCase, :) = pumpingPower(:); 
    
    [...
        crankAngle, ...
        cylPresTheta, ...
        cylVol, ...
        cylPresVol, ...        
        regime ...
    ] = readGasdynCylPresFiles(gasdynCase{1});
    confCrankAngle(iCase, :, :) = crankAngle(:,:);
    confCylPresTheta(iCase, :, :) = cylPresTheta(:,:);
    confCylVol(iCase, :, :) = cylVol(:,:);
    confCylPresVol(iCase, :, :) = cylPresVol(:,:);
    
    [...
        crankAngle, ...
        cylBurntMass, ...
        cylBMF, ...
        cylAHRR ...        
    ] = readGasdynCombustionFiles(gasdynCase{1});

    for n = 1:nRegimes
        for i = 1:360
            confCylCombCrankAngle(iCase, i, n) = i-180;
        end
        startDeg = crankAngle(1, n);
        
        for i = 1:360
            k = i-(180+startDeg)+1;
            if (confCylCombCrankAngle(iCase, i, n) < min(crankAngle(:, n)))
                confCylBurntMass(iCase, i, n) = min(cylBurntMass(:, n));    
                confCylBMF(iCase, i, n) = min(cylBMF(:, n));   
                confCylAHRR(iCase, i, n) = 0;
            elseif (confCylCombCrankAngle(iCase, i, n) > max(crankAngle(:, n)))
                confCylBurntMass(iCase, i, n) = max(cylBurntMass(:, n));    
                confCylBMF(iCase, i, n) = max(cylBMF(:, n));    
                confCylAHRR(iCase, i, n) = 0;
            else 
                confCylBurntMass(iCase, i, n) = cylBurntMass(k, n);    
                confCylBMF(iCase, i, n) = cylBMF(k, n);    
                confCylAHRR(iCase, i, n) = cylAHRR(k, n);
            end
        end
            
%             if (crankAngle(i,n) == confCylCombCrankAngle(iCase, i, n))
%                 confCylCombCrankAngle(iCase, i, n) = crankAngle(i, n);
%                 confCylBurntMass(iCase, i, n) = cylBurntMass(i, n);    
%                 confCylBMF(iCase, i, n) = cylBMF(i, n);
%                 confCylAHRR(iCase, i, n) = cylAHRR(i, n);
%             else
%                 confCylCombCrankAngle(iCase, i, n) = i;
%                 if i < cylBurntMass(1, n)
%                     confCylBurntMass(iCase, i, n) = cylBurntMass(1, n);    
%                     confCylBMF(iCase, i, n) = cylBMF(1, n);
%                 elseif i > cylBurntMass(end, n)
%                     confCylBurntMass(iCase, i, n) = cylBurntMass(end, n);    
%                     confCylBMF(iCase, i, n) = cylBMF(end, n);                    
%                 end
%                 confCylAHRR(iCase, i, n) = 0;
%             end
%         end
    end


%     [...
%         crankAngle, ...
%         pres, ...
%         mass, ...
%         regime ...
%     ] = readGasdynFiles(gasdynCase{1});
%     confCrankAngle(iCase, :, :) = crankAngle(:,:);
%     confPres(iCase, :, :, :) = pres(:,:,:);
%     confMass(iCase, :, :, :) = mass(:,:,:);

    waitbar(iCase/length(confCase));

end
close(h);
    
%%

if exist('./graphs', 'dir')
   rmdir ./graphs s
end
mkdir ./graphs

nFig = 0;

% % pressure aspirazione
% nFig = nFig + 1;
% figure(nFig)
% title ('pressure')
% subPlot = 0;
% for n=1:nRegimes
% 
%     if (find(plotRegimes == regime(n)) > 0)
%         subPlot = subPlot + 1;
% 
%         hold on
%         subplot(ceil(0.5*length(plotRegimes)),ceil(0.5*length(plotRegimes)),subPlot);
% 
%         hold on
%         ymin = 10e10;
%         ymax = -10e10;
%         for iCase=1:length(confCase)
%             plot(confCrankAngle(iCase,:,n), confPres(iCase,:,probe(1),n), color(iCase));
%             ymin = min([ymin, confPres(iCase,:,probe(1),n)]);
%             ymax = max([ymax, confPres(iCase,:,probe(1),n)]);
%         end
%         legend(confCase, 'Location', 'northwest');
%         legend('boxoff');
%         title(['intake - ', num2str(regime(n))]);
%         xlabel('crank angle [°]');
%         ylabel('pressure [Pa]');
%         amp = ymax - ymin;
%         ymin = ymin - 0.25*amp;
%         ymax = ymax + 0.25*amp;
%         axis([0 720 ymin ymax]);
%         plot([intakeTiming(1), intakeTiming(1)], [ymin, ymax], 'b');
%         plot([intakeTiming(2), intakeTiming(2)], [ymin, ymax], 'b');
%     end
% 
% end
% saveas(nFig, './graphs/intake-pressure.png')
% 
% % pressure scarico
% nFig = nFig + 1;
% figure(nFig)
% title ('pressure')
% subPlot = 0;
% for n=1:nRegimes
% 
%     if (find(plotRegimes == regime(n)) > 0)
%         subPlot = subPlot + 1;
% 
%         hold on
%         subplot(ceil(0.5*length(plotRegimes)),ceil(0.5*length(plotRegimes)),subPlot);
% 
%         hold on
%         ymin = 10e10;
%         ymax = -10e10;
%         for iCase=1:length(confCase)
%             plot(confCrankAngle(iCase,:,n), confPres(iCase,:,probe(2),n), color(iCase));
%             ymin = min([ymin, confPres(iCase,:,probe(2),n)]);
%             ymax = max([ymax, confPres(iCase,:,probe(2),n)]);
%         end
%         legend(confCase, 'Location', 'northwest');
%         legend('boxoff');
%         title(['exhaust - ', num2str(regime(n))]);
%         xlabel('crank angle [°]');
%         ylabel('pressure [Pa]');
%         amp = ymax - ymin;
%         ymin = ymin - 0.25*amp;
%         ymax = ymax + 0.25*amp;
%         axis([0 720 ymin ymax]);
%         plot([exhaustTiming(1), exhaustTiming(1)], [ymin, ymax], 'r');
%         plot([exhaustTiming(2), exhaustTiming(2)], [ymin, ymax], 'r');
%     end
% 
% end
% 
% saveas(nFig, './graphs/exhaust-pressure.png')
% 
% % mass flow aspirazione
% nFig = nFig + 1;
% figure(nFig)
% title ('mass flow')
% subPlot = 0;
% for n=1:nRegimes
% 
%     if (find(plotRegimes == regime(n)) > 0)
%         subPlot = subPlot + 1;
% 
%         hold on
%         subplot(ceil(0.5*length(plotRegimes)),ceil(0.5*length(plotRegimes)),subPlot);
% 
%         hold on
%         ymin = 10e10;
%         ymax = -10e10;
%         for iCase=1:length(confCase)
%             plot(confCrankAngle(iCase,:,n), confMass(iCase,:,probe(1),n), color(iCase));
%             ymin = min([ymin, confMass(iCase,:,probe(1),n)]);
%             ymax = max([ymax, confMass(iCase,:,probe(1),n)]);
%         end
%         legend(confCase, 'Location', 'northwest');
%         legend('boxoff');
%         title(['intake - ', num2str(regime(n))]);
%         xlabel('crank angle [°]');
%         ylabel('mass flow [kg/s]');
%         amp = ymax - ymin;
%         ymin = ymin - 0.25*amp;
%         ymax = ymax + 0.25*amp;
%         axis([0 720 ymin ymax]);
%         plot([intakeTiming(1), intakeTiming(1)], [ymin, ymax], 'b');
%         plot([intakeTiming(2), intakeTiming(2)], [ymin, ymax], 'b');
%     end
% end
% 
% saveas(nFig, './graphs/intake-massflow.png')
% 
% % mass flow scarico
% nFig = nFig + 1;
% figure(nFig)
% title ('mass flow')
% subPlot = 0;
% for n=1:nRegimes
% 
%     if (find(plotRegimes == regime(n)) > 0)
%         subPlot = subPlot + 1;
% 
%         hold on
%         subplot(ceil(0.5*length(plotRegimes)),ceil(0.5*length(plotRegimes)),subPlot);
% 
%         hold on
%         ymin = 10e10;
%         ymax = -10e10;
%         for iCase=1:length(confCase)
%             plot(confCrankAngle(iCase,:,n), confMass(iCase,:,probe(2),n), color(iCase));
%             ymin = min([ymin, confMass(iCase,:,probe(2),n)]);
%             ymax = max([ymax, confMass(iCase,:,probe(2),n)]);
%         end
%         legend(confCase, 'Location', 'northwest');
%         legend('boxoff');
%         title(['exhaust - ', num2str(regime(n))]);
%         xlabel('crank angle [°]');
%         ylabel('mass flow [kg/s]');
%         amp = ymax - ymin;
%         ymin = ymin - 0.25*amp;
%         ymax = ymax + 0.25*amp;
%         axis([0 720 ymin ymax]);
%         plot([exhaustTiming(1), exhaustTiming(1)], [ymin, ymax], 'r');
%         plot([exhaustTiming(2), exhaustTiming(2)], [ymin, ymax], 'r');
%     end
% 
% end
% 
% saveas(nFig, './graphs/exhaust-massflow.png')
% 
% % mass flow aspirazione - time
% nFig = nFig + 1;
% figure(nFig)
% title ('mass flow')
% subPlot = 0;
% for n=1:nRegimes
% 
%     if (find(plotRegimes == regime(n)) > 0)
%         subPlot = subPlot + 1;
% 
%         hold on
%         subplot(ceil(0.5*length(plotRegimes)),ceil(0.5*length(plotRegimes)),subPlot);
% 
%         hold on
%         ymin = 10e10;
%         ymax = -10e10;
%         for iCase=1:length(confCase)
%             time = confCrankAngle(iCase,:,n)/regime(n)/60;
%             plot(time, confMass(iCase,:,probe(1),n), color(iCase));
%             ymin = min([ymin, confMass(iCase,:,probe(1),end)]);
%             ymax = max([ymax, confMass(iCase,:,probe(1),end)]);
%             xmax =  720/regime(1)/60;
%         end
%         legend(confCase, 'Location', 'northwest');
%         legend('boxoff');
%         title(['intake - ', num2str(regime(n))]);
%         xlabel('time [s]');
%         ylabel('mass flow [kg/s]');
%         amp = ymax - ymin;
%         ymin = ymin - 0.25*amp;
%         ymax = ymax + 0.25*amp;
%         axis([0 xmax ymin ymax]);
%         plot([intakeTiming(1), intakeTiming(1)], [ymin, ymax], 'b');
%         plot([intakeTiming(2), intakeTiming(2)], [ymin, ymax], 'b');
%     end
% end
% 
% saveas(nFig, './graphs/intake-massflow-time.png')
% 
% % mass flow scarico - time
% nFig = nFig + 1;
% figure(nFig)
% title ('mass flow')
% subPlot = 0;
% for n=1:nRegimes
% 
%     if (find(plotRegimes == regime(n)) > 0)
%         subPlot = subPlot + 1;
% 
%         hold on
%         subplot(ceil(0.5*length(plotRegimes)),ceil(0.5*length(plotRegimes)),subPlot);
% 
%         hold on
%         ymin = 10e10;
%         ymax = -10e10;
%         for iCase=1:length(confCase)
%             time = confCrankAngle(iCase,:,n)/regime(n)/60;
%             plot(time, confMass(iCase,:,probe(2),n), color(iCase));
%             ymin = min([ymin, confMass(iCase,:,probe(2),end)]);
%             ymax = max([ymax, confMass(iCase,:,probe(2),end)]);
%             xmax =  720/regime(1)/60;
%         end
%         legend(confCase, 'Location', 'northwest');
%         legend('boxoff');
%         title(['exhaust - ', num2str(regime(n))]);
%         xlabel('time [s]');
%         ylabel('mass flow [kg/s]');
%         amp = ymax - ymin;
%         ymin = ymin - 0.25*amp;
%         ymax = ymax + 0.25*amp;
%         axis([0 xmax ymin ymax]);
%         plot([exhaustTiming(1), exhaustTiming(1)], [ymin, ymax], 'r');
%         plot([exhaustTiming(2), exhaustTiming(2)], [ymin, ymax], 'r');
%     end
% end
% 
% saveas(nFig, './graphs/exhaust-massflow-time.png')


% cylinder pressure
nFig = nFig + 1;
figure(nFig)
title ('cylinder pressure')
subPlot = 0;
for n=1:nRegimes

    if (find(plotRegimes == regime(n)) > 0)
        subPlot = subPlot + 1;

        hold on
        subplot(ceil(0.5*length(plotRegimes)),ceil(0.5*length(plotRegimes)),subPlot);

        hold on
        ymin = 10e10;
        ymax = -10e10;
        for iCase=1:length(confCase)
            plot(confCrankAngle(iCase,:,n), confCylPresTheta(iCase,:,n), color(iCase));
            ymin = min([ymin, confCylPresTheta(iCase,:,n)]);
            ymax = max([ymax, confCylPresTheta(iCase,:,n)]);
        end
        legend(confCase, 'Location', 'northwest');
        legend('boxoff');
        title([num2str(regime(n)), ' rpm']);
        xlabel('crank angle [°]');
        ylabel('pressure [Pa]');
        amp = ymax - ymin;
        ymin = ymin - 0.25*amp;
        ymax = ymax + 0.25*amp;
        axis([0 720 ymin ymax]);
%         plot([intakeTiming(1), intakeTiming(1)], [ymin, ymax], 'b');
%         plot([intakeTiming(2), intakeTiming(2)], [ymin, ymax], 'b');
    end

end
saveas(nFig, './graphs/cylinder-pressure.png')

% cylinder pressure - volume
nFig = nFig + 1;
figure(nFig)
title ('cylinder pressure vs volume')
subPlot = 0;
for n=1:nRegimes

    if (find(plotRegimes == regime(n)) > 0)
        subPlot = subPlot + 1;

        hold on
        subplot(ceil(0.5*length(plotRegimes)),ceil(0.5*length(plotRegimes)),subPlot);

        hold on
        ymin = 10e10;
        ymax = -10e10;
        for iCase=1:length(confCase)
            plot(confCylVol(iCase,:,n), confCylPresVol(iCase,:,n), color(iCase));
            ymin = min([ymin, confCylPresVol(iCase,:,n)]);
            ymax = max([ymax, confCylPresVol(iCase,:,n)]);
        end
        legend(confCase, 'Location', 'northwest');
        legend('boxoff');
        title([num2str(regime(n)), ' rpm']);
        xlabel('volume [m3]');
        ylabel('pressure [Pa]');
        amp = ymax - ymin;
        ymin = ymin - 0.25*amp;
        ymax = ymax + 0.25*amp;
    end

end
saveas(nFig, './graphs/cylinder-pressure-volume.png')


% cylinder pressure - volume (log10)
nFig = nFig + 1;
figure(nFig)
title ('cylinder pressure vs volume (log10)')
subPlot = 0;
for n=1:nRegimes

    if (find(plotRegimes == regime(n)) > 0)
        subPlot = subPlot + 1;

        hold on
        subplot(ceil(0.5*length(plotRegimes)),ceil(0.5*length(plotRegimes)),subPlot);

        hold on
        ymin = 10e10;
        ymax = -10e10;
        for iCase=1:length(confCase)
            plot(confCylVol(iCase,:,n), confCylPresVol(iCase,:,n), color(iCase));
            %loglog(confCylPresVol(iCase,:,n));
            ymin = min([ymin, confCylPresVol(iCase,:,n)]);
            ymax = max([ymax, confCylPresVol(iCase,:,n)]);
        end
        legend(confCase, 'Location', 'northwest');
        legend('boxoff');
        title([num2str(regime(n)), ' rpm']);
        xlabel('volume [m3]');
        ylabel('pressure [Pa]');
        set(gca, 'YScale', 'log')
        amp = ymax - ymin;
        ymin = ymin - 0.01;
        ymax = ymax + 0.01;
        axis([0.9*min(confCylVol(1,:,1)) 1.05*max(confCylVol(1,:,1)) ymin ymax]);
    end

end
saveas(nFig, './graphs/cylinder-pressure-volume-log10.png')



% volumetric efficiency
nFig = nFig + 1;
figure(nFig)
title('volumetric efficiency')
for iCase=1:length(confCase)
    hold on
    plot(confReg(iCase,:), confVolEfficiency(iCase,:), color(iCase));
end
legend(confCase, 'Location', 'southwest');
legend('boxoff');
xlabel('regimes [gir/min]');
ylabel('volumetric efficiency [%]');
axis([regime(1)-500 regime(end)+500 0 140]);

for iCase=1:length(confCase)
    for n=1:length(regime)
        if (find(plotRegimes == regime(n)) > 0)
            hold on
            plot(confReg(iCase,n), confVolEfficiency(iCase,n), 'o', 'MarkerEdgeColor', color(iCase));
        end
    end
end

saveas(nFig, './graphs/volumetric-efficiency.png')


% BMEP
nFig = nFig + 1;
figure(nFig)
title('BMEP')
ymin = 10e10;
ymax = -10e10;
for iCase=1:length(confCase)
    hold on
    plot(confReg(iCase,:), confBMEP(iCase,:), color(iCase));
    ymin = min([ymin, confBMEP(iCase,:)]);
    ymax = max([ymax, confBMEP(iCase,:)]);
end
legend(confCase, 'Location', 'southwest');
legend('boxoff');
xlabel('regimes [gir/min]');
ylabel('BMEP [bar]');
axis([regime(1)-500 regime(end)+500 0 ymax+5]);

for iCase=1:length(confCase)
    for n=1:length(regime)
        if (find(plotRegimes == regime(n)) > 0)
            hold on
            plot(confReg(iCase,n), confBMEP(iCase,n), 'o', 'MarkerEdgeColor', color(iCase));
        end
    end
end

saveas(nFig, './graphs/BMEP.png')


% BSFC
nFig = nFig + 1;
figure(nFig)
title('BSFC')
ymin = 10e10;
ymax = -10e10;
for iCase=1:length(confCase)
    hold on
    plot(confReg(iCase,:), confBSFC(iCase,:), color(iCase));
    ymin = min([ymin, confBSFC(iCase,:)]);
    ymax = max([ymax, confBSFC(iCase,:)]);
end
legend(confCase, 'Location', 'southwest');
legend('boxoff');
xlabel('regimes [gir/min]');
ylabel('BSFC [g/MJ]');
axis([regime(1)-500 regime(end)+500 ymin-5 ymax+5]);

for iCase=1:length(confCase)
    for n=1:length(regime)
        if (find(plotRegimes == regime(n)) > 0)
            hold on
            plot(confReg(iCase,n), confBSFC(iCase,n), 'o', 'MarkerEdgeColor', color(iCase));
        end
    end
end

saveas(nFig, './graphs/BSFC.png')

% indEfficiency
nFig = nFig + 1;
figure(nFig)
title('indicated efficiency')
ymin = 10e10;
ymax = -10e10;
for iCase=1:length(confCase)
    hold on
    plot(confReg(iCase,:), confIndEfficiency(iCase,:), color(iCase));
    ymin = min([ymin, confIndEfficiency(iCase,:)]);
    ymax = max([ymax, confIndEfficiency(iCase,:)]);
end
legend(confCase, 'Location', 'southwest');
legend('boxoff');
xlabel('regimes [gir/min]');
ylabel('Indicated efficiency [-]');
axis([regime(1)-500 regime(end)+500 ymin-0.05 ymax+0.05]);

for iCase=1:length(confCase)
    for n=1:length(regime)
        if (find(plotRegimes == regime(n)) > 0)
            hold on
            plot(confReg(iCase,n), confIndEfficiency(iCase,n), 'o', 'MarkerEdgeColor', color(iCase));
        end
    end
end

saveas(nFig, './graphs/indEfficiency.png')

% brakeTorque
nFig = nFig + 1;
figure(nFig)
title('brake torque')
ymin = 10e10;
ymax = -10e10;
for iCase=1:length(confCase)
    hold on
    plot(confReg(iCase,:), confBrakeTorque(iCase,:), color(iCase));
    ymin = min([ymin, confBrakeTorque(iCase,:)]);
    ymax = max([ymax, confBrakeTorque(iCase,:)]);
end
legend(confCase, 'Location', 'southwest');
legend('boxoff');
xlabel('regimes [gir/min]');
ylabel('Brake torque [Nm]');
axis([regime(1)-500 regime(end)+500 0 ymax+5]);

for iCase=1:length(confCase)
    for n=1:length(regime)
        if (find(plotRegimes == regime(n)) > 0)
            hold on
            plot(confReg(iCase,n), confBrakeTorque(iCase,n), 'o', 'MarkerEdgeColor', color(iCase));
        end
    end
end

saveas(nFig, './graphs/brakeTorque.png')


% brakePower
nFig = nFig + 1;
figure(nFig)
title('brake power')
ymin = 10e10;
ymax = -10e10;
for iCase=1:length(confCase)
    hold on
    plot(confReg(iCase,:), confBrakePower(iCase,:), color(iCase));
    ymin = min([ymin, confBrakePower(iCase,:)]);
    ymax = max([ymax, confBrakePower(iCase,:)]);
end
legend(confCase, 'Location', 'southwest');
legend('boxoff');
xlabel('regimes [gir/min]');
ylabel('Brake torque [Nm]');
axis([regime(1)-500 regime(end)+500 0 ymax+5]);

for iCase=1:length(confCase)
    for n=1:length(regime)
        if (find(plotRegimes == regime(n)) > 0)
            hold on
            plot(confReg(iCase,n), confBrakePower(iCase,n), 'o', 'MarkerEdgeColor', color(iCase));
        end
    end
end

saveas(nFig, './graphs/brakePower.png')


% pumpingPower
nFig = nFig + 1;
figure(nFig)
title('pumping power')
ymin = 10e10;
ymax = -10e10;
for iCase=1:length(confCase)
    hold on
    plot(confReg(iCase,:), confPumpingPower(iCase,:), color(iCase));
    ymin = min([ymin, confPumpingPower(iCase,:)]);
    ymax = max([ymax, confPumpingPower(iCase,:)]);
end
legend(confCase, 'Location', 'southwest');
legend('boxoff');
xlabel('regimes [gir/min]');
ylabel('Brake torque [Nm]');
axis([regime(1)-500 regime(end)+500 ymin-1 0]);

for iCase=1:length(confCase)
    for n=1:length(regime)
        if (find(plotRegimes == regime(n)) > 0)
            hold on
            plot(confReg(iCase,n), confPumpingPower(iCase,n), 'o', 'MarkerEdgeColor', color(iCase));
        end
    end
end

saveas(nFig, './graphs/pumpingPower.png')


% burnt mass
nFig = nFig + 1;
figure(nFig)
title ('burnt mass')
subPlot = 0;
for n=1:nRegimes

    if (find(plotRegimes == regime(n)) > 0)
        subPlot = subPlot + 1;

        hold on
        subplot(ceil(0.5*length(plotRegimes)),ceil(0.5*length(plotRegimes)),subPlot);

        hold on
        ymin = 10e10;
        ymax = -10e10;
        for iCase=1:length(confCase)
            plot(confCylCombCrankAngle(iCase,:,n), confCylBurntMass(iCase,:,n), color(iCase));
            %loglog(confCylPresVol(iCase,:,n));
            ymin = min([ymin, confCylBurntMass(iCase,:,n)]);
            ymax = max([ymax, confCylBurntMass(iCase,:,n)]);
        end
        legend(confCase, 'Location', 'northwest');
        legend('boxoff');
        title([num2str(regime(n)), ' rpm']);
        xlabel('crank angle [°]');
        ylabel('burnt mass [kg]');
        amp = ymax - ymin;
        ymin = 0.95*ymin;
        ymax = 1.05*ymax;
        axis([-90 90 ymin ymax]);
    end

end
saveas(nFig, './graphs/burntMass.png')

% BMF
nFig = nFig + 1;
figure(nFig)
title ('Apparent Heat Release Rate')
subPlot = 0;
for n=1:nRegimes

    if (find(plotRegimes == regime(n)) > 0)
        subPlot = subPlot + 1;

        hold on
        subplot(ceil(0.5*length(plotRegimes)),ceil(0.5*length(plotRegimes)),subPlot);

        hold on
        ymin = 10e10;
        ymax = -10e10;
        for iCase=1:length(confCase)
            plot(confCylCombCrankAngle(iCase,:,n), confCylBMF(iCase,:,n), color(iCase));
            ymin = min([ymin, confCylBMF(iCase,:,n)]);
            ymax = max([ymax, confCylBMF(iCase,:,n)]);
        end
        legend(confCase, 'Location', 'northwest');
        legend('boxoff');
        title([num2str(regime(n)), ' rpm']);
        xlabel('crank angle [°]');
        ylabel('BMF [-]');
        amp = ymax - ymin;
        ymin = 0.95*ymin;
        ymax = 1.05*ymax;
        axis([-90 90 ymin ymax]);
    end

end
saveas(nFig, './graphs/BMF.png')

% AHRR
nFig = nFig + 1;
figure(nFig)
title ('Apparent Heat Release Rate')
subPlot = 0;
for n=1:nRegimes

    if (find(plotRegimes == regime(n)) > 0)
        subPlot = subPlot + 1;

        hold on
        subplot(ceil(0.5*length(plotRegimes)),ceil(0.5*length(plotRegimes)),subPlot);

        hold on
        ymin = 10e10;
        ymax = -10e10;
        for iCase=1:length(confCase)
            plot(confCylCombCrankAngle(iCase,:,n), confCylAHRR(iCase,:,n), color(iCase));
            ymin = min([ymin, confCylAHRR(iCase,:,n)]);
            ymax = max([ymax, confCylAHRR(iCase,:,n)]);
        end
        legend(confCase, 'Location', 'northwest');
        legend('boxoff');
        title([num2str(regime(n)), ' rpm']);
        xlabel('crank angle [°]');
        ylabel('AHRR [J/deg]');
        amp = ymax - ymin;
        ymin = 0.95*ymin;
        ymax = 1.05*ymax;
        axis([-90 90 ymin ymax]);
    end

end
saveas(nFig, './graphs/AHRR.png')

% confCylBurntMass(iCase, i, n) 
% confCylBMF(iCase, i, n) 
% confCylAHRR(iCase, i, n) 

end

%%
% /////////////////////////////////////////////////////////////////////////

function [...
            crankAngle, ...
            pres, ...
            mass, ...
            regime ...
    ] = ...
    readGasdynFiles(gasdynCase)

global nRegimes

file_{1}=[gasdynCase '/pressure.axi'];
file_{2}=[gasdynCase '/mass.axi'];

for i=1:2

    j=1;

    if(exist(file_{i},'file'))
        fp=fopen(file_{i});
    else
        error(['File: ' file_{i} ' not found. ']);
        return
    end
    
    % read data
    n = 1;
%     data = zeros(721, 5, 2, 7);
    while(~feof(fp))        
        try
            % read data one line at once
            data_ = textscan(fp,' %f, %f, %f , %f , %f',1,'collectOutput',1);
            data(j,:,i, n) = data_{1};
            
            if (data(j,1,i,n) == 720) 
                n=n+1;
                j=1;
            else
                j=j+1;
            end
      
        catch
            
            try
                regime_ = textscan(fp,'--------------------Engine speed[rpm] =  %f',1,'returnOnError',0);
                regime(n) = regime_{1};
            catch
                tmp =textscan(fp,' %s ',1,'MultipleDelimsAsOne',1,'delimiter','\n');
            end 
        end
    end
    fclose(fp);
end
nRegimes = n-1;

for n=1:nRegimes
    crankAngle(:,n) = data(:,1,1,n);
    pres(:,:,n) = data(:,2:end,1,n)*1e5;
    mass(:,:,n) = data(:,2:end,2,n);
end


end

%%
% /////////////////////////////////////////////////////////////////////////

function [...
            crankAngle, ...
            cylPresTheta, ...
            cylVol, ...            
            cylPresVol, ...
            reg ...
    ] = ...
    readGasdynCylPresFiles(gasdynCase)

global nRegimes

file_{1}=[gasdynCase '/cylPres.axi'];
file_{2}=[gasdynCase '/cylPresVol.axi'];

for i=1:2

    j=1;

    if(exist(file_{i},'file'))
        fp=fopen(file_{i});
    else
        error(['File: ' file_{i} ' not found. ']);
        return
    end
    
    % read data
    n = 1;
%     data = zeros(721, 5, 2, 7);
    while(~feof(fp))        
        try
            % read data one line at once
            data_ = textscan(fp,' %f, %f',1,'collectOutput',1);
            data(j,:,i, n) = data_{1};
            
            if (i == 1 && data(j,1,i,n) == 720)               
                n=n+1;
                j=1;
            elseif (i == 2 && data(j,1,1,n) == 720) 
                n=n+1;
                j=1;                    
            else
                j=j+1;
            end
      
        catch
            
            try
                regime_ = textscan(fp,'--------------------Engine speed[rpm] =  %f',1,'returnOnError',0);
                regime(i, n) = regime_{1};
            catch
                tmp =textscan(fp,' %s ',1,'MultipleDelimsAsOne',1,'delimiter','\n');
            end 
        end
    end
    fclose(fp);
end
nRegimes = n-1;

for n=1:nRegimes
    crankAngle(:,n) = data(:,1,1,n);
    cylPresTheta(:,n) = data(:,2,1,n)*1e5;
    cylVol(:,n) = data(:,1,2,n);    
    cylPresVol(:,n) = data(:,2,2,n);
end
reg = regime(1,:)


end

% /////////////////////////////////////////////////////////////////////////

function [...
            crankAngle, ...
            cylBurntMass, ...
            cylBMF, ...
            cylAHRR ...        
    ] = ...
    readGasdynCombustionFiles(gasdynCase)

    global nRegimes

    file_=[gasdynCase '/burntMass.axi'];

    j=0;

    if(exist(file_,'file'))
        fp=fopen(file_);
    else
        error(['File: ' file_ ' not found. ']);
        return
    end
    
    % read data
    n = 0;
%     data = zeros(721, 5, 2, 7);
    while(~feof(fp))        
        try
            % read data one line at once
            data_ = textscan(fp,' %f, %f, %f, %f',1,'collectOutput',1);
            j=j+1;
            data(j,:, n) = data_{1};
            
%             if (data(j,1,n) == 720)               
%                 n=n+1;
%                 j=1;                  
%             else
%                 j=j+1;
%             end
      
        catch
            
            try
                regime_ = textscan(fp,'--------------------Engine speed[rpm] =  %f',1,'returnOnError',0);
                n=n+1;
                j=0;
                regime(n) = regime_{1};
            catch
                tmp =textscan(fp,' %s ',1,'MultipleDelimsAsOne',1,'delimiter','\n');
                j=0;
            end 
        end
    end
    fclose(fp);

    nRegimes = n-1;

    for n=1:nRegimes
        crankAngle(:,n) = data(:,1,n);
        cylBurntMass(:,n) = data(:,2,n);
        cylBMF(:,n) = data(:,3,n);    
        cylAHRR(:,n) = data(:,4,n);
    end
    reg = regime(1,:)
end

% /////////////////////////////////////////////////////////////////////////

function [...
            regime, ...
            volEfficiency, ...
            BMEP, ...
            BSFC, ...            
            indEfficiency, ...
            brakeTorque, ...
            brakePower, ...
            pumpingPower ...
    ] = ...
    readGasdynRegimesFile(gasdynCase)


file_=[gasdynCase '/regimes.axi'];

for i=1:2

    j=1;

    if(exist(file_,'file'))
        fp=fopen(file_);
    else
        error(['File: ' file_ ' not found. ']);
        return
    end
    
    
    % read data
    while(~feof(fp))        
        try
            % read data one line at once
            data_ = textscan(fp,' %f, %f, %f , %f , %f, %f, %f, %f , %f , %f ,%f, %f, %f , %f , %f, %f, %f, %f',1,'collectOutput',1);
            data(j,:) = data_{1};
            j=j+1;     
        catch
            tmp =textscan(fp,' %s ',1,'MultipleDelimsAsOne',1,'delimiter','\n');
        end
    end
    fclose(fp);
end

for n = 1:j-1
    regime(n) = data(n,1);
    volEfficiency(n) = data(n,3);
    BMEP(n) = data(n,13);
    BSFC(n) = data(n,10);    
    indEfficiency(n) = data(n,18);
    brakeTorque(n) = data(n,6);
    brakePower(n) = data(n,8);
    pumpingPower(n) = data(n,9);
end

nRegimes = length(regime);

end


% /////////////////////////////////////////////////////////////////////////

