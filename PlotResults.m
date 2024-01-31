clear, close all

%ImportData
data = importdata('Results.mat');
CV = importdata('CVs.mat');
Ures = data.Ures;
Vres = data.Vres;
Lres = data.Lres;
MaxIter = size(Ures,2);

gfcPosition = [300, 300, 1200, 300];

% Deformed Mesh Results
figure;
for i = 1:MaxIter
    getR('ShowPlot', 'mesh', Ures(:,i));
    hold on
    xline((1+CV.umax), 'r', 'linewidth', 2)
    xline((1+CV.umax)*1.01, 'r', 'linewidth', 2)
    hold off
    title(num2str(CV.dt*i));
    set(gcf, 'Position',gfcPosition)
    
    xlimit = [-1 ,1.5+CV.umax];
    xlim(xlimit);
    ym = xlimit(2) - xlimit(1);
    ylim([-ym/2+0.5,ym/2+0.5]/5)
    pause(0.02);
end

% Stress Plot
figure;
for i = 1:MaxIter
    getR('ShowPlot', 'stre',Ures(:,i), 1, 1);
    titlestr = sprintf('\\sigma_{11} at t = %2.4f',i*CV.dt);
    title(num2str(titlestr));
  %  set(gcf, 'Position', gfcPosition)
    a = colorbar;
    a.Label.String = 'Stress (MPa)';
    clim([-16,10]);        %Change to clim/caxis for newer version of MATLAB
    view(2);
    xlimit = [-1 ,1.5+CV.umax];
    xlim(xlimit);
    ym = xlimit(2) - xlimit(1);
    ylim([-ym/2+0.5,ym/2+0.5]/5)
    pause(0.02);
end

% Total Internal Energy
figure;
Ener = zeros(MaxIter, 1);
for i = 1:MaxIter
    Ener(i) = getR('InteEner',Ures(:,i),Vres(:,i))*1e3;
end
xps = linspace(0,CV.dt*MaxIter*1e3,MaxIter);
plot(xps,Ener, 'linewidth', 2);
title('Total Internal Energy');
xlabel('Time (ms)');
ylabel('Internal Energy (kJ)')
ylim([0, max(Ener)*1.2])

% Lamda values
figure;
xps = linspace(0,CV.dt*MaxIter*1e3,MaxIter);
plot(xps, Lres, 'linewidth', 2);
title('\Lambda vs time');
xlabel('Time (ms)');
ylabel('\Lambda')
ylim([min(Lres)*1.2, max(Lres)*1.2])


% Contraint values at middle node
figure;
g = Ures(CV.CN(3),:) - CV.umax;
xps = linspace(0,CV.dt*MaxIter*1e3,MaxIter);
plot(xps, g, 'linewidth', 2);
title('Constraint Value vs time');
xlabel('Time (ms)');
ylabel('Constraint Value (g)')
ylim([min(g)*1.2, -min(g)*0.8])

