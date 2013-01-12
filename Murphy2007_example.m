% We first make the implementation of the model
gillespie_maker('Murphy2007');

% Now, to average say 10 trajectories
nruns = 10; % number of simulations
maxt = 5000; % maximum time
dist = cell(1,nruns); % to save all trajectories
tbins= 0:50:maxt; % time bins to calculate averages
colors = hsv(nruns); % for plotting

% This is simply to illustrate the use of setvalues.
% setvalues.('N') = 2; % This would be for two copies of N
% setvalues.('p') = 2; % This would increase protein translation twice

figure; hold on;
for i = 1:nruns
    time_serie = Murphy2007(1000000,maxt,{'P'},[]); % Simulation
    subplot(1,2,1); hold on;
    plot(time_serie(:,1),time_serie(:,2),'Color', colors(i,:)); % Plotting
    title('Single Trajectories');
    xlabel('time');ylabel('P (a.u.)');
    dist{i} = time_serie; % Saving the trayectory
end

[tmeans,tstds] = average_gillespie(dist,tbins); % Averaging
subplot(1,2,2);
plot(tbins,tmeans,'Color','red'); % Ploting average
title('Average trajectory');
xlabel('time');ylabel('P (a.u.)');