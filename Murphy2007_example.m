% We first make the implementation of the model
gillespie_maker('Murphy2007');
% This will create or overwrite a file Murphy2007.m

% Now, to average say nruns trajectories, we set
nruns = 5; % number of simulations
maxt = 5000; % maximum time
dist = cell(1,nruns); % to save all trajectories
tbins= 0:50:maxt; % time bins to calculate averages
colors = hsv(nruns); % for plotting

% This is simply to illustrate the use of setvalues.
% setvalues.('N') = 2; 
% setvalues.('p') = 2; 
% If used, it would set initial N to 2, and translation constant p
% to 2. You can try replacing [] by setvalus below

figure; hold on;
for i = 1:nruns
    time_serie = Murphy2007(1000000,maxt,{'P'},[]); % Simulation
    %% time_serie = Murphy2007(1000000,maxt,{'P'},setvalues); 
    subplot(1,2,1); hold on;
    plot(time_serie(:,1),time_serie(:,2),'Color', colors(i,:)); % Plotting
    axis([0 maxt 0 1200]);
    title('Single Trajectories');
    xlabel('time');ylabel('P (a.u.)');
    dist{i} = time_serie; % Saving the trayectory
end

% Now we average all our trajectories using the time intervals from
% tbins
[tmeans,tstds] = average_gillespie(dist,tbins); % Averaging
subplot(1,2,2);
plot(tbins,tmeans,'Color','black'); % Ploting average
axis([0 maxt 0 1200]);
title('Average trajectory');
xlabel('time');ylabel('P (a.u.)');