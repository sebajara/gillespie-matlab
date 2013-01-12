function [tmeans,tstds] = averate_gillespie(dist,tbins)
    
% Average the values from DIST given TBINS It expects DIST to be a
% cell array dimensions (NDIST,1). Each value DIST{j} (j=1:NDIST) is
% expected to be a matrix, where the first columns are the time
% values, and the rest are the number of instances for each specie, at
% each time in the simulation. For each specie, it will collect all
% the values for all j that are within the time intervals specified in
% TBINS. Finally it will return TMEANS and TSTDS, being the mean and
% standard deviation for each specie at each time interval in TBINS.
     
[ndist,~] = size(dist);
[~,nvals] = size(dist{1});
nvals = nvals - 1; % we expect the first column to be time
ntbins = length(tbins);
nvtdist = cell(ndist,nvals,ntbins);

% group values in the time dimention
for j = 1:ndist
    serie = dist{j};
    [~,timebin] = histc(serie(:,1),tbins);
    for v = 1:nvals
        vals = serie(:,v+1);
        for t = 1:ntbins
            nvtdist{j,v,t} = vals(find(timebin==t));
        end
    end
end

tmeans = zeros(ntbins,nvals);
tstds  = zeros(ntbins,nvals);
% now average for each time interval
for t = 1:ntbins
    for v = 1:nvals
        vals = cat(1,nvtdist{:,v,t});
        tmeans(t,v) = mean(vals);
        tstds(t,v) = std(vals);
    end
end

end