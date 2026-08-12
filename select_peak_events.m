function selected = select_peak_events(candidates, series, minSep)
% SELECT_PEAK_EVENTS Cluster candidate peak indices that fall within
% minSep samples of one another, keeping only the strongest peak (by
% value in series) from each cluster.
%
%   selected = SELECT_PEAK_EVENTS(candidates, series, minSep)

if isempty(candidates)
    selected = [];
    return
end

clusterID = cumsum([1; diff(candidates) > minSep]);
selected = [];

for k = 1:max(clusterID)
    clusterMembers = candidates(clusterID == k);
    [~, idxMax] = max(series(clusterMembers));
    selected(end+1) = clusterMembers(idxMax); %#ok<AGROW>
end

selected = selected(:);

end
