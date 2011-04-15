function [results] = plotResults(query_set,method,fignum)

% [results] = plotResults(query_set,method,fignum)
% 
% DESCRIPTION
%   This functions 

if method.fuzzy
    results_file = ['.\',method.decision,'\query',num2str(query_set), ...
        'fuzzy',num2str(dRound(method.cell_dist,0)),'_results.mat'];
else
    results_file = ['.\',method.decision,'\query',num2str(query_set), ...
        'exact',num2str(dRound(method.cell_dist,0)),'_results.mat'];
end
load(results_file)

wv = method.wv;
ws = method.ws;

% Check for incomplete  results
todo = [];
for k=1:length(wv)
    idxv = find(wv(k)==results.wv);
    idxs = find(ws(k)==results.ws);
    idx = intersect(idxv,idxs);
    if isempty(idx)
        todo = [todo,k];
    end
end
if ~isempty(todo)
    disp('Some results not complete. Computing these now...')
    method_todo = method;
    method_todo.wv = wv(todo);
    method_todo.ws = ws(todo);
    results = post_process(query_set,method_todo);
end

for k=1:length(wv)
    
    idxv = find(wv(k)==results.wv);
    idxs = find(ws(k)==results.ws);
    idx = intersect(idxv,idxs);
    match = results.match(idx,:);
    poss = results.poss(idx,:);
    total = results.total(idx,:);
    
    figure(fignum+k)
    hold off
    bar(1:10,poss./total,'FaceColor','y')
    hold on
    bar(1:10,match./total,'FaceColor','g')
    axis([0.5 10.5 0 1])
    title(['Post-processing performance on top N for v = ', ...
        num2str(wv(k)),', r = ',num2str(ws(k))])
    xlabel('Top N')
    ylabel('Fraction of successful matches; Y=possible, G=actual')
    
end