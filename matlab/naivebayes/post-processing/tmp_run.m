% method.set='query2';
% method.cell_dist=271.6;
% method.canddiv=[];
% method.decision='bayes-dv';
% method.distribution='exact-'
% classifier = trainQueryClassifier(method,1);
% 
% % -----------------------------------------------------
% 
% method.set='query1';
% method.cell_dist=271.6;
% method.canddiv=[];
% method.distribution='exact-';
% method.decision='bayes-v'
% results = post_process(method,0);
% method.decision='bayes-dv'
% results = post_process(method,0);
% 
% method.set='query3';
% method.cell_dist=271.6;
% method.canddiv=[];
% method.distribution='exact-';
% method.decision='bayes-v'
% results = post_process(method,0);
% method.decision='bayes-dv'
% results = post_process(method,0);
% 
% method.set='query4';
% method.cell_dist=271.6;
% method.canddiv=[];
% method.distribution='exact-';
% method.decision='bayes-v'
% results = post_process(method,0);
% method.decision='bayes-dv'
% results = post_process(method,0);
% 
% method.set='query5horizontal';
% method.cell_dist=271.6;
% method.canddiv=[];
% method.distribution='exact-';
% method.decision='bayes-v'
% results = post_process(method,0);
% method.decision='bayes-dv'
% results = post_process(method,0);
% 
method.set='query5vertical';
method.cell_dist=336.6;
method.canddiv=[];
method.distribution='exact-';
method.decision='bayes-v'
results = post_process(method,1);
method.decision='bayes-dv'
results = post_process(method,0);

% ------------------------------------------------
% 
% method.set='query1';
% method.cell_dist=271.6;
% method.canddiv=[];
% method.distribution='exact-';
% method.decision='vote-'
% results = post_process(method,0);
% 
% method.set='query3';
% method.cell_dist=271.6;
% method.canddiv=[];
% method.distribution='exact-';
% method.decision='vote-'
% results = post_process(method,0);
% 
% method.set='query4';
% method.cell_dist=271.6;
% method.canddiv=[];
% method.distribution='exact-';
% method.decision='vote-'
% results = post_process(method,0);
% 
% method.set='query5horizontal';
% method.cell_dist=271.6;
% method.canddiv=[];
% method.distribution='exact-';
% method.decision='vote-'
% results = post_process(method,0);
% 
% method.set='query5vertical';
% method.cell_dist=271.6;
% method.canddiv=[];
% method.distribution='exact-';
% method.decision='vote-'
% results = post_process(method,0);