% method.cell_dist=336.6;
% method.canddiv=[];
% method.decision='bayes-v';
% method.distribution='unif-75';
% method.set=2
% classifier=trainQueryClassifier(method,1);
% method.set=3
% classifier=trainQueryClassifier(method,0);
% 
% method.cell_dist=336.6;
% method.canddiv=[];
% method.decision='bayes-v';
% method.distribution='expo-50';
% method.set=2
% classifier=trainQueryClassifier(method,1);
% method.set=3
% classifier=trainQueryClassifier(method,0);
% 
% method.set=1;
% method.cell_dist=336.6;
% method.canddiv=[];
% method.distribution='unif-75';
% method.decision='bayes-v'
% results = post_process(method);
% method.decision='bayes-dv'
% results = post_process(method);
% 
% method.set=1;
% method.cell_dist=336.6;
% method.canddiv=[];
% method.distribution='expo-50';
% method.decision='bayes-v'
% results = post_process(method);
% method.decision='bayes-dv'
% results = post_process(method);
% 
% method.set=4;
% method.cell_dist=336.6;
% method.canddiv=[];
% method.distribution='unif-75';
% method.decision='bayes-v'
% results = post_process(method);
% method.decision='bayes-dv'
% results = post_process(method);
% 
% method.set=4;
% method.cell_dist=336.6;
% method.canddiv=[];
% method.distribution='expo-50';
% method.decision='bayes-v'
% results = post_process(method);
% method.decision='bayes-dv'
% results = post_process(method);
% 
% % -----------------------------------------
% 
% method.cell_dist=336.6;
% method.canddiv=[10 100];
% method.decision='bayes-v';
% method.distribution='unif-75';
% method.set=2
% classifier=trainQueryClassifier(method,1);
% method.set=3
% classifier=trainQueryClassifier(method,0);

method.cell_dist=336.6;
method.canddiv=[10 100];
method.decision='bayes-v';
method.distribution='expo-50';
method.set=2
classifier=trainQueryClassifier(method,1);
method.set=3
classifier=trainQueryClassifier(method,0);

method.set=1;
method.cell_dist=336.6;
method.canddiv=[10 100];
method.distribution='unif-75';
method.decision='bayes-v'
results = post_process(method);
method.decision='bayes-dv'
results = post_process(method);

method.set=1;
method.cell_dist=336.6;
method.canddiv=[10 100];
method.distribution='expo-50';
method.decision='bayes-v'
results = post_process(method);
method.decision='bayes-dv'
results = post_process(method);

method.set=4;
method.cell_dist=336.6;
method.canddiv=[10 100];
method.distribution='unif-75';
method.decision='bayes-v'
results = post_process(method);
method.decision='bayes-dv'
results = post_process(method);

method.set=4;
method.cell_dist=336.6;
method.canddiv=[10 100];
method.distribution='expo-50';
method.decision='bayes-v'
results = post_process(method);
method.decision='bayes-dv'
results = post_process(method);