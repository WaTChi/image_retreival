% testTrain

classifier = struct;
classifier.dists = {'Normal','Normal','Normal'};

sampling = -3:1e-3:+3;
ntrain = length(sampling);

f = zeros(ntrain,3);
c = zeros(ntrain,1);
w = zeros(ntrain,1);
for k = 1:ntrain
    
    % features
    f(k,1) = sampling(k);
    f(k,2) = randn;
    f(k,3) = randn;
    
    % classes
    c(k) = 1;
    
    % weights
%     w(k) = exp(-f(k,1)^2/2);
%     w(k) = rand;
    w(k) = 1;
    
end

classifier = trainbayes(f,c,w,classifier)