% testTrain

classifier = struct;
classifier.spacing = 1e-1;

ntrain = 10000;

f = zeros(ntrain,1);
c = zeros(ntrain,1);
w = zeros(ntrain,1);
for k = 1:ntrain
    
    % features
    f(k) = random('gam',5,1);
    
    % classes
    c(k) = 1;
    
    % weights
    w(k) = 1;
    
end

classifier = trainbayes(f,c,classifier,w)