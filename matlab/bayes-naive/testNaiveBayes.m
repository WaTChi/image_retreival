% testNaiveBayes

feat_chance = .5;
ntrain = 1000;
ntest = 10000;
C = 3; priors = [.6;.3;.1]; cdf = cumsum(priors);
M = 5; spacing = .1 * ones(1,M);

% Generate class parameters
prmA = zeros(C,M);
prmB = zeros(C,M);
for m = 1:M
    prmA(:,m) = random('gam',10,1/5,C,1);
    prmB(:,m) = 0.5 + 1.5*rand;
end

% generate training and testing samples
train_class = zeros(ntrain,1);
train_features = zeros(ntrain,M);
for k=1:ntrain
    seed = rand; c = 1;
    while seed > cdf(c)
        c = c+1;
    end
    train_class(k) = c;
    for m=1:M
        if rand < feat_chance
            train_features(k,m) = random('gam',prmA(c,m),prmB(c,m));
        else
            train_features(k,m) = NaN;
        end
    end
end
test_class = zeros(ntest,1);
test_features = zeros(ntest,M);
for k=1:ntest
    seed = rand; c = 1;
    while seed > cdf(c)
        c = c+1;
    end
    test_class(k) = c;
    for m=1:M
        if rand < feat_chance
            test_features(k,m) = random(dists{m},prmA(c,m),prmB(c,m));
        else
            test_features(k,m) = NaN;
        end
    end
end

% train classifier
classifier = trainbayes(train_features,train_class,spacing);

% test classifier on training and testing samples
[train_class_pred,train_condP] = classifybayes(train_features,classifier);
[test_class_pred,test_condP] = classifybayes(test_features,classifier);