% testNaiveBayes

feat_chance = .5;
ntrain = 1000;
ntest = 10000;
C = 3; priors = [.6;.3;.1]; cdf = cumsum(priors);
M = 5; dists = cellstr(strvcat('Gamma','Normal','Gamma','Gamma','Normal'))';

% Generate class parameters
prmA = zeros(C,M);
prmB = zeros(C,M);
for m = 1:M
    if strcmp(dists{m},'Gamma')
        prmA(:,m) = random('gam',2,2,C,1);
        prmB(:,m) = random('gam',1,5,C,1);
    else % if strcmp(dists{c},'Normal')
        prmA(:,m) = random('gam',3,2,C,1);
        prmB(:,m) = random('gam',2,1,C,1);
    end
end

% generate training and testing samples
train_class = zeros(ntrain,1);
train_features = zeros(ntrain,M);
train_weights = zeros(ntrain,1);
for k=1:ntrain
    seed = rand; c = 1;
    while seed > cdf(c)
        c = c+1;
    end
    train_class(k) = c;
    train_weights(k) = rand;
    for m=1:M
        if rand < feat_chance
            train_features(k,m) = random(dists{m},prmA(c,m),prmB(c,m));
        else
            train_features(k,m) = NaN;
        end
    end
end
test_class = zeros(ntest,1);
test_features = zeros(ntest,M);
test_weights = zeros(ntrain,1);
for k=1:ntest
    seed = rand; c = 1;
    while seed > cdf(c)
        c = c+1;
    end
    test_class(k) = c;
    test_weights(k) = rand;
    for m=1:M
        if rand < feat_chance
            test_features(k,m) = random(dists{m},prmA(c,m),prmB(c,m));
        else
            test_features(k,m) = NaN;
        end
    end
end

% train classifier
classifier = trainbayes(train_features,train_class,train_weights,dists);

% test classifier on training and testing samples
[train_class_pred,train_condP] = classifybayes(train_features,classifier);
[test_class_pred,test_condP] = classifybayes(test_features,classifier);