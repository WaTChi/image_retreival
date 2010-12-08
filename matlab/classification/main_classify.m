% Main script for running classification analysis

% Variables
new_train = true;
query_dir = '..........';
method_dir = 'E:\Research\app\code\matlab\classification\NaiveBayes-Gaussian\';


addpath 'E:\Research\app\code\matlab\classification\NaiveBayes-Gaussian\'


%------------------------------------------------------
% 
% CODE TO PARSE AND GET FEATURES AND CLASSES FROM FILES
% 
%------------------------------------------------------

distributions = parseVotes(query_dir);
features = getFeatures(distributions);
classes = getClasses();

if new_train
    classifier = train(features,classes);
else
    load([method_dir,'classifier.mat'])
end

class_predictions = classify(features,classifier);