function [classifier] = collapsebayes(file,comp,naive)

file = ['./bayes/classifier/',file];
load(file)

if naive
    
    M = length(classifier.spacing);
    
    for m=1:M
        
        classifier.spacing(m) = classifier.spacing(m) * comp(m);
        K = ceil( size(classifier.bins{m},1) / comp(m) );
        
        bins = zeros(K,2);
        for k=1:K
            if k==K
                bins(k,:) = sum( classifier.bins{m}(comp(m)*(k-1)+1:end,:) , 1 );
            else
                bins(k,:) = sum( classifier.bins{m}(comp(m)*(k-1)+1:comp(m)*k,:) , 1 );
            end
        end
        
        classifier.bins{m} = bins;
        
    end      
    
else % 2D distribution
    
end

save(file,'classifier')