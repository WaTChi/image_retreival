function [classifier] = collapsebayes(file,comp,naive)

if naive
    
    file = ['./bayes/classifier/',file];
    load(file)
    
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
    
    file = ['./bayes_future/classifier/',file];
    load(file)
    
    classifier.spacing = classifier.spacing .* comp;
    
    J = ceil( size(classifier.bins,1) / comp(1) );
    K = ceil( size(classifier.bins,2) / comp(2) );
    bins = zeros(J,K,2);
    for j=1:J
        for k=1:K
            if j==J && k==K
                bins(j,k,:) = sum( sum( classifier.bins(comp(1)*(j-1)+1:end,comp(2)*(k-1)+1:end,:) , 1 ) , 2 );
            elseif j==J
                bins(j,k,:) = sum( sum( classifier.bins(comp(1)*(j-1)+1:end,comp(2)*(k-1)+1:comp(2)*k,:) , 1 ) , 2 );
            elseif k==K
                bins(j,k,:) = sum( sum( classifier.bins(comp(1)*(j-1)+1:comp(1)*j,comp(2)*(k-1)+1:end,:) , 1 ) , 2 );
            else
                bins(j,k,:) = sum( sum( classifier.bins(comp(1)*(j-1)+1:comp(1)*j,comp(2)*(k-1)+1:comp(2)*k,:) , 1 ) , 2 );
            end
        end
    end
    
    classifier.bins = bins;
    
end

save(file,'classifier')