function saveClassifier(method,cls)

% Parse distribution parameter
idx = strfind(method.distribution,'-');
distr = method.distribution(1:idx-1);
distr_prm = str2double(method.distribution(idx+1:end));

% Parse the decision parameter
idx = strfind(method.decision,'-');
decision = method.decision(1:idx-1);
decis_prm = method.decision(idx+1:end);

D = 1 + length(method.featdiv);
for d=1:D
    if d==1
        num1=0;
    else
        num1=method.featdiv(d-1);
    end
    if d==D
        num2=Inf;
    else
        num2=method.featdiv(d);
    end
    file = ['./bayes/classifier/',distr,num2str(dRound(method.cell_dist,0)), ...
        '_',num2str(num1),',',num2str(num2),'_',decision,'.mat'];
    classifier = cls{d};
    save(file,'classifier');
end