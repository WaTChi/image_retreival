% fix ground truth 3

gt3 = fopen('gt3.txt','w');

g = textread('gt3g.txt','%s');
y = textread('gt3y.txt','%s');
r = textread('gt3r.txt','%s');
b = textread('gt3b.txt','%s');
o = textread('gt3o.txt','%s');

idxG = strmatch('DSC_',strvcat(g));
idxY = strmatch('DSC_',strvcat(y));
idxR = strmatch('DSC_',strvcat(r));
idxB = strmatch('DSC_',strvcat(b));
idxO = strmatch('DSC_',strvcat(o));

for k=1:length(idxG)
    
    if k==length(idxG)
        
        for j=idxG(k):length(g)
            fprintf(gt3,'%s\n',g{j});
        end
        for j=idxY(k)+1:length(y)
            fprintf(gt3,'%s\n',y{j});
        end
        for j=idxR(k)+1:length(r)
            fprintf(gt3,'%s\n',r{j});
        end
        for j=idxB(k)+1:length(b)
            fprintf(gt3,'%s\n',b{j});
        end
        for j=idxO(k)+1:length(o)
            fprintf(gt3,'%s\n',o{j});
        end
        
    else
        
        for j=idxG(k):idxG(k+1)-1
            fprintf(gt3,'%s\n',g{j});
        end
        for j=idxY(k)+1:idxY(k+1)-1
            fprintf(gt3,'%s\n',y{j});
        end
        for j=idxR(k)+1:idxR(k+1)-1
            fprintf(gt3,'%s\n',r{j});
        end
        for j=idxB(k)+1:idxB(k+1)-1
            fprintf(gt3,'%s\n',b{j});
        end
        for j=idxO(k)+1:idxO(k+1)-1
            fprintf(gt3,'%s\n',o{j});
        end
        
    end
    
end

fclose(gt3);