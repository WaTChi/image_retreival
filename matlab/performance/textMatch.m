function isMatch = textMatch(listA,listB)

A = length(listA);
if A == 1
    isMatch = ~isempty(find(strcmp(listA,listB),1));
else
    isMatch = false;
    for k=1:A
        isMatch = isMatch || ~isempty(find(strcmp(listA(k),listB),1));
    end
end