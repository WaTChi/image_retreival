function isMatch = textMatch(listA,listB)

isMatch = false;
for k=1:length(listA)
    isMatch = isMatch || ~isempty(find(strcmp(listA(k),listB),1));
end