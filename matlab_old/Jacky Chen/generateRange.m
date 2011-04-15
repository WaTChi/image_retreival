function [R] = generateRange(p, range)

p = p{1}; % to extract char from cell
R = repmat(p, 2*range+1, 1);

places = 6;
extL = 5; %.feat

str = p(end-places-extL+1: end-extL);
n = str2double(str);
k = n-range;
for i=0:2*range
    new = num2str(k+i, ['%0', num2str(places), 'd']);
    R(i+1,:) = strrep(p, str, new);
end