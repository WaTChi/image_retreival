function num = dRound(num,k)

if nargin < 2
    k = 0;
end
num = 10^k * round(10^(-k)*num);

end