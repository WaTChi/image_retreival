function num = dRound(num,k)

num = 10^k * round(10^(-k)*num);

end