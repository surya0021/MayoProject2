function d = getDPrime(x1,x2)
stdVal = sqrt((var(x1)+var(x2))/2);
d = (mean(x1)- mean(x2))/stdVal;
end