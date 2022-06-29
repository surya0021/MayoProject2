function d = getDPrime(x1,x2)
n1 = length(x1);    n2 = length(x2);
stdVal = sqrt(((n1-1)*var(x1,'omitnan')+(n2-1)*var(x2,'omitnan'))/(n1+n2-2));
d = (mean(x1,'omitnan')- mean(x2,'omitnan'))/stdVal;
end