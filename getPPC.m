function ppc = getPPC(thetaList)
input = exp(1i*thetaList(:));
siz = size(input);
n = siz(1);
if n>1
    outsum        = nansum(input);      % compute the sum; this is 1 x size(2:end)
    ppc  = (outsum.*conj(outsum) - n)./(n*(n-1)); % do the pairwise thing in a handy way
else
    error('computation of PPC requires >1 trial, please feed all trial dataset into computeCoherencyFromSpectrum program')
end

end