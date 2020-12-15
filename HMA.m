function MAOnBuffer = HMA(Close, HMAPeriod, InpMAMethod)

p = floor(sqrt(HMAPeriod));
medp = floor(HMAPeriod/2);

if strcmp(InpMAMethod,'MODE_LWMA')
    HMABuffer = 2*LWMA(Close,medp)-LWMA(Close,HMAPeriod);
    MAOnBuffer = LWMA(HMABuffer, p);
elseif strcmp(InpMAMethod,'MODE_SMA')
    HMABuffer = 2*SMA(Close,medp)-SMA(Close,HMAPeriod);
elseif strcmp(InpMAMethod,'MODE_EMA')
    HMABuffer = 2*EMA(Close,medp)-EMA(Close,HMAPeriod);
elseif strcmp(InpMAMethod,'MODE_SMMA')
    HMABuffer = 2*SMMA(Close,medp)-SMMA(Close,HMAPeriod);
elseif strcmp(InpMAMethod,'MODE_JMA')
    HMABuffer = 2*JMA(Close,medp,0)-JMA(Close,HMAPeriod,0);
    MAOnBuffer = JMA(HMABuffer, p, 0);
end



end