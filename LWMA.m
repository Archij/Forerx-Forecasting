function lwma = LWMA(Y, Period)

     s = size(Y);
     
     if (s(1) < s(2)) && ((s(1) ~= 1) || s(2)<1)
         error('Y is not a 1D time series!')
     end 
     
     if (s(1) > s(2)) && ((s(2) ~= 1) || s(1)<1)
         error('Y is not a 1D time series!')
     end
     
     if (Period > s(1)) && (s(1) > s(2))
         error('Period is higher than length of time series!')
     end
     
     if (Period > s(2)) && (s(2) > s(1))
         error('Period is higher than length of time series!')
     end
     
     lwma = [];
     N = length(Y);
     if (s(1) > s(2))
         lwma = zeros(N,1);
     elseif (s(1) < s(2))
         lwma = zeros(1,N);
     end
     
     lwma(Period) = sum([Y(1:Period)].*[Period:-1:1]')/sum([1:Period]);
     
     
     for i=Period+1:N
         
         lwma(i) = sum([Y(i-Period+1:i)].*[Period:-1:1]')/sum([1:Period]);
         
     end
      
end