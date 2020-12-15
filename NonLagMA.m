function MABuffer = NonLagMA(Price, Length, Displace, PctFilter, Deviation)

%---- indicator buffers
Bars = length(Price);
MABuffer = zeros(Bars,1);
UpBuffer = zeros(Bars,1);
DnBuffer = zeros(Bars,1);
trend = [];
Del = [];
AvgDel = [];


alfa = [];
i = []; 
Phase = []; 
Len = [];
Cycle=4;
Coeff = [];
beta = [];
t = [];
Sum = [];
Weight = [];
g = [];
UpTrendAlert=false;
DownTrendAlert=false;
 
Coeff =  3*pi;
Phase = Length-1;
Len = Length*4 + Phase;  
alfa = zeros(Len,1);
Weight=0;    
      
for i=0:Len-1

    if (i<=Phase-1) 
         t = 1.0*i/(Phase-1);
    else
         t = 1.0 + (i-Phase+1)*(2.0*Cycle-1.0)/(Cycle*Length-1.0);  
    end
    beta = cos(pi*t);
    g = 1.0/(Coeff*t+1);   
    if t <= 0.5
          g = 1;
    end
    alfa(i+1) = g * beta;
    Weight = Weight +  alfa(i+1);
end
 
%+------------------------------------------------------------------+
%| NonLagMA_v7.1                                                      |
%+------------------------------------------------------------------+

price = Price;      
   
for shift=Len:Bars
	
    Sum = 0;
    for i=0:Len-1
	      
      Sum = Sum + alfa(i+1)*price(shift-i);
      
    end
   
	if (Weight > 0) 
        MABuffer(shift) = (1.0+Deviation/100)*Sum/Weight;
    end
   
    if (PctFilter>0)

      Del(shift) = abs(MABuffer(shift) - MABuffer(shift-1));
   
      sumdel=0;
      for i=0:Length-1
          sumdel = sumdel + Del(shift-i);
      end
      AvgDel(shift) = sumdel/Length;
    
      sumpow = 0;
      for i=0:Length-1
          sumpow = sumpow + pow(Del(shifti-i)-AvgDel(shift-i),2);
      end
      StdDev = sqrt(sumpow/Length); 
     
      
      Filter = PctFilter * StdDev;

     
    end
         


end