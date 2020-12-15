function ExtMapBuffer1 = ZeroLagJMA(Close, Period)

ExtMapBuffer1  = zeros(length(Close),1);

Lag = floor((Period-1)/1);

Data = Close(Lag:end) + (Close(Lag:end) - Close(1:end-Lag+1));

ExtMapBuffer1(Lag:end) = JMA(Data,Period,0);

end