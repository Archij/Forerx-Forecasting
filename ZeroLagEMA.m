function ExtMapBuffer1 = ZeroLagEMA(Close, Period)

N = length(Close);

ExtMapBuffer1 = zeros(N,1);

Lag = floor((Period-1)/2);

Data = Close(Lag:end) + (Close(Lag:end) - Close(1:end-Lag+1));

ExtMapBuffer1(Lag:end) = EMA(Data,Period);

end