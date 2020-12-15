function [z, peaks, peaks_indices] = ZigZag(y, percent)


    s = size(y);
    N = s(1);
    
    if s(1) < s(2)
        N = s(2);
    end

    z = zeros(N,1);
    z(1) = y(1);
    
    peaks = [];
    peaks(1) = y(1);
    
    peaks_indices = [];
    peaks_indices(1) = 1;
    
    peak_index = 1;
    number_of_candles = 0;
    candle_index = 1;
    
    trend = 0;
    
    for i=1:N
        

        if abs((y(i)-peaks(peak_index))/peaks(peak_index)*100)>=percent
            z(i) = y(i);
            number_of_candles = number_of_candles+1;
            peak_index = peak_index+1;
            peaks(peak_index) = y(i);
            peaks_indices(peak_index) = i;
            difference = abs(peaks(peak_index)-peaks(peak_index-1))/number_of_candles;
            if peaks(peak_index)>peaks(peak_index-1)
                trend = 1;
                for j=candle_index+1:candle_index+number_of_candles
                    z(j) = z(j-1)+difference;
                end
            elseif peaks(peak_index)<peaks(peak_index-1)
                trend = -1;
                 for j=candle_index+1:candle_index+number_of_candles
                    z(j) = z(j-1)-difference;
                end
            end
            
            candle_index = i;
            number_of_candles = 0;
        else
            number_of_candles = number_of_candles+1;
        end
        
        
    end


end