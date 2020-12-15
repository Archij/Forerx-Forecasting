
warning('off')
format long

%Input Data
EURUSD_Ask = csvread('EURUSD_Renko_10_PIPS_Ticks_Ask_2003.05.05_2020.06.26.csv',1,2);
EURUSD_Bid = csvread('EURUSD_Renko_10_PIPS_Ticks_Bid_2003.05.05_2020.06.26.csv',1,2);
Open = EURUSD_Ask(:,1);
High = EURUSD_Ask(:,2);
Low = EURUSD_Ask(:,3);
Close = EURUSD_Ask(:,4);
Volume_Ask = EURUSD_Ask(:,5);
Volume_Bid = EURUSD_Bid(:,5);


%JMA

JMAs = zeros(length(Close),24);
for i=7:30
    JMAs(:,i-6) = JMA(Close,i,0);
end

JMAs = JMAs(30:end,:);

Target = Close(30:end);
diff=Target(2:end)-Target(1:end-1);

Open = Open(30:end);
High = High(30:end);
Low = Low(30:end);
Close = Close(30:end);
Volume_Ask = Volume_Ask(30:end);
Volume_Bid = Volume_Bid(30:end);

for t=100:length(Close)-30

   if diff(t+1) < 0
       Class = 0;
   elseif diff(t+1)> 0
       Class = 1;
   end
    
  Data = zeros(1,3000);
  Data(1,1:100) = Open(t-100+1:t)';
  Data(1,101:200) = High(t-100+1:t)';
  Data(1,201:300) = Low(t-100+1:t)';
  Data(1,301:400) = Close(t-100+1:t)';
  Data(1,401:500) = Volume_Ask(t-100+1:t)';
  Data(1,501:600) = Volume_Bid(t-100+1:t)';
  temp = JMAs(t-100+1:t,:);
  Data(1,601:end) = temp(:)';
 

  
  dlmwrite('eurusd_train.csv',[Data Class],'delimiter',',','precision', 16, '-append');

end


