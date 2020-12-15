F = csvread('C:\Users\Artur\Dropbox\Input.csv',1,0);

N = length(F);

Tr = ceil(0.7*N);

XTr = F(1:Tr,3:end);
YTr = F(1:Tr,2);

XTe = F(Tr+1:end,3:end);
YTe = F(Tr+1:end,2);

Best_OA = 0;
Best_k = 0;
for k=1:20
    Mdl = fitcknn(XTr,YTr,'CategoricalPredictors',[],'Distance','euclidean','DistanceWeight','equal','NumNeighbors',k);
    Yp = predict(Mdl,XTe);
    Results = assessment(YTe,Yp,'class');
    
    if Results.OA > Best_OA
        Best_OA = Results.OA;
        Best_k = k;
    end
end 

 Mdl = fitcknn(XTr,YTr,'CategoricalPredictors',[],'Distance','euclidean','DistanceWeight','equal','NumNeighbors',Best_k);
 Yp = predict(Mdl,XTe);
 RESULTS = assessment(YTe,Yp,'class');
 RESULTS.OA