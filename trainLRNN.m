function [model, trbest, RMSE] = trainLRNN(XTrain, YTrain, XVal, YVal, XTest, YTest)

Train = length(YTrain);
Val = length(YVal);
Test = length(YTest);


[C R] = size(XTrain);

hidden_neurons = 1:C;

redes = cell(C,numel(hidden_neurons)*10);
RMSE = zeros(C,numel(hidden_neurons)*10);
trs =  cell(C,numel(hidden_neurons)*10);

MinRMSE = Inf;
Best_net = [];
Best_tr = [];
   
for hn=hidden_neurons
   

            net = layrecnet(1,hn);
            net.trainFcn = 'trainbr';
            net.adaptFcn = 'adaptwb';
            net.trainParam.min_grad = 0;
            net.trainParam.max_fail = 1000;
            net.trainParam.goal = 0;

            
            net.initFcn = 'initlay';
            
            for l=1:length(net.layers)
                net.layers{l}.initFcn = 'initwb';
            end
            
            %Weight initialization
            
            s = size(net.inputWeights);
            for s1=1:s(1)
             
                 for s2=1:s(2)
                 
                     net.inputWeights{s1,s2}.initFcn = 'randsmall';
                     net.inputWeights{s1,s2}.learnFcn = 'learnwh';
                 
                 end
             
            end
         
            s = size(net.layerWeights);
            for s1=1:s(1)
             
                 for s2=1:s(2)
                 
                      net.layerWeights{s1,s2}.initFcn = 'randsmall';
                      net.layerWeights{s1,s2}.learnFcn = 'learnwh';
                 
                 end
             
            end
            
            %Bias initialization
         
            s = size(net.biases);
            for s1=1:s(1)
             
                 for s2=1:s(2)
                 
                      net.biases{s1,s2}.initFcn = 'randsmall';
                      net.biases{s1,s2}.learnFcn = 'learnwh';
                 
                 end
             
            end
            
             net.divideFcn = 'divideind';
             net.divideParam.trainInd = 1:Train;
             net.divideParam.valInd = Train+1:Train+Val+1;
             net.divideParam.testInd = Train+Val+2:Train+Val+Test;

             net.divideMode = 'time';
             net.trainParam.epochs = 10000;


             net.plotFcns = {'plotperform', 'ploterrcorr', 'ploterrhist', 'plotfit', 'plotinerrcorr', 'plotregression', 'plotresponse', 'plotwb', 'plottrainstate'};
    
            % Train
            
            XtrainSeq = con2seq([XTrain XVal XTest]);
            YtrainSeq = con2seq([YTrain YVal YTest]);
            
            net.trainParam.showWindow = false;
            net.trainParam.showCommandLine = false; 

            % Do not display anything
            net.trainParam.show = NaN;
            
            [Xs,Xi,Ai,Ts] = preparets(net,XtrainSeq,YtrainSeq);
            [net, tr] = train(net,Xs,Ts,Xi,Ai,'useParallel','yes');

            Ypred=net(XTest);
            
            p = size(XTest,1);
            resultsRLR     = assessment(YTest', Ypred','regress',p);
         
                
     
                 if resultsRLR.RMSE < MinRMSE
                
                          MinRMSE = resultsRLR.RMSE;
                          Best_net = net;
                          Best_tr = tr;
                 end
                          

end

model = Best_net;
trbest = Best_tr;
RMSE = MinRMSE;

end
