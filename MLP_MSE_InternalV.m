%MLP model for 2 random cracks

clear

% Load data


load("Flextural_Crack_Train.mat")
load("Flextural_Crack_Validation.mat")
load("Node_injection_Crack.mat")

% build a model

layers = [...
    sequenceInputLayer(3024)
    fullyConnectedLayer(3000,'WeightsInitializer','he')
    reluLayer
    dropoutLayer(0.5)
    fullyConnectedLayer(3000,'WeightsInitializer','he')
    reluLayer
    dropoutLayer(0.5)
    fullyConnectedLayer(3000,'WeightsInitializer','he')
    reluLayer
    dropoutLayer(0.5)
    fullyConnectedLayer(1406,'WeightsInitializer','he')
    internalV('mse+intV')];


options = trainingOptions("adam", ...
    MaxEpochs=1000, ...
    InitialLearnRate=0.001, ...
    MiniBatchSize=20, ...
    ValidationData={XValidation,YValidation}, ...
    ValidationFrequency=10, ...
    L2Regularization=0, ...
    ExecutionEnvironment="cpu");

[net, trainInfo] = trainNetwork(XTrain, YTrain, layers, options);
% save('D:\matlabwork\crack-liang\neural networks\adam_rd\net_MSEonly.mat','net_MSEonly')

% save('D:\matlabwork\crack-liang\neural networks\adam_rd\trainInfo_MAE.mat',"trainInfo")
% 
% new_sigma=predict(net, XValidation(:,11));

% figure(1),clf, PlotSolScaled(g,H,new_sigma,min(new_sigma),max(new_sigma)),
% colormap(hot),
% axis image,colorbar,drawnow
% plot(trainInfo.TrainingLoss)