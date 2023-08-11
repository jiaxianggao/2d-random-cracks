classdef internalV < nnet.layer.RegressionLayer % ...
        % & nnet.layer.Acceleratable % (Optional)
        
    properties
        % (Optional) Layer properties.

        % Layer properties go here.
    end
 
    methods
        function layer = internalV(name)           
            % Create a physics-based loss function
            % that includes the internal voltage 

            % Set layer name.
            layer.Name=name;
            % Set layer description.
            layer.Description = 'MSE + Internal Voltage loss';
        end

        function loss = forwardLoss(~,Y,T)
            % Return the loss between the predictions Y and the training
            % targets T.
            % Calculate MSE.

            R = size(Y,1);
            meanSquaredError = sum((Y-T).^2, 1)/R;
    
            % Take mean over mini-batch.
            N = size(Y,2);
            MSEloss = sum(meanSquaredError)/N;
            
            % Internal voltage loss goes here.

            pred_internalV = calculateInternalV(Y); 
            true_internalV = calculateInternalV(T);
             
            MSEintV = mean(mean((true_internalV-pred_internalV).^2, 1), 2);
            internal_loss = sum(MSEintV)/N;
           
            % Calculate MAE.

            % meanAbsoluteError = sum(abs(Y-T),3)/R;
    
            % Take mean over mini-batch.
            % MAEloss = sum(meanAbsoluteError)/N;

            % Add a random number
            % rand_num = rand*0;



            lambda=1;

            loss =  MSEloss+lambda*internal_loss;

        end
    end
end