function [testxnscores, magenta_model, sigma_delta_scores] = ...
    magenta_rf_3(traindrugs, trainchemgen, train_interactions, ...
    trainxnscores, testdrugs, testchemgen, test_interactions, ...
    indigo_mode, magenta_model)

    % DESCRIPTION 
    % This function implements the Random Forests algorithm to either train
    % or test a machine learning model predictive of drug combination 
    % effect. 
    % 
    % STEPS 
    % 1. Input processing based on mode (training or testing)
    % 2. Run RF algorithm
    % 
    % Author:   Murat Cokol
    % Created:  October 23, 2018
    % Updates:  October 21, 2020 (Carolina H. Chung)
    %           September 17, 2020 (Carolina H. Chung)
    %           August 27, 2020 (Carolina H. Chung)
    
    % I/O
    %{
    REQUIRED INPUTS: 
        1. traindrugs:          drug list for training set
        2. trainchemgen:        chemogenomic data for training set
        3. train_interactions:  interaction names for training set
        4. trainxnscores:       interaction scores for training set
        5. testdrugs:           drug list for testing set
        6. testchemgen:         chemogenomic data for testing set
        7. test_interactions:   interaction names for testing set
        8. indigo_mode:         numeric value specifying mode for model
                                development 
                                (training:  indigo_mode = 1, 
                                 testing:   indigo_mode = 2)
        9. magenta_model:       MAGENTA model (specify if mode = 2)
    
    OUTPUTS:
        1. testxnscores:        predicted interaction scores for test set
        2. magenta_model:       trained MAGENTA model 
        3. sigma_delta_scores:  numerical matrix of combination profiles
    %}

%% INPUT PROCESSING BASED ON MODE (TRAINING OR TESTING) %%%%%%%%%%%%%%%%%%%
   
    % Determine sigma and delta scores for interaction combinations 
    if (indigo_mode == 1) || (indigo_mode == 0)     % training mode
        chemgen = trainchemgen;             % chemogenomic data
        alldrugs = traindrugs;              % drug list
        interactions = train_interactions;  % interaction names
        for i = 1:size(interactions,1)
            ix1 = find(ismember(alldrugs, interactions(i,:)));
            if numel(ix1) == 1
                ix1 = repelem(ix1, size(interactions, 2)); 
            end
            t1 = chemgen(:,ix1);
            t2 = sum(t1,2) * 2/length(ix1); % sigma scores
            traindiffdat1xxz2(:,i) = [t2; [(sum(logical(t1')) ==1)]'];
        end
    elseif (indigo_mode == 2) || (indigo_mode == 0) % testing mode
        chemgen = testchemgen;              % chemogenomic data
        alldrugs = testdrugs;               % drug list
        interactions = test_interactions;   % interaction names
        for i = 1:size(interactions,1)
            ix1 = find(ismember(alldrugs, interactions(i,:)));
            if numel(ix1) == 1
                ix1 = repelem(ix1, size(interactions, 2)); 
            end
            t1 = chemgen(:,ix1);
            t2 = sum(t1,2) * 2/length(ix1); % sigma scores
            testdiffdat1xxz2(:,i) = [t2;[(sum(logical(t1')) ==1)]'];
        end
    end

%% RUN RF ALGORITHM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if indigo_mode == 1
        magenta_model = regRF_train(traindiffdat1xxz2', trainxnscores);
        sigma_delta_scores = traindiffdat1xxz2;
        testxnscores = [];
    elseif indigo_mode == 2
        testxnscores = regRF_predict(testdiffdat1xxz2', magenta_model);
        sigma_delta_scores = testdiffdat1xxz2;
    else
        error("Value for 'indigo_mode' is invalid.")
    end

end