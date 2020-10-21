function [train_interactions, trainxnscores, phenotype_labels, ...
    magenta_model, sigma_delta_scores, conditions] = ...
    MAGENTA_train(interaction_filename, annotation_filename, ...
    chemogenomics_filename, z, phenotype_data, phenotype_labels, ...
    conditions, interaction_scores, interaction_pairs)

    % DESCRIPTION 
    % This function constructs a MAGENTA model using regression-based
    % random forests. 
    % 
    % STEPS 
    % 1. Input processing
    % 2. Convert and match interaction labels with chemogenomic data
    % 3. Filter out interactions without chemogenomic data
    % 4. Define inputs to MAGENTA and train model
    % 
    % Author:   Murat Cokol
    % Created:  October 23, 2018
    % Updates:  August 27, 2020 (Carolina H. Chung)
    
    % I/O
    %{
    REQUIRED INPUTS: 
        1. interaction_filename:    filename for drug interactions
        2. annotation_filename:     filename for matching drug names to 
                                    chemogenomic condition names
        3. chemogenomics_filename:  filename for chemogenomic data
    OPTIONAL INPUTS: 
        1. z:                   threshold value to define significant 
                                effect on fitness (default: z = 2)
        2. phenotype_data:      binary matrix defining sensitive and 
                                resistant phenotypes
        3. phenotype_labels:    labels (i.e. genes) for phenotype data
        4. conditions:          list of conditions (i.e. treatments)
        5. interaction_scores:  numerical array of drug interaction scores
        6. interaction_pairs:   cell array of interaction pair names
    
    OUTPUTS:
        1. train_interactions:  interaction pairs used for model training
                                (i.e. chemogenomic data available)
        2. trainxnscores:       interaction scores corresponding to 
                                train_interactions
        3. phenotype_labels:    labels (i.e. genes) for phenotype data
        4. magenta_model:       trained MAGENTA model 
        5. sigma_delta_scores:  numerical matrix of combination profiles
        6. conditions:          list of conditions (i.e. treatments) in 
                                chemogenomic data
    %}

%% INPUT PROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~exist('interaction_scores','var') || isempty(interaction_scores)
        [interaction_scores, interaction_pairs] = ...
            xlsread(interaction_filename);
    end
    drugs_all = unique(interaction_pairs);
    if ~exist('z','var') || isempty(z)
        z = 2;
    end
    if ~exist('phenotype_data','var') || isempty(phenotype_data)
        [phenotype_data, phenotype_labels, conditions] = ...
            process_chemgen_v2(chemogenomics_filename, z);
    end
    
%% CONVERT AND MATCH INTERACTION LABELS WITH CHEMOGENOMIC DATA %%%%%%%%%%%%
    [~, txt] = xlsread(annotation_filename);
    [drugxn_id, chemgen_id] = deal(txt(:,1), txt(:,2));
    drugnames_cell = drugs_all;
    for i = 1:length(drugxn_id)
        drugnames_cell(ismember(drugnames_cell, ...
            drugxn_id(i))) = chemgen_id(i);
    end
    drugpairsname_cell = interaction_pairs;
    for i = 1:length(drugxn_id)
        drugpairsname_cell(ismember(drugpairsname_cell, ...
            drugxn_id(i))) = chemgen_id(i);
    end

%% FILTER OUT INTERACTIONS WITHOUT CHEMOGENOMIC DATA %%%%%%%%%%%%%%%%%%%%%%
    ix = ismember(drugpairsname_cell, conditions); 
    ix = (sum(~cellfun(@isempty, drugpairsname_cell), 2)) == sum(ix, 2);
    train_interactions = drugpairsname_cell(ix,:); 
    trainxnscores = interaction_scores(ix);

%% DEFINE INPUTS FOR MAGENTA AND TRAIN MODEL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    traindrugs = unique(train_interactions(:));
    traindrugs(cellfun(@isempty, traindrugs)) = []; 
    [~, pos] = ismember(traindrugs,conditions); 
    trainchemgen = phenotype_data(:,pos);
    [~, magenta_model, sigma_delta_scores] = ...
        magenta_rf_3(traindrugs, trainchemgen, train_interactions, ...
        trainxnscores, [], [], [], 1, []);
 
end