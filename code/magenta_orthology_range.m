function deviations = ...
    magenta_orthology_range(phenotype_labels, ecoli_staph_orth, ...
    sigma_delta_input, magenta_model, ml_method)

    % DESCRIPTION 
    % This function calculates interaction score deviations due to 
    % species difference (when orthology information is used)
    % 
    % STEPS 
    % 1. Input processing
    % 2. Get E. coli model and orthologs
    % 3. Determine predicted variable interactions b/t species
    % 
    % Author:   Sriram Chandrasekaran
    % Created:  2018-10-23
    % Updated:  2021-05-21 (Carolina H. Chung)
    
    % I/O
    %{
    REQUIRED INPUTS: 
        1. phenotype_labels:    labels for phenotype_data (i.e. genes) 
        2. ecoli_staph_orth:    list of ortholog genes b/t species
        3. sigma_delta_input:   numeric matrix of combination profiles
        4. magenta_model:       trained MAGENTA model
    
    OPTIONAL INPUTS: 
        1. ml_method:           ML method to use for training/testing
    
    OUTPUTS:
        1. deviations:          deviations in predicted scores
    %}
  
%% INPUT PROCESSING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ADDED NEW VARIABLE TO ACCOUNT FOR MULTIPLE RF FUNCTIONS (2021-05-21)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~exist('ml_method', 'var')
        ml_method = []; 
    end

%% GET E. COLI MODEL AND ORTHOLOGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Find non-orthologs
    nic_row = [phenotype_labels; phenotype_labels];
    ix = ismember(nic_row, ecoli_staph_orth);
    nonorthtop = nic_row(~ix);
    teststaphdat = sigma_delta_input; 
    teststaphdat(~ix,:) = 0;
    % Set the sigma scores to be 2 or 0
    ix2 = ismember(phenotype_labels, nonorthtop);
    ix2 = find(ix2);
    teststaphdat1 = sigma_delta_input;
    teststaphdat1(ix2,:) = 2;
    teststaphdat1(ix2 + length(phenotype_labels),:) = 0;
    
%% DETERMINE PREDICTED VARIABLE INTERACTIONS B/T SPECIES %%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ADDED FUNCTIONALITY TO USE regRF OR fitrensemble (2021-05-21)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmpi(ml_method, 'regRF')
        testpredictions_staphchem_ecolixns2 = ...
            regRF_predict(teststaphdat1', magenta_model);
        testpredictions_staphchem_ecolixns20 = ...
            regRF_predict(sigma_delta_input', magenta_model);
        testpredictions_staphchem_ecolixns21 = ...
            regRF_predict(teststaphdat', magenta_model);
    elseif strcmpi(ml_method, 'fitrensemble')
        testpredictions_staphchem_ecolixns2 = ...
            predict(magenta_model, teststaphdat1');
        testpredictions_staphchem_ecolixns20 = ...
            predict(magenta_model, sigma_delta_input');
        testpredictions_staphchem_ecolixns21 = ...
            predict(magenta_model, teststaphdat');
    else
        error('Provide valid ML method.')
    end
    deviations_0 = testpredictions_staphchem_ecolixns20(:) - ...
        testpredictions_staphchem_ecolixns21(:); 
    deviations_1 = testpredictions_staphchem_ecolixns20(:) - ...
        testpredictions_staphchem_ecolixns2(:); 
    deviations = [deviations_0(:) deviations_1(:)];

end