function deviations = ...
    magenta_orthology_range(phenotype_labels, ecoli_staph_orth, ...
    sigma_delta_input, indigo_model)

    % DESCRIPTION 
    % This function calculates interaction score deviations due to 
    % species difference (when orthology information is used)
    % 
    % STEPS 
    % 1. Get E. coli model and orthologs
    % 2. Determine predicted variable interactions b/t species
    % 
    % Author:   Murat Cokol
    % Created:  October 23, 2018
    % Updates:  August 27, 2020 (Carolina H. Chung)
    
    % I/O
    %{
    REQUIRED INPUTS: 
        1. phenotype_labels:    labels for phenotype_data (i.e. genes) 
        2. ecoli_staph_orth:    list of ortholog genes b/t species
        3. sigma_delta_input:   numeric matrix of combination profiles
        4. magenta_model:       trained MAGENTA model
    
    OUTPUTS:
        1. deviations:          deviations in predicted scores
    %}
    
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
    testpredictions_staphchem_ecolixns2 = ...
        regRF_predict(teststaphdat1', indigo_model);
    testpredictions_staphchem_ecolixns20 = ...
        regRF_predict(sigma_delta_input', indigo_model);
    testpredictions_staphchem_ecolixns21 = ...
        regRF_predict(teststaphdat', indigo_model);
    deviations_0 = testpredictions_staphchem_ecolixns20(:) - ...
        testpredictions_staphchem_ecolixns21(:); 
    deviations_1 = testpredictions_staphchem_ecolixns20(:) - ...
        testpredictions_staphchem_ecolixns2(:); 
    deviations = [deviations_0(:) deviations_1(:)];

end