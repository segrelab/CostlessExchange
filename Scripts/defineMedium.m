function newModels  = defineMedium(definedMedium,models)

% Takes in a defined medium and tests all models for growth. If the models
% do not grow, calls getMinimalMedium and tests for minimal growth
% conditions. Adds those metabolites to returned medium. Returns models
% with bounds altered with new medium conditions.
%
% Inputs:
%    definedMedium: A cell array of external metabolites in the starting
%        medium
%    models: A struct of metabolic models
%
% Outputs:
%    newModels: A struct of the metabolic models with their LBs altered to
%        match the new medium
%
% Alan Pacheco 11/5/16, updated 8/3/17

modelNames = fieldnames(models);

for i = 1:length(modelNames)
    model = models.(modelNames{i});
    
    modelOrig = model;
    
    model.lb(find(strncmp('EX_',model.rxns,3))) = 0; % First constrain all uptake of exchange reactions
    medRxns = intersect(find(strncmp('EX_',model.rxns,3)),find(ismember(model.rxns,findRxnsFromMets(model,definedMedium))));

    for j = 1:length(medRxns)
        if modelOrig.lb(medRxns(j)) >= 0 % if a reaction was initially bounded or forced to secrete, unconstrain it
            model.lb(medRxns(j)) = -10;
        else % if a reaction already had a lower bound, set it to that
            model.lb(medRxns(j)) = modelOrig.lb(medRxns(j));
        end
    end
    
    % Upper bound correction
    model.ub(find(strncmp('EX_',model.rxns,3))) = 1000; % Unconstrain upper bounds on all exchange reactions
    noSecreteMets = {'fe2','fe3','fe','n2'};
    model.ub(intersect(find(strncmp('EX_',model.rxns,3)),find(ismember(model.rxns,findRxnsFromMets(model,noSecreteMets))))) = 0;
    
    [newModels.(modelNames{i})] = model;
end