function [newGrowthRates,speciesMetList,speciesSecMatMerge,speciesAbsMatMerge,netSummary,absMetsAll,secMetsAll,absFluxesAll,secFluxesAll] = getSecAbsMetsMSAV(models,growthRates,aerobic)
% Takes metabolic models and solves FBA by minimizing sum of absolute
% values of fluxes, returns two matrices defining which metabolites are
% secreted and absorbed by each model. Returns complete list of metabolites
% secreted or absorbed.
%
% Inputs:
%    models: Struct of metabolic models
%    growthRates: Growth rates of models
% Outputs:
%    speciesMetList: List of all metabolites secreted or absorbed
%    speciesSecMatMerge: m by n matrix of m models and n metabolites
%       secreted. Entry is 1 if metabolite is secreted by model.
%    speciesAbsMatMerge: m by n matrix of m models and n metabolites
%       absorbed. Entry is 1 if metabolite is absorbed by model.
%    netSummary: A table of all metabolites secreted and absorbed by any
%       species
%   absMetsAll/secMetsAll - lists of secreted and
%       absorbed metabolites by species in columns
%
% Alan Pacheco 11/6/16, MSAV added 4/4/17

modelNames = fieldnames(models);
modelAlphabet = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','X'}; %up to 26 species...

newGrowthRates = zeros(1,length(modelNames));
speciesAbsMat = zeros(0,length(modelNames));
speciesAbsList = cell(0,0);
speciesSecMat = zeros(0,length(modelNames));
speciesSecList = cell(0,0);

[absMetsAll,secMetsAll,absFluxesAll,secFluxesAll] = deal(struct());

for i = 1:length(modelNames)
    secMets = [];
    absMets = [];
    secFluxes = [];
    absFluxes = [];
    model = models.(modelNames{i});
    
    model.lb(find(model.c)) = growthRates(i); %Place lower bound on growth rate
    
    % Add trace amount of o2 for yeast to grow anaerobically
    if strcmp(model.description,'iAZ900_noCycles_03_25_2015') && ~aerobic
        model.lb(find(ismember(model.rxns,'EX_o2(e)'))) = -0.01;
    end
    % change R. sphaeroides objective function to phototroph in anaerobic
    if strcmp(model.description,'iRsp1095') && ~aerobic
        model.c(:) = 0;
        model.c(1188) = 1;
    end
    
    FBAsoln = optimizeCbModel(model,'max','one'); %minimize taxicab norm: min |v| s.t.: S*v = b, c'v = f, lb <= v <= ub
    newGrowthRates(i) = FBAsoln.f;
    
    if ~isempty(FBAsoln.f) && FBAsoln.f > 0
    
        outRxns = intersect(find(FBAsoln.x > 1e-6), find(strncmp('EX_',model.rxns,3)));
        inRxns = intersect(find(FBAsoln.x < -1e-6), find(strncmp('EX_',model.rxns,3)));

        for j = 1:length(outRxns)
            secMet = model.mets(find(model.S(:,outRxns(j))));
            secFlux = FBAsoln.x(outRxns(j));
            if length(secMet) == 1
                if ~strcmp(secMet,'h[e]') && isExtMetab(char(secMet),'[') %don't include proton or non-external mets
                    secMets = [secMets secMet];
                    secFluxes = [secFluxes secFlux];
                end
            end
        end
        for j = 1:length(inRxns)
            absMet = model.mets(find(model.S(:,inRxns(j))));
            absFlux = FBAsoln.x(inRxns(j));
            if length(absMet) == 1
                if ~strcmp(absMet,'h[e]') && ~strcmp(absMet,'h2o[e]') && isExtMetab(char(absMet),'[') %don't include proton or non-external mets
                    absMets = [absMets absMet];
                    absFluxes = [absFluxes absFlux];
                end
            end
        end

        absMetsAll.(modelAlphabet{i}) = absMets';
        secMetsAll.(modelAlphabet{i}) = secMets';
        absFluxesAll.(modelAlphabet{i}) = absFluxes';
        secFluxesAll.(modelAlphabet{i}) = secFluxes';

        % Build species-secreted and species-absorbed matrices
        for j = 1:length(secMets)
            if sum(ismember(speciesSecList,secMets{j})) > 0
                speciesSecMat(find(ismember(speciesSecList,secMets{j})),i) = 1;
            else
                speciesSecMat(length(speciesSecList) + 1,i) = 1;
                speciesSecList{length(speciesSecList) + 1} = secMets{j};
            end
        end
        for j = 1:length(absMets)
            if sum(ismember(speciesAbsList,absMets{j})) > 0
                speciesAbsMat(find(ismember(speciesAbsList,absMets{j})),i) = 1;
            else
                speciesAbsMat(length(speciesAbsList) + 1,i) = 1;
                speciesAbsList{length(speciesAbsList) + 1} = absMets{j};
            end
        end
    end
end

%% Condense secreted and absorbed metabolite lists
speciesMetList = union(speciesSecList,speciesAbsList)';
speciesSecMatMerge = zeros(length(speciesMetList),length(modelNames));
for i = 1:length(speciesSecList)
    speciesSecMatMerge(find(ismember(speciesMetList,speciesSecList{i})),:) = speciesSecMat(i,:);
end
speciesAbsMatMerge = zeros(length(speciesMetList),length(modelNames));
for i = 1:length(speciesAbsList)
    speciesAbsMatMerge(find(ismember(speciesMetList,speciesAbsList{i})),:) = speciesAbsMat(i,:);
end

%% Summarize metabolite network in table
netSummary = cell(1,3);
netSummary(1,:) = {'Metabolite','Produced by','Consumed by'};
for i = 1:length(speciesMetList)
    if sum(speciesAbsMatMerge(i,:)) > 0 && sum(speciesSecMatMerge(i,:)) > 0 %if a metabolite is both produced and consumed
        netSummary{end+1,1} = speciesMetList{i};
        netSummary{end,2} = modelNames(find(speciesSecMatMerge(i,:)));
        netSummary{end,3} = modelNames(find(speciesAbsMatMerge(i,:)));
    end
end