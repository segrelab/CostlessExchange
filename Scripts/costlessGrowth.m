function costlessGrowth(saveFileName,modelsFile,CSourceFile,minMedFile,numSpecies,numCSources,isAerobic,simStart,simEnd)

% Takes in a set of models and outputs n-way crossfeeding results in
% aerobic or anaerobic conditions.
%
% INPUTS:
%    saveFileName: The name of the .mat file to which to save the data
%    modelsFile: A .mat file with genome-scale metabolic models
%    numSpecies: The number of species to grow together (ex. 2 = pairwise)
%    numCSources: The number of carbon sources to use for growth
%    isAerobic: A logical indicating whether or not to include o2
%    simStart,simEnd (optional): The simulation number at which to begin
%    and end the run (for parallelization)
%
% OUTPUTS:
%    A .mat data file of name 'saveFileName' with the following contents:
%    modelNames: List of models used
%    speciesPairCombos: Enumerated combinations of species used
%    CSourcePairCombos: Enumerated combinations of medium conditions used
%    growthRatesAlone: Growth rates of the models grown alone for each
%    scenario
%    growthRatesCross: Growth rates of the models grown together for each
%    scenario
%    Crossfeed: Vector whose values are 1 if crossfeeding occured in
%    scenario 
%    K: Carbon source consumption matrix at all expansions. An SxKxCxE
%    matrix where S is the number of species in each simulation, K is the
%    number of carbon sources provided, C is the number of simulations, and
%    E is the number of expansions. Values (s,k) in K are 1 if primary
%    carbon source k is consumed by species s.
%    I: Interaction matrix at final expansion. An SxMxC matrix where S is
%    the number of species in each simulation, M is the number of all
%    exchange metabolites for all models, and C is the number of
%    simulations. -1 denotes that a metabolite was produced, and 1 denotes
%    A: Absorption matrix at all expansions An SxMxCxE matrix where S is
%    the number of species in each simulation, M is the number of all
%    exchange metabolites for all models, C is the number of simulations,
%    and E is the number of expansions.
%    S: Secretion matrix at all expansions An SxMxCxE matrix where S is
%    the number of species in each simulation, M is the number of all
%    exchange metabolites for all models, C is the number of simulations,
%    and E is the number of expansions.
%    G: Growth matrix at all expansions An SxCxE matrix where S is the
%    number of species in each simulation, C is the number of simulations,
%    and E is the number of expansions.
%    fullMetList: Enumerated all exchange metabolites for all models
%    expansions: Number of medium network expansions undertaken
%
% Example, runs simulations for all pairs of species in the threeModels
% file, with the carbon sources in CSources, the minimal medium in minMed,
% one carbon source at a time without oxygen:
%
% costlessGrowth('Results/pairGrowthDemoNoO2','Models/threeModels.mat','Medium/CSources.mat','Medium/minMed.mat',2,1,0)
%
% Notes:
%    Requires COBRA toolbox for MATLAB. Downloadable at
%    https://opencobra.github.io/cobratoolbox/
%
% Alan R. Pacheco 03/08/2017, last modified 09/20/18

%% Load files and initialize COBRA toolbox
initCobraToolbox;
changeCobraSolver('gurobi')

f=load(modelsFile);
fn=fieldnames(f);
models = f.(fn{1});
modelNames = fieldnames(models);

%% Prepare simulation parameters
% Define how many species to interact
modelAlphabet = {'A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','X'}; %up to 26 species...
nS = numSpecies;

% List of all of the exchange metabolites in all of the models
fullMetList = {};
for i = 1:length(modelNames)
    model = models.(modelNames{i});
    excRxns = find(strncmp('EX_',model.rxns,3));
    for j = 1:length(excRxns)
        fullMetList = vertcat(fullMetList,model.mets(find(model.S(:,excRxns(j)))));
    end
end
fullMetList = unique(fullMetList);

% Define a list of possible carbon sources for pairs to use
load(CSourceFile);
nCS = numCSources;
CSourceComboLength = length(nchoosek(CSources,nCS));
speciesPairCombos = repelem(nchoosek(modelNames,nS),CSourceComboLength,1);
CSourcePairCombos = nchoosek(CSources,nCS);
CSourceNamePairCombos = nchoosek(CSourceNames,nCS);
CSourcePairCombos = repmat(CSourcePairCombos,length(nchoosek(modelNames,nS)),1);
CSourceNamePairCombos = repmat(CSourceNamePairCombos,length(nchoosek(modelNames,nS)),1);
totalLength = length(speciesPairCombos);

% Load absolute minimal medium with required mets
load(minMedFile)
sharedMedium = minMed;

% Initialize data structures
[growthRatesAlone,growthRatesCross] = deal(zeros(totalLength,nS));
Crossfeed = zeros(totalLength,1);
A = zeros(nS,length(fullMetList),totalLength,10);
S = zeros(nS,length(fullMetList),totalLength,10);
K = zeros(nS,length(CSources),totalLength,10);
G = zeros(nS,totalLength,10);
I=zeros(nS,length(fullMetList),totalLength);
FS=zeros(nS,length(fullMetList),totalLength);
FA=zeros(nS,length(fullMetList),totalLength);
expansions = zeros(totalLength,1);
secMetCountsAlone = zeros(totalLength,1);
secMetCountsCross = zeros(totalLength,1);
expansionFail = zeros(totalLength,1);

% Specify simulation range
if nargin > 7
    simStart = simStart;
    simEnd = simEnd;
else
    simStart = 1;
    simEnd = totalLength;
end

%% Begin simulations
for i = simStart:simEnd
    currentModels= cell(1,nS);
    [currentCSources, currentCSourceNames] = deal(cell(1,nCS));
    modelsTest=struct();
    for j = 1:nS
        currentModels{j} = speciesPairCombos{i,j};
        modelsTest.(modelAlphabet{j}) = models.(speciesPairCombos{i,j});
    end
    for j = 1:nCS
        currentCSources{j} = CSourcePairCombos{i,j};
        currentCSourceNames{j} = CSourceNamePairCombos{i,j};
    end
    
    if nS == 1
        if nCS == 1
            disp(strcat({'Simulation '},num2str(i),{' of '},num2str(totalLength),{': '},currentModels{1},{' with '},currentCSourceNames{1}))
        else
            disp(strcat({'Simulation '},num2str(i),{' of '},num2str(totalLength),{': '},currentModels{1},{' with '},currentCSourceNames{1},{' and '},currentCSourceNames{2}))
        end
    else
        if nCS == 1
            disp(strcat({'Simulation '},num2str(i),{' of '},num2str(totalLength),{': '},currentModels{1},{' and '},currentModels{2},{' with '},currentCSourceNames{1}))
        else
            disp(strcat({'Simulation '},num2str(i),{' of '},num2str(totalLength),{': '},currentModels{1},{' and '},currentModels{2},{' with '},currentCSourceNames{1},{' and '},currentCSourceNames{2}))
        end
    end
    
    modelsTestNative = modelsTest;
    
    % Define base medium as shared medium
    baseMedium = sharedMedium;
    
    if ~isAerobic
        baseMedium(find(strcmp(baseMedium,'o2[e]'))) = [];
    end
    
    % Reset models to native medium
    modelsTest = modelsTestNative;

    % Fully unconstrain oxygen reaction if aerobic
    if isAerobic
        for j = 1:nS
            modelsTest.(modelAlphabet{j}).lb(find(ismember(modelsTest.(modelAlphabet{j}).rxns,'EX_o2(e)'))) = -1000;
        end
    end

    % Modify models with base medium and C source pairs
    testMediumAlone = vertcat(baseMedium,CSourcePairCombos{i,:});
    modelsTestAlone = defineMedium(testMediumAlone,modelsTest);
    
    % Grow and get secreted and absorbed metabolites
    [growthRatesAlone(i,:),metListAlone,secMatAlone,absMatAlone,~,absMetsAlone,secMetsAlone,absFluxesAlone,secFluxesAlone] = getSecAbsMetsMSAV(modelsTestAlone,zeros(1,nS),isAerobic);
    newMetsAlone = setdiff(metListAlone,testMediumAlone);
    
    % Record secreted and absorbed metabolites and carbon sources
    secMetCountsAlone(i) = length(newMetsAlone);
    
    for s = 1:nS
        if isfield(secMetsAlone,modelAlphabet{s})
            if ~isempty(secMetsAlone.(modelAlphabet{s}))
                S(s,find(ismember(fullMetList,secMetsAlone.(modelAlphabet{s}))),i,1) = 1;
            end
        end
        if isfield(absMetsAlone,modelAlphabet{s})
            if ~isempty(absMetsAlone.(modelAlphabet{s}))
            A(s,find(ismember(fullMetList,absMetsAlone.(modelAlphabet{s}))),i,1) = 1;
            K(s,intersect(find(ismember(CSources,currentCSources)),find(ismember(CSources,absMetsAlone.(modelAlphabet{s})))),i,1) = 1;
            end
        end
        G(s,i,1) = growthRatesAlone(i,s);
    end
    
    % If at least one of the species in the pair grows and new metabolites are secreted
    if ~isempty(find(growthRatesAlone(i,:), 1)) && ~isempty(newMetsAlone)
       
        % Perform network expansion of secreted metabolites
        expansion = 0;
        testMediumOld = testMediumAlone;
        growthRatesOld = growthRatesAlone(i,:);
        modelsTestOld = modelsTestAlone;
        metListOld = metListAlone;
        secMatOld = secMatAlone;
        absMatOld = absMatAlone;
        absMetsOld = absMetsAlone;
        secMetsOld = secMetsAlone;
        absFluxesOld = absFluxesAlone;
        secFluxesOld = secFluxesAlone;
        newMetsOld = newMetsAlone;

        while ~isempty(newMetsOld)

            expansion = expansion + 1;
            
            % Redefine the medium set for the models
            if ~iscolumn(newMetsOld)
                newMetsOld = newMetsOld';
            end
            testMediumNew = vertcat(testMediumOld,newMetsOld);
            modelsTestNew = defineMedium(testMediumNew,modelsTestOld);
            
            % Regrow and get new secreted metabolites
            [growthRatesNew,metListNew,secMatNew,absMatNew,~,absMetsNew,secMetsNew,absFluxesNew,secFluxesNew] = getSecAbsMetsMSAV(modelsTestNew,growthRatesOld,isAerobic);
            newMetsNew = setdiff(metListNew,testMediumNew);

            % If a model stops growing after the expansion, go back
            if any(growthRatesNew-growthRatesOld < -0.01) %some tolerance for stochasticity in solutions
                expansionFail(i) = 1;
                expansion = expansion - 1;
                growthRatesNew = growthRatesOld;
                modelsTestNew = modelsTestOld;
                metListNew = metListOld;
                secMatNew = secMatOld;
                absMatNew = absMatOld;
                absMetsNew = absMetsOld;
                secMetsNew = secMetsOld;
                absFluxesNew = absFluxesOld;
                secFluxesNew = secFluxesOld;
                break
            end
            
            % Record secreted and absorbed metabolites and carbon sources
            for s = 1:nS
                if isfield(secMetsNew,modelAlphabet{s})
                    if ~isempty(secMetsNew.(modelAlphabet{s}))
                        S(s,find(ismember(fullMetList,secMetsNew.(modelAlphabet{s}))),i,expansion+1) = 1;
                    end
                end
                if isfield(absMetsNew,modelAlphabet{s})
                    if ~isempty(absMetsNew.(modelAlphabet{s}))
                        A(s,find(ismember(fullMetList,absMetsNew.(modelAlphabet{s}))),i,expansion+1) = 1;
                        K(s,intersect(find(ismember(CSources,currentCSources)),find(ismember(CSources,absMetsNew.(modelAlphabet{s})))),i,expansion+1) = 1;
                
                    end
                end
                G(s,i,expansion+1) = growthRatesNew(s);
            end
            
            metListOld = metListNew;
            testMediumOld = testMediumNew;
            modelsTestOld = modelsTestNew;
            growthRatesOld = growthRatesNew;
            secMatOld = secMatNew;
            absMatOld = absMatNew;
            absMetsOld = absMetsNew;
            secMetsOld = secMetsNew;
            absFluxesOld = absFluxesNew;
            secFluxesOld = secFluxesNew;
            newMetsOld = newMetsNew; 
        end % end expansion

        expansions(i) = expansion;
        modelsPairsCross = modelsTestNew;
        growthRatesCross(i,:) = growthRatesNew;
        metListCross = metListNew;
        secMatCross = secMatNew;
        absMatCross = absMatNew;
        absMetsCross = absMetsNew;
        secMetsCross = secMetsNew;
        absFluxesCross = absFluxesNew;
        secFluxesCross = secFluxesNew;

        secMetCountCross = 0;
        if isfield(secMetsCross,'A') || isfield(secMetsCross,'B')
            if isfield(secMetsCross,'A')
                secMetCountCross = secMetCountCross + length(secMetsCross.A);
            end
            if isfield(secMetsCross,'B')
                secMetCountCross = secMetCountCross + length(secMetsCross.B);
            end
        end
        secMetCountsCross(i) = secMetCountCross;
        
        if ~isempty(find(absMatCross+fliplr(secMatCross) > 1, 1)) %if they crossfeed on combined medium
            Crossfeed(i) = 1;
        end

        % make I matrix by collapsing expansions from S and A
        if nS > 1
            ICurr = zeros(nS,length(fullMetList));
            if all(growthRatesCross(i,:)) %If they both grew
                for e = 1:expansions(i)
                    for s = 1:nS
                        secM = find(S(s,:,i,e)); %secreted metabolites by this species in this expansion
                        if ~isempty(secM)
                            tempA = A(:,secM,i,e+1:end);
                            tempA(s,:,:,:) = []; %the absorption matrix corresponding to the secM metabolites in all other species in future expansions
                            tempA = squeeze(tempA); %rows: secreted metabolites, columns: expansions. sum to get consumption

                            if size(tempA,1) ~= length(secM) %correct for wrong assignment of dimensions to make secreted metabolites rows and columns expansions
                                tempA = tempA';
                            end

                            tempA = sum(tempA,2);
                            ICurr(s,secM(find(tempA))) = 1; % Record absorption
                            ICurr(setdiff([1:nS],s),secM(find(tempA))) = -1; % Record secretion
                        end
                    end
                end
            end
            I(:,:,i) = ICurr;
        end
    end
end

% Truncate 4D matrices
e = max(expansions);
A(:,:,:,e+1:end) = [];
S(:,:,:,e+1:end) = [];
K(:,:,:,e+1:end) = [];
G(:,:,e+1:end) = [];

disp(['saving to ',saveFileName])
save(saveFileName,'Crossfeed','CSourcePairCombos','fullMetList','growthRatesAlone','growthRatesCross','A','S','G','I','K','FS','FA','modelNames','speciesPairCombos','expansions','expansionFail','secMetCountsAlone','secMetCountsCross','sharedMedium','-v7.3');

end