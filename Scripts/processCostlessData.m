function processCostlessData(saveFileNameMMatrix,resultsFileNames,conditions)

% Takes in raw results from costlessGrowth and formats information
% into 'M' matrix.
%
% INPUTS:
%    saveFileNameMMatrix: The name of the .mat file to which to save the
%    formatted data
%    resultsFileNames: The name(s) of the .mat file output by
%    costlessGrowth (e.g. resultsAerobic resultsAnaerobic)
%    conditions: The conditions used in the simulation set (e.g.
%    aer/anaer)
%
% OUTPUTS:
%    A .mat data structure containing the following fields:
%
%         Data structures containing the following results for the
%         simulations for each condition:
% 
%             allGrewAlone: Binary vectors denoting whether or not both
%             organisms grew in the first iteration of the algorithm.
% 
%             allGrewCross: Binary vectors denoting whether or not both
%             organisms grew in the last iteration of the algorithm (after
%             cross-feeding).
% 
%             expansions: Numerical vector detailing the number of medium
%             expansions that occured in each simulation
% 
%             secMetCountsAlone: Numerical vector containing the number of
%             metabolites secreted at the first iteration of each
%             simulation.
% 
%             secMetCountsCross: Numerical vector containing the number of
%             metabolites secreted at the last iteration of each
%             simulation.
% 
%             sharedMedium: Cell array detailing the minimal medium
%             composition.
% 
%             sub: Sparse binary double matrix detailing phenotypic
%             information for the first (subject) organism in each
%             simulation. Rows are simulations, column identities are
%             contained in M.col.
% 
%             par: Sparse binary double matrix detailing phenotypic
%             information for the second (partner) organism in each
%             simulation. Rows are simulations, column identities are
%             contained in M.col.
% 
%             Motifs: Data structure with two fields:
% 
%                 MotifsGen: Cell array detailing the general interaction
%                 type in each simulation (C: commensal, M: mutualistic, N:
%                 non-interacting).
% 
%                 MotifsSp: Cell array detailing the specific interaction
%                 motifs in each simulation (detailed in main text).
% 
%         col: Data structure defining which columns in the 'sub' and 'par'
%         matrices correspond to phenotypic information (Species:
%         organisms, CSources: carbon sources, absMets: absorbed
%         metabolites corresponding to M.fullMetList, secMets: secreted
%         metabolites corresponding to M.fullMetList, O2: oxygen
%         availability, growthAlone: whether the organism grew alone in
%         minimal medium, growthCross: whether the organism grew after
%         cross-feeding).
% 
%         fullMetList: Cell array of all extracellular metabolites in all
%         organism models.
% 
%         speciesPairCombos: Cell array of organisms used in each of the
%         simulations.
% 
%         CSource pair combos: Cell array of carbon sources used in each of
%         the simulations.
%
% Example 
%    processCostlessData('Results/MDataDemo',{'Results/pairGrowthDemoWithO2','Results/pairGrowthDemoNoO2'},{'aer','anaer'})
%
% Alan R. Pacheco 03/08/2017, last modified 09/20/18

%% Get general and specific motifs and add to data structure
for f = 1:length(resultsFileNames)
    resultsFileName = resultsFileNames{f};
    condition = conditions{f};
    disp(['Loading file ',resultsFileName,'...'])

    load(resultsFileName)
    disp(['Formatting data for ',condition,'...'])

    allGrewAlone = find(all(growthRatesAlone,2));
    allGrewCross = find(all(growthRatesCross,2));
    growthPercentIncrease = (length(allGrewCross)-length(allGrewAlone))/length(allGrewAlone)*100;

    % Sort by general interaction motif 
    % No Growth 'X', Non-Interacting 'N', Commensal 'C', Mutualistic 'M'

    MotifsGen = cell(length(growthRatesCross),1);
    MotifsGen(setdiff(1:length(growthRatesCross),allGrewCross)) = {'No Growth'};

    % Search for interactions by comparing secretions in S to absorption in A
    for i = 1:length(allGrewCross)
        ind = allGrewCross(i);
        currI = squeeze(I(:,:,ind));
        
        if ~isempty(find((currI), 1))
            [row,col]=find(currI);
            for m = 1:length(col) %go by m column in I matrix and find 1 (giving)
                if currI(row(m),col(m)) == 1
                    givingS = row(m); %the species that gives the shared met
                    receivingS = find(currI(:,col(m))==-1); %the receiving species

                    receivingGive = find(currI(:,currI(receivingS,:)==1)==-1);%find the species to which the receiving species gives
                    if ~isempty(intersect(givingS,receivingGive)) %If the receiving species gives back in any other way to the original giver -> mutualism
                        MotifsGen{ind} = 'M';
                    else
                        MotifsGen{ind} = 'C';
                    end
                end
            end
        else
            MotifsGen{ind} = 'N';
        end
    end

    NG = find(strncmp('No Growth',MotifsGen,2)); 
    NI = find(strncmp('N',MotifsGen,2));
    CO = find(strncmp('C',MotifsGen,2)); 
    MU = find(strncmp('M',MotifsGen,2));

    % Establish competition status. If each species only consumes its own
    % carbon source, then status is 'a'. Otherwise, status is 'b'. Also count
    % the total number of primary carbon sources they are consuming

    Competition = cell(length(growthRatesCross),1);
    nCS = zeros(length(growthRatesCross),1);

    for i = 1:length(allGrewCross)
        ind = allGrewCross(i);
        currK = squeeze(K(:,:,ind));
        if ~isempty(find(sum(currK)>1, 1))
           Competition{ind} = 'b';
        else
           Competition{ind} = 'a';
        end

        nCS(ind) = length(find(sum(currK)));
    end

    % Assign specific motifs
    MotifsSp = cell(length(growthRatesCross),1);
    MotifsSp(setdiff(1:length(growthRatesCross),allGrewCross)) = {'No Growth'};

    for i = 1:length(allGrewCross)
        ind = allGrewCross(i);
        MotifsSp{ind} = [MotifsGen{ind},num2str(nCS(ind)),Competition{ind}];
    end

    specMotifNames = unique(MotifsSp);
    specMotifCounts = zeros(length(specMotifNames),1);
    for s = 1:length(specMotifNames)
        specMotifCounts(s) = length(find(strcmp(MotifsSp,specMotifNames{s})));
    end

    % Make a data structure
    data.(condition).A = A;
    data.(condition).allGrewAlone = allGrewAlone;
    data.(condition).allGrewCross = allGrewCross;
    data.(condition).Competition = Competition;
    data.(condition).CSourcePairCombos = CSourcePairCombos;
    data.(condition).expansionFail = expansionFail;
    data.(condition).expansions = expansions;
    data.(condition).fullMetList = fullMetList;
    data.(condition).G = G;
    data.(condition).growthPercentIncrease = growthPercentIncrease;
    data.(condition).growthRatesAlone = growthRatesAlone;
    data.(condition).growthRatesCross = growthRatesCross;
    data.(condition).I = I;
    data.(condition).K = K;
    data.(condition).modelNames = modelNames;
    data.(condition).MotifsGen = MotifsGen;
    data.(condition).MotifsSp = MotifsSp;
    motifs.(condition).MotifsGen = MotifsGen;
    motifs.(condition).MotifsSp = MotifsSp;
    data.(condition).S = S;
    data.(condition).speciesPairCombos = speciesPairCombos;
    data.(condition).sharedMedium = sharedMedium;
    data.(condition).secMetCountsAlone = secMetCountsAlone;
    data.(condition).secMetCountsCross = secMetCountsCross;
    data.(condition).growthPercentIncrease = growthPercentIncrease;
    clearvars -except saveFileName data motifs resultsFileNames conditions f saveFileNameData saveFileNameMMatrix
    
end

%% Build full multidimensional "M" matrices with species, carbon sources, absmets, secmets, o2, and growth

%Define feature column sections
MHeight = length(data.(conditions{1}).speciesPairCombos);
speciesCol = [1:length(unique(data.(conditions{1}).speciesPairCombos))];
CSourceCol = [speciesCol(end)+1:speciesCol(end)+length(unique(data.(conditions{1}).CSourcePairCombos))];
absMetCol = [CSourceCol(end)+1:CSourceCol(end)+length(data.(conditions{1}).fullMetList)];
secMetCol = [absMetCol(end)+1:absMetCol(end)+length(data.(conditions{1}).fullMetList)];
O2Col = secMetCol(end)+1; %Presence or absence of oxygen
growthAloneCol = O2Col+1; %Growth of species alone
growthCrossCol = growthAloneCol+1; %Growth of species in coculture
MWidth = 1+length(speciesCol)+length(CSourceCol)+length(absMetCol)+length(secMetCol)+3;

d=data.(conditions{1}); % Use one condition just to get dimensions and indices
uniqueSpecies = unique(d.speciesPairCombos);
uniqueSpeciesLength = length(uniqueSpecies);
[speciesMatSub,speciesMatPar] = deal(zeros(length(d.speciesPairCombos),uniqueSpeciesLength));
for i = 1:uniqueSpeciesLength
    speciesMatSub(find(ismember(d.speciesPairCombos(:,1),uniqueSpecies{i})),i) = 1; %subject species
    if size(d.speciesPairCombos,2) > 1
        speciesMatPar(find(ismember(d.speciesPairCombos(:,2),uniqueSpecies{i})),i) = 1; %partner species
    end
end

for c = 1:length(conditions)
    condition = conditions{c};
    d = data.(condition);
    
    disp(['Preparing M matrix for ',condition,'...'])

    % Binary carbon source consumption
    K = squeeze(sum(d.K,4));
    if size(d.speciesPairCombos,2) > 1
        [KNewSub,KNewPar] = deal(zeros(size(K,3),size(K,2)));
        for i = 1:size(K,3)
            KNewSub(i,:) = squeeze(K(1,:,i));
            KNewPar(i,:) = squeeze(K(2,:,i));
        end
    else
        [KNewSub,KNewPar] = deal(zeros(size(K,2),size(K,1)));
        for i = 1:size(K,2)
            KNewSub(i,:) = K(:,i)';
        end
    end
    KSub = KNewSub;
    KSub(find(KSub)) = 1;
    KPar = KNewPar;
    KPar(find(KPar)) = 1;
    

    % Binary metabolite absorption
    A = squeeze(sum(d.A,4));
    if size(d.speciesPairCombos,2) > 1
        [ANewSub,ANewPar] = deal(zeros(size(A,3),size(A,2)));
        for i = 1:size(A,3)
            ANewSub(i,:) = squeeze(A(1,:,i));
            ANewPar(i,:) = squeeze(A(2,:,i));
        end
    else
        [ANewSub,ANewPar] = deal(zeros(size(A,2),size(A,1)));
        for i = 1:size(A,2)
            ANewSub(i,:) = A(:,i)';
        end
    end
    ASub = ANewSub;
    ASub(find(ASub)) = 1;
    APar = ANewPar;
    APar(find(APar)) = 1;

    % Binary metabolite secretion
    S = squeeze(sum(d.S,4));
    if size(d.speciesPairCombos,2) > 1
        [SNewSub,SNewPar] = deal(zeros(size(S,3),size(S,2)));
        for i = 1:size(S,3)
            SNewSub(i,:) = squeeze(S(1,:,i));
            SNewPar(i,:) = squeeze(S(2,:,i));
        end
    else
        [SNewSub,SNewPar] = deal(zeros(size(S,2),size(S,1)));
        for i = 1:size(S,2)
            SNewSub(i,:) = S(:,i)';
        end
    end
    SSub = SNewSub;
    SSub(find(SSub)) = 1;
    SPar = SNewPar;
    SPar(find(SPar)) = 1;

    % Binary growth: growth subject alone, growth partner alone, growth
    % subject cross, growth partner cross
    [GSub,GPar] = deal(zeros(length(d.growthRatesCross),2));
    for i = 1:length(d.growthRatesCross)
        if d.growthRatesAlone(i,1) > 0 %growth subject alone
            GSub(i,1) = 1;
        end
        if size(d.speciesPairCombos,2) > 1
            if d.growthRatesAlone(i,2) > 0 %growth partner alone
                GPar(i,1) = 1;
            end
        end
        if d.growthRatesCross(i,1) > 0 %growth subject cross
            GSub(i,2) = 1;
        end
        if size(d.speciesPairCombos,2) > 1
            if d.growthRatesCross(i,2) > 0 %growth partner cross
                GPar(i,2) = 1;
            end
        end
    end

    if strcmp(condition,'aer') % If aerobic
        if size(d.speciesPairCombos,2) > 1
            M.(condition).sub = sparse(horzcat(sparse(speciesMatSub),sparse(KSub),sparse(ASub),sparse(SSub),ones(size(S,3)),sparse(GSub)));
            M.(condition).par = sparse(horzcat(sparse(speciesMatPar),sparse(KPar),sparse(APar),sparse(SPar),ones(size(S,3)),sparse(GPar)));
        else
            M.(condition).sub = sparse(horzcat(sparse(speciesMatSub),sparse(KSub),sparse(ASub),sparse(SSub),ones(size(S,3)),sparse(GSub)));
        end
    else
        if size(d.speciesPairCombos,2) > 1
            M.(condition).sub = sparse(horzcat(sparse(speciesMatSub),sparse(KSub),sparse(ASub),sparse(SSub),sparse(zeros(size(S,3),1)),sparse(GSub)));
            M.(condition).par = sparse(horzcat(sparse(speciesMatPar),sparse(KPar),sparse(APar),sparse(SPar),sparse(zeros(size(S,3),1)),sparse(GPar)));
        else
            M.(condition).sub = sparse(horzcat(sparse(speciesMatSub),sparse(KSub),sparse(ASub),sparse(SSub),sparse(zeros(size(S,3),1)),sparse(GSub)));
        end
    end
    M.(condition).expansions = d.expansions;
    M.(condition).sharedMedium = d.sharedMedium;
    M.(condition).secMetCountsAlone = d.secMetCountsAlone;
    M.(condition).secMetCountsCross = d.secMetCountsCross;
    M.(condition).Motifs.motifsGen = d.MotifsGen;
    M.(condition).Motifs.motifsSp = d.MotifsSp;
    M.(condition).allGrewAlone = d.allGrewAlone;
    M.(condition).allGrewCross = d.allGrewCross;
end

%% Make supporting data structures for M matrix
M.col.Species = speciesCol;
M.col.CSources = CSourceCol;
M.col.absMets = absMetCol;
M.col.secMets = secMetCol;
M.col.O2 = O2Col;
M.col.growthAlone = growthAloneCol;
M.col.growthCross = growthCrossCol;

M.sta.MWidth = MWidth;
M.sta.MHeight = MHeight;

M.motifs = motifs;
M.fullMetList = d.fullMetList;

M.speciesPairCombos = d.speciesPairCombos;
M.CSourcePairCombos = d.CSourcePairCombos;

%% Save Data
disp(['Saving M Matrix data to ',saveFileNameMMatrix])
save(saveFileNameMMatrix,'M','-v7.3');
disp('Done.')