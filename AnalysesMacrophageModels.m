clear all
close all
clc
%% Load Generated Models
initCobraToolbox(false) %don't update the toolbox 
fprintf(1, 'Now reading Recon3D \n')
load('Recon3d_Thermocurated.mat')
recon3d = GEMmodel;
% Specify the folder where the models are.
myFolder = ; % '/Users/patrickgelbach/Dropbox/MacrophageModelingProject/20220618_AllTechniques';
% Check to make sure that folder actually exists.  Warn user if it doesn't.
if ~isfolder(myFolder)
    errorMessage = sprintf('Error: The following folder does not exist:\n%s\nPlease specify a new folder.', myFolder);
    uiwait(warndlg(errorMessage));
    myFolder = uigetdir(); % Ask for a new one.
    if myFolder == 0
         % User clicked Cancel
         return;
    end
end
% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(myFolder, '*.mat'); % Change to whatever pattern you need.
theFiles = dir(filePattern);
for k = 1 : length(theFiles)
    baseFileName = theFiles(k).name;
    fullFileName = fullfile(theFiles(k).folder, baseFileName);
    fprintf(1, 'Now reading %s\n', baseFileName);
    % Now do whatever you want with this file name,
    % such as reading it in as an image array with imread()
    load(fullFileName);
end


%% Read in Consolidated M1 and M2 models
fprintf(1, 'Now reading M1 Model \n')
%Load M1 Model
fprintf(1, 'Now reading M2 Model \n')
% Load M2 model
% Change subsystem names back from Redgem approach to standard version
for i = 1:size(m1_model.rxns)
    ind = find(ismember(recon3d.rxns, m1_model.rxns(i)));
    m1_model.subSystems(i) = recon3d.subSystems(ind);
    m1_model.rxnNames(i) = recon3d.rxnNames(ind);
    
end
for i = 1:size(m2_model.rxns)
    ind = find(ismember(recon3d.rxns, m2_model.rxns(i)));
    m2_model.subSystems(i) = recon3d.subSystems(ind);
    m2_model.rxnNames(i) = recon3d.rxnNames(ind);
end

%% Clean up generated models
% remove unused genes
% make sure stochiometry is all correct
fprintf(1, 'Now adjusting M1 iMAT \n')
M1_iMAT = removeUnusedGenes(M1_iMAT);
M1_iMAT = removeTrivialStoichiometry(M1_iMAT);

fprintf(1, 'Now adjusting M1 GIMME, biomass \n')
M1_GIMME_biomassmmain = removeUnusedGenes(M1_GIMME_biomassmmain);
M1_GIMME_biomassmmain = removeTrivialStoichiometry(M1_GIMME_biomassmmain);

fprintf(1, 'Now adjusting M1 GIMME, atp maintenance \n')
M1_GIMME_atpmain = removeUnusedGenes(M1_GIMME_atpmain);
M1_GIMME_atpmain = removeTrivialStoichiometry(M1_GIMME_atpmain);

fprintf(1, 'Now adjusting M1 INIT \n')
M1_INIT = removeUnusedGenes(M1_INIT);
M1_INIT = removeTrivialStoichiometry(M1_INIT);

fprintf(1, 'Now adjusting M1 CORDA \n')
M1_CORDA = removeUnusedGenes(M1_CORDA);
M1_CORDA = removeTrivialStoichiometry(M1_CORDA);

fprintf(1, 'Now adjusting M1 Fastcore \n')
M1_FastCore = removeUnusedGenes(M1_FastCore);
M1_FastCore = removeTrivialStoichiometry(M1_FastCore);

fprintf(1, 'Now adjusting M2 iMAT \n')
M2_iMAT = removeUnusedGenes(M2_iMAT);
M2_iMAT = removeTrivialStoichiometry(M2_iMAT);

fprintf(1, 'Now adjusting M2 GIMME, biomass \n')
M2_GIMME_biomassmmain = removeUnusedGenes(M2_GIMME_biomassmmain);
M2_GIMME_biomassmmain = removeTrivialStoichiometry(M2_GIMME_biomassmmain);

fprintf(1, 'Now adjusting M2 GIMME, atp maintenance \n')
M2_GIMME_atpmain = removeUnusedGenes(M2_GIMME_atpmain);
M2_GIMME_atpmain = removeTrivialStoichiometry(M2_GIMME_atpmain);

fprintf(1, 'Now adjusting M2 INIT \n')
M2_INIT = removeUnusedGenes(M2_INIT);
M2_INIT = removeTrivialStoichiometry(M2_INIT);

fprintf(1, 'Now adjusting M2 CORDA \n')
M2_CORDA = removeUnusedGenes(M2_CORDA);
M2_CORDA = removeTrivialStoichiometry(M2_CORDA);

fprintf(1, 'Now adjusting M2 Fastcore \n')
M2_FastCore = removeUnusedGenes(M2_FastCore);
M2_FastCore = removeTrivialStoichiometry(M2_FastCore);


%% Load M1 and M2 sampled flux distributions
load('samplesM1.mat')
load('samplesM2.mat')

%% Clear unneeded data  
clearvars -except recon3d M1_iMAT M1_GIMME_biomassmmain M1_GIMME_atpmain M1_INIT M1_CORDA M1_FastCore M2_iMAT M2_GIMME_biomassmmain M2_GIMME_atpmain M2_INIT ...
    M2_CORDA M2_FastCore samplesM1 samplesM2 m1_model m2_model

%% Make structure to analyze all models together
models = [recon3d, M1_GIMME_atpmain, M1_iMAT, M1_INIT, M1_CORDA, M1_FastCore, M2_GIMME_atpmain, M2_iMAT, M2_INIT, M2_CORDA, M2_FastCore];
modelNames = {'recon3d'; 'M1_GIMME_atpmain'; 'M1_iMAT'; 'M1_INIT';'M1_CORDA'; 'M1_FastCore';'M2_GIMME_atpmain';'M2_iMAT'; 'M2_INIT'; 'M2_CORDA'; 'M2_FastCore'};
%% Compare Reaction Composition
% Get Reactions 
allRxns = recon3d.rxns;
clear results
results = [];
for j = 2:size(models, 2)
    for i = 1:size(allRxns, 1)
        results(i,j-1) = ismember(allRxns(i), models(j).rxns); 
    end
    j
end
% M1 reactions
M1 = results(:, 1:5);
M2 = results(:, 6:10);
cutoff = [1 2 3 4 5];
for i = 1:size(cutoff,2)
    cutoff(i)
    M1_Unique_Rxn(i)= size(setdiff(recon3d.rxns(find(sum(M1, 2) >= cutoff(i))), recon3d.rxns(find(sum(M2, 2) >= cutoff(i)))),1);
    M2_Unique_Rxn(i)= size(setdiff(recon3d.rxns(find(sum(M2, 2) >= cutoff(i))), recon3d.rxns(find(sum(M1, 2) >= cutoff(i)))),1);
    Shared_Rxn(i)= size(intersect(recon3d.rxns(find(sum(M2, 2) >= cutoff(i))), recon3d.rxns(find(sum(M1, 2) >= cutoff(i)))),1);
    M1_allRxns(i) = numel(recon3d.rxns(find(sum(M1, 2) >= cutoff(i))));
    M2_allRxns(i) = numel(recon3d.rxns(find(sum(M2, 2) >= cutoff(i))));
end
%% Compare Metabolite Composition
% Get Mets 
allMets = recon3d.mets;
clear results
results = [];
for j = 2:size(models, 2)
    for i = 1:size(allMets, 1)
        results(i,j-1) = ismember(allMets(i), models(j).mets); 
    end
    j
end
M1 = results(:, 1:5);
M2 = results(:, 6:10);
cutoff = [1 2 3 4 5];
for i = 1:size(cutoff,2)
    cutoff(i)
    M1_Unique_Mets(i)= size(setdiff(recon3d.mets(find(sum(M1, 2) >= cutoff(i))), recon3d.mets(find(sum(M2, 2) >= cutoff(i)))),1);
    M2_Unique_Mets(i)= size(setdiff(recon3d.mets(find(sum(M2, 2) >= cutoff(i))), recon3d.mets(find(sum(M1, 2) >= cutoff(i)))),1);
    Shared_Mets(i)= size(intersect(recon3d.mets(find(sum(M2, 2) >= cutoff(i))), recon3d.mets(find(sum(M1, 2) >= cutoff(i)))),1);
    M1_allMets(i) = size(recon3d.mets(find(sum(M1, 2) >= cutoff(i))),1);
    M2_allMets(i) = size(recon3d.mets(find(sum(M2, 2) >= cutoff(i))),1);
end

%% Compare Gene Composition
% Get Reactions 
allGenes = recon3d.genes;
clear results
results = [];
for j = 2:size(models, 2)
    for i = 1:size(allGenes, 1)
        results(i,j-1) = ismember(allGenes(i), models(j).genes); 
    end
    j
end
M1 = results(:, 1:5);
M2 = results(:, 6:10);
for i = 1:size(cutoff,2)
    cutoff(i)
    M1_Unique_Genes(i)= size(setdiff(recon3d.genes(find(sum(M1, 2) >= cutoff(i))), recon3d.genes(find(sum(M2, 2) >= cutoff(i)))),1);
    M2_Unique_Genes(i)= size(setdiff(recon3d.genes(find(sum(M2, 2) >= cutoff(i))), recon3d.genes(find(sum(M1, 2) >= cutoff(i)))),1);
    Shared_Genes(i)= size(intersect(recon3d.genes(find(sum(M2, 2) >= cutoff(i))), recon3d.genes(find(sum(M1, 2) >= cutoff(i)))),1);
    M1_allGenes(i) = size(recon3d.genes(find(sum(M1, 2) >= cutoff(i))),1);
    M2_allGenes(i) = size(recon3d.genes(find(sum(M2, 2) >= cutoff(i))),1);
end

%% M1 and M2 models have been generated with REDGEM. Basic details of models
% Cutoff of 3 models

%% Subsystem analysis
subsystems = unique(recon3d.subSystems);
subsystem_Models_M1 =[];
subsystem_Models_M2 =[];
M1_notPresentSubs = setdiff(subsystems,unique(m1.subSystems));
M2_notPresentSubs = setdiff(subsystems,unique(m2.subSystems));
M1_PresentSubs = intersect(subsystems,unique(m1.subSystems));
M2_PresentSubs = intersect(subsystems,unique(m2.subSystems));
% How many subsystems, and how many reactions in recon3d
[uniqueXX, ~, J]=unique(recon3d.subSystems);
occ = histc(J, 1:numel(uniqueXX));
% 104 subsystems
% For figure with % bar charts
for i = 1:size(uniqueXX,1)
    i
M1_rxns = m1_model.rxns(find(contains(m1_model.subSystems, uniqueXX(i))))
M2_rxns = m2_model.rxns(find(contains(m2_model.subSystems, uniqueXX(i))))
recon3dRxns = recon3d.rxns(find(contains(recon3d.subSystems, uniqueXX(i))))
count = length(recon3dRxns);
SubSystem_int(i) = length(intersect(M1_rxns,M2_rxns)) /count
Subsystem_M1spec(i) = length(setdiff(M1_rxns,M2_rxns)) /count
Subsystem_M2spec(i) = length(setdiff(M2_rxns,M1_rxns)) /count
Neither(i) = 1 -(SubSystem_int(i)+Subsystem_M1spec(i)+Subsystem_M2spec(i)) 
end

% FEA 
% M1 vs Recon3d
a = find(contains(recon3d.rxns, m1_model.rxns));
resultCellM1 = FEA(recon3d, a, 'subSystems');
a = find(contains(recon3d.rxns, m2_model.rxns));
resultCellM2 = FEA(recon3d, a, 'subSystems');

%% Sampled fluxes
% M1
allRxns = recon3d.rxns;
M1_samples_recon3d = [];
for i = 1:size(allRxns, 1)
    if ismember(allRxns(i), m1_model.rxns)
        M1_samples_recon3d(i,:) = samplesM1(find(strcmp(m1_model.rxns, allRxns(i))), :);
    else
         M1_samples_recon3d(i,:) = NaN(1,50000,'single');
    end
    if mod(i,100) == 0
        fprintf('M1 samples: At iteration %d...\n',i);
    end
end
%M2
M2_samples_recon3d = [];
for i = 1:size(allRxns, 1)
    if ismember(allRxns(i), m2_model.rxns)
        M2_samples_recon3d(i,:) = samplesM2(find(strcmp(m2_model.rxns, allRxns(i))), :);
    else
         M2_samples_recon3d(i,:) = NaN(1,50000,'single');
    end
    if mod(i,100) == 0
        fprintf('M1 samples: At iteration %d...\n',i);
    end
end

%% Plot all Sampled Fluxes
aa = find(contains(ihuman.rxns,Shared_Rxn(:)));
input = aa;
nrows = 6;
ncols = 6;
nplts=nrows*ncols;
nfigs=ceil(length(input)/nplts);
h=zeros(nplts,nfigs);              % hang onto the handles...
k=1;
for j=1:nfigs
    figure(j)
    for i = 1:nplts
        h(i,j)=subplot(nrows, ncols, i);
        x = samplesM1(k,:);
        edges=linspace(-100,100,10);
        histogram(x, edges, 'FaceColor', '#0000FF')
        hold on
        y = samplesM2(k,:);
        histogram(y, edges, 'FaceColor', '#FF0000')
        title(gca,string(m1_model.rxns((aa(k)))));
        hold off
        k=k+1;
    end
end

%% Sampled fluxes: KL Divergence and Discrepency Scores
a = intersect(m1_model.rxns, m2_model.rxns);
normalized_sampleM1 = [];
normalized_sampleM2 = [];
edges=linspace(-1000,1000,100);
for i = 1:50000
   normalized_sampleM1(:,i) = samplesM1(:,i) ./ (sum(abs(samplesM1(:,i)),1, 'omitnan'));
   normalized_sampleM2(:,i) = samplesM2(:,i) ./ (sum(abs(samplesM2(:,i)),1, 'omitnan'));  
end
AveDist = [];
for i = 1:length(a)
    if mod(i,100) == 0
        fprintf('Samples: At iteration %d...\n',i);
    end
    m1Rxni = find(ismember(m1_model.rxns, a(i)));
    m2Rxni = find(ismember(m2_model.rxns, a(i)));
    [f1,x1]= ksdensity(samplesM1(m1Rxni,:));
    [f2,x2]= ksdensity(samplesM2(m2Rxni,:));
    dist1 = KLDis(f1,f2);
    dist2 = KLDis(f2,f1);
    AveDist(i) = (dist1+dist2) /2;
    
    
    testM1 = histcounts(samplesM1(m1Rxni,:), edges, 'Normalization', 'probability');
    testM2 = histcounts(samplesM2(m2Rxni,:), edges, 'Normalization', 'probability');
    discrepency(i) = sum(abs(testM1 - testM2), 'omitnan');
    
end

m1Rxns = m1_model.subSystems(find(ismember(m1_model.rxns, a)));
[aa,~,c] = unique(m1Rxns(:,1));
b = [aa, accumarray(c,AveDist(:),[],@mean)];
    
%% Flux sum for each metabolite: From Samples
Shared_Mets = intersect(m1_model.mets, m2_model.mets);
M1_samplingResult = samplesM1;
M2_samplingResult = samplesM2;
for i = 1:size(Shared_Mets)
    met = Shared_Mets(i);
    met_location = find(ismember(m1_model.mets,Shared_Mets(i)));
    indices = find(m1_model.S(met_location,:));
    for d = 1:500 %size(M1_samplingResult,2)
        total1(i, d) = 0.5 * sum( abs(M1_samplingResult(indices, d)), 'omitnan');
    end
    if mod(i,100) == 0
        fprintf('M1 samples: At iteration %d...\n',i);
    end
end
  M1mean =  mean(total1,2);
  M1stdev = std(total1, 0, 2);
  
for i = 1:size(Shared_Mets)
    met = Shared_Mets(i);
    met_location = find(ismember(m2_model.mets,Shared_Mets(i)));
    indices = find(m2_model.S(met_location,:));
    for d = 1:500 %size(M2_samplingResult,2)
        total2(i, d) = sum( abs(M2_samplingResult(indices, d)), 'omitnan');
    end
    if mod(i,100) == 0
        fprintf('M2 samples: At iteration %d...\n',i);
    end
end
  M2mean =  mean(total2,2);
  M2stdev = std(total2, 0, 2);
  
a = intersect(m1_model.mets, m2_model.mets);
  for i = 1:50000

      [Pi,Ci,vPi,vCi,si] = computeFluxSplits(m1_model,a,M1_samplingResult(:, i));
      P1(:,i) = Pi;
      C1(:,i) = Ci;
      vP1(:,i) = vPi;
      vC1(:,i) = vCi;
      s1(:,i) = si;
      if mod(i,100) == 0
          fprintf('M1 samples: At iteration %d...\n',i);
      end
  end
  for i = 1:50000

      [Pi,Ci,vPi,vCi,si] = computeFluxSplits(m2_model,a,M2_samplingResult(:, i));
      P2(:,i) = Pi;
      C2(:,i) = Ci;
      vP2(:,i) = vPi;
      vC2(:,i) = vCi;
      s2(:,i) = si;
      if mod(i,100) == 0
          fprintf('M2 samples: At iteration %d...\n',i);
      end
  end
 
  

%% Flux through each pathway- normalized
% M1
normalized_sampleM1 = [];
normalized_sampleM2 = [];
for i = 1:50000
    normalized_sampleM1(:,i) = samplesM1(:,i) ./ (sum(abs(samplesM1(:,i)),1, 'omitnan'));
    normalized_sampleM2(:,i) = samplesM2(:,i) ./ (sum(abs(samplesM2(:,i)),1, 'omitnan'));
end
for i = 1:size(unique(m1_model.subSystems),1)
    pathway = unique(m1_model.subSystems(i));
    indices = find(ismember(m1_model.subSystems, pathway));
    for j = 1:size(samplesM1,2)
        pathA = samplesM1(indices, j);
        a = size(pathA,1);
        fluxRxn_aveM1(i,j) = sum(pathA, 'omitnan');
    end
    %fluxM1n(i,1) = mean(fluxRxn_aveM1, 2);
    %fluxM1n(i,2) = std(fluxRxn_aveM1);
end
%M2
for i = 1:size(unique(m2_model.subSystems),1)
    pathway = unique(m2_model.subSystems(i));
    indices = find(ismember(m2_model.subSystems, pathway));
    for j = 1:size(samplesM2,2)
        pathA = samplesM2(indices, j);
        a = size(pathA,1);
        fluxRxn_aveM2(i,j) = sum(pathA, 'omitnan');
    end
    %fluxM2n(i,1) = mean(fluxRxn_aveM2,2);
    %fluxM2n(i,2) = std(fluxRxn_aveM2);
end
intSubs = intersect(m1_model.subSystems, m2_model.subSystems);
indM1 = find(ismember(unique(m1_model.subSystems), intSubs));
indM2 = find(ismember(unique(m2_model.subSystems), intSubs));
m1sSub = mean(fluxRxn_aveM1(indM1,:),2);
m2sSub = mean(fluxRxn_aveM2(indM2,:),2);
m1stSub = std(fluxRxn_aveM1(indM1,:),[],2);
m2stSub = std(fluxRxn_aveM2(indM2,:),[],2);

%% 
intRxns = intersect(m1_model.rxns, m2_model.rxns);
indM1 = find(ismember(m1_model.rxns, intRxns));
indM2 = find(ismember(m2_model.rxns, intRxns));
%pat
intSubs = intersect(m1_model.subSystems, m2_model.subSystems);
indM1 = find(ismember(unique(m1_model.subSystems), intSubs));
indM2 = find(ismember(unique(m2_model.subSystems), intSubs));
D = zeros(size(fluxRxn_aveM1, 2), size(fluxRxn_aveM2, 2));
for i = 1:size(fluxRxn_aveM1, 2)
       if mod(i,100) == 0
        fprintf('At M1 iteration %d...\n',i);
       end
    for j = 1:size(fluxRxn_aveM2, 2)
        G = fluxRxn_aveM1(indM1, i);
        G2 = fluxRxn_aveM2(indM2, j);
        D(i,j) = norm(G-G2);
    end
end


D = zeros(size(M1_samples, 2), size(M2_samples, 2));
for i = 1:size(M1_samples, 2)
    for j = 1:size(M2_samples, 2)
        G = M1_samples(indM1, i);
        G2 = M2_samples(indM2, j);
        D(i,j) = norm(G-G2);
    end
end
baseDistance = mean(D, 'all');
save('baseDistance.mat', 'baseDistance')

%% PageRank percentile
pagerankM1 = [];
allPercents = zeros(2*size(recon3d.rxns,1),1);
for j = 1:size(M1models, 2)
    model1 = M1models(j);
    try
        model1 =  changeObjective(model1,'MAR03993');
        b1 = optimizeCbModel(model1);
        size(b1.v)
        M1 = MFG_cancer(model1, b1.v);
        v = b1.v;
        %pg1 = pagerank(M1);
        indc = find(ismember( model1.rxns, Shared_Rxn));
        M1 = M1(indc, indc);
        pga = pagerank(M1);
        %pga = pg1(indc);
        v_plus = (0.5.*(abs(v) + v))';
        v_minus = (0.5.*(abs(v) - v))';
        v2m = [v_plus v_minus]';
        [percentile_test] = calculate_percentile(pga);
%         indices = find(ismember(ihuman.rxns, model1.rxns(ind)));
%         indices = [indices; indices+size(indices, 1)];
        for i = 1:size(percentile_test,1)
            if indc(i) <= size(model1.rxns,1)
               indexIhuman =  find(ismember(recon3d.rxns, model1.rxns(indc(i))));
               allPercents(indexIhuman,j) =  percentile_test(i);
               pagerankM1(indexIhuman,j) = pga(i);
            else
               indexIhuman =  find(ismember(recon3d.rxns, model1.rxns(indc(i-(size(model1.rxns,1))))));
               allPercents(indexIhuman+13078,j) =  percentile_test(i);
               pagerankM1(indexIhuman+13078,j) =  pga(i);
            end
        end
    end
end

pagerankM2 = [];
allPercentsM2 = zeros(2*size(ihuman.rxns,1),1);
for j = 1:size(M2models, 2)
    model1 = M2models(j);
    try
        model1 =  changeObjective(model1,'biomass_maintenance');
        b1 = optimizeCbModel(model1);
        size(b1.v)
        M1 = MFG_cancer(model1, b1.v);
        v = b1.v;
        %pg1 = pagerank(M1);
        indc = find(ismember( model1.rxns, Shared_Rxn));
        %pgb = pg1(indc);
        M1 = M1(indc, indc);
        pgb = pagerank(M1);
        [percentile_test] = calculate_percentile(pgb);
%         indices = find(ismember(ihuman.rxns, model1.rxns(ind)));
%         indices = [indices; indices+size(indices, 1)];
        for i = 1:size(percentile_test,1)
            if indc(i) <= size(model1.rxns,1)
               indexIhuman =  find(ismember(recon3d.rxns, model1.rxns(indc(i))));
               allPercentsM2(indexIhuman,j) =  percentile_test(i);
               pagerankM2(indexIhuman,j) = pgb(i);
            else
               indexIhuman =  find(ismember(recon3d.rxns, model1.rxns(indc(i-(size(model1.rxns,1))))));
               allPercentsM2(indexIhuman+13078,j) =  percentile_test(i);
               pagerankM2(indexIhuman+13078,j) =  pgb(i);
            end
        end
    end
end


allPercentsM2 = zeros(2*size(recon3d.rxns,1),1);
M2models = models(7:11);
for j = 1:size(M2models, 2)
    model1 = M2models(j);
    try
        model1 =  changeObjective(model1,'biomass_maintenance');
        b1 = optimizeCbModel(model1);
        size(b1.v)
        M1 = MFG_cancer(model1, b1.v);
        pg1 = pagerank(M1);
        [percentile_test] = calculate_percentile(pg1)
        indices = find(ismember(recon3d.rxns, model1.rxns));
        indices = [indices; indices+size(indices, 1)];
        for i = 1:size(percentile_test,1)
            allPercentsM2(indices(i),j) =  percentile_test(i);
        end
    end
end
ind = find(contains(ihuman.rxns, Shared_Rxn));
ind = [ind;ind+size(ihuman.rxns, 1)];
p1 = allPercents(ind, :);
p2 = allPercentsM2(ind, :);



pM2 = pagerankM2;
pM1 = pagerankM1;
pM2(pM2==0) = NaN;
pM1(pM1==0) = NaN;
med1 = mean(pM1,2, 'omitnan');
med2 = mean(pM2,2, 'omitnan');
ploop = horzcat(med1, med2);
ploop(any(isnan(ploop), 2), :) = [];
med1 = ploop(:,1);
med2 = ploop(:,2);
[percentile_test1] = calculate_percentile(med1);
[percentile_test2] = calculate_percentile(med2);


pM2 = allPercentsM2;
pM1 = allPercents;
pM2(pM2==0) = NaN;
pM1(pM1==0) = NaN;
med1 = nan(pM1,2, 'omitnan');
med2 = mean(pM2,2, 'omitnan');

tic
for i = 1:5000
    n = randi([ 1 50000]);
    [i n]
    m1s = M1_allSamples(:,n);
    m2s = M2_allSamples(:,n);
    m1s(isnan(m1s)) = 0;
    m2s(isnan(m2s)) = 0;
    M1 = MFG_cancer(ihuman, m1s);
    M2 = MFG_cancer(ihuman, m2s);
    pg1 = pagerank(M1);
    pg2 = pagerank(M2);
    [percentile_M1] = calculate_percentile(pg1);
    [percentile_M2] = calculate_percentile(pg2);
    percM1(:,i) = percentile_M1;
    percM2(:,i) = percentile_M2;
    pageM1(:,i) = pg1;
    pageM2(:,i) = pg2;      
end
toc


for i = 1:5000 
    aPage1(:,i) = pageM1(:, i);
    aPage2(:,i) = pageM2(:, i);
    [percentile_M1] = calculate_percentile(aPage1(:,i));
    [percentile_M2] = calculate_percentile( aPage2(:,i));
    apercM1(:,i) = percentile_M1;
    apercM2(:,i) = percentile_M2;
end

