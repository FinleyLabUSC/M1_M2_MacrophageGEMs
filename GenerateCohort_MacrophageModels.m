%% Generate M1 and M2 models
%% Clean workspace
clear all
close all
clc
runName = %Define folder to save to, eg ['MacrophageModelingProject/20220618_AllTechniques'];
mkdir(runName);
ProtectedRxns = {'3DPHBH2';'3DSPHR';'3SALATAi';'3SPYRSP';'4HBZtm';'ACCOAC';'ACCOAC';'ACGAM2E';'ACGAM6PSi';'ACGAMPM';'ACITL';'ACNAM9PL';'ACNAMPH';'ACONT';'ACONT';'ACONTm';'ACS';'ACSm';'ACYP';'ADK1';'ADK3';'ADMDC';'ADNK1';'ADPT';'ADSK';'ADSL2';'ADSL2';'ADSS';'AHC';'AICART';'AIRCr';'AKGDm';'AKGDm';'ALASm';'ALATA_L';'ALCD2if';'ALCD2yf';'ALDD2xm';'ALDD2y';'AMANK';'AMPDA';'arachidonates';'ARGN';'ARGSL';'ARGSS';'ASP1DC';'ASPCTr';'ASPTA';'ASPTA';'ASPTAm';'B3GALT5g';'biomass_components';'biomass_components';'C14STRr';'C3STDH1r';'C3STKR2r';'C4STMO2Pr';'CAT2p';'CATm';'CBPPer';'CBPS';'CBPSam';'CBPSam';'CDIPTr';'CDIPTr';'CDS';'CEPTC;HMR_0625';'CEPTE;HMR_0614';'CHLPCTD';'CHOLK';'CITL';'CITtbm';'CLS_hs';'CMPSAS';'cofactors_vitamins';'CPPPGO';'CSm';'CTPS1';'CYSO';'CYSTGL';'DCMPDA';'DESAT16_2';'DESAT18_4';'DESAT18_5';'DESAT18_5';'DESAT18_9';'DESAT20_1';'DGULND';'DHCR243r';'DHCR71r';'DHCR72r';'DHDPBMTm';'DHORD9';'DHORTS';'DHPM1';'DM_atp_c_';'DM_atp_c_';'DMATT';'DPCOAK';'DPGase';'DPGM';'DPHMBDCm';'DPMVDc';'DURAD';'EBP2r';'ENO';'ETHAK';'ETOHMO';'F1PGT';'F1PGT';'FACOAL160i';'FACOAL160i';'FACOAL1812;HMR_0259';'FACOAL1813';'FACOAL1821';'FACOAL1831';'FACOAL204';'FAEL183';'FAOXC16080m';'FAOXC161802m';'FAOXC1811601m';'FAOXC182806m';'FAOXC183806m';'FAOXC204';'FAOXC80';'FAOXC80';'FAS180COA';'FAS180COA';'FBA';'FBA4';'FBA4';'FBP';'FDH';'FMNAT';'FORt2m';'FRDPtr';'FTHFLm';'FUM';'FUMm';'FUMm';'G3PD2m';'G5SADrm';'G5SADrm';'G5SDym';'G6PPer';'GACMTRc';'GALK';'GALT';'GALU';'GALU';'GAPD';'GARFT';'GCALDD';'GCALDD';'GCALDDm';'GF6PTA';'GHMT2r';'GK1';'GLBRAN';'GLGNS1';'GLPASE1';'GLUDxm';'GluForTx';'GLUNm';'GLUPRT';'GLUt2m';'GLXO1';'GLXO2p';'GLYAMDTRc';'GLYCK2';'GLYCLTDym';'GLYCTO1p';'GLYK';'GMPS2';'GND';'GND';'GNDer';'GNMT';'GRTT';'GUAD';'GUAPRT';'GULN3D';'GULNter';'HBZOPT10m';'HEX1';'HEX4';'HEX7';'HISD';'HISGLUr';'HMBS';'HMGCOARc;r0488';'HMGCOARc;r0488';'HMGCOASim';'HMGCOASim';'HMGLm';'HMR_0410';'HMR_0482';'HMR_0517';'HMR_0588';'HMR_0589';'HMR_0597';'HMR_0750';'HMR_2287';'HMR_2293';'HMR_2963';'HMR_4284';'HMR_5395;GGNG';'HMR_6782';'HMR_6784';'HMR_6785';'HMR_6786';'HMR_6907';'HMR_6908';'HMR_6909';'HMR_6910';'HMR_7745';'HMR_7747';'HMR_9718';'HMR_9725';'HMR_9802';'HPYRDC';'HPYRDCm';'HPYRR2x';'HPYRRy';'HXPRT';'ICDHxm';'ICDHxm';'ICDHy';'ICDHyrm';'IMPC';'IMPD';'INOSTO';'IPDDI';'IZPN';'KAS8';'KAS8';'KHK2';'KHK2';'LDH_L';'LINOFATPtc';'LNSTLSr';'LSTO1r';'MACACI';'MAN1PT2';'MAN6PI';'MDH';'MDH';'MDHm';'MDHm';'MDHx';'ME1m';'METAT';'MEVK1c';'MGSA';'MI1PP';'MI1PS';'MMEm';'MMMm';'MTHFCm';'MTHFD2m';'NDP6';'NDPK1';'NDPK2';'NDPK3';'NDPK5';'NDPK8';'NNATr';'NOS1';'NOS2;r0024';'NTD11';'OBDHc';'OCBTm';'OCBTm';'OCOAT1m';'OIVD2m';'OMPDC';'ORNDC';'ORNTArm';'ORNTArm';'ORPT';'P5CDm';'P5CDm';'P5CRm';'P5CRm';'PALFATPtc';'PDHm';'PDHm';'PDHm';'PETHCT';'PFK';'PGCD';'PGI';'PGK';'PGL';'PGLc';'PGLer';'PGM';'PGMT';'PGMT';'PGMT';'PMANM';'PMEVKc';'PNTK';'PPBNGS';'PPCDC';'PPNCL3';'PPPG9tm';'PPPGOm';'PRAGSr';'PRASCS';'PRFGS';'PRPPS';'PSERT';'PSP_L';'PSSA1_hs';'PTPAT';'PUNP5';'PYK';'r0055';'r0060';'r0068';'r0082';'r0083';'r0084';'r0097';'r0163';'r0173';'r0191';'r0210';'r0249';'r0309';'r0317';'r0383';'r0384';'r0391';'r0395';'r0407';'r0422';'r0423';'r0424';'r0425';'r0426';'r0502';'r0504';'r0509';'r0555';'r0556';'r0571';'r0571';'r0620';'r0666';'r0780';'r0781';'r0783';'r0783';'r1109';'r1380';'r1393';'r1433';'r1514';'RBFK';'RE0383C';'RE1514M';'RE2675C2';'RE2677G';'RE3103C';'RNDR1';'RNDR2';'RNDR3';'RNDR4';'RPI';'RPI';'SADT';'SARCStm';'SARDHm';'SERPT';'SMS';'SPMS';'SQLEr';'SQLSr';'steroids';'SUCD1m';'SUCD1m';'SUCD1m';'SUCOAS1m';'SUCOAS1m';'SUCOASm';'THRS';'TKT2';'TPI';'TRDR3';'TYRTA';'TYRTAm';'UAG2EMAi';'UAG4E';'UAGDP';'UDPGD';'UDPGD';'UDPGNP';'UGLT';'UMPK';'UPP3S';'UPPDC1';'UPPN';'URCN';'URIDK3';'VALTA';'VALTAm';'vitaminA';'vitaminD';'vitaminE';'xenobiotics';'XYLK';'XYLTD_Dr;r0784';'XYLTD_Dr;r0784'};

%% Perform M1 macrophage model with Imat
load('M1_expression_recon3d.mat')
% Load Model
initCobraToolbox(false) %don't update the toolbox
human = readCbModel('Recon3D.mat')
% Map genomic data to reactions
[expressionRxns parsedGPR] = mapExpressionToReactions(human, M1_Expression_recon3D)
% Make cutoffs: top 1/4, bottom 1/4, ignore -1's (where there is no
% measurement)
indices = find(expressionRxns< 0);
expressionRxns_Thresh = expressionRxns;
expressionRxns_Thresh(indices) = NaN;
threshold_lb = prctile(expressionRxns_Thresh,25);
threshold_ub = prctile(expressionRxns_Thresh,75);

%% iMAT Model
options = 'options_iMAT';
load(['options_methods' filesep options]);
expressionRxns(isnan(expressionRxns))=threshold_lb + 1;
options.expressionRxns = expressionRxns;
options.threshold_lb = threshold_lb;
options.threshold_ub = threshold_ub;
% Include core reactions (esp biomass and atp maintenance)
options.core = ProtectedRxns;
%options.runtime = 21600;
% Make Model
M1_iMAT = createTissueSpecificModel(human, options)
save([runName '/M1_iMAT.mat'], 'M1_iMAT');

%% INIT Model
options = 'options_INIT';
load(['options_methods' filesep options]);
weights = 5 .* log(expressionRxns ./ threshold_ub);
weights(find(ismember(human.rxns, ProtectedRxns))) = max(weights);
weights(isinf(weights))= -2;
weights(isnan(weights))=-2;
options.weights = weights;
options.runtime = 18000;
M1_INIT = createTissueSpecificModel(human, options)
save([runName '/M1_INIT.mat'], 'M1_INIT');

%% GIMME
human = changeObjective(human,'ATPM')
options = 'options_GIMME';
load(['options_methods' filesep options]);
options.obj_frac = 0.75;
options.threshold = threshold_lb;
options.expressionRxns = expressionRxns;
M1_GIMME_atpmain = createTissueSpecificModel(human, options);
save([runName '/M1_GIMME_atpmain.mat'], 'M1_GIMME_atpmain');

human = changeObjective(human,'BIOMASS_maintenance')
options = 'options_GIMME';
load(['options_methods' filesep options]);
options.obj_frac = 0.75;
options.threshold = threshold_lb;
options.expressionRxns = expressionRxns;
M1_GIMME_biomassmmain = createTissueSpecificModel(human, options);
save([runName '/M1_GIMME_biomassmmain.mat'], 'M1_GIMME_biomassmmain');

%% Fastcore
changeCobraSolverParams('LP', 'feasTol', 1e-8)
epsilon = 1e-6; % smallest flux considered non-zero
printLevel = 2; % summary print level
[A,modelFlipped,V] = fastcc(human, epsilon, printLevel);
remove = setdiff(1:numel(human.rxns), A);
rxnRemoveList = human.rxns(remove);
options = 'options_fastCore';
load(['options_methods' filesep options]);
base_model_FluxConsistent = removeRxns(human, rxnRemoveList);


[expressionRxns parsedGPR] = mapExpressionToReactions(base_model_FluxConsistent, M1_Expression_recon3D)
indices = find(expressionRxns< 0);
expressionRxns_Thresh = expressionRxns;
expressionRxns_Thresh(indices) = NaN;
threshold_ub = prctile(expressionRxns_Thresh,85);
%options.core = (find(ismember(base_model_FluxConsistent.rxns, ProtectedRxns)));
high_set = vertcat(base_model_FluxConsistent.rxns(find(expressionRxns > threshold_ub)), ProtectedRxns);
options.core =(find(ismember(base_model_FluxConsistent.rxns, high_set)));
M1_FastCore = createTissueSpecificModel(base_model_FluxConsistent, options);
save([runName '/M1_FastCore.mat'], 'M1_FastCore')

% CORDA
epsilon = 1e-6; % smallest flux considered non-zero
printLevel = 2; % summary print level
[A,modelFlipped,V] = fastcc(human, epsilon, printLevel);
remove = setdiff(1:numel(human.rxns), A);
rxnRemoveList = human.rxns(remove);
base_model_FluxConsistent = removeRxns(human, rxnRemoveList);
[expressionRxns parsedGPR] = mapExpressionToReactions(base_model_FluxConsistent, M1_Expression_recon3D)
indices = find(expressionRxns< 0);
expressionRxns_Thresh = expressionRxns;
expressionRxns_Thresh(indices) = NaN;
threshold_ub = prctile(expressionRxns_Thresh,75);
threshold_lb = prctile(expressionRxns_Thresh,25);
ES = {};
PR = {};
NP = {};
for i = 1:length(expressionRxns)
    temp = base_model_FluxConsistent.rxns(i);
    if expressionRxns(i) > threshold_ub || ismember(human.rxns(i), ProtectedRxns)
        ES{end+1,1} = temp{1,1};  
    elseif expressionRxns(i) < threshold_lb
        NP{end+1,1} = temp{1,1};
    else
        PR{end+1,1} = temp{1,1};
    end
end
constraint = 1e-4;

[M1_CORDA, rescue, HCtoMC, HCtoNC, MCtoNC] = CORDA(base_model_FluxConsistent, {}, ES,PR,NP, 2, constraint);
% [M1_CORDA, rescue, HCtoMC, HCtoNC, MCtoNC] = CORDA(base_model_FluxConsistent,[],high_set,medium_set,low_set, 2, constraint)
save([runName '/M1_CORDA.mat'], 'M1_CORDA');


%% M2 Model

load('M2_expression_recon3d.mat')

% Load Model
human = readCbModel('Recon3D.mat');


% Map genomic data to reactions
[expressionRxns parsedGPR] = mapExpressionToReactions(human, M2_Expression_recon3D);

% Make cutoffs
%top 1/4, bottom 1/4
indices = find(expressionRxns< 0);
expressionRxns_Thresh = expressionRxns;
expressionRxns_Thresh(indices) = NaN;
threshold_lb = prctile(expressionRxns_Thresh,25);
threshold_ub = prctile(expressionRxns_Thresh,75);

%% iMAT Model
options = 'options_iMAT';
load(['options_methods' filesep options]);
expressionRxns(isnan(expressionRxns))=threshold_lb + 1;
options.expressionRxns = expressionRxns;
options.threshold_lb = threshold_lb;
options.threshold_ub = threshold_ub;
% Include core reactions (esp biomass and atp maintenance)
options.core = ProtectedRxns;
%options.runtime = 21600;
% Make Model
M2_iMAT = createTissueSpecificModel(human, options)
save([runName '/M2_iMAT.mat'], 'M2_iMAT');
% 
%% INIT Model
options = 'options_INIT';
load(['options_methods' filesep options]);
weights = 5 .* log(expressionRxns ./ threshold_ub);
weights(find(ismember(human.rxns, ProtectedRxns))) = max(weights);
weights(isinf(weights))= -2;
weights(isnan(weights))=-2;
options.weights = weights;
options.runtime = 18000;
M2_INIT = createTissueSpecificModel(human, options);
save([runName '/M2_INIT.mat'], 'M2_INIT');
% 
%% GIMME
human = changeObjective(human,'ATPM')
options = 'options_GIMME';
load(['options_methods' filesep options]);options.obj_frac = 0.75;
options.threshold = threshold_lb;
options.expressionRxns = expressionRxns;
M2_GIMME_atpmain = createTissueSpecificModel(human, options);
save([runName '/M2_GIMME_atpmain.mat'], 'M2_GIMME_atpmain');

human = changeObjective(human,'BIOMASS_maintenance')
options = 'options_GIMME';
load(['options_methods' filesep options]);options.obj_frac = 0.75;
options.threshold = threshold_lb;
options.expressionRxns = expressionRxns;
M2_GIMME_biomassmmain = createTissueSpecificModel(human, options);
save([runName '/M2_GIMME_biomassmmain.mat'], 'M2_GIMME_biomassmmain');

%% Fastcore
epsilon = 1e-6; % smallest flux considered non-zero
printLevel = 2; % summary print level
[A,modelFlipped,V] = fastcc(human, epsilon, printLevel);
remove = setdiff(1:numel(human.rxns), A);
rxnRemoveList = human.rxns(remove);
options = 'options_fastCore';
load(['options_methods' filesep options]);

base_model_FluxConsistent = removeRxns(human, rxnRemoveList);
[expressionRxns parsedGPR] = mapExpressionToReactions(base_model_FluxConsistent, M2_Expression_recon3D)
indices = find(expressionRxns< 0);
expressionRxns_Thresh = expressionRxns;
expressionRxns_Thresh(indices) = NaN;
threshold_ub = prctile(expressionRxns_Thresh,75);
%options.core = (find(ismember(base_model_FluxConsistent.rxns, ProtectedRxns)));
high_set = vertcat(base_model_FluxConsistent.rxns(find(expressionRxns > threshold_ub)), ProtectedRxns);
options.core =(find(ismember(base_model_FluxConsistent.rxns, high_set)))
FastCore_model = createTissueSpecificModel(base_model_FluxConsistent, options);
save([runName '/M2_FastCore.mat'], 'M2_FastCore')
%% Corda
% CORDA
epsilon = 1e-6; % smallest flux considered non-zero
printLevel = 2; % summary print level
[A,modelFlipped,V] = fastcc(human, epsilon, printLevel);

remove = setdiff(1:numel(human.rxns), A);
rxnRemoveList = human.rxns(remove);
base_model_FluxConsistent = removeRxns(human, rxnRemoveList);

ES = {};
PR = {};
NP = {};
for i = 1:length(expressionRxns)
    temp = base_model_FluxConsistent.rxns(i);
    if expressionRxns(i) > threshold_ub || ismember(human.rxns(i), ProtectedRxns)
        ES{end+1,1} = temp{1,1};  
    elseif expressionRxns(i) < threshold_lb
        NP{end+1,1} = temp{1,1};
    else
        PR{end+1,1} = temp{1,1};
    end
end
constraint = 1e-4;
[M2_CORDA, rescue, HCtoMC, HCtoNC, MCtoNC] = CORDA(base_model_FluxConsistent, {}, ES,PR,NP, 2, constraint);
%[M1_CORDA, rescue, HCtoMC, HCtoNC, MCtoNC] = CORDA(base_model_FluxConsistent,[],high_set,medium_set,low_set, 2, constraint)
save([runName '/M2_CORDA.mat'], 'M2_CORDA');


 



