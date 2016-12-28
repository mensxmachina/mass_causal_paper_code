% Type the project directory here
mainDir = 'C:\Users\striant\Documents\MATLAB\MC';

fcsDir = [mainDir filesep '_fcs' filesep 'BM'];
matDir =  [dataDir filesep '_mat' filsep 'BM'];
if ~isdir(matDir);mkdir(matDir);end
codeDir = [mainDir filesep '_code'];
addpath(genpath(codeDir));
cd(mainDir);

scrTic = tic;

% Donors: fsc file donor ids and mat file donor ids
donorIDs = {'Marrow1', 'Marrow2'};
donors = {'M2', 'M2'};

% Activators: fsc file activator ids and mat file activator ids
activatorIDs = {'BCR', 'Flt3L', 'GCSF', 'GMCSF', 'IFNa', 'IL3', 'IL7', 'LPS', 'PMAiono', 'PVO4','SCF', 'TNFa', 'TPO'};
activators = {'BCR', 'Flt3L', 'GCSF', 'GMCSF', 'IFNa', 'IL3', 'IL7', 'LPS', 'PMA-IONO', 'PVO4','SCF', 'TNFa', 'TPO'};

referenceIDs = {'Basal1', 'Basal2', 'Basal3', 'Basal4', 'Basal5'};
reference = {'Ref1', 'Ref2', 'Ref3', 'Ref4', 'Ref5'};

% Populations: fsc file population ids and mat file population ids
populationIDs= {'CD11b- Monocyte', 'CD11bhi Monocyte', 'CD11bmid Monocyte', 'NK',  'Mature CD4+ T', 'Naive CD4+ T', 'Naive CD8+ T', ...
    'Mature CD8+ T', 'Pre-B II', 'Mature CD38lo B', 'Pre-B I', 'Mature CD38mid B', 'Immature B', 'Plasmacytoid DC', 'HSC', 'MPP', 'CMP',...,
    'GMP', 'MEP', 'Erythroblast', 'Megakaryocyte', 'Platelet', 'Myelocyte', 'Plasma Cell'};

populations = {'monocyte_cd11b-', 'monocyte_cd11bhi', 'monocyte_cd11bmid', 'nk',  'mature_cd4+t', 'naive_cd4+t', 'mature_cd8+t', ...
    'naive_cd8+t', 'pre-bII', 'pre-bI', 'mature_CD38lo_b', 'mature_CD38mid_b', 'immature_b', 'dc', 'hsc', 'mpp', 'cmp',...,
    'gmp', 'mep', 'erythroblast', 'megakaryocyte', 'platelet', 'myelocyte', 'plasma_cell'};

originalHeaders = {'Time', 'Cell Length', '191-DNA','193-DNA','103-Viability', '115-CD45', '110-CD3', '111-CD3','112-CD3',...
    '114-CD3','139-CD45RA','141-pPLCgamma2', '142-CD19', '144-CD11b', '145-CD4',  '146-CD8', '148-CD34' , '150-pSTAT5', ...
    '147-CD20', '152-Ki67' ,'154-pSHP2','151-pERK1/2', '153-pMAPKAPK2', '156-pZAP70/Syk',  '158-CD33', '160-CD123' , '159-pSTAT3', ...
    '164-pSLP-76', '165-pNFkB', '166-IkBalpha', '167-CD38','168-pH3', '170-CD90', '169-pP38', '171-pBtk/Itk',  '172-pS6','174-pSrcFK', ...
    '176-pCREB', '175-pCrkL',  '110_114-CD3', 'EventNum'};
headers = {'Time', 'Cell_length', 'DNA-1','DNA-2','103-Viability', 'CD45', '110-CD3', '111-CD3','112-CD3',...
    '114-CD3','CD45RA','pPlcg2', 'CD19', 'CD11b', 'CD4',  'CD8', 'CD34' , 'pStat5', 'CD20', 'Ki67' ,'pSHP2','pErk', 'pMAPKAPK2', ...
    'pZap70',  'CD33', 'CD123' , 'pStat3', 'pSlp76', 'pNFkB', 'IkBa', 'CD38','pH3', 'CD90', 'pp38', 'pBtk',  'pS6','pSrcFK', ...
    'pCREB', 'pCrkL',  'CD3', 'EventNum'};

nPopulations = length(populations);
nActivators = length(activators);
nDonors = length(donors);
nVariables = length(headers);

% General populations: monocytes, T-cells, dendritic cells, B-cell, natural killer cells.
isMonocyte = false(1, nPopulations); isMonocyte (1:3)=true;
isTcell = false(1, nPopulations); isTcell(5:8)= true;
isDendritic = false(1, nPopulations); isDendritic(14) = true;
isBcell= false(1, nPopulations); isBcell(9:13)= true;
isNk = false(1, nPopulations); isNk(4)=true;

% Syk is measured instead of pZap70 in monocytes/dendritic/b-cells. BLNK is
% measured instead of pSlp76 in monocytes/b-cells. This matrix is used to
% correct the headers.
isSyk = isMonocyte|isDendritic|isBcell;
isBLNK = isMonocyte|isBcell;

%CD_, HLA-DR, Igm are surface markers used for gating.
isSurfMarker = ~cellfun('isempty', strfind(headers, 'CD'));
%DNA-1/2 are markers used for distinghuishing if an event is a cell
isDNAMarker = ~cellfun('isempty', strfind(headers, 'DNA'));
% the rest are functional markers - Time, Cell_length, EventNum
isFuncMarker = ~(isSurfMarker|isDNAMarker|~cellfun('isempty', strfind(headers, 'Time'))| ~cellfun('isempty', strfind(headers, 'Cell'))| ~cellfun('isempty', strfind(headers, 'Event')));

save([dataDir filesep 'info.mat']);

logFile = 'log.txt';
logFid = fopen(logFile, 'a+');
fprintf(logFid, '\n=======================================================================================\n');
fprintf(logFid, '                 Finished script_BND_01_data_info, time elapsed %.3f                   \n', toc(scrTic));
fprintf(logFid, '\n=======================================================================================\n');
