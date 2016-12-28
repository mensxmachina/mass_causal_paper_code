strTic = tic;

mainDir = 'C:\Users\striant\Documents\MATLAB\MC';
codeDir = 'C:\Users\striant\Dropbox\MATLAB\MC\_code';
rmpath([codeDir filesep '_trash']);
dataDir = [mainDir filesep 'BM'];
matDir =  [dataDir filesep '_mat'];
fcsDir = [dataDir filesep '_fcs'];
if ~isdir(matDir);mkdir(matDir);end
cd(mainDir)

inhibitors = {'Akti', 'BTKi', 'Crassin' ,'Dasatinib', 'GDC-0941', 'Go69', 'H89', 'IKKi', 'Imatinib', 'Jak1i', 'Jak2i',...
    'Jak3i', 'Lcki', 'Lestaurtinib', 'PP2', 'Rapamycin', 'Ruxolitinib', 'SB202', 'SP6', 'Sorafenib', 'Staurosporine', ... 
    'Streptonigrin', 'Sunitinib', 'Syki',  'Tofacitinib', 'U0126', 'VX680'};

% inhibitorFolders = inhibitors;
% % For three inhbitors, the folder name is different than the filename.
% [inhibitorFolders{[1, 2, 5, 24]}] = deal('AKTi', 'BTK', 'GDC', 'SykInh');

inhibitorFolders = dir(fcsDir);
inhibitorFolders = {inhibitorFolders(3:end).name};

% Populations and 8 donor populations differ in monocytes: cd14+/-hladrhigh
% and cd14+/-hladrmid are missing, each replaced by a single population cd14+/-hladr+
populations = {'cd14+hladr-', 'cd14+hladrhigh', 'cd14+hladrmid', 'cd14+surf-', 'cd14-hladr-', 'cd14-hladrhigh', 'cd14-hladrmid' ...
     'cd14-surf-', 'cd4+', 'cd8+', 'dendritic', 'igm+', 'igm-', 'nk'};
eightDonorPopulationIDs = {'CD14+HLADR-Monocytes', 'CD14+HLADR+Monocytes', 'Surf-CD14+',  'CD14-HLADR-Monocytes', 'CD14-HLADR+Monocytes', 'Surf-CD14-', 'CD4+Tcells', 'CD8+Tcells', 'Dendritic', 'IgM+Bcells', 'IgM-Bcells', 'NKcells'};
eightDonorPopulations ={'cd14+hladr-', 'cd14+hladr+', 'cd14+surf-', 'cd14-hladr-', 'cd14-hladr+', 'cd14-surf-', 'cd4+', 'cd8+', 'dendritic', 'igm+', 'igm-', 'nk'};


%Dosages: A is the highest dosage, H is zero dosage.
dosageIDs = { 'A' , 'B', 'C', 'D', 'E', 'F', 'G', 'H'};
dosages ={'D7', 'D6', 'D5', 'D4', 'D3', 'D2', 'D1', 'D0'};

% Activators: fcs File ids and activator names
activatorIDs  = {'01','02', '03', '04', '06', '07', '08', '09', '10', '11', '12'};
activators = {'pVO4', 'IL3', 'IL2',  'IL12', 'GCSF', 'GMCSF', 'BCR', 'IFNg', 'IFNa', 'LPS', 'PMA-IONO'};

referenceID = '05';
reference = 'Ref';


% donor IDs and donor names
donorIDs = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'};
donors = {'donor1', 'donor2', 'donor3', 'donor4', 'donor5', 'donor6', 'donor7', 'donor8'};

% time point IDs and time points
timePointIDs = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'};
timePoints = {0, 1, 5, 15, 30, 60, 120, 240};


% load a file to save experiment headers.
loadFile =  [fcsDir filesep inhibitorFolders{1}  filesep inhibitors{1} '_' populations{9} '_A01.fcs'];
[~, fcshdr, ~] = fca_readfcs(loadFile);
headers =  {fcshdr.par.name2};
clear  loadFile fcshdr
% molecular weight of the langthanide tag of each header antibody. Can be
% used to explore channel bleeding.
mw = [nan, nan, 114, 115, 139,141,142,144,145,146,147,148,150:154,156,158:160,164:172,174:176, nan, nan];



nInhibitors = length(inhibitors);
nPopulations = length(populations);
nActivators = length(activators);
nVariables = length(headers);
nDosages = length(dosageIDs);


% General populations: monocytes, T-cells, dendritic cells, B-cell, natural killer cells.
isMonocyte = false(1, 14); isMonocyte (1:8)=true;
isTcell = false(1, 14); isTcell(9:10)= true;
isDendritic = false(1, 14); isDendritic(11) = true;
isBcell= false(1, 14); isBcell(12:13)= true;
isNk = false(1, 14); isNk(14)=true;


% Syk is measured instead of pZap70 in monocytes/dendritic/b-cells. BLNK is
% measured instead of pSlp76 in monocytes/b-cells. However, this is not
% reported in the original fcs files. These matrices are used to correct
% the headers.

isSyk = isMonocyte|isDendritic|isBcell;
isBLNK = isMonocyte|isBcell;

%BC1-BC7 are markers denoting barcoding channels.
isBCMarker = ~cellfun('isempty', strfind(headers, 'BC'));
%CD_, HLA-DR, Igm are surface markers used for gating.
isSurfMarker = ~cellfun('isempty', strfind(headers, 'CD'))| ~cellfun('isempty', strfind(headers, 'Ig'))| ~cellfun('isempty', strfind(headers, 'HLA'));
%DNA-1/2 are markers used for distinghuishing if an event is a cell
isDNAMarker = ~cellfun('isempty', strfind(headers, 'DNA'));
% the rest are functional markers - Time, Cell_length
isFuncMarker = ~(isSurfMarker|isBCMarker|isDNAMarker|~cellfun('isempty', strfind(headers, 'Time'))| ~cellfun('isempty', strfind(headers, 'Cell')));

save([dataDir filesep 'info.mat']);

logFile = 'log.txt';
logFid = fopen(logFile, 'a+');
fprintf(logFid, '\n=======================================================================================\n');
fprintf(logFid, '                 Finished script_01_data_info, time elapsed %.3f                   \n', toc(strTic));
fprintf(logFid, '\n=======================================================================================\n');





