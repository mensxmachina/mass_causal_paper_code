scrT = tic;
% Set and cd to the main directory
mainDir = 'C:\Users\striant\Documents\MATLAB\MC';
%addpath(genpath(curPath));

% Where the code is
codeDir = 'C:\Users\striant\Dropbox\MATLAB\MC\_code';
addpath(genpath(codeDir));


% Where you want to save all data
dataDir = [mainDir filesep 'BM'];
if ~isdir(dataDir);mkdir(dataDir);end

% Fcs file folder
fcsDir = [dataDir filesep '_fcs'];
if ~isdir(fcsDir);mkdir(fcsDir);end

% open log file
logFile = 'log.txt';
logFid = fopen(logFile, 'wt');

% Download zip files
urlFile = [dataDir filesep 'BMreport.html'];
url = 'http://supplemental.cytobank.org/report_data/report_105/fcs_file_downloads.html';
urlwrite(url,urlFile);

inhibitors = {'Akti', 'BTKi', 'Crassin' ,'Dasatinib', 'GDC-0941', 'Go69', 'H89', 'IKKi', 'Imatinib', 'Jak1i', 'Jak2i',...
    'Jak3i', 'Lcki', 'Lestaurtinib', 'PP2', 'Rapamycin', 'Ruxolitinib', 'SB202', 'SP6', 'Sorafenib', 'Staurosporine', ... 
    'Streptonigrin', 'Sunitinib', 'Syki',  'Tofacitinib', 'U0126', 'VX680'};
inhibitorUrls = cell(27,1);
iInh =0;

urlFid = fopen(urlFile, 'r');
while ~feof(urlFid)
    tline = fgetl(urlFid);
    if strfind(tline, '.zip')
        iInh =iInh+1;
        inhibitorUrls{iInh} = tline(strfind(tline, 'http'):strfind(tline, 'zip')+2);
        zipFile = [dataDir filesep inhibitors{iInh}, '.zip'];
        dT = tic;
        fprintf('Downloading %s\n', inhibitorUrls{iInh});
        urlwrite(inhibitorUrls{iInh}, zipFile);
        fprintf(logFid, 'Downloaded data set from %s, %.3f sec\n', inhibitorUrls{iInh}, toc(dT));
        fprintf('\t Unzipping %s\n', zipFile);
        zT =tic;
        unzip(zipFile, fcsDir);
        fprintf(logFid, 'Unziped %s in %s, %.3f sec\n', zipFile, fcsDir, toc(zT));
    end
end
fprintf(logFid, '\n=======================================================================================\n');
fprintf(logFid, '            Finished script_00_download_data, time elapsed: %.3f sec                  \n', toc(scrT));
fprintf(logFid, '\n=======================================================================================\n');


fclose(logFid);
fclose(urlFid);

clearvars -except *Dir inhibitors logFile


