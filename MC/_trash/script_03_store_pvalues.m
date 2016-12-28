% script for calculating and storing p-values and such.
% info.mat must be loaded
% info.mat must be loaded
mainDir = 'C:\Users\striant\Documents\MATLAB\MC';
codeDir = 'C:\Users\striant\Dropbox\MATLAB\MC\_code';
dataDir = [mainDir filesep 'BM'];
matDir =  [dataDir filesep '_mat'];
cd(mainDir)
load([dataDir filesep 'info.mat']);

% open log file
logFile = 'log_nbc.txt';
logFid = fopen(logFile, 'a');


% do not include barcoding markers, time and DNA-2 variable.
isIncluded  = ~(isBCMarker|~cellfun('isempty', strfind(headers, 'Time'))| ~cellfun('isempty', strfind(headers, 'DNA-2')));

% choose which inhibitors, populations, activators,  you care about.
inhibitorInds =1:2%1:nInhibitors;
populationInds = 1:length(populations);%length(populations);
activatorInds = 1:length(activators);

% each triplet (activator --> k --> m) is represented by a row with the following meaning: 
% [k m pvalue(activator,k) pvalue(k,m) pvalue(activator,m|k) pvalue(activator,m) pvalue(activator,k|m) pvalue(k,m|activator) activation(activator,k) activation(activator,m) R^2(k,m) sampleSize]
nStatistics = 9;


% If you have run an experiment before, load statitstics, else create new
% matrix to store statistics and sample sizes.
% statFile = [mainDir filesep 'statistics_nbc.mat'];
% % if you have already performed some of the calculation
% if ~exist(statFile, 'file')
%     statistics = nan(nInhibitors, nPopulations, nActivators, nVariables, nVariables, nStatistics);
% else 
%     load(statFile);
% end
% count = 1;
statFile = [mainDir filesep 'statistics.mat'];
statistics = nan(nInhibitors, nPopulations, nActivators, nVariables, nVariables, nStatistics);
% open pool for parallel calculations
matlabpool open 12;
% for Matlab 2014+ use the following:
% myPool = parpool(11);

% timing
scrTic =  tic;
% let's calculate the pvalues!   
% for each inhibitor
for iInh = inhibitorInds 
    tInh = tic;
    curInhibitor = inhibitors{iInh};
    fprintf('Inhibitor %s (%d of %d),\n', curInhibitor, iInh, length(inhibitorInds));
    fprintf(logFid, 'Inhibitor %s (%d of %d),\n', curInhibitor, iInh, length(inhibitorInds));
    clear plate*;
    loadFile = [matDir filesep curInhibitor '.mat'];
    load(loadFile);
    % for each population
    parfor iPop = populationInds
        curPopulation = populations{iPop};
        %fprintf(logFid, '\t Population %s (%d of %d) ...', curPopulation, iPop, length(populationInds));
        %for each activator
        populationData = plateData(ismember(plateData(:, [36 38]), [iPop 8],  'rows'), :);
        for iAct = 1:nActivators
            tmpStats = nan(nVariables, nVariables, 9);
            %curActivator = activators{iAct};
            %fprintf(logFid, '\t\t Activator %s\n', curActivator);
            %loading the data
            actData =  populationData(ismember(populationData(:, 37), iAct), :);
            if size(actData,1) < 10
                %fprintf(logFid, 'Less than 10 activator samples, skipping\n')
                continue;
            end
            refData =  populationData(ismember(populationData(:, 37), 12), :);
            if size(refData,1) < 10
                %fprintf(logFid, 'Less than 10 reference samples, skipping\n')
                continue;
            end

            %creating the data
            data = [[actData 1*ones(size(actData,1),1)];...
                   [refData 0*ones(size(refData,1),1)]];
            activatorId = size(data,2);

            %calculating the pvalues for each couple of measurements

            for k = 1:nVariables
                if ~isIncluded(k)
                    continue;
                end

                for m = 1:nVariables

                    if ~isIncluded(m)
                        continue;
                    end

                    %skip if same variable
                    if m == k
                        continue;
                    end

                    %calculating the pvalues
                    %filling the triplet
                    %each triplet (activator --> k --> m) is represented by a row with the following meaning: 
                    %[1.k 2.m 3.pvalue(activator,k) 4.pvalue(activator, m) 5.pvalue(k,m) 6.pvalue(activator,m|k) 7.pvalue(activator,k|m) 8.pvalue(k,m|activator) ...
                    % 9. sampleSize]
                    tmpStats(k, m,1) = k;
                    tmpStats(k, m,2) = m;
                    idx1 = data(:,activatorId) == 1;
                    idx0 = data(:,activatorId) == 0;   
                    % 3.pvalue(activator,k)
                    [~, tmpStats(k, m, 3)] = ttest2(data(idx1,k), data(idx0,k));
                    % 4.pvalue(activator,m)
                    [~, tmpStats(k, m, 4)] = ttest2(data(idx1,m), data(idx0,m));
                    % 5.pvalue(k,m)
                    tmpStats(k, m, 5) = fisher(k, m, [], data);
                    % 6.pvalue(activator,m|k) IND
                    tmpStats(k, m, 6) = logisticTest(activatorId, m, k, data);
                    % 7.pvalue(activator,k|m)
                    tmpStats(k, m, 7) = logisticTest(activatorId, k, m, data);
                    % 8.pvalue(k,m|activator)
                    tmpStats(k, m, 8) = fisher(k, m, activatorId, data);
                    % 9. correlation(k, m)
                    tmpStats(k, m, 9) = corr(data(:, k), data(:, m));                    
                end
            end
            statistics(iInh, iPop, iAct,:, :, :) = tmpStats;            
            %fprintf(logFid, '... done\n');
        end % end activators
        %fprintf(logFid, 'done\n');

    end % end populations
   fprintf(logFid, '------------Finished Inhbitor %s, time elapsed %.3f----------------\n\n', curInhibitor,  toc(tInh));

end % end inhibitors
    

fprintf(logFid, 'Saving statistics\n');
%saving the data...
save(statFile, 'statistics');

matlabpool close;
% for Matlab 2014+ use the following:
% delete(myPool);


% 
fprintf(logFid, '\n=======================================================================================\n');
fprintf(logFid, '                 Finished script_03_store_pvalues, time elapsed %.3f                     \n', toc(scrTic));
fprintf(logFid, '\n=======================================================================================\n');
fclose(logFid);
   