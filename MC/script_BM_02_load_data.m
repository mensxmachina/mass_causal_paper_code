% Script for loading .fcs files and saving as .mat files. 
% info.mat must be loaded.
infoFile = [dataDir filesep 'info.mat'];
load(infoFile);


% Choose arcsin cofactor for data transformation 
arcsinh_cofactor = 5;

% choose which inhibitors, populations, activators, dosages you care about.
inhibitorInds =1:length(inhibitors);
populationInds = 1:length(populations);
activatorInds = 1:length(activators);
dosageInds=1:length(dosageIDs);

strTic =  tic;

% Inhibitor Data
% variable for caching inhibitor sample sizes
sampleSizes = nan(nInhibitors, nPopulations, 12, nDosages);

fprintf(logFid, 'Saving .mat files for inhibitor data\n');
for iInh=inhibitorInds
    curInhibitor = inhibitors{iInh};
    curInhibitorFolder= inhibitorFolders{iInh};    
    fprintf('\n Processing inhibitor %s,----------------\n',  curInhibitor);
    fprintf(logFid, '\n Processing inhibitor %s,----------------\n',  curInhibitor);
    tInh =tic;
    % which populations
    for iPop = populationInds
        curPopulation = populations{iPop};
        fprintf(logFid, '\t Population %s,----------------\n',  curPopulation);
        % which activators
        for iAct = 1:12
            if iAct<12
                curActivatorID = activatorIDs{iAct};
                curActivator = activators{iAct};
            else 
                curActivatorID = referenceID;
                curActivator = reference;
            end
            fprintf(logFid, '\t \t Activator %s,----------------\n',  curActivator);

            % which dosages
            for iDos =dosageInds
                curDosageID = dosageIDs{iDos};
                curDosage = dosages{iDos};          
                fprintf(logFid, '\t \t \t Dosage %s: ...',  curDosage);
                fcsFile = [fcsDir, filesep, curInhibitorFolder, filesep, curInhibitor, '_' , curPopulation, '_', curDosageID, curActivatorID, '.fcs'];
                matFile = [matDir, filesep, curInhibitor, '_'  curDosage, '_' , curPopulation, '_', curActivator,  '.mat'];
                if exist(fcsFile, 'file');
                    [fcsdat, fcshdr, ~] = fca_readfcs(fcsFile);
                    fprintf(logFid, 'Loading file %s ...', fcsFile);

                    data = flow_arcsinh(fcsdat',arcsinh_cofactor);
                    %data = data(:,1:min(end,500000));
                    dataset.data =  data';
                    dataset.headers =  {fcshdr.par.name2};
                    if isSyk(iPop)
                        dataset.headers{ismember(dataset.headers, 'pZap70')} = 'Syk';
                    end
                    if isBLNK(iPop)
                        dataset.headers{ismember(dataset.headers, 'pSlp76')} = 'BLNK';
                    end
                    dataset.arcsinh_cofactor = arcsinh_cofactor;
                    sampleSizes(iInh, iPop, iAct, iDos)= size(dataset.data, 1);
                    fprintf(logFid, 'saving file %s\n', matFile);
                else
                    fprintf('File %s does not exist ...\n', fcsFile);
                    dataset.data =  [];
                    dataset.headers =  {};
                    sampleSizes(iInh, iPop, iAct, iDos)= 0;
                    fprintf(logFid, '\t saving empty file %s\n', matFile);                
                end % end if
               % save(matFile, 'dataset');
                clear dataset;
            end %end for iDos
            
        end% end activator
        
    end % end population    
   fprintf(logFid, '\n ------------Done with inhibitor %s, time elapsed %.3f----------------\n',  curInhibitor, toc(tInh));
end % end inhibitor 

%% Now load 8donor data
fprintf(logFid, 'Saving .mat files for  8donor data\n');
 for iDonor = 1:nDonors
    tDon = tic;
    curTimePointID = donorsIDs{iDonor};
    curDonor = donors{iDonor};
    fprintf(logFid, 'Donor %s,----------------\n',  curDonor);
    for iPop = 1:nEightDonorPopulations
        curPopulationID = eightDonorPopulationsIDs{iPop};
        curPopulation = eightDonorPopulations{iPop};
        if isempty(curPopulation)
            continue;
        end
        fprintf(logFid, '\t Population %s,----------------\n',  curPopulation);
   
        for iAct =1:nActivators
            curActivator = activators{iAct};
            curActivatorID = activatorsIDs{iAct};
            fcsFile = [fscDir, filesep, '8_donor', filesep,  'EightDonor_30min_', curPopulationID, '_', curTimePointID, curActivatorID, '.fcs'];
            matFile = [matDir, filesep,  curDonor, '_',  curPopulation, '_', curActivator '.mat'];
            if exist(fcsFile, 'file');
                [fcsdat, fcshdr, ~] = fca_readfcs(fcsFile);
                fprintf(logFid, 'Loading file %s ...', fcsFile);
                data = flow_arcsinh(fcsdat',arcsinh_cofactor);
                dataset.data =  data';
                dataset.headers =  {fcshdr.par.name2};
                if isSyk(iPop)
                    dataset.headers{ismember(dataset.headers, 'pZap70')} = 'Syk';
                end
                if isBLNK(iPop)
                    dataset.headers{ismember(dataset.headers, 'pSlp76')} = 'BLNK';
                end
                dataset.arcsinh_cofactor = arcsinh_cofactor;
                fprintf(logFid, 'saving file %s\n', matFile);
            else
                fprintf('File %s does not exist ...\n', fcsFile);
                dataset.data =  [];
                dataset.headers =  {};
                fprintf(logFid, '\t saving empty file %s\n', matFile);  
            end % end if
        %    save(matFile, 'dataset');
            clear dataset;
        end
    end
   fprintf(logFid, '\n ------------Done with donor %s, time elapsed %.3f----------------\n',  curDonor, toc(tDon));

end

%% Now load time course data

fprintf(logFid, 'Saving .mat files for time course data\n');
 for iTP =1:nTimePoints
    tTP =tic;
    curTimePointID = timePointIDs{iTP};
    curTimePoint = timePoints{iTP};
    fprintf(logFid, '\t Time Point %s,----------------\n',  curTimePoint);
    for iPop = 1:nPopulations
        curPopulation = eightDonorePopulations{iPop};
        fprintf(logFid, '\t Population %s,----------------\n',  curPopulation);
        for iAct =1:nActivators
            curActivator = activators{iAct};
            curActivatorID = activatorsIDs{iAct};
            fcsFile = [fscDir, filesep, '8_TimePoints', filesep, curPopulation, '_', curTimePointID, curActivatorID, '.fcs'];
            matFile = [matDir, filesep,  curTimePoint , '_',  curPopulation, '_', curActivator '.mat'];
            if exist(fcsFile, 'file');
                [fcsdat, fcshdr, ~] = fca_readfcs(fcsFile);
                fprintf(logFid, 'Loaded file %s ...', fcsFile);
                data = flow_arcsinh(fcsdat',arcsinh_cofactor);
                dataset.data =  data';
                dataset.headers =  {fcshdr.par.name2};
                if isSyk(iPop)
                    dataset.headers{ismember(dataset.headers, 'pZap70')} = 'Syk';
                end
                if isBLNK(iPop)
                    dataset.headers{ismember(dataset.headers, 'pSlp76')} = 'BLNK';
                end
                dataset.arcsinh_cofactor = arcsinh_cofactor;
                fprintf(logFid, 'saved file %s\n', matFile);
            else
                fprintf('File %s does not exist ...\n', fcsFile);
                dataset.data =  [];
                dataset.headers =  {};
                fprintf(logFid, '\t saved empty file %s\n', matFile);                
            end % end if
            %save(matFile, 'dataset');
            clear dataset;
        end
    end
   fprintf(logFid, '\n ------------Done with TimePoint %s, time elapsed %.3f----------------\n',  curTimePoint, toc(tTP));

end

fprintf(logFid, '\n=======================================================================================\n');
fprintf(logFid, '                 Finished script_02_load_data, time elapsed %.3f                     \n', toc(strTic));
fprintf(logFid, '\n=======================================================================================\n');

% clear directory variables
clearvars *Dir;
