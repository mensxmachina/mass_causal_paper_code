% Script for loading .fcs files and saving as .mat files. 
% info.mat must be loaded.
infoFile = [dataDir filesep 'info.mat'];
load(infoFile);
% variable for caching sample sizes
sampleSizes = nan(nInhibitors, nPopulations,12, nDosages);

% Choose arcsin cofactor for data transformation. XXXX
arcsinh_cofactor = 5;


% choose which inhibitors, populations, activators, dosages you care about.
inhibitorInds =1:length(inhibitors);
populationInds = 1:length(populations);
activatorInds = 1:length(activators);
dosageInds=1:length(dosageIDs);

scrTic =  tic;
% which inhibitors
for iInh=inhibitorInds
    curInhibitor = inhibitors{iInh};
    curInhibitorFolder= inhibitorFolders{iInh};    
    fprintf('\n Processing inhibitor %s,----------------\n',  curInhibitor);
    fprintf(logFid, '\n Processing inhibitor %s,----------------\n',  curInhibitor);

    if ~isdir([matDir filesep curInhibitor])
        mkdir([matDir filesep curInhibitor]);
        fprintf('\t Creating Folder %s\n',  [matDir filesep curInhibitor]);
        fprintf(logFid, '\t Creating Folder %s\n',  [matDir filesep curInhibitor]);
    end
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
                if exist(fcsFile, 'file');
                    [fcsdat, fcshdr, ~] = fca_readfcs(fcsFile);
                    fprintf(logFid, 'Loaded file %s ...', fcsFile);

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
                    matFile = [matDir, filesep, curInhibitor, filesep, curInhibitor, '_' , curPopulation, '_', curActivator, '_'  curDosage,  '.mat'];
                    fprintf(logFid, 'saved file %s\n', matFile);
                    save(matFile, 'dataset') 
                else
                    fprintf('File %s does not exist ...\n', fcsFile);
                    dataset.data =  [];
                    dataset.headers =  {};
                    sampleSizes(iInh, iPop, iAct, iDos)= 0;
                    matFile = [matDir, filesep, curInhibitor, filesep, curInhibitor, '_' , curPopulation, '_', curActivator, '_'  curDosage,  '.mat'];
                    fprintf(logFid, '\t saved empty file %s\n', matFile);                
                    save(matFile, 'dataset');
                end % end if
                clear dataset;
            end %end for iDos
            
        end% end activator
        
    end % end population    
   fprintf(logFid, '\n ------------Done with inhibitor %s, time elapsed %.3f----------------\n',  curInhibitor, toc(tInh));
end % end inhibitor 
fprintf(logFid, '\n=======================================================================================\n');
fprintf(logFid, '                 Finished script_02_load_data, time elapsed %.3f                     \n', toc(scrTic));
fprintf(logFid, '\n=======================================================================================\n');

   
