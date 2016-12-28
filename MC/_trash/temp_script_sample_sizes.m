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
    tInh =tic;
    % which populations
    for iPop = populationInds
        curPopulation = populations{iPop};
        fprintf( '\t Population %s,----------------\n',  curPopulation);
        % which activators
        for iAct = 1:12
            if iAct<12
                curActivatorID = activatorIDs{iAct};
                curActivator = activators{iAct};
            else 
                curActivatorID = referenceID;
                curActivator = reference;
            end
            fprintf( '\t \t Activator %s,----------------\n',  curActivator);

            % which dosages
            for iDos =dosageInds
                curDosageID = dosageIDs{iDos};
                curDosage = dosages{iDos};          
                fprintf( '\t \t \t Dosage %s: ...',  curDosage);
                matFile = [matDir, filesep, curInhibitor, filesep, curInhibitor, '_' , curPopulation, '_', curActivator, '_'  curDosage,  '.mat'];
                load(matFile);
                sampleSizes(iInh, iPop, iAct, iDos)= size(dataset.data, 1);
            end %end for iDos
            
        end% end activator
        
    end % end population    
   fprintf('\n ------------Done with inhibitor %s, time elapsed %.3f----------------\n',  curInhibitor, toc(tInh));
end % end inhibitor 
fprintf('\n=======================================================================================\n');
fprintf('                 Finished script_02_load_data, time elapsed %.3f                     \n', toc(scrTic));
fprintf('\n=======================================================================================\n');
save(infoFile, 'sampleSizes', '-append');

