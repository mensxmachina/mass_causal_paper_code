% Pulling all data together for removing BC bleeding:
% BC bleeding affects the entire tube processed in the machine. The tube
% includes one plate (?), so we pull all data from each plate together to
% do the removal.
% i.e. we pull together all populations, all activators and all dosages for each
% inhibitor.
mainDir = 'C:\Users\striant\Documents\MATLAB\MC';
codeDir = 'C:\Users\striant\Dropbox\MATLAB\MC\_code';
dataDir = [mainDir filesep 'BM'];
matDir =  [dataDir filesep '_mat'];
cd(mainDir)
load([dataDir filesep 'info.mat']);


for iInh =1:nInhibitors
    curInhibitor = inhibitors{iInh};
    curSampleSizes = squeeze(sampleSizes(iInh, :, :, :));
    plateSampleSize = sum(sum(sum(curSampleSizes)));
    plateData  = nan(plateSampleSize, 38);
    saveFile = [matDir filesep curInhibitor '.mat'];
    sampleInd =0;
    for iPop =1:nPopulations
        curPopulation = populations{iPop};
        fprintf('\nPopulation %s\n', curPopulation);
        % load activator data
        for iAct =1:nActivators            
            curActivator = activators{iAct};                
            fprintf('\n Activator %s\n', curActivator);
            for iDos = 1:nDosages
                nSamples = curSampleSizes(iPop, iAct, iDos);
                repCol = ones(nSamples, 1);
                curDosage = dosages{iDos};
                fprintf('\t dosage %s\t', curDosage);
                loadFile = [matDir, filesep, curInhibitor, filesep, curInhibitor, '_' , curPopulation, '_', curActivator, '_'  curDosage,  '.mat'];
                ds =load(loadFile);
                plateData(sampleInd+1:sampleInd+nSamples, :) = [ds.dataset.data iPop*repCol iAct*repCol iDos*repCol];
                sampleInd = sampleInd+nSamples;
            end % end for iDos  
        end %end for iAct
        % load reference data
        fprintf('\n Reference \n');
        for iDos = 1:nDosages
            nSamples = curSampleSizes(iPop, 12, iDos);
            repCol = ones(nSamples, 1);
            curDosage = dosages{iDos};
            fprintf('\t dosage %s\t', curDosage);
            loadFile = [matDir, filesep, curInhibitor, filesep, curInhibitor, '_' , curPopulation, '_Ref_'  curDosage,  '.mat'];
            ds =load(loadFile);
            plateData(sampleInd+1:sampleInd+nSamples, :)= [ds.dataset.data iPop*repCol 12*repCol iDos*repCol];
            sampleInd = sampleInd+nSamples;
        end % end for iDos                 
    end %end for 
    plateHeaders = cell(1, 38);
    plateHeaders(1:35) = headers;
    [plateHeaders{36:38}] = deal('Population', 'Activator', 'Dosage');
    save(saveFile, 'plateData', 'plateHeaders');
    clc
end %end for iInh







