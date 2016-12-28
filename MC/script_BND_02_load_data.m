% Script for loading .fcs files and saving as .mat files. 
% info.mat must be loaded.
infoFile = [dataDir filesep 'info.mat'];
load(infoFile);


allConditionIDs = [referenceIDs, activatorIDs];
allConditions = [reference, activators];
nConditions = length(allConditionIDs);
arcsinh_cofactor =5;

for iPop = 1:nPopulations
    curPopulation = populations{iPop};
    curPopulationID = populationIDs{iPop};
    fprintf(logFid, '\t Population %s,----------------\n',  curPopulation);
    for iDonor = 1:nDonors 
        curDonorID = donorIDs{iDonor};
        curDonor = donors{iDonor};
        for iCond= nConditions        
            clear dataset
            curConditionID = allConditionIDs{iCond};
            curCondition = allConditions{iCond};
            curNum = num2str(iCond, '%02d');
            fcsFile = [fcsDir filesep curDonorID, '_', curNum, '_',  curConditionID, '_', curDonorID, '_', curConditionID, '_', curPopulationID, '.fcs'];
            matFile = [matDir filesep curPopulation, '_', curCondition, '_', curDonor, '.mat']; 
            [fcsdat, fcshdr, fcsdatscaled] = fca_readfcs(fcsFile);
            data = flow_arcsinh(fcsdat',arcsinh_cofactor);
            dataset.data =data';
            h ={fcshdr.par.name2};
            if ~isequal(originalHeaders,  h)
                fprintf(logFid, 'Wrong headers in file %s\n\n', fcsFile);
                return;
            end      
            dataset.headers = headers;
            if isSyk(iPop)
                dataset.headers{ismember(dataset.headers, 'pZap70')} = 'Syk';
            end
            if isBLNK(iPop)
                dataset.headers{ismember(dataset.headers, 'pSlp76')} = 'BLNK';
            end
            dataset.arcsinh = arcsinh_cofactor;
            fprintf(logFid, 'Saving %s to \n\t %s with ss %d\n', fcsFile, matFile, size(dataset.data, 1));        
            save(matFile, 'dataset');
        end
    end
end

% clear directory variables
clearvars *Dir;