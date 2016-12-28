% Pulling all data together for removing BC bleeding:
% BC bleeding affects the entire tube processed in the machine. The tube
% includes one plate (?), so we pull all data from each plate together to
% do the removal.
% i.e. we pull together all populations, all activators and all dosages for each
% inhibitor.

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
    save(saveFile, 'plateData');
    clc
end %end for iInh

%% 2 Check overall BC bleeding:
% Set a threshold for molecular weight interval. Let mw(X) be the molecular
% weight of the lanthanide attached to the antibody used to tag X. We then
% check the correlation between the measured levels each variable P and any
% barcoding channel BC_i where |mw(P)-mw(BC_i)|< mwWindow. 

% do not include barcoding markers, time and DNA-2 variable.
isIncluded  = ~(isBCMarker|~cellfun('isempty', strfind(headers, 'Time'))| ~cellfun('isempty', strfind(headers, 'DNA-2')));


plotsDir = [mainDir filesep '_plots' filesep 'BC_plots'] ;
if ~isdir(plotsDir)
    mkdir(plotsDir)
end
mwWindow = 10;
iPair =0;
pairsIds = nan(nVariables*7, 2);
% matrix to hold nPairs x nPlates correlations
correlations = nan(nVariables*7,  27);
pvalues = nan(nVariables*7,  27);
pairHeaders = cell(nVariables*7, 1);
mwDif = nan(nVariables*7, 1);
% load colormap for scatter plot
cc = colormap(colorcube);
close all;

for iInh =1:27
    curInhibitor = inhibitors{iInh};
    loadFile = [matDir filesep curInhibitor '.mat'];
    load(loadFile);
    fprintf('Inhibitor %d of %d\n', iInh, nInhibitors);
    %saveFile = [matDir filesep curInhibitor '_BCrem' num2str(mwWindow) '.mat'];
    %nSamples = size(plateData, 1);
    for iVar = 1:nVariables
        % find which BCs to regress
        if ~isIncluded(iVar)
            continue;
        end
        curMw = mw(iVar);
        checkMws = mw<curMw;
        whichBCs = find(isBCMarker & mw>=curMw-mwWindow & mw<=curMw+mwWindow);
        if isempty(whichBCs)
            continue;
        end
       % fprintf('Checking BC markers [%s] with mw [%s] on variable %s with mw %d\n', [headers{whichBCs}], num2str(mw(whichBCs)), headers{iVar}, mw(iVar));
        % for iBC = whichBCs
        % plot BC vs variable
%         varData= plateData(:, iVar);
%         bcData = [plateData(:, whichBCs) ones(nSamples, 1)];
%         [b, bint, r, rint, stats] = regress(varData, bcData);
        % XXX Print some stuff about which regressors turned out important;
        % Replace variable data with residuals.
        % plateData(:, iVar) = r;
        % downsample to 100000 bit for the plots
        % randInds = randsample(1:nSamples, 100000);   
        for iBC = 1:length(whichBCs)
            curBC = whichBCs(iBC);
            if iInh ==1
                iPair =iPair+1;
                pairsIds(iPair,:) = [iVar, curBC];
                pairHeaders{iPair} = sprintf('%s-%s', headers{curBC}, headers{iVar});
                mwDif(iPair) = abs(mw(iVar)-mw(curBC));
            else
                iPair = find(ismember(pairsIds, [iVar curBC], 'rows'));
            end
            [correlations(iPair, iInh), pvalues(iPair, iInh)] = corr(plateData(:, iVar), plateData(:, curBC));      
%             scFig = figure;
%             scatplot2(plateData(randInds, curBC), plateData(randInds, iVar)); 
%             hold all
%             ylabel(headers{iVar});
%             xlabel(headers{curBC});
%             title(sprintf('%s, reg.coef.: %.5f', curInhibitor, b(iBC)));
%             plotFileName = [plotsDir filesep scatterplots filesep curInhibitor '_' headers{iVar} '_' headers{curBC}];
%             saveas(gcf, plotFileName, 'fig');
%             saveas(gcf, plotFileName, 'png');
%             close(scFig);
        end
    end% end for iVar
    if iInh ==1
        nPairs = iPair;
        pairsIds = pairsIds(1:nPairs, :);
        pairHeaders = pairHeaders(1:nPairs);
        correlations = correlations(1:nPairs, :);
        pvalues = pvalues(1:nPairs, :);
        mwDif = mwDif(1:nPairs, :);
        [~, ord] = sort(mwDif);

    end
    %x= 1:nPairs;
    %isSignificant = pvalues(:, iInh) <0.05;
    %scatter(corFigAx, x(isSignificant), correlations(isSignificant, iInh), 60, cc(iInh, :), 's',  'filled');   

    %scatter(corFigAx, x(~isSignificant), correlations(~isSignificant, iInh), 30, cc(iInh, :), 'filled'); 
    


  %  hold all;
end
%%  Figure: Correlations per BC channel

[unBC, unBCLoc] = unique(orderedPairsIds(:, 2));
cc2 =zeros(nPairs, 3);
cc2(unBC, :) =[0 0 1; 0 1 1; 0 1 0; 1 0 0; 1 1 0; 1 0 1;0.5 0.5 0.5];
figure('units','normalized','outerposition',[0 0 1 1]);
corFig2Ax = gca;
corFig2 = gcf;
orderedCorrelations = correlations(ord, :);
orderedPairsIds = pairsIds(ord, :);

for iPair =1:nPairs
    h(iPair) = scatter(corFig2Ax, iPair*ones(1, 27),  orderedCorrelations(iPair, :), 60, cc2(orderedPairsIds(iPair,2), :), 's',  'filled');   
    w(iPair) = orderedPairsIds(iPair, 2);
    hold all;
end

% x axis
set(corFigAx, 'xtick', 1:nPairs);
set(corFi2Ax, 'xticklabel',pairHeaders(ord));
rotateXLabels(corFig2Ax, 90);
% For matlab 2014b+ :
% set(gca, 'tickLabelInterpreter', 'none');
% set(gca, 'xtickLabelRotation', 90);

%legend
legend(corFigAx, h(unBCLoc), headers(isBCMarker));

%y axis
ylim([-0.5 0.5]);
ylabel('correlation coefficient');

% double x, y axes
[~, dXtick] = unique(mwDif(ord));
corFig2Ax2 = axes('Position', get(corFig2Ax, 'Position'), 'xAxisLocation', 'top', 'yAxisLocation', 'right', 'Color', 'none',...
    'xlim', get(corFig2Ax, 'xlim'), 'xTick', dXtick, 'xTicklabel', 1:mwWindow,  ...
    'ylim', get(corFig2Ax, 'ylim'), 'yTick', [], 'yTickLabel', {});
xlabel(corFig2Ax2, 'molecular weight difference');

% line in correlation = 0 
line(get(corFig2Ax, 'xlim'), [0 0], 'LineStyle', '--', 'color', 'k')

% title
title(sprintf('Correlation of variables with BC channels within +-%d mw window', mwWindow));

% save file
figFile = [plotsDir filesep 'correlations_per_plate_mw_' sprintf('%d', mwWindow)];
saveas(corFig2, figFile, 'png');
%% Figure: Correlations per plate (inhibitor)
figure('units','normalized','outerposition',[0 0 1 1]);
corFigAx = gca;
corFig = gcf;
hold all;
for iInh = 1:nInhibitors
    scatter(corFigAx, 1:nPairs, correlations(ord, iInh), 60, cc(iInh, :), 's',  'filled');   
end


% x axis
set(corFigAx, 'xtick', 1:nPairs);
set(corFigAx, 'xticklabel',pairHeaders(ord));
% enlarge a bit the x axis for legend
set(corFigAx, 'xlim', [0 nPairs+10]);
rotateXLabels(corFigAx, 90);
% For matlab 2014b+ :
% set(gca, 'tickLabelInterpreter', 'none');
% set(gca, 'xtickLabelRotation', 90);

%legend
legend(corFigAx, inhibitors);

%y axis
ylim([-0.5 0.5]);
ylabel('correlation coefficient');

% double x, y axes
[~, dXtick] = unique(mwDif(ord));
corFig2Ax2 = axes('Position', get(corFigAx, 'Position'), 'xAxisLocation', 'top', 'yAxisLocation', 'right', 'Color', 'none',...
    'xlim', get(corFigAx, 'xlim'), 'xTick', dXtick, 'xTicklabel', 1:mwWindow,  ...
    'ylim', get(corFigAx, 'ylim'), 'yTick', [], 'yTickLabel', {});
xlabel(corFig2Ax2, 'molecular weight difference')
% title
title(sprintf('Correlation of variables with BC channels within +-%d mw window', mwWindow));


% line in correlation = 0 
line([0 nPairs+1], [0 0], 'LineStyle', '--', 'color', 'k')

% save file
figFile = [plotsDir filesep 'correlations_per_plate_mw_' sprintf('%d', mwWindow)];
saveas(corFig, figFile, 'png');

%%



%%
set(corFigAx, 'xtick', 1:nPairs);
%set(corFigAx, 'xtickLabels', pairHeaders);
xticklabel_rotate(1:nPairs, 60, pairHeaders(ord), 'interpreter', 'none');
% For matlab 2014b+ :
% set(gca, 'tickLabelInterpreter', 'none');
% set(gca, 'xtickLabelRotation', 60);
ylim([-0.5 0.5]);
ylabel('correlation coefficient');
%legend(corFigAx, inhibitors, 'Location', 'SouthEastOutside');
legend(corFigAx, headers(isBCMarker), 'Location', 'SouthEastOutside');
figFile = [plotsDir filesep 'correlations_bc_per_plate_mw_' sprintf('%d', mwWindow)];

[~, dXtick] = unique(mwDif);
corFigAx2 = axes('Position', get(corFigAx, 'Position'), 'xAxisLocation', 'top', 'yAxisLocation', 'right', 'Color', 'none',...
    'xlim', get(corFigAx, 'xlim'), 'xTick', 1:nPairs, 'xTicklabel', mwDif(ord), ...
    'yTick', [], 'yTickLabel', {});
xlabel(corFigAx2, 'mw difference');
title(sprintf('Correlation of variables with BC channels within +-%d mw window', mwWindow));

set(corFig,'PaperUnits', 'inches')
set(corFig, 'PaperPosition', [0 0 15 15*6/9]);
set(corFig, 'PaperPositionMode', 'manual');
saveas(corFig, figFile, 'png');

% For Matlab 2014a+ 
% corFig.PaperUnits = 'inches';
% corFig.PaperPosition = [0 0 13 13*6/9];
% corFig.PaperPositionMode = 'manual';
% saveas(corFig, figFile, 'png');

%%
% For each variable, remove bleeding from  BC channels with antibody tag
% within \pm 5 molecular weight.
mainDir = 'BM';

mwWindow = 5;
indsBC = find(isBCMarker);

iPair =0;
for iInh =inhibitorInds
    curInhibitor = inhibitors{iInh};
    loadFile = [matDir filesep curInhibitor '.mat'];
   % saveFile = [matDir filesep curInhibitor '_BCrem' num2str(mwWindow) '.mat'];
    nSamples = size(plateData, 1);   
    for iVar = 1:nVariables
        if ~isIncluded(iVar)
            continue;
        end
        curMw = mw(iVar);             
        checkMws = mw<curMw;
        whichBCs = find(isBCMarker & mw>=curMw-mwWindow & mw<=curMw+mwWindow);
        if isempty(whichBCs)
            continue;
        end
        % pick a specific experimental condition
        for iAct =1:12
            for iDos = 1:nDosages
                experimentInds = plateData(:,37)==iAct & plateData(:, 38)==iDos;    
                % find which BCs to regress
               
                fprintf('Regressing BC markers [%s] with mw [%s] on variable %s with mw %d\n', [headers{whichBCs}], num2str(mw(whichBCs)), headers{iVar}, mw(iVar));
                % for iBC = whichBCs
                % plot BC vs variable
                varData= plateData(experimentInds, iVar);
                bcData = [plateData(experimentInds, whichBCs) ones(sum(experimentInds), 1)];
                [b, bint, r, rint, stats] = regress(varData, bcData);
                % XXX Print some stuff about which regressors turned out important;
                % Replace variable data with residuals.
             %   plateData(:, iVar) = r;

                for iBC = 1:length(whichBCs)    
                    close all;  
                    curBC = whichBCs(iBC);
                    scatplot2(plateData(experimentInds, curBC), plateData(experimentInds, iVar)); 
                    hold all
                    ylabel(headers{iVar});
                    xlabel(headers{curBC});
                    title(sprintf('%s, reg.coef.: %.5f', curInhibitor, b(iBC))); plotsDir = [mainDir filesep '_plots' filesep 'BC_plots' filesep headers{iVar} '_' headers{curBC}];
                    if ~isdir(plotsDir)
                        mkdir(plotsDir);
                    end

                    if iAct ==12
                        plotFileName = [plotsDir filesep curInhibitor '_Ref_' dosages{iDos} '_'  headers{iVar} '_' headers{curBC}];
                    else
                        plotFileName = [plotsDir filesep curInhibitor '_' activators{iAct} '_' dosages{iDos} '_'  headers{iVar} '_' headers{curBC}];
                    end
                    saveas(gcf, plotFileName, 'fig');
                    saveas(gcf, plotFileName, 'png');
                end
            end % end for dosgage
        end % end for activator
    end% end for iVar

    %save(saveFile, 'plateData');
end





