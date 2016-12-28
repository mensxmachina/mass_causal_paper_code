%% 2 Check BC bleeding:
% Set a threshold for molecular weight interval. Let mw(X) be the molecular
% weight of the lanthanide attached to the antibody used to tag X. We then
% check the correlation between the measured levels each variable P and any
% barcoding channel BC_i where |mw(P)-mw(BC_i)|< mwWindow. 
% mainDir = 'C:\Users\striant\Documents\MATLAB\MC';
% codeDir = 'C:\Users\striant\Dropbox\MATLAB\MC\_code';
% addpath(genpath(codeDir))
% 
% rmpath([codeDir filesep '_trash']);
% dataDir = [mainDir filesep 'BM'];
% matDir =  [dataDir filesep '_mat'];
% cd(mainDir)
load([dataDir filesep 'info.mat']);

% do not include barcoding markers, time and DNA-2 variable.
isIncluded  = ~(isBCMarker|~cellfun('isempty', strfind(headers, 'Time'))|isDNAMarker|~cellfun('isempty', strfind(headers, 'Cell')));
nIncluded = sum(isIncluded);
includedVariables = find(isIncluded);

plotsDir = [mainDir filesep '_plots' filesep 'BC_plots'] ;
if ~isdir(plotsDir)
    mkdir(plotsDir)
end
isBCMarker = isBCMarker|isDNAMarker|~cellfun('isempty', strfind(headers, 'Cell'));
iPair =0;
nPairs = nIncluded*sum(isBCMarker);
pairsIds = nan(nPairs, 2);
% matrix to hold nPairs x nPlates correlations
nWells = (nActivators+1)*nDosages;
wells = zeros(nWells, 2);
correlations = nan(nPairs, nWells, nInhibitors);
pairHeaders = cell(nPairs, 1);
mwDif = nan(nPairs, 1);
% load colormap for scatter plot
cc = colormap(colorcube);
close all;
bcInds = find(isBCMarker);
barCode = nan(nInhibitors, nActivators, nDosages, length(bcInds));
barCodeWell = nan(nInhibitors, 96, length(bcInds));
for iInh =1:nInhibitors
    curInhibitor = inhibitors{iInh};
    loadFile = [matDir filesep curInhibitor '.mat'];
    load(loadFile);
    fprintf('Inhibitor %d of %d\n', iInh, nInhibitors);
    iWell = 0;
    for iAct =1:nActivators+1
        for iDosage =1:nDosages
            iWell =iWell+1;
            wells(iWell, :) =[iAct iDosage];
            wellDataInds = ismember(plateData(:, [37 38]), [iAct iDosage], 'rows');
            wellData = plateData(wellDataInds, :);           
            for iBC = 1:length(bcInds)
                curBC = bcInds(iBC); 
                if any(wellData(:, curBC)<0)
                    %fprintf('BC %d is 0 in well %d, continuing\n', iBC, iWell);
                    barCode(iInh, iAct, iDosage, iBC) = false;
                    barCodeWell(iInh, iWell, iBC) = false;
                    continue;
                else
                    barCode(iInh, iAct, iDosage, iBC) = true;
                    barCodeWell(iInh, iWell, iBC) = true;
                    for iVar =1:nIncluded
                        curVar = includedVariables(iVar);
                        iPair = (iBC-1)*nIncluded+iVar;

                        if isnan(pairsIds(iPair,1))
                            pairsIds(iPair,:) = [curVar, curBC];
                            pairHeaders{iPair} = sprintf('%s-%s', headers{curBC}, headers{curVar});
                            mwDif(iPair) = abs(mw(curVar)-mw(curBC));
                        end
                        correlations(iPair, iWell, iInh) = corr(wellData(:, curVar) , wellData(:, curBC));
                    end %iVar
                end %if BC =1;
            end % iBC
        end% iDos
    end% iAct
end% iInh

save('bc_bleeding.mat', 'correlations', 'barCode*');
%%  Figure: histogram of a BC channel
close all;
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'PaperUnits','centimeters','PaperPosition',[0 0 40 30])
bcFigAx= gca;
bcFig = gcf;
%box on; hold on;

hist(plateData(:, bcInds(3)), 1000);

xlabel(headers{bcInds(3)}, 'FontSize', 14);
title([headers{bcInds(3)} 'in plate '  curInhibitor], 'FontSize', 18);

% save file
figFile = [plotsDir filesep 'hist_BC'];
saveas(corFig2, figFile, 'png');
%%  Figure: Correlations per variable channel
cc2 = colormap(jet);
cStep =1:floor(length(cc2)/length(bcInds)):length(cc2);
close all;
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'PaperUnits','centimeters','PaperPosition',[0 0 40 30])
corFig2Ax = gca;
corFig2 = gcf;
box on; hold on;

medianCorsPerInh = squeeze(nanmedian(correlations, 2));
medianCorsPerPair = nanmedian(medianCorsPerInh,2);
uqCorsPerPair = quantile(medianCorsPerInh, .90, 2);
lqCorsPerPair = quantile(medianCorsPerInh, .10, 2);
hold all;
clear h
for iBC =1:length(bcInds)-3
    curBC = bcInds(iBC+1);
    curPairsIds = find(ismember(pairsIds(:, 2), curBC));
    curVar = includedVariables(iVar);
    h(iBC) = scatter(corFig2Ax, [1:nIncluded]-0.15+0.3*rand(1, length(curPairsIds)),  medianCorsPerPair(curPairsIds), 30, cc2(cStep(iBC), :), 'o', 'filled');   
end
%
% x axis
set(corFig2Ax, 'xtick', 1:nIncluded);
set(corFig2Ax, 'xticklabel',headers(isIncluded));
rotateXLabels(corFig2Ax, 45);
% For matlab 2014b+ :
% set(gca, 'tickLabelInterpreter', 'none');
% set(gca, 'xtickLabelRotation', 90);

%legend
legend(corFig2Ax, h, headers(bcInds(2:end-2)), 'FontSize', 12, 'FontName', 'Helvetica');

%y axis
ylim([-0.6 0.6]);


ylabel('corr(P_i, BC_j)', 'FontSize', 14);
xlabel(corFig2Ax, 'protein P_i', 'FontSize', 14, 'FontName', 'Helvetica');

% line in correlation = 0 
line(get(corFig2Ax, 'xlim'), [0 0], 'LineStyle', '--', 'color', 'k')

% title
title('Correlation of variables with BC channels', 'FontSize', 18, 'FontName', 'Helvetica');

% save file
figFile = [plotsDir filesep 'correlations_VARS_BC'];
saveas(corFig2, figFile, 'png');

%% now plot mean over all BCs, Cell Lenght, DNA-1, DNA-2
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'PaperUnits','centimeters','PaperPosition',[0 0 40 30])
corFigAx = gca;
corFig22 = gcf;
box on; hold on;


clear h
corsPerVarBC = nan(length(bcInds), nIncluded);

for iVar =1:nIncluded
    curVar = includedVariables(iVar);
    curPairInds= find(ismember(pairsIds(:, 1), curVar));
    corsPerVarBC(:, iVar) = medianCorsPerPair(curPairInds);
end
hold all;
   
h(1) = scatter(corFigAx, [1:nIncluded], corsPerVarBC(1, :), 30, cc2(1, :), 's', 'filled');   
h(2) = scatter(corFigAx, [1:nIncluded], mean(corsPerVarBC(2:8, :)), 50, cc2(20, :), 'o', 'filled');  
h(3) = scatter(corFigAx, [1:nIncluded], corsPerVarBC(end-1, :), 50, cc2(60, :), 'v', 'filled');
h(4) = scatter(corFigAx, [1:nIncluded], corsPerVarBC(end, :), 50, cc2(62, :), 'v', 'filled');   



% x axis
set(corFigAx, 'xtick', 1:nIncluded);
set(corFigAx, 'xticklabel', headers(isIncluded));
rotateXLabels(corFigAx, 45);
% For matlab 2014b+ :
% set(gca, 'tickLabelInterpreter', 'none');
% set(gca, 'xtickLabelRotation', 90);

%legend
legend(corFigAx, h,{'Cell Length', 'mean BC', 'DNA-1', 'DNA-2'}, 'FontSize', 12, 'FontName', 'Helvetica');

%y axis
ylim([-0.6 0.6]);


ylabel(corFigAx, 'correlation coefficient', 'FontSize', 14, 'FontName', 'Helvetica');
xlabel(corFigAx, 'protein P_i', 'FontSize', 14, 'FontName', 'Helvetica');


% line in correlation = 0 
line(get(corFigAx, 'xlim'), [0 0], 'LineStyle', '--', 'color', 'k')

% title
title(['Correlation of variables with BC channels, Cell Length, DNA'], 'FontSize', 18, 'FontName', 'Helvetica');

% save file
figFile = [plotsDir filesep 'correlations_VARS_BC_CL_DNA'];
saveas(corFig22, figFile, 'png');
print('-dpng','-r150', figFile)

%%  Figure: Correlations per mw difference
close all;
[tmp, ord] = sort(mwDif);
ord = ord(~isnan(tmp));
orderedPairsIds = pairsIds(ord, :);
orderedPairHeaders = pairHeaders(ord);
orderedMwDif = mwDif(ord);
nPlotPairs =length(orderedPairsIds);

cc2 = colormap(jet);
close all;
figure('Units', 'normalized', 'outerposition', [0 0 1 1], 'PaperUnits','centimeters','PaperPosition',[0 0 40 30]);
corFig2Ax = gca;
corFig2 = gcf;
box on; hold on;

orderedMedianCorsPerPair = medianCorsPerPair(ord);

h = gscatter(orderedMwDif, orderedMedianCorsPerPair, reshape(headers(orderedPairsIds(:, 2)),nPlotPairs,1), cc2(1:floor(length(cc2)/7):length(cc2), :), '.', 25);%, '50*ones(1,nPlotPairs));%50*ones(1,nPlotPairs), , 'filled');
[a, b] = unique(sort(mwDif(~isnan(mwDif))), 'last');mc =[];
for i =1:length(a)
    mc(i) =nanmean(orderedMedianCorsPerPair(orderedMwDif==a(i)));
end
% tmp = b';
% x axis
% x
% [a, b] = unique(sort(mwDif(~isnan(mwDif))), 'last');
% tmp = b';
% skipLabel= [tmp(2:end)-tmp(1:end-1)==1 0]'|a>35;
% dXTickLabel = cellfun(@num2str, num2cell(a), 'UniformOutput', false);
% [dXTickLabel{skipLabel}] =deal('');
% set(corFig2Ax, 'xtick', b);
% set(corFig2Ax, 'xticklabel',dXTickLabel);
set(corFig2Ax, 'xtick', 1:5:max(orderedMWDif));
% set(corFig2Ax, 'xticklabel',dXTickLabel);
xlabel('Absolute pair molecular weight differrence |mwDif_{ij}|', 'FontSize', 14);  

%y axis
ylim([-0.6 0.6]);
ylabel('correlation corr(P_i, BC_j)', 'FontSize', 14);

%legend
[a,b,c,text]=legend;
ord2 =[7 1:6];
legend(c(ord2), text(ord2), 'FontSize', 12);

% line in correlation = 0 
line(get(corFig2Ax, 'xlim'), [0 0], 'LineStyle', '--', 'color', 'k')

% title
title('Correlation of variables with BC markers', 'FontSize', 18);

% save file
figFile = [plotsDir filesep 'correlations_BC_MWDIF'];
saveas(corFig2, figFile, 'png');
