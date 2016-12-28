%% 2 Remove BC bleeding:
% Regress all BC markers against each variable, and replace values with
% residuals.

% do not include barcoding markers, time and DNA-2 variable.
isIncluded  = ~(isBCMarker|~cellfun('isempty', strfind(headers, 'Time'))| ~cellfun('isempty', strfind(headers, 'DNA-2')));
% regress all BC markers
bcInds = find(isBCMarker);
r_stats = cell(nInhibitors, nVariables);
beta_coeffs = nan(nInhibitors, nVariables, 8);

for iInh =1:27
    clear plate*
    curInhibitor = inhibitors{iInh};
    loadFile = [matDir filesep curInhibitor '.mat'];
    load(loadFile);
    fprintf('Inhibitor %d of %d\n', iInh, nInhibitors);
    plateDataNBC = plateData;
    %saveFile = [matDir filesep curInhibitor '_BCrem' num2str(mwWindow) '.mat'];
    nSamples = size(plateData, 1);
    bcData = [plateData(:, bcInds) ones(nSamples, 1)];
    parfor iVar = 1:nVariables
        % find which BCs to regress
        if ~isIncluded(iVar)
            continue;
        end
       % fprintf('Checking BC markers [%s] with mw [%s] on variable %s with mw %d\n', [headers{whichBCs}], num2str(mw(whichBCs)), headers{iVar}, mw(iVar));
        % for iBC = whichBCs
        % plot BC vs variable
         varData= plateData(:, iVar);
         [beta_coeffs(iInh, iVar, :), bint, r, rint, r_stats{iInh, iVar}] = regress(varData, bcData);
        % XXX Print some stuff about which regressors turned out important;
        % Replace variable data with residuals.
        plateDataNBC(:, iVar) = r;
    end% end for iVar
    saveFile = [matDir filesep curInhibitor '_nbc.mat'];
    save(saveFile, 'plateDataNBC');
end