function [p, x2, exitflag] = logisticTest(yIdx, xIdx, condvarset, data)

    x = data(:, xIdx);
    y = data(:, yIdx);
    cs = data(:, condvarset); 
    
    if ~isempty(cs)     
        [~, dev1] = glmfit(cs,y,'binomial');
        [~, dev2] = glmfit([x cs],y,'binomial');
    else
        [~, dev1] = glmfit(ones(size(data, 1), 1), y, 'binomial', 'constant', 'off');
        [~, dev2] = glmfit(x,y,'binomial');
    end
    
    %calculate the p value and stat.
    x2 = abs(dev1 - dev2);
    p = 1 - chi2cdf(x2, 1);
    exitflag = 1;

end