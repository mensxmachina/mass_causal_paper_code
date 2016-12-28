function [p, r, exitflag] = fisher(var1, var2, condvarset,data, domain_counts)
%
% Fisher Z-test
% =============
%
% p-value = Probability of accepting H0 (null hypothesis, i.e.
%           that var1 and var2 are independent, given condvarset).
%           If p-value > 0.05 => Accept H0, i.e. vars are independent (rho=0)
%           If p-value <= 0.05 => Reject H0, i.e. vars are dependent (rho~=0)
%
%**********************************************************************
% Copyright (C), DSL 2002
%**********************************************************************

% Compute partial correlation coefficient
S=corrcoef(data(:,[var1 var2 condvarset]));

condvarset=3:length(condvarset)+2;
vars = [1 2];
S2 = S(vars,vars) - S(vars,condvarset)*inv(S(condvarset,condvarset))*S(condvarset,vars);
r = abs(S2(1,2) / sqrt(S2(1,1) * S2(2,2)));

% Sample size
N=size(data,1);

% Compute Z and P-value
z = 0.5*log( (1+r)/(1-r));
df=N - length(condvarset) - 3;
W = sqrt(df)*z; 
p = 2*tcdf(-abs(W),df);

exitflag=1; % Dummy variable indicated reliablity of computation of p-value



