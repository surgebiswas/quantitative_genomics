function [ D, pval ] = tajimasD( d, n, race )
% tajimasD - Calculates Tajima's D for the PS2 datset.
%
% d = dataset structure representing PS2_data_modified.dat or a subset of
% it.
% n = number of DNA sequences in the sample.
% race = 'AF' or 'EF'

a1 = sum( 1./(1:n-1) );
a2 = sum( 1 ./ ((1:n-1).^2) );
b1 = (n + 1)/(3*(n-1));
b2 = 2*(n^2 + n + 3)/(9*n*(n-1));
c1 = b1 - 1/a1;
c2 = b2 - (n+2)/(a1*n) + a2/(a1^2);
e1 = c1/a1;
e2 = c2/(a1^2 + a2);

% pi = number of sites that differ. 
switch race
    case 'AF'
        freq = d.AF;
    case 'EF'
        freq = d.EF;
end
pi = 2*n*sum( freq.*(1 - freq) )/(n - 1);
S = sum(freq > 0 & freq < 1);

D = (pi - S/a1)/sqrt( e1*S + e2*S*(S-1) );
pval = 2*min(normcdf(D, 0, 1), 1 - normcdf(D,0,1));









end

