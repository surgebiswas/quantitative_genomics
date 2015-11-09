% Read the data
d = dataset('file', 'PS2_data_modified.dat', 'ReadObsNames', false, 'ReadVarNames', true);
d.LocationType = strtrim(d.LocationType);

%%% Problem 1
% Part (a)
[D_af,p_af] = tajimasD(d, 24, 'AF')
[D_ef,p_ef] = tajimasD(d, 23, 'EF')

% Part (b)
% Subset by synonymous vs. nonsynonymous
syn = strcmpi(d.LocationType, 'SYNON');
non = strcmpi(d.LocationType, 'NON-SYN');

[D_af_syn, p_af_syn] = tajimasD(d(syn,:), 24, 'AF')
[D_ef_syn, p_ef_syn] = tajimasD(d(syn,:), 23, 'EF')
[D_af_non, p_af_non] = tajimasD(d(non,:), 24, 'AF')
[D_ef_non, p_ef_non] = tajimasD(d(non,:), 23, 'EF')

% Part (c)
noncode = strcmpi(d.LocationType, '------------');
[D_af_nonc, p_af_nonc] = tajimasD(d(noncode,:), 24, 'AF')
[D_ef_nonc, p_ef_nonc] = tajimasD(d(noncode,:), 23, 'EF')


% Part (d)
% It might be interesting to see if there are gene-specific differences in
% Tajima's D.
genes = unique(d.Gene_Name);
Ds_AF = zeros(length(genes),1);
Ds_EF = zeros(length(genes),1);
for i = 1 : length(genes)
    Ds_AF(i) = tajimasD(d(strcmpi(d.Gene_Name, genes{i}),:), 24, 'AF');
    Ds_EF(i) = tajimasD(d(strcmpi(d.Gene_Name, genes{i}),:), 23, 'EF');
end

figure;
plot(Ds_AF, Ds_EF, 'ok');
hold on
plot([-1.96 1.96], [1.96 1.96], '-r');
plot([-1.96 1.96], -[1.96 1.96], '-r');
plot([1.96 1.96], [-1.96 1.96], '-r');
plot(-[1.96 1.96], [-1.96 1.96], '-r');
axis([-3 3 -3 3]);
xlabel('African', 'FontSize', 18);
ylabel('European', 'FontSize', 18);
title('Tajima''s D, gene-by-gene', 'FontSize', 20);
set(gca, 'FontSize', 16);
axis square
