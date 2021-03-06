clear; 

% Read the data
d = dataset('file', 'PS2_data_modified.txt', 'ReadObsNames', false, 'ReadVarNames', true);
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

[~, p_AF] = ttest(Ds_AF)
[~, p_EF] = ttest(Ds_EF)



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
plotSave('TajimasD_genebygene.png');
close


%%% Problem 2
% Get the BLOSUM62 matrix
B = blosum62;

code = strcmpi(d.LocationType, 'SYNON') | strcmpi(d.LocationType, 'NON-SYN');
dc = d(code,:);
O = zeros(size(B));
for i = 1 : size(dc,1)
    if ~strcmpi(dc.aa1{i}, '***')
        from = aminolookup('Abbreviation', dc.aa1{i});
        to = aminolookup('Abbreviation', dc.aa2{i});

        i1 = aa2int(to(1));
        i2 = aa2int(from(1));

        O(i1,i2) = O(i1,i2) + 1;
    end
    
end

figure;
hold on
plot(reshape(B,24*24,1), reshape(O, 24*24,1), 'ok')
axis square
xlabel('BLOSUM62 score', 'FontSize', 18);
ylabel('Observed number of transitions', 'FontSize', 18);
set(gca, 'FontSize', 16);
box on
plotSave('BLOSUM62_vs_observed.png');
close




