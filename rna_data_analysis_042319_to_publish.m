%% MCB-111 project 

%% change the current folder 

cd 'directory'

%% create a table with only data from specific Rs: 

% name of M and Rs:
Rs_names = {'variable_1','variable_2', 'variable_3', 'variable_4', 'variable_5', 'variable_6', 'variable_7', 'variable_8','variable_9' , 'variable_10'};
M_names = {'M_2', 'M_1','F1_M_2','F1_M_1'};

% get the genes names:
genes_names = M_2M_1F1split(:,1);
genes_names = lower(genes_names{:,:});

% create a new table M_2M_1F1split with only the parts we want:
M_2M_1F1split = M_2M_1F1split(:,2:241);
M_2M_1F1split = table2array(M_2M_1F1split);

% find ps genes index:

% ps genes are genes we expect have not been under selection.

p_s = pswithoutproduct(:,1);
p_s = lower(p_s{:,:});

p_s_ind = zeros(length(p_s),1);

for ite=1:length(p_s)
    p_s_ind(ite) = find(contains(genes_names,num2str(p_s(ite))));
end 

% here I also get the index for all genes that are not ps genes:
non_p_s_ind = setdiff(1:length(genes_names),p_s_ind);

% create a table for ps genes and one for non ps
M_2M_1F1split_np = M_2M_1F1split(non_p_s_ind,:);
M_2M_1F1split_p = M_2M_1F1split(p_s_ind,:);

% remove all rows that contain NAs:
M_2M_1F1split_np(any(isnan(M_2M_1F1split_np),2),:) = [];
M_2M_1F1split_p(any(isnan(M_2M_1F1split_p),2),:) = [];

% create a table for each R:
variable_1_np = M_2M_1F1split_np(:,1:24);
variable_2_np = M_2M_1F1split_np(:,25:48);
variable_3_np = M_2M_1F1split_np(:,49:72);
variable_4_np = M_2M_1F1split_np(:,73:96);
variable_5_np = M_2M_1F1split_np(:,97:120);
variable_6_np = M_2M_1F1split_np(:,121:144);
variable_7_np = M_2M_1F1split_np(:,145:168);
variable_8_np = M_2M_1F1split_np(:,169:192);
variable_9_np = M_2M_1F1split_np(:,193:216);
variable_10_np = M_2M_1F1split_np(:,217:240);

variable_1_p = M_2M_1F1split_p(:,1:24);
variable_2_p = M_2M_1F1split_p(:,25:48);
variable_3_p = M_2M_1F1split_p(:,49:72);
variable_4_p = M_2M_1F1split_p(:,73:96);
variable_5_p = M_2M_1F1split_p(:,97:120);
variable_6_p = M_2M_1F1split_p(:,121:144);
variable_7_p = M_2M_1F1split_p(:,145:168);
variable_8_p = M_2M_1F1split_p(:,169:192);
variable_9_p = M_2M_1F1split_p(:,193:216);
variable_10_p = M_2M_1F1split_p(:,217:240);

% make table for each mouse for each R

my_beacons = [1 7 13 19];

for regi=1:length(Rs_names)
    for mou=1:length(M_names)   
        my_name = strcat(Rs_names(regi),'_np');
        my_new_name = strcat(Rs_names(regi),'_',M_names(mou),'_np');
        my_variable = eval(char(my_name));
        my_new_variable = my_variable(:,my_beacons(mou):(my_beacons(mou)+5));
        assignin('base', char(my_new_name),my_new_variable);
    end 
end 

for regi=1:length(Rs_names)
    for mou=1:length(M_names)   
        my_name = strcat(Rs_names(regi),'_p');
        my_new_name = strcat(Rs_names(regi),'_',M_names(mou),'_p');
        my_variable = eval(char(my_name));
        my_new_variable = my_variable(:,my_beacons(mou):(my_beacons(mou)+5));
        assignin('base', char(my_new_name),my_new_variable);
    end 
end

% Here multiply the values of F1 by two as there were not
% normalized

for regi=1:length(Rs_names)
    for mou=3:length(M_names)  
        my_name = strcat(Rs_names(regi),'_',M_names(mou),'_p');
        my_variable = eval(char(my_name));
        my_new_variable = my_variable*2;
        assignin('base', char(my_name),my_new_variable);
    end 
end 
 
for regi=1:length(Rs_names)
    for mou=3:length(M_names)  
        my_name = strcat(Rs_names(regi),'_',M_names(mou),'_np');
        my_variable = eval(char(my_name));
        my_new_variable = my_variable*2;
        assignin('base', char(my_name),my_new_variable);
    end 
end 

%% calculate the means for each R for all M

% by looking at the data we see that the 6 different M replicates have
% very similar results, so I will caluculate the mean and use this as the
% data I will work with:

for regi=1:length(Rs_names)
    for mou=1:length(M_names)
        my_name = strcat(Rs_names(regi),'_',M_names(mou),'_np');
        my_new_name = strcat(Rs_names(regi),'_',M_names(mou), '_mean_np');
        my_variable = eval(char(my_name));
        my_variable_mean = mean(my_variable,2);
        assignin('base', char(my_new_name),my_variable_mean);
    end 
end 

for regi=1:length(Rs_names)
    for mou=1:length(M_names)
        my_name = strcat(Rs_names(regi),'_',M_names(mou),'_p');
        my_new_name = strcat(Rs_names(regi),'_',M_names(mou), '_mean_p');
        my_variable = eval(char(my_name));
        my_variable_mean = mean(my_variable,2);
        assignin('base', char(my_new_name),my_variable_mean);
    end 
end 


%% Calculate fold change 

% One main thing I want to look at is the change in expression between the
% two strains M_1 and M_2

for regi=1:length(Rs_names)
    % for non ps genes
    M_2_name = strcat(Rs_names(regi),'_M_2_mean_np');
    my_M_2 = eval(char(M_2_name));
    % here I add 'ps counts' to prevent inf with 0. I choose 0.0017
    % which is the smallest value in both matrices
    my_M_2 = my_M_2 + 0.0017;
    M_1_name = strcat(Rs_names(regi),'_M_1_mean_np');
    my_M_1 = eval(char(M_1_name));
    my_M_1 = my_M_1 + 0.0017;
    my_new_name_P = strcat(Rs_names(regi),'_P_ratio_np');
    my_new_variable_1 = my_M_2./my_M_1;
    assignin('base', char(my_new_name_P),my_new_variable_1);
    
    % now for ps genes
    M_2_name = strcat(Rs_names(regi),'_M_2_mean_p');
    my_M_2 = eval(char(M_2_name));
    % here I add 'ps counts'
    my_M_2 = my_M_2 + 0.0017;
    M_1_name = strcat(Rs_names(regi),'_M_1_mean_p');
    my_M_1 = eval(char(M_1_name));
    my_M_1 = my_M_1 + 0.0017;
    my_new_name_P = strcat(Rs_names(regi),'_P_ratio_p');
    my_new_variable_1 = my_M_2./my_M_1;
    assignin('base', char(my_new_name_P),my_new_variable_1);
end 

%% Create a table with Pal fold-change for all Rs:

% the table "P_diff_all" will have diff expression for genes in rows
% and colums will corresM_1nd to different brain Rs. 

P_ratio_all_np = zeros(size(M_2M_1F1split_np,1), length(Rs_names));

for regi=1:length(Rs_names)
    my_name = strcat(Rs_names(regi),'_P_ratio_np');
    my_data = eval(char(my_name));
    P_ratio_all_np(:,regi) = my_data;
end 

P_ratio_all_p = zeros(size(M_2M_1F1split_p,1), length(Rs_names));

for regi=1:length(Rs_names)
    my_name = strcat(Rs_names(regi),'_P_ratio_p');
    my_data = eval(char(my_name));
    P_ratio_all_p(:,regi) = my_data;
end 

% I also create here a table with the mean expression for each genes
% between the two species:

P_mean_all_np = zeros(size(M_2M_1F1split_np,1), length(Rs_names));

for regi=1:length(Rs_names)
    my_name = strcat(Rs_names(regi),'_P_mean_np');
    my_data = eval(char(my_name));
    P_mean_all_np(:,regi) = my_data;
end 

P_mean_all_p = zeros(size(M_2M_1F1split_p,1), length(Rs_names));

for regi=1:length(Rs_names)
    my_name = strcat(Rs_names(regi),'_P_mean_p');
    my_data = eval(char(my_name));
    P_mean_all_p(:,regi) = my_data;
end 

% create table with mean expression for all Rs for each M. The
% table will have genes as rows and Rs as colums:

M_1_mean_all_np = zeros(size(M_2M_1F1split_np,1), length(Rs_names));

for regi=1:length(Rs_names)
    my_name = strcat(Rs_names(regi),'_M_1_mean_np');
    my_data = eval(char(my_name));
    M_1_mean_all_np(:,regi) = my_data;
end 

M_2_mean_all_np = zeros(size(M_2M_1F1split_np,1), length(Rs_names));

for regi=1:length(Rs_names)
    my_name = strcat(Rs_names(regi),'_M_2_mean_np');
    my_data = eval(char(my_name));
    M_2_mean_all_np(:,regi) = my_data;
end 

M_1_mean_all_p = zeros(size(M_2M_1F1split_p,1), length(Rs_names));

for regi=1:length(Rs_names)
    my_name = strcat(Rs_names(regi),'_M_1_mean_p');
    my_data = eval(char(my_name));
    M_1_mean_all_p(:,regi) = my_data;
end 

M_2_mean_all_p = zeros(size(M_2M_1F1split_p,1), length(Rs_names));

for regi=1:length(Rs_names)
    my_name = strcat(Rs_names(regi),'_M_2_mean_p');
    my_data = eval(char(my_name));
    M_2_mean_all_p(:,regi) = my_data;
end 
%% Look at overall expression:

% here I want to look at the overall expression. We should expect both
% Ps M_1 and M_2 to have about the same overall expression (this is
% a way to check that the normalization was made correclty

% looks like right

figure
scatter(log(M_1_mean_all_np(:,1)),log(M_2_mean_all_np(:,1)),[],[0.4940, 0.1840, 0.5560])
xlabel('Log of M_1 expression')
ylabel('Log of M_2 expression')
title('Expression')
hold on
scatter(log(M_1_mean_all_p(:,1)),log(M_2_mean_all_p(:,1)),[], [0.9290, 0.6940, 0.1250])
plot([-8 10],[-8 10],'k--', 'LineWidth',2)
legend('np Genes','ps Genes', 'One-one Line' ,'Location','northwest') 
hold off
saveas(gcf,'M_2_vs_M_1_overall_expression.png')


%% plot diff expression vs mean expression  

% here I plot the mean expression against the expresion fold change

% I use a cut off to remove 

my_mean = P_mean_all_np(:); 
my_ratio = P_ratio_all_np(:);
my_ratio = my_ratio(my_mean>1);
my_mean = my_mean(my_mean>1);

my_mean_p = P_mean_all_p(:); 
my_ratio_p = P_ratio_all_p(:);
my_ratio_p = my_ratio_p(my_mean_p>1);
my_mean_p = my_mean_p(my_mean_p>1);

figure
subplot(1,2,1)
dscatter(log(my_mean),log(my_ratio))
xlabel('Mean expression level / threshold = 1TPM')
ylabel('Fold change')
title('MA plot np genes')
subplot(1,2,2)
dscatter(log(my_mean_p),log(my_ratio_p))
xlabel('Mean expression level')
ylabel('Fold change')
title('MA plot ps genes / threshold = 1TPM')
%saveas(gcf,'MA_plot.png')

% look at mean expression between ps and np genes:
figure
histogram((log2(my_mean)),'BinEdges',(0:0.2:15),'Normalization','Probability','FaceColor',[0.4940, 0.1840, 0.5560] )
hold on 
histogram((log2(my_mean_p)),'BinEdges',(0:0.2:15),'Normalization','Probability','FaceColor',[0.9290, 0.6940, 0.1250])
xlabel('Log2 mean expression')
ylabel('Probability')
title('Log2 Mean Expression')
legend('Non ps Genes', 'ps Genes')
hold off
saveas(gcf,'mean_expression_all.png')

% now look at the fold change distribution 

figure
histogram((log2(my_ratio)),'BinEdges',(-10:0.2:10),'Normalization','Probability','FaceColor',[0.4940, 0.1840, 0.5560] )
hold on 
histogram((log2(my_ratio_p)),'BinEdges',(-10:0.2:10),'Normalization','Probability','FaceColor',[0.9290, 0.6940, 0.1250])
xlabel('Log2 of fold change TPM')
ylabel('Probability')
title('Fold change TPM distribution')
legend('Non ps Genes', 'ps Genes')
hold off

% abs value

figure
histogram(abs(log2(my_ratio)),'BinEdges',(0:0.2:10),'Normalization','Probability','FaceColor',[0.4940, 0.1840, 0.5560] )
hold on 
histogram(abs(log2(my_ratio_p)),'BinEdges',(0:0.2:10),'Normalization','Probability','FaceColor',[0.9290, 0.6940, 0.1250])
xlabel('Abs Log2 of fold change TPM')
ylabel('Probability')
title('Absolute fold change TPM distribution')
legend('Non ps Genes', 'ps Genes')
hold off

%% Calculate Empirical p-value: 

% Now I will use the ps genes distribution as my 'random' distribution
% and then look at the np genes ditribution. The idea is to look at
% a fold change and see what proportion of the "random" distribution is
% smaller than this value. The smaller the Empirical p-value, the more
% likelly it is that the gene I look at is under stabilizing selection. 

my_pval_stabi = zeros(length(b),1);

for ite=1:length(b)
    my_pval_stabi(ite) = sum(d<=b(ite))/length(d);
end 

figure
histogram(my_pval_stabi,'FaceColor',[0.8500, 0.3250, 0.0980])
alpha(1)
ylabel('Probability')
xlabel('Empirical p value')
title('Empirical p-value for stabilizing selection')

% now I can look at the percentage of genes form the np
% distribution that we can say are under stabilizing selection:
n_genes_stabi = length(my_pval_stabi(my_pval_stabi < 0.05));
p_genes_stabi = (n_genes_stabi-10000)/(length(my_pval_stabi))*100;
p_genes_stabi;

a = mafdr(my_pval_stabi);
b = length(a(a < 0.05));
c = b/(length(my_pval_stabi))*100;
c;

% plot the p-values distribution after FDR correction 
figure
histogram(a,'Normalization','Probability' )
alpha(1)
ylabel('Probability')
xlabel('Empirical p value')
title('Empirical p-value for stabilizing selection')
%% Plot the expression for each M for different Rs 

% as I have the data for the different Rs, and I have so far merge all
% the data, I will now look to see if there are differences between the
% Rs. 

% first for M_1
figure
for regi=1:length(Rs_names)
    my_name = char(Rs_names(regi));
    my_name_np = strcat(Rs_names(regi),'_M_1_mean_np');
    my_data_np = eval(char(my_name_np));
    my_name_p = strcat(Rs_names(regi),'_M_1_mean_p');
    my_data_p = eval(char(my_name_p));
    subplot(2,5,regi)
    histogram(log(my_data_np),'BinEdges',(-10:0.2:10),'Normalization','Probability','FaceColor',[0.4940, 0.1840, 0.5560])
    hold on
    histogram(log(my_data_p),'BinEdges',(-10:0.2:10),'Normalization','Probability','FaceColor',[0.9290, 0.6940, 0.1250])
    legend('np','ps','Location','northwest')
    title(strrep(my_name,'_',' '))
    hold off
end
%saveas(gcf,'mean_expression_per_Rs_M_1.png')

% first for M_2
figure
for regi=1:length(Rs_names)
    my_name = char(Rs_names(regi));
    my_name_np = strcat(Rs_names(regi),'_M_2_mean_np');
    my_data_np = eval(char(my_name_np));
    my_name_p = strcat(Rs_names(regi),'_M_2_mean_p');
    my_data_p = eval(char(my_name_p));
    subplot(2,5,regi)
    histogram(log(my_data_np),'BinEdges',(-10:0.2:10),'Normalization','Probability','FaceColor',[0.4940, 0.1840, 0.5560])
    hold on
    histogram(log(my_data_p),'BinEdges',(-10:0.2:10),'Normalization','Probability','FaceColor',[0.9290, 0.6940, 0.1250])
    legend('np','ps','Location','northwest')
    title(strrep(my_name,'_',' '))
    hold off
end
%saveas(gcf,'mean_expression_per_Rs_M_2.png')

%% look at sum of diff expression between different Rs

my_sum_diff = zeros(1,10);
for ite=1:10
    my_sum_diff(ite) = sum(abs(P_diff_all(:,ite)));
end 

% we can see that the sum is similar for different Rs
figure
bar(my_sum_diff)

%% Plot heatmap:

% non ps genes:
data_np = P_ratio_all_np;
data_np_log2 = log2(data_np);
data_np_log2(isinf(data_np_log2)) = NaN;
data_np_log2(any(isnan(data_np_log2),2),:) = [];

% heatmap before doing anything:
figure
heatmap(data_np_log2, 'GridVisible','off', 'CellLabelColor','none')
colormap(redbluecmap)
ylabel('Gene')
xlabel('Brain Rs')
ax = gca;
ax.XData = ["variable_1" "variable_2" "variable_3" "variable_4" "variable_5" "variable_6" "variable_7" "variable_8" "variable_9" "variable_10"];


% now heatmap once organized according to the first column
figure
heatmap(sortrows(data_np_log2,1), 'GridVisible','off', 'CellLabelColor','none')
colormap(redbluecmap)
ylabel('Gene')
xlabel('Brain Rs')
ax = gca;
ax.XData = ["variable_1" "variable_2" "variable_3" "variable_4" "variable_5" "variable_6" "variable_7" "variable_8" "variable_9" "variable_10"];


% ps genes:
data_p = P_ratio_all_p;
data_p_log2 = log2(data_p);
data_p_log2(isinf(data_p_log2)) = NaN;
data_p_log2(any(isnan(data_p_log2),2),:) = [];

% heatmap before doing anything:
figure
heatmap(data_p_log2, 'GridVisible','off', 'CellLabelColor','none')
colormap(redbluecmap)
ylabel('Gene')
xlabel('Brain Rs')
ax = gca;
ax.XData = ["variable_1" "variable_2" "variable_3" "variable_4" "variable_5" "variable_6" "variable_7" "variable_8" "variable_9" "variable_10"];

% now heatmap once organized according to the first column
figure
heatmap(sortrows(data_p_log2,1), 'GridVisible','off', 'CellLabelColor','none')
colormap(redbluecmap)
ylabel('Gene')
xlabel('Brain Rs')
ax = gca;
ax.XData = ["variable_1" "variable_2" "variable_3" "variable_4" "variable_5" "variable_6" "variable_7" "variable_8" "variable_9" "variable_10"];

%% PCA on Rs 

% I will do PCA on the differential expression table (rows are diff
% expression and colums are different Rs)

% data I will work with:
% data_np_log2 and data_p_log2


[coeff_np,score_np,latent_np,tsquared_np,explained_np,mu_np] = pca(data_np_log2);
[coeff_p,score_p,latent_p,tsquared_p,explained_p,mu_p] = pca(data_p_log2);

d = zeros(10,2);
d(:,1) = explained_np;
d(:,2) = explained_p;

% plot the variance explained by each PCs

figure
b = bar(d);
legend('Non ps genes','ps genes')
b(1).FaceColor = [0.4940, 0.1840, 0.5560];
b(2).FaceColor = [0.9290, 0.6940, 0.1250];
xlabel('PCs')
ylabel('Variance explained')

figure
scatter(score_np(:,1),score_np(:,2),[],[0.4940, 0.1840, 0.5560]);
xlabel('PC-1')
ylabel('PC-2')

figure
scatter(score_p(:,1),score_p(:,2),[],[0.4940, 0.1840, 0.5560]);
xlabel('PC-1')
ylabel('PC-2')

%% PCA on genes 

% I will do PCA on the differential expression table (rows are diff
% expression and colums are different Rs)

% data I will work with:
% data_np_log2 and data_p_log2


[coeff_np,score_np,latent_np,tsquared_np,explained_np,mu_np] = pca(data_np_log2');
[coeff_p,score_p,latent_p,tsquared_p,explained_p,mu_p] = pca(data_p_log2');

d = zeros(9,2);
d(:,1) = explained_np;
d(:,2) = explained_p;

% plot the variance explained by each PCs

figure
b = bar(d);
legend('Non ps genes','ps genes')
b(1).FaceColor = [0.4940, 0.1840, 0.5560];
b(2).FaceColor = [0.9290, 0.6940, 0.1250];
xlabel('PCs')
ylabel('Variance explained')


figure
scatter(score_np(:,1),score_np(:,2),80,[0.4940, 0.1840, 0.5560], 'filled');
xlabel('PC-1')
ylabel('PC-2')
title('PCA against genes - np genes')

figure
scatter(score_p(:,1),score_p(:,2),80,[0.9290, 0.6940, 0.1250], 'filled');
xlabel('PC-1')
ylabel('PC-2')
title('PCA against genes - ps genes')

%%

% everything after here is very exploratory, it's just things I tried but I
% didn't keep for the presentation.

%% Cluster genes 

% here I try to see of I can culster the data where it would form clusters
% relating to different fold change in different Rs. 

% it didn't seem to work much either because I did something wrong either
% beacause as we saw that most expression is similar between Rs, then
% clustering would not reveal anything.

cgo_p = clustergram(data_p_log2,'Standardize','none''Linkage','centroid','Colormap',redbluecmap, 'ColumnLabels',Rs_names);

cgo_all = clustergram(data_p_log2,'Colormap',redbluecmap,'Standardize','none')

cgo_np = clustergram(data_np_log2,'Standardize','none');
set(cgo_np,'Linkage','ward','Colormap',redbluecmap, 'ColumnLabels',Rs_names)

% other clustering I tried:
Z_np = linkage(data_np,'average');
c_np = cluster(Z_np,'Maxclust',20);
[sorted_np, order_np] = sort(c_np);
new_data_np = data_np(order_np,:);

figure
heatmap(log(new_data_np), 'GridVisible','off', 'CellLabelColor','none')
colormap(gray)
ylabel('Gene')
xlabel('Brain Rs')

Z_p = linkage(data_p,'complete');
c_p = cluster(Z_p,'Maxclust',10);
[sorted_p, order_p] = sort(c_p);
new_data_p = data_np(order_p,:);
new_data_p_log = log(new_data_p);
new_data_p_log(any(isnan(new_data_p_log),2),:) = [];

figure
heatmap(new_data_p_log, 'GridVisible','off', 'CellLabelColor','none')
colormap(redbluecmap)
ylabel('Gene')
xlabel('Brain Rs')

figure
dendrogram(Z_p)

data_p = P_ratio_all_norm_p;
data_p_log2 = log2(data_p);
data_p_log2(isinf(data_p_log2)) = NaN;
data_p_log2(any(isnan(data_p_log2),2),:) = [];

cgo = clustergram(data_p_log2,'Standardize','none');
set(cgo,'Linkage','ward','Colormap',redbluecmap, 'ColumnLabels',Rs_names)


Z = linkage(data_p,'ward');
c = cluster(Z,'Maxclust',20);
[a_sorted, a_order] = sort(c);
new_data_p = data_p(a_order,:);

figure
heatmap(log(new_data_p), 'GridVisible','off', 'CellLabelColor','none')
colormap(gray)
ylabel('Gene')
xlabel('Brain Rs')


[idx_np,C_np] = kmeans(data_np,10);
[a_sorted_np, a_order_np] = sort(idx_np);
new_data_np = data_np(a_order_np,:);

[idx,C] = kmeans(data_p,10);
[a_sorted_p, a_order_p] = sort(idx);
new_data_p = data_p(a_order,:);
T = clusterdata(data_np,'cutoff',1.154695);
unique(T)
[a_sorted, a_order] = sort(T);
new_data_np = data_np(a_order,:);


figure
heatmap(log(new_data_np), 'GridVisible','off', 'CellLabelColor','none')
colormap(gray)
ylabel('Gene')
xlabel('Brain Rs')

%% Look at specific genes
drd1_ind = find(contains(genes_names,'drd1'));
drd2_ind = find(contains(genes_names,'drd2'));
drd3_ind = find(contains(genes_names,'drd3'));
drd4_ind = find(contains(genes_names,'drd4'));

%find(contains(genes_names,'slc3a6'))
drd_ind = [transM_1se(drd1_ind), transM_1se(drd2_ind), transM_1se(drd3_ind), transM_1se(drd4_ind)]

% find avp genes
avp_ind = find(contains(genes_names,'gene7200_avp'));
%% plot the DRD genes expression in variable_9 and variable_1

figure
errorbar(1:length(drd_ind),variable_9_M_2_mean(drd_ind,:),variable_9_M_2_std(drd_ind,:), ".", 'MarkerSize',20)
hold on
errorbar(1:length(drd_ind),variable_9_M_1_mean(drd_ind,:),variable_9_M_1_std(drd_ind,:), ".", 'MarkerSize',20)
hold on
errorbar(1:length(drd_ind),variable_9_F1_M_2_mean(drd_ind,:),variable_9_F1_M_2_std(drd_ind,:), ".", 'MarkerSize',20)
hold on
errorbar(1:length(drd_ind),variable_9_F1_M_1_mean(drd_ind,:),variable_9_F1_M_1_std(drd_ind,:), ".", 'MarkerSize',20)
hold off

figure
errorbar(1:length(drd_ind),variable_1_M_2_mean(drd_ind,:),variable_1_M_2_std(drd_ind,:), ".", 'MarkerSize',20)
hold on
errorbar(1:length(drd_ind),variable_1_M_1_mean(drd_ind,:),variable_1_M_1_std(drd_ind,:), ".", 'MarkerSize',20)
hold on
errorbar(1:length(drd_ind),variable_1_F1_M_2_mean(drd_ind,:),variable_1_F1_M_2_std(drd_ind,:), ".", 'MarkerSize',20)
hold on
errorbar(1:length(drd_ind),variable_1_F1_M_1_mean(drd_ind,:),variable_1_F1_M_1_std(drd_ind,:), ".", 'MarkerSize',20)
hold off

%% plot avp exprssion

figure
plot(1,log(variable_5_M_2(avp_ind,:)), ".", 'MarkerSize',20)
hold on
plot(2,log(variable_5_F1_M_2(avp_ind,:)), ".", 'MarkerSize',20)
hold on
plot(3,log(variable_5_F1_M_1(avp_ind,:)), ".", 'MarkerSize',20)
hold on
plot(4,log(variable_5_M_1(avp_ind,:)), ".", 'MarkerSize',20)

%%
errorbar(1:length(drd_ind),variable_1_M_1_mean(drd_ind,:),variable_1_M_1_std(drd_ind,:), ".", 'MarkerSize',20)
hold on
errorbar(1:length(drd_ind),variable_1_F1_M_2_mean(drd_ind,:),variable_1_F1_M_2_std(drd_ind,:), ".", 'MarkerSize',20)
hold on
errorbar(1:length(drd_ind),variable_1_F1_M_1_mean(drd_ind,:),variable_1_F1_M_1_std(drd_ind,:), ".", 'MarkerSize',20)
hold off

%%

diff_M_2_M_1 = variable_9_M_2_mean - variable_9_M_1_mean;

%% plot of the DDC genes: 

% info about ddc:  http://www.informatics.jax.org/marker/MGI:94876

ddc_ind = find(contains(genes_names,'ddc'));
ddc_ind = transM_1se(ddc_ind);

% plot in variable_9
figure
errorbar(1:length(ddc_ind),variable_9_M_2_mean(ddc_ind,:),variable_9_M_2_std(ddc_ind,:), ".", 'MarkerSize',20)
hold on
errorbar(1:length(ddc_ind),variable_9_M_1_mean(ddc_ind,:),variable_9_M_1_std(ddc_ind,:), ".", 'MarkerSize',20)
hold on
errorbar(1:length(ddc_ind),variable_9_F1_M_2_mean(ddc_ind,:),variable_9_F1_M_2_std(ddc_ind,:), ".", 'MarkerSize',20)
hold on
errorbar(1:length(ddc_ind),variable_9_F1_M_1_mean(ddc_ind,:),variable_9_F1_M_1_std(ddc_ind,:), ".", 'MarkerSize',20)
hold off

% plot in variable_1
figure
errorbar(1:length(ddc_ind),variable_1_M_2_mean(ddc_ind,:),variable_1_M_2_std(ddc_ind,:), ".", 'MarkerSize',20)
hold on
errorbar(1:length(ddc_ind),variable_1_M_1_mean(ddc_ind,:),variable_1_M_1_std(ddc_ind,:), ".", 'MarkerSize',20)
hold on
errorbar(1:length(ddc_ind),variable_1_F1_M_2_mean(ddc_ind,:),variable_1_F1_M_2_std(ddc_ind,:), ".", 'MarkerSize',20)
hold on
errorbar(1:length(ddc_ind),variable_1_F1_M_1_mean(ddc_ind,:),variable_1_F1_M_1_std(ddc_ind,:), ".", 'MarkerSize',20)
hold off
%%
figure
plot(sort(diff_M_2_M_1), 'o')

%% try to see if significant difference

% I will try with boostrapping:

for ite=1:1000;
    my_boot_M_2(ite) = mean(datasample(variable_9_M_2(ddc_ind(1),:),6));
end; 
   
for ite=1:1000;
    my_boot_M_1(ite) = mean(datasample(variable_9_M_1(ddc_ind(1),:),6));
end; 
   

figure 
histogram(my_boot_M_2)
hold on 
histogram(my_boot_M_1)
%hold on;
%line([mean(variable_9_M_2(ddc_ind(1),:)), mean(variable_9_M_2(ddc_ind(1),:))], ylim, 'LineWidth', 2, 'Color', 'r');
%hold off

[h,p] = ttest2(variable_9_M_2(ddc_ind(1),:),variable_9_M_1(ddc_ind(1),:))


%% Look at differential expression by brain R.  

a = variable_9_M_2_mean - variable_9_M_1_mean;

figure
plot((log(a)), 'o')


% get the standard deviation 

variable_9_M_2_std = std(variable_9_M_2,[],2)/sqrt(size(variable_9_M_2,2));
variable_9_M_1_std = std(variable_9_M_1,[],2)/sqrt(size(variable_9_M_1,2));
variable_9_F1_M_2_std = std(variable_9_F1_M_2,[],2)/sqrt(size(variable_9_F1_M_2,2));
variable_9_F1_M_1_std = std(variable_9_F1_M_1,[],2)/sqrt(size(variable_9_F1_M_1,2));

variable_1_M_2_std = std(variable_1_M_2,[],2)/sqrt(size(variable_1_M_2,2));
variable_1_M_1_std = std(variable_1_M_1,[],2)/sqrt(size(variable_1_M_1,2));
variable_1_F1_M_2_std = std(variable_1_F1_M_2,[],2)/sqrt(size(variable_1_F1_M_2,2));
variable_1_F1_M_1_std = std(variable_1_F1_M_1,[],2)/sqrt(size(variable_1_F1_M_1,2));


%% create tables with log differential expression by R between the two Ps and between the two F1

for regi=1:length(Rs_names)
    % start with the Ps:
    M_2_name = strcat(Rs_names(regi),'_M_2_mean_np');
    my_M_2 = eval(char(M_2_name));
    M_1_name = strcat(Rs_names(regi),'_M_1_mean_np');
    my_M_1 = eval(char(M_1_name));
    my_new_name_P = strcat(Rs_names(regi),'_log_P_diff_np');
    my_new_variable_1 = log(my_M_2) - log(my_M_1);
    assignin('base', char(my_new_name_P),my_new_variable_1);
    
    M_2_name = strcat(Rs_names(regi),'_M_2_mean_p');
    my_M_2 = eval(char(M_2_name));
    M_1_name = strcat(Rs_names(regi),'_M_1_mean_p');
    my_M_1 = eval(char(M_1_name));
    my_new_name_P = strcat(Rs_names(regi),'_log_P_diff_p');
    my_new_variable_1 = log(my_M_2) - log(my_M_1);
    assignin('base', char(my_new_name_P),my_new_variable_1);
 
end 

% here I just look at the difference

for regi=1:length(Rs_names)
    M_2_name = strcat(Rs_names(regi),'_M_2_mean_np');
    my_M_2 = eval(char(M_2_name));
    my_M_2 = my_M_2 + min([min(my_M_2(my_M_2>0)) min(my_M_1(my_M_1>0))]);
    M_1_name = strcat(Rs_names(regi),'_M_1_mean_np');
    my_M_1 = eval(char(M_1_name));
    my_M_1 = my_M_1 + min([min(my_M_2(my_M_2>0)) min(my_M_1(my_M_1>0))]);
    my_new_name_P = strcat(Rs_names(regi),'_P_diff_np');
    my_new_variable_1 = my_M_2 - my_M_1;
    assignin('base', char(my_new_name_P),my_new_variable_1);
    
    M_2_name = strcat(Rs_names(regi),'_M_2_mean_p');
    my_M_2 = eval(char(M_2_name));
    my_M_2 = my_M_2 + min([min(my_M_2(my_M_2>0)) min(my_M_1(my_M_1>0))]);
    M_1_name = strcat(Rs_names(regi),'_M_1_mean_p');
    my_M_1 = eval(char(M_1_name));
    my_M_1 = my_M_1 + min([min(my_M_2(my_M_2>0)) min(my_M_1(my_M_1>0))]);
    my_new_name_P = strcat(Rs_names(regi),'_P_diff_p');
    my_new_variable_1 = my_M_2 - my_M_1;
    assignin('base', char(my_new_name_P),my_new_variable_1);
end 

% here I want to look at the mean between the two, I will use this to
% weight the fold growth

for regi=1:length(Rs_names)
    M_2_name = strcat(Rs_names(regi),'_M_2_mean_np');
    my_M_2 = eval(char(M_2_name));
    my_M_2 = my_M_2 + min([min(my_M_2(my_M_2>0)) min(my_M_1(my_M_1>0))]);
    M_1_name = strcat(Rs_names(regi),'_M_1_mean_np');
    my_M_1 = eval(char(M_1_name));
    my_M_1 = my_M_1 + min([min(my_M_2(my_M_2>0)) min(my_M_1(my_M_1>0))]);
    my_new_name_P = strcat(Rs_names(regi),'_P_mean_np');
    my_new_variable_1 = (my_M_2+my_M_1)/2;
    assignin('base', char(my_new_name_P),my_new_variable_1);
    
    M_2_name = strcat(Rs_names(regi),'_M_2_mean_p');
    my_M_2 = eval(char(M_2_name));
    my_M_2 = my_M_2 + min([min(my_M_2(my_M_2>0)) min(my_M_1(my_M_1>0))]);
    M_1_name = strcat(Rs_names(regi),'_M_1_mean_p');
    my_M_1 = eval(char(M_1_name));
    my_M_1 = my_M_1 + min([min(my_M_2(my_M_2>0)) min(my_M_1(my_M_1>0))]);
    my_new_name_P = strcat(Rs_names(regi),'_P_mean_p');
    my_new_variable_1 = (my_M_2+my_M_1)/2;
    assignin('base', char(my_new_name_P),my_new_variable_1);
end

%% create a table with Pal diff for all Rs:

% the table "P_diff_all" will have diff expression for genes in rows
% and colums will corresM_1nd to different brain Rs. 

P_diff_all_np = zeros(size(M_2M_1F1split_np,1), length(Rs_names));

for regi=1:length(Rs_names)
    my_name = strcat(Rs_names(regi),'_P_diff_np');
    my_data = eval(char(my_name));
    P_diff_all_np(:,regi) = my_data;
end 

P_diff_all_p = zeros(size(M_2M_1F1split_p,1), length(Rs_names));

for regi=1:length(Rs_names)
    my_name = strcat(Rs_names(regi),'_P_diff_p');
    my_data = eval(char(my_name));
    P_diff_all_p(:,regi) = my_data;
end 

log_P_diff_all_np = zeros(size(M_2M_1F1split_np,1), length(Rs_names));

for regi=1:length(Rs_names)
    my_name = strcat(Rs_names(regi),'_log_P_diff_np');
    my_data = eval(char(my_name));
    log_P_diff_all_np(:,regi) = my_data;
end 

log_P_diff_all_p = zeros(size(M_2M_1F1split_p,1), length(Rs_names));

for regi=1:length(Rs_names)
    my_name = strcat(Rs_names(regi),'_log_P_diff_p');
    my_data = eval(char(my_name));
    log_P_diff_all_p(:,regi) = my_data;
end 

P_ratio_all_np = zeros(size(M_2M_1F1split_np,1), length(Rs_names));

for regi=1:length(Rs_names)
    my_name = strcat(Rs_names(regi),'_P_ratio_np');
    my_data = eval(char(my_name));
    P_ratio_all_np(:,regi) = my_data;
end 

P_ratio_all_p = zeros(size(M_2M_1F1split_p,1), length(Rs_names));

for regi=1:length(Rs_names)
    my_name = strcat(Rs_names(regi),'_P_ratio_p');
    my_data = eval(char(my_name));
    P_ratio_all_p(:,regi) = my_data;
end 

P_mean_all_np = zeros(size(M_2M_1F1split_np,1), length(Rs_names));

for regi=1:length(Rs_names)
    my_name = strcat(Rs_names(regi),'_P_mean_np');
    my_data = eval(char(my_name));
    P_mean_all_np(:,regi) = my_data;
end 

P_mean_all_p = zeros(size(M_2M_1F1split_p,1), length(Rs_names));

for regi=1:length(Rs_names)
    my_name = strcat(Rs_names(regi),'_P_mean_p');
    my_data = eval(char(my_name));
    P_mean_all_p(:,regi) = my_data;
end 

P_ratio_all_norm_np = P_ratio_all_np.*P_mean_all_np;
P_ratio_all_norm_p = P_ratio_all_p.*P_mean_all_p;

P_ratio_all_norm_np(isinf(P_ratio_all_norm_np)) = NaN; 
P_ratio_all_norm_p(isinf(P_ratio_all_norm_p)) = NaN;

P_diff_all_norm_np = P_diff_all_np.*P_mean_all_np;
P_diff_all_norm_p = P_diff_all_p.*P_mean_all_p;

P_diff_all_norm_np(isinf(P_diff_all_norm_np)) = NaN; 
P_diff_all_norm_p(isinf(P_diff_all_norm_p)) = NaN;

log_P_diff_all_norm_np = log_P_diff_all_np.*P_mean_all_np;
log_P_diff_all_norm_p = log_P_diff_all_p.*P_mean_all_p;

log_P_diff_all_norm_np(isinf(log_P_diff_all_norm_np)) = NaN; 
log_P_diff_all_norm_p(isinf(log_P_diff_all_norm_p)) = NaN;