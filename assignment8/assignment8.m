% read data
vals = textread('data.txt', '%f');

% plot histogram of data
hist(vals)
xlabel('data values')
ylabel('frequency')
title('histogram of data values')
saveas(gcf, 'hist.png')

% qq plot of data
qqplot(vals)
saveas(gcf, 'qq.png')

% sort data for cdf
val_sort = sort(vals);

% calc p values
pd = makedist('Normal', 0, 1);
pvals = 2*cdf(pd, -abs(vals));

% plot histogram of p values
hist(pvals,20)
xlabel('data values')
ylabel('frequency')
title('histogram of p-values')
saveas(gcf, 'hist-p.png')

% calculate False Discovery Rate using Benjamini & Hochberg method
fdr = 0.05;
[pvals_sort, sortIndex] = sort(pvals);
vals_sort2 = val_sort(sortIndex);
for i=1:length(pvals_sort)
    newvals(i,1) = length(pvals_sort)*pvals_sort(i,1)*(1/fdr);
end

r=0;
for i=1:length(newvals)
   if i < newvals(i)
       r = i;
       break;
   end
end

r = r - 1;
p_cutoff = pvals_sort(r);

totalp = 0;
% # pvals less than 0.05
for i=1:length(pvals_sort)
   if pvals_sort(i) < 0.05
       totalp = totalp + 1;
   end
end

% plot raw data as points
% indicate which are below FDR and which have p-value below 0.05
for i=1:length(pvals_sort)
   if i < r
       fdr_cut(i, 1) = vals_sort2(i);
       fdr_cut(i, 2) = pvals_sort(i);
   elseif pvals_sort(i) < 0.05
       pval_cut(i,1) = vals_sort2(i);
       pval_cut(i,2) = pvals_sort(i);
   else
       rest(i,1) = vals_sort2(i);
       rest(i,2) = pvals_sort(i);
   end
       
end

scatter(fdr_cut(:,1), fdr_cut(:,2), 'bo')
hold on
scatter(pval_cut(:,1), pval_cut(:,2), 'rx')
scatter(rest(:,1), rest(:,2), 'g+')
legend('below FDR cutoff', 'p-value below 0.05', 'all other data points')
xlabel('data value')
ylabel('p-value')
title('Data points plotted with respect to their p-values')
saveas(gcf, 'all-data.png')

