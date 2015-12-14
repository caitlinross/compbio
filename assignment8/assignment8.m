% Author: Caitlin Ross
% Due date: 12/15/15
% CSCI 6971
% Assignment 8

clear all
% read data
vals = textread('data.txt', '%f');

% plot histogram of data
figure(1)
hist(vals, 20)
xlabel('data values')
ylabel('frequency')
title('histogram of data values')
saveas(gcf, 'hist.png')

% qq plot of data
figure(2)
qqplot(vals)
saveas(gcf, 'qq.png')

% calc p values
pd = makedist('Normal', 0, 1);
pvals = 2*cdf(pd, -abs(vals));

% plot histogram of p values
figure(3)
hist(pvals,20)
xlabel('data values')
ylabel('frequency')
title('histogram of p-values')
saveas(gcf, 'hist-p.png')

% calculate False Discovery Rate using Benjamini & Hochberg method
fdr = 0.05;
[pvals_sort, sortIndex] = sort(pvals);
vals_sort = vals(sortIndex);
for i=1:length(pvals_sort)
    newvals(i,1) = length(pvals_sort)*pvals_sort(i,1)*(1/fdr);
end

% find r, the cutoff point
r=0;
for i=1:length(newvals)
   if i < newvals(i)
       r = i;
       break;
   end
end

r = r - 1;
p_cutoff = pvals_sort(r);
fprintf('FDR cutoff = 0.05\n');
fprintf('r = %d.\n', r);
fprintf('p-value at FDR cutoff = %f\n', p_cutoff);

totalp = 0;
% find # of pvals less than 0.05
for i=1:length(pvals_sort)
   if pvals_sort(i) < 0.05
       totalp = totalp + 1;
   end
end

fprintf('number of p-values less than 0.05 is %d\n', totalp);
fprintf('number of p-values less than %f is %d\n', p_cutoff, r);

% plot raw data as points
% indicate which are below FDR and which have p-value below 0.05
f_i = 1;
p_i = 1;
r_i = 1;
for i=1:length(pvals_sort)
   if i <= r
       fdr_cut(f_i, 1) = vals_sort(i);
       fdr_cut(f_i, 2) = pvals_sort(i);
       f_i = f_i + 1;
   elseif pvals_sort(i) < 0.05
       pval_cut(p_i,1) = vals_sort(i);
       pval_cut(p_i,2) = pvals_sort(i);
       p_i = p_i + 1;
   else
       rest(r_i,1) = vals_sort(i);
       rest(r_i,2) = pvals_sort(i);
       r_i = r_i + 1;
   end
       
end

figure(4)
scatter(fdr_cut(:,1), fdr_cut(:,2), 'bo')
hold on
scatter(pval_cut(:,1), pval_cut(:,2), 'rx')
scatter(rest(:,1), rest(:,2), 'g+')
legend('below FDR cutoff', 'p-value below 0.05', 'all other data points')
xlabel('data value')
ylabel('p-value')
title('Data points plotted with respect to their p-values')
saveas(gcf, 'all-data.png')

