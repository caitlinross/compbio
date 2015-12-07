% read data
vals = textread('data.txt', '%f');

% plot histogram of data
hist(vals)

% qq plot of data
qqplot(vals)

% sort data for cdf
val_sort = sort(vals);

% calc p values
pd = makedist('Normal', 0, 1);
pvals = 2*cdf(pd, -abs(vals));

% plot histogram of p values
hist(pvals,20)

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
