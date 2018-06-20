aefig = figure;
srfig = figure;
nd = size(dsols, 1);

for i=1:nd
    figure(aefig)
    subplot(nd, 1, i);
    scatter(aesols, dsols(i, :), 'o');
    ylabel(['Offset ', num2str(i)]);
    xlabel('Truncation Error');
    title(['Offsets over Truncation Error for n = ', num2str(n)]);
    
    figure(srfig)
    subplot(nd, 1, i);
    scatter(srsols, dsols(i, :), 'o');
    ylabel(['Offset ', num2str(i)]);
    xlabel('Spectral Radius');
    title(['Offsets over Spectral Radius for n = ', num2str(n)]);
end