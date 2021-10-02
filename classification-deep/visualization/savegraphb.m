% Paweł Antoniuk 2021
% Bialystok University of Technology

function savegraph(h, filename)
    destpath = fullfile(filename);
%     saveas(gcf, destpath, 'PaperSize', [21 29.7]);
    set(h,'Units','Inches');
    pos = get(h,'Position');
    set(h,'PaperPositionMode','Auto','PaperUnits','Inches','PaperSize',[pos(3), pos(4)])

    print(h, destpath, '-dpdf', '-r0');
    print(h, destpath, '-dpng', '-r300');
end
