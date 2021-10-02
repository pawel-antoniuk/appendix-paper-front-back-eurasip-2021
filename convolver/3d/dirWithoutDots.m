% Paweł Antoniuk 2021
% Bialystok University of Technology

function items = dirWithoutDots(path)
items = dir(path);
items = items(~ismember({items.name},{'.','..'}));
