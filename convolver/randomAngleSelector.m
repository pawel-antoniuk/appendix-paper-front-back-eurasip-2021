% Sławomir Zieliński 2021
% Bialystok University of Technology

function azi = randomAngleSelector(NofSources,lim1,lim2)
% Calculate random angles between lim1 and lim2
a = lim1;
b = lim2;
azi = zeros(1,NofSources);
for ss = 1:NofSources
    azi(ss) = a + (b-a).*rand(1,1);
end
azi = wrapTo180(azi);