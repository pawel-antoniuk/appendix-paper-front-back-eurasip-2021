% Paweł Antoniuk 2021
% Bialystok University of Technology

function [bestFitAngles,bestFitAnglesI] = findBestFitAnglePairs(...
    anglePairsToFit, desiredAnglePairs, directionAngles, ...
    capWidth, capWidthEps)

desDirAngSub = angleDistance(desiredAnglePairs,directionAngles);
anglesInAreaIdx = find(all(desDirAngSub <= capWidth + capWidthEps, 2));
desiredAnglePairs = desiredAnglePairs(anglesInAreaIdx, :);

% bestFitAngles = zeros(size(anglePairsToFit));

a = repelem(anglePairsToFit, size(desiredAnglePairs, 1), 1);
b = repmat(desiredAnglePairs, size(anglePairsToFit, 1), 1);
dist1 = angleDistance(a, b);
dist2 = reshape(dist1, size(desiredAnglePairs, 1), size(anglePairsToFit, 1));
[~,minDistI] = min(dist2, [], 1);
bestFitAngles = desiredAnglePairs(minDistI, :);
bestFitAnglesI = anglesInAreaIdx(minDistI);


% for iAnglePairsToFit = 1:size(anglePairsToFit, 1)
%     anglePairToFit = anglePairsToFit(iAnglePairsToFit, :);
%     minAngDiff = Inf;
%     
%     for iDesiredAnglePairs = 1:size(desiredAnglePairs, 1)
%         desiredAnglePair = desiredAnglePairs(iDesiredAnglePairs, :);
%         
% %         angDiff = angSub(anglePairToFit, desiredAnglePair);
% %         if angDiff < minAngDiff
% %             bestFitAngles(iAnglePairsToFit, :) = desiredAnglePair;
% %             minAngDiff = angDiff;
% %         end
%         
% %         anglePairToFitRad = deg2rad(anglePairToFit);
% %         desiredAnglePairRad = deg2rad(desiredAnglePair);
% %         [ax,ay,az] = sph2cart(anglePairToFitRad(1), anglePairToFitRad(2), 1);
% %         [dx,dy,dz] = sph2cart(desiredAnglePairRad(1), desiredAnglePairRad(2), 1);
% %         dist = pdist([ax ay az; dx dy dz]);       
%         dist = distance('gc',flip(anglePairToFit),flip(desiredAnglePair));
%         
%         if dist < minAngDiff
%             bestFitAngles(iAnglePairsToFit, :) = desiredAnglePair;
%             minAngDiff = dist;
%         end
%     end
% end
% 
% if ~all(bestfit2 == bestFitAngles, 'all')
%     error('test');
% end

end

% for iAngle = 1:length(anglesToFit)
%     minAngDiff = 720;
%     angleToFit = anglesToFit(iAngle);
%     for desiredAngle = desiredAngles
%         angDiff = abs(rad2deg(angdiff(deg2rad(angleToFit), deg2rad(desiredAngle))));
%         if angDiff < minAngDiff
%             bestFitAngles(iAngle) = desiredAngle;
%             minAngDiff = angDiff;
%         end
%     end
% end

function d = angleDistance(a, b)
d = distance('gc', flip(a, 2), flip(b, 2));
end
