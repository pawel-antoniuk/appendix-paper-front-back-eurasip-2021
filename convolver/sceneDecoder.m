% Sławomir Zieliński 2021
% Bialystok University of Technology

function scene = sceneDecoder(n) 
scene = [];
if n == 1
    scene = 'front';    
end
if n == 2
    scene = 'back';   
end
if n == 3
    scene = 'up';   
end
if n == 4
    scene = 'down';    
end
end
