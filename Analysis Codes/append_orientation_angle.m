function [assay] = append_orientation_angle(assay)
% function [assay] = orientation_angle(assay)
% 
% Dinesh Natesan
% 5th Oct, 2017

treatments = fieldnames(assay);

for i=1:length(treatments)
    
    for j=1:length(assay.(treatments{i}).trajData)
        
        assay.(treatments{i}).trajData{j,1}.globalpara.OrientationAngle = ...
            orientation_angle(assay.(treatments{i}).trajData{j,1}.expdata.pts);
        
    end    
    
end

% Done. Exit

end

function [angle] = orientation_angle(pts)
% function [] = orientation_angle(pts)
% 
% 

% Obtain 2D vectors
vectors_2D = diff(pts);
vectors_2D(:,3) = 0;

% xaxis reference 
reference = zeros(size(vectors_2D));
reference(:,1) = 1;

% Find the angle w.r.t the x-axis (with sign)
dot_product = dot(unitizeVect(vectors_2D), reference, 2);
cross_product = cross(unitizeVect(vectors_2D), reference, 2);

angle = sign(cross_product(:,3)).* ...
    atan2(matrixNorm(cross_product), dot_product) .* 180/pi;


end
