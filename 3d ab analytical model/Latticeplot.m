function [] = Latticeplot(coords,IEN)



conVec = [1,9,2,10,3,11,4,12,1,...
         1,17,5,16,8,20,4,...
         4,20,8,15,7,19,3,...
         3,19,7,14,6,18,2,...
         2,18,6,13,5,17,1];

% Connection vector of unit cube
    figure()
    ne = size(IEN,1);   % number of elements
    for e = 1:ne
        coord = coords(IEN(e,conVec),:);
        plot3(coord(:,1),coord(:,2),coord(:,3),'-g')
        hold on
%         % plot element skeleton
%         coordD = zeros(size(conVec,2),3);
        
    end
 end

