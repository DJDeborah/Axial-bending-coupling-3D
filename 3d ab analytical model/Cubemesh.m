 function [coords,IEN] = Cubemesh()
global nelx nely nelz macroLength macroWidth macroHeight nn 
%% 3d cubic mesh, creating undeformed coordinates and IEN
    % ne = nelx*nely*nelz number of elements xyz
    fullrow = 2*nelx+1;
    halfrow = nelx+1;
    fullplane = fullrow*(nely+1)+halfrow*nely;
    halfplane = halfrow*(nely+1);
    nn = fullplane*(nelz+1)+halfplane*nelz;
%         nn2 = (2*nelx+1)*(nely+1)*(nelz+1)+(nelx+1)*nely*(nelz+1)+...
%             (nelx+1)*(nely+1)*nelz; % total number of nodes
    el = 1; % loop index of element number
    IEN = zeros(nelx*nely*nelz,20);
    for i = 1:nelz
        planeInitialnode = (i-1)*(fullplane+halfplane);
        midplaneInitial = planeInitialnode + fullplane;
        topplaneInitial = midplaneInitial + halfplane;
        for j = 1:nely
            % bottom layer
            brownode = (j-1)*(fullrow+halfrow);
            browInitialnode = planeInitialnode + brownode;
            
            % mid layer
            mrownode = (j-1)*halfrow;
            mrowInitialnode = midplaneInitial + mrownode;

            % top layer
            trownode = (j-1)*(fullrow+halfrow);
            trowInitialnode = topplaneInitial + trownode;
            for k = 1:nelx
                % bottom layer
                IEN(el,[1,9,2]) = browInitialnode + (k-1)*2 +(1:3);
                IEN(el,[12,10]) = browInitialnode + fullrow + k + [0,1];
                IEN(el,[4,11,3]) = browInitialnode + fullrow + halfrow + (k-1)*2 +(1:3);
                
                % middle layer
                IEN(el,[17,18]) = mrowInitialnode + k + [0,1];
                IEN(el,[20,19]) = mrowInitialnode + halfrow+ k + [0,1];
                % top layer
                IEN(el,[5,13,6]) = trowInitialnode + (k-1)*2 +(1:3);
                IEN(el,[16,14]) = trowInitialnode + fullrow + k + [0,1];
                IEN(el,[8,15,7]) = trowInitialnode + fullrow + halfrow + (k-1)*2 +(1:3);
                el = el+1;
            end    
        end
    end

%% Cooords
    % Global undeformed coordinates (x,y,z)
    coords = zeros(nn,3);
    node = 0;
    for h = 1:nelz          % h height
        % full plane
        coords(node+(1:fullplane),3) = (h-1)*macroHeight/nelz;
        for w = 1:nely    % w width
            coords(node+(1:fullrow),1) = linspace(0,macroLength,fullrow);
            coords(node+(1:fullrow),2) = (w-1)*macroWidth/nely;
            node = node+fullrow;
            
            coords(node+(1:halfrow),1) = linspace(0,macroLength,halfrow);
            coords(node+(1:halfrow),2) = (w-0.5)*macroWidth/nely;
            node = node+halfrow;
        end
            % end row of fullplane
                coords(node+(1:fullrow),1) = linspace(0,macroLength,fullrow);
                coords(node+(1:fullrow),2) = macroWidth;
                node = node+fullrow;

        % half plane
        for w = 1:nely
            coords(node+(1:halfrow),1) = linspace(0,macroLength,halfrow);
            coords(node+(1:halfrow),2) = (w-1)*macroWidth/nely;
            coords(node+(1:halfrow),3) = (h-.5)*macroHeight/nelz;
            node = node+halfrow;
        end
            % end row of halfplane
                coords(node+(1:halfrow),1) = linspace(0,macroLength,halfrow);
                coords(node+(1:halfrow),2) = macroWidth;
                coords(node+(1:halfrow),3) = (h-.5)*macroHeight/nelz;       
                node = node + halfrow;
    end

%     Top nodes coordinates
        coords(node+(1:fullplane),3) = macroHeight;
        for w = 1:nely    % w width
            coords(node+(1:fullrow),1) = linspace(0,macroLength,fullrow);
            coords(node+(1:fullrow),2) = (w-1)*macroWidth/nely;
            node = node+fullrow;
            
            coords(node+(1:halfrow),1) = linspace(0,macroLength,halfrow);
            coords(node+(1:halfrow),2) = (w-0.5)*macroWidth/nely;
            node = node+halfrow;
        end
            % end row of fullplane
                coords(node+(1:fullrow),1) = linspace(0,macroLength,fullrow);
                coords(node+(1:fullrow),2) = macroWidth;
%                 node = node+fullrow;

end