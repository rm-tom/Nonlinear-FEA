function [Nodes, Els] = MakeNodes(pts, GridxIP, GridyIP)

    persistent xmin
    persistent xmax
    persistent ymin
    persistent ymax
    persistent Gridx
    persistent Gridy

    if nargin == 0
        if isempty(xmin)
            error('Error: Nodes not Initialized');
        end
        
    else
        xmin = pts(1);
        ymin = pts(2);
        xmax = pts(3);
        ymax = pts(4);
        
        Gridx = GridxIP;
        Gridy = GridyIP;
    end
    xv = linspace(xmin, xmax, Gridx+1);
    yv = linspace(ymin, ymax, Gridy+1);

    NumNodes = (Gridx + 1)*(Gridy+1);
    Nodes = zeros(NumNodes, 2);

    cnt = 0;
    for i = 1:Gridy+1
        for j = 1:Gridx+1
            cnt = cnt+1;
            Nodes(cnt,1) = xv(j);
            Nodes(cnt,2) = yv(i);
        end
    end

    NumEl = Gridx*Gridy;
    Els = zeros(NumEl,4);

    cnt = 0;
    for i = 1:Gridy
        for j = 1:Gridx
            cnt = cnt+1;
            Els(cnt,1) = (i-1)*(Gridx+1) + j;
            Els(cnt,2) = (i-1)*(Gridx+1) + j + 1;
            Els(cnt,3) = i*(Gridx+1) + j + 1;
            Els(cnt,4) = i*(Gridx+1) + j;

        end
    end 
end