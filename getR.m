function Res = getR(flag, varargin)

    %Function to get results at certain locations and points
    
    switch flag
        case 'StressState'
            Res = getStressState(varargin{1}, varargin{2}, varargin{3});
        case 'ShowPlot'
            flag2 = varargin{1};
            switch flag2
                case 'mesh'
                    Res = 0;
                    ShowPlot(flag2, varargin{2});
                case 'stre'
                    Res = 0;
                    ShowPlot(flag2, varargin{2}, varargin{3}, varargin{4});
                otherwise
                    error('Error: Unkown showplot flag');
            end
        case 'InteEner'
            Res = getInternalEnergy(varargin{1}, varargin{2});
        otherwise
            error('Error: Unknown Result flag');
    end
    
end

function TotalEner = getInternalEnergy(u,v)

    % Get Mesh
    [Nodes, Els] = MakeNodes();
    NumEls = size(Els, 1);
    
    StrEner = 0;
    KEEner = 0;

    for i = 1:NumEls
        %Area
        CurNodes = Els(i,:);
        xLen = Nodes(CurNodes(2),1) - Nodes(CurNodes(1),1);
        yLen = Nodes(CurNodes(4),2) - Nodes(CurNodes(1),2);
        Ar = xLen*yLen;
        DT = [1/xLen, 1/yLen];

        %Assemble Matrix
        Gdofs = getD('SortDof', CurNodes);
        ElE = getMV('StrEner', u(Gdofs), DT);
        ELK = getMV('KE',v(Gdofs));
        StrEner = StrEner + Ar*ElE;
        KEEner = KEEner + Ar*ELK;
    end
    
    TotalEner = StrEner+KEEner;

end

function ShowPlot(flag, u, varargin)

    [Nodes, ~] = MakeNodes();
    NumNodes = size(Nodes, 1);
    
    if strcmp(flag, 'mesh')
        x = zeros(NumNodes, 1);
        y = zeros(NumNodes, 1);

        for i = 1:NumNodes
            x(i) = Nodes(i,1) + u(2*i-1);
            y(i) = Nodes(i,2) + u(2*i);
        end
        scatter(x,y, '.');

    elseif strcmp(flag, 'stre')
        
        Nums = 10;

        x = linspace(0,1,5*Nums);
        y = linspace(0,0.2,Nums);
        [X,Y]=meshgrid(x,y);
        Z = zeros(5*Nums, Nums);
        
        xp = zeros(5*Nums*Nums, 1);
        yp = zeros(5*Nums*Nums, 1);
        
        cnt = 0;
        for i = 1:length(x)
            for j = 1:length(y)
                cnt = cnt+1;
                xp(cnt) = x(i);
                yp(cnt) = y(j);
            end
        end
        Sn = getR('StressState',u,xp, yp);
        
        cnt = 0;
        for i = 1:length(x)
            for j = 1:length(y)
                cnt = cnt+1;
                Z(i,j) = Sn(cnt, (varargin{1}-1)*2 + varargin{2});
            end
        end
        surf(X,Y,Z');
   
    end
end

% Get stress state at a particular point in the body
function res = getStressState(u, x, y)

    [Nodes, Els] = MakeNodes();
    NumEls = size(Els, 1);
    

    NumVals  = length(x);
    res = zeros(NumVals, 4);
    
    for cn = 1:NumVals
        for i = 1:NumEls
            CurNodes = Els(i,:);
            xLen = Nodes(CurNodes(2),1) - Nodes(CurNodes(1),1);
            yLen = Nodes(CurNodes(4),2) - Nodes(CurNodes(1),2);
            DT = [1/xLen, 1/yLen];

            xn = Nodes(CurNodes,1);
            yn = Nodes(CurNodes,2);

            xmin = min(xn);
            xmax = max(xn);
            ymin = min(yn);
            ymax = max(yn);

            if x(cn) >= xmin && x(cn) <= xmax && y(cn) >= ymin && y(cn) <= ymax

                Gdofs = getD('SortDof', CurNodes);
                us = getD('ConvV', u(Gdofs));

                xs = (x(cn) - xmin)/(xmax - xmin);
                ys = (y(cn) - ymin)/(ymax - ymin);

                Sp = getD('S',us, xs,ys,DT);
                res(cn, :) = reshape(Sp,4,1);
            end
        end
    end
    
end