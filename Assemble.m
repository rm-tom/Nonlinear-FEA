function GlobalMatrix = Assemble(flag, varargin)

    % The assembly function does what it is called, it assembles the global
    % matrix or vector from the constituent element matrices. 
    
    persistent Nodes
    persistent Els
    persistent NumNodes
    persistent NumEls
    persistent dof
    
    if isempty(NumEls)
        [Nodes, Els] = MakeNodes();
        NumNodes = size(Nodes, 1);
        NumEls = size(Els, 1);
        dof = 2;
    end
    
    % Initialize Dimensions
    switch flag
        case 'Mass'
            Ys = 1;
            u = zeros(dof*NumNodes, 1);
            Const = 1;
        case 'StrDiv'
            Ys = 0;
            u = varargin{1};
            Const = 1;
        case 'Force'
            Ys = 0;
            u = varargin{1};
            Const = 1;
        case 'Stiff'
            Ys = 1;
            u = varargin{1};
            Const = varargin{2};
        otherwise
            error('Error: Unknown Assemble ID');
    end
    
    % Initialize Matrices
    switch Ys
        case 0
            GlobalMatrix = zeros(dof*NumNodes,1);
        case 1
            GlobalMatrix = zeros(dof*NumNodes,dof*NumNodes);
        otherwise
            error('Error: Unknown size of Ys');
    end
    
    for i = 1:NumEls
    
        %Area
        CurNodes = Els(i,:);
        xLen = Nodes(CurNodes(2),1) - Nodes(CurNodes(1),1);
        yLen = Nodes(CurNodes(4),2) - Nodes(CurNodes(1),2);
        Ar = xLen*yLen;
        DT = [1/xLen, 1/yLen];
        
        % Assemble Matrix
        Gdofs = getD('SortDof', CurNodes);
        ElemMat = getMV(flag, u(Gdofs), DT, Const);
        
        if Ys
            GlobalMatrix(Gdofs, Gdofs) = GlobalMatrix(Gdofs, Gdofs) + Ar*ElemMat;
        else
            GlobalMatrix(Gdofs) = GlobalMatrix(Gdofs) + Ar*ElemMat;
        end
    end
end