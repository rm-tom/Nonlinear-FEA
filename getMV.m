function ElMat = getMV(flag, varargin)

    %This function calculates the elemental matrices and vectors. This
    %function is typically called by the assembly function. There are size
    %kinds of matrices that can be called. THe elemental mass matrix, force
    %and stress divergence vectors. The stiffness matrix for the Newton
    %Raphson method. It also calculates the strain and kinetic energy of an
    %element.

    switch flag
        case 'Mass'
            ElMat = getELM();
        case 'Force'
            ElMat = getElemF(varargin{1});
        case 'StrDiv'
            ElMat = getElemDiv(varargin{1}, varargin{2});
        case 'Stiff'
            ElMat = getElemK(varargin{1}, varargin{2}, varargin{3});
        case 'KE'
            ElMat = getElemKE(varargin{1});
        case 'StrEner'
            ElMat = getElemStr(varargin{1}, varargin{2});
        otherwise
            error('Error: Unknown getMV flag');
    end

end

%Elemental Strain Energy
function res = getElemStr(u, DT)
    
    u = getD('ConvV', u);
    QPts = getD('Q', 2);
        
    res = 0;
    for k = 1:4
        F = getD('F', u, QPts(k,1), QPts(k,2), DT);
        E = 0.5*(F'*F - eye(2));
        S = getD('S', u, QPts(k,1), QPts(k,2), DT);
        
        E = reshape(E,4, []);
        S = reshape(S,4, []);
        
        res = res + QPts(k,3)*(E'*S)/2;
        
    end
end

%Elemental Kinetic Energy
function res = getElemKE(us)
    
    QPts = getD('Q', 2);
    rho = getD('Constants', 'rho');    
    res = 0;
    for k = 1:4
        Ne1 = getD('Ne1',QPts(k,1),QPts(k,2));
        res = res + rho*QPts(k,3)*norm(Ne1*us)^2/2;
    end
    
end

function ElemK = getElemK(u, DT, Fac)

    %Initialising
    K1 = zeros(8,8);
    K2 = zeros(8,8);
    Fac = getD('Constants', 'rho')/Fac;
    u = getD('ConvV', u);
    QPts = getD('Q', 1);
    NumQ = size(QPts, 1);

    % Forming Matrices
    for k = 1:NumQ
        Ne1 = getD('Ne1',QPts(k,1),QPts(k,2));
        K1 = K1 + Fac*QPts(k,3)*(Ne1'*Ne1);
    end

    QPts = getD('Q', 2);
   
    for k = 1:NumQ
        Be = getD('Be', QPts(k,1), QPts(k,2), DT);
        dPdF = getD('dPdF', u, QPts(k,1), QPts(k,2), DT);
        K2 = K2 + QPts(k,3)*(Be'*dPdF*Be);
    end
    
    ElemK = K1 + K2;
    
end

% Term for the Stress Divergence
function F = getElemDiv(u, DT)
    u = getD('ConvV', u);
    QPts = getD('Q', 2);
    F = zeros(8,1);
    
    for k = 1:4
        P = getD('P', u, QPts(k,1), QPts(k,2), DT);
        P2 = [P(1,1), P(2,2),P(1,2),P(2,1)]';
        B = getD('Be', QPts(k,1), QPts(k,2), DT);
        F = F + QPts(k,3)*(B'*P2);
    end

end

% Elemental Mass matrix
function M = getELM()

    M = zeros(8,8);
    QPts = getD('Q', 1);
    rho = getD('Constants', 'rho');

    for k = 1:4
        Ne1 = getD('Ne1',QPts(k,1),QPts(k,2));
        M = M + rho*QPts(k,3)*(Ne1'*Ne1);
    end  
end

function F = getElemF(t)
    % Nodal Forces
    F1 = zeros(8,1);
    F2 = zeros(8,1);
    t = getD('ConvV', t);
    rho = getD('Constants', 'rho');
    QPts = getD('Q', 2);
   
    for k = 1:4
        Ne1 = getD('Ne1',QPts(k,1),QPts(k,2));
        b =  getD('B',QPts(k,1),QPts(k,2),1);
        F1 = F1 + rho*QPts(k,3)*QPts(k,3)*(Ne1'*b);
    end

    for k = 1:2
        Ne1 = getD('Ne1',1,QPts(k,2));
        pb =  getD('T',1,QPts(k,2),t);
        F2 = F2 + 2*QPts(k,3)*(Ne1'*pb);
    end
    F = F1 + F2;
       
end