function OParg = getD(flag, varargin)

    %This function makes getting the data easy, for example, it calculates
    %the F/S/P at a point within the element. It also contains the function
    %to calculate the shape functions and the material constants.

    switch flag
        case 'N'
            OParg = getN(varargin{1}, varargin{2}, varargin{3});
        case 'Q'
            OParg = getQPts(varargin{1});
        case 'DN'
            OParg = getDN(varargin{1}, varargin{2}, varargin{3}, varargin{4}, varargin{5});
        case 'B'
            OParg = getB(varargin{1}, varargin{2}, varargin{3});
        case 'T'
            OParg = getT(varargin{1}, varargin{2}, varargin{3});
        case 'ConvV'
            OParg = ConvV(varargin{1});
        case 'Constants'
            OParg = Constants(varargin{1});
        case 'H'
            OParg = getH(varargin{1}, varargin{2}, varargin{3}, varargin{4});
        case 'F'
            OParg = getF(varargin{1}, varargin{2}, varargin{3}, varargin{4});
        case 'S'
            OParg = getS(varargin{1}, varargin{2}, varargin{3}, varargin{4});
        case 'P'
            OParg = getP(varargin{1}, varargin{2}, varargin{3}, varargin{4});
        case 'Be'
            OParg = getBe(varargin{1}, varargin{2}, varargin{3});
        case 'dPdF'
            OParg = getdPdF(varargin{1}, varargin{2}, varargin{3}, varargin{4});
        case 'Ne1'
            OParg = getNe1(varargin{1}, varargin{2});
        case 'SortDof'
            OParg = Sdof(varargin{1});
        otherwise
            error('Error: Unknown getD flag');
    end
    
end

function Ne1 = getNe1(xp,yp)
    N = zeros(4,1);
    for m = 1:4
        N(m) =  getD('N',m,xp,yp);
    end
    Ne1 = [N(1)*eye(2),N(2)*eye(2),N(3)*eye(2),N(4)*eye(2)];
end

function dPdF = getdPdF(u, x,y, DT)
    
    dPdF = zeros(4,4);
    F = getF(u, x, y, DT);
    lam = Constants('lam');
    mu = Constants('mu');
    
    del = @(a,b) eq(a,b);
    FjCFjC = F(1,1)^2 + F(1,2)^2 + F(2,1)^2 + F(2,2)^2;
    Fj = @(a,b) F(1, a)*F(1, b) + F(2, a)*F(2, b);
    FA = @(a,b) F(a, 1)*F(b, 1) + F(a, 2)*F(b, 2);
    
    for i = 1:2
        for B = 1:2
            for k = 1:2
                for D = 1:2
                    xp = (i-1)*2 + B;
                    yp = (k-1)*2 + D;
                    
                    T1 = lam/2*(2*F(k,D)*F(i,B) + FjCFjC*del(i,k)*del(B,D));
                    T2 = mu*(Fj(D,B)*del(i,k) + F(i,D)*F(k,B) + FA(i,k)*del(B,D));
                    T3 = -(lam + mu)*del(i,k)*del(B,D);

                    dPdF(xp,yp) = T1 + T2 + T3;
                    
                end
            end
        end
    end
    
    % This rearrangement is to make sure that the ordering of dPdF is
    % according to the convention P = [p11, p22, p12, p21]
    rev = [1, 4, 2, 3];
    dPdF = dPdF(rev, rev);

end

% sB is the Be Matrix differentiated along a single direction. 
function Be = getBe(x,y, DT)
    Be = zeros(4,8);
    sB = zeros(4,2); 
    for i = 1:4
        sB(1,1) = getDN(i, x,y,1, DT);
        sB(2,2) = getDN(i,x,y,2, DT);
        sB(3,1) = getDN(i,x,y,2, DT);
        sB(4,2) = getDN(i,x,y,1, DT); 
        Be(:,2*i-1:2*i) = sB;
    end
end

function P = getP(u, x, y, DT)
    F = getF(u, x, y, DT);
    S = getS(u, x, y, DT);
    P = F*S;
end

function S = getS(u, x, y, DT)
    lam = Constants('lam');
    mu = Constants('mu');
    
    F = getF(u, x, y, DT);
    E = 0.5*(F'*F - eye(2));
    S = lam*trace(E)*eye(2) + 2*mu*E;
    S(abs(S)<1e-15) = 0;
end

function F = getF(u, x, y, DT)
    H = getH(u, x, y, DT);
    F = eye(2) + H;
end

function H = getH(u, x, y, DT)
    
    H = zeros(2,2);
    for i = 1:2
        for j = 1:2
            for k = 1:4
               H(i,j) = H(i,j) + getDN(k,x,y,j, DT)*u(k,i);
            end
        end
    end
    
    H(abs(H)<1e-10) = 0;
    
end

function res = Constants(what)
    persistent vals
    if isempty(vals)
        vals = importdata('Mats.mat');
    end
    switch what
        case 'rho'
            res = vals.rho;
        case 'lam'
            res = vals.lam;
        case 'mu'
            res = vals.mu;
        otherwise
            error('Error: Unknown constant name');
    end
end

%Converts Values to another format
function t = ConvV(ts)
    t = zeros(4,2);
    t(1, 1) = ts(1);
    t(1, 2) = ts(2);
    t(2, 1) = ts(3);
    t(2, 2) = ts(4);
    t(3, 1) = ts(5);
    t(3, 2) = ts(6);
    t(4, 1) = ts(7);
    t(4, 2) = ts(8);
end

% Traction at a certain point
function res = getT(x,y,t)
    res = zeros(2,1);
    for i = 1:4
        res(1) = res(1) + getN(i,x,y)*t(i,1);
        res(2) = res(2) + getN(i,x,y)*t(i,2);
    end   
end

% Body Force at a point
function res = getB(x,y,t)
    res = zeros(2,1);
    res(2) = -1e7*x*y*t;
    res(2) = 0;
end

function QPts = getQPts(flag)
    switch flag
        case 1
            QPts = [0, 0, 0.25;
                    0, 1, 0.25;
                    1, 0, 0.25;
                    1, 1, 0.25];
        case 2
            QPts = [0.788675134594813, 0.788675134594813, 0.25;
                    0.211324865405187, 0.211324865405187, 0.25;
                    0.211324865405187, 0.788675134594813, 0.25;
                    0.788675134594813, 0.211324865405187, 0.25];
        otherwise
            disp('Error: Unknown QPts flag');
    end

end

function res = getDN(i,x,y,d, DT)
    
    % d = order of differentiation
    res = 0;
    switch d
        case 1
            switch i
                case 1
                    res = -(1-y);
                case 2
                    res = (1-y);
                case 3
                    res = y;
                case 4
                    res = -y;
            end
        case 2
            switch i
                case 1
                    res = -(1-x);
                case 2
                    res = -x;
                case 3
                    res = x;
                case 4
                    res = (1-x);
            end   
    end
    res = res*DT(d);
end

function res = getN(i,x,y)
    res = 0;
    switch i
        case 1
            res = (1-x)*(1-y);
        case 2
            res = x*(1-y);
        case 3
            res = x*y;
        case 4
            res = (1-x)*y;
    end   
end

%Function to sort dofs
function res = Sdof(CurNodes)
    dof = [2*CurNodes,2*CurNodes-1];
    dof = sort(dof);
    res(1:4) = dof(1:4);
    res(5:6) = dof(7:8);
    res(7:8) = dof(5:6);

end