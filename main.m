clc, clear, close all

% Initialize Nodes
[Nodes, ~] = MakeNodes([0 0 1 0.2],20,4);
NumNodes = size(Nodes, 1);
dof = 2;

% Newmark Method factors
gam = 0.5;
bet = 1/4;

% Materials Constants
rho = 1e-2;
lam = 100;
mu = 40;

%Time Conditions
MaxIter = 200;
umax = 0.15;
dt = 0.001;
MaxTol = 1e-5;

%Contact Nodes | Constrain Function definition
CN = [];
for i = 1:NumNodes
    if abs(Nodes(i,1) - 1) < 1e-10
        CN = [CN; i];
    end 
end
CN = 2*CN-1;
NumC = length(CN);
duk = zeros(dof*NumNodes, NumC);
Gf = max(1, lam/1e6);
for i = 1:NumC
   duk(CN(i), i) = Gf; 
end

% Initial Conditions
Un1 = zeros(dof*NumNodes, 1);
Vn1 = 10*ones(dof*NumNodes, 1);
Vn1(2:2:end) = 0;
An1 = zeros(dof*NumNodes, 1);
Lam = zeros(NumC,1);
Exflag = 0;

%GlobalResult Vectors
Ures = zeros(length(Un1), MaxIter);
Vres = zeros(length(Un1), MaxIter);
Lres = zeros(1, MaxIter);

% Save material constants and Iteration Condiitons
clear NRIteration
clear Assemble
clear getMV
clear getD
save('CVs.mat', 'dt', 'bet', 'MaxTol', 'CN', 'duk', 'umax', 'Gf');
save('Mats.mat', 'rho', 'lam', 'mu');

% Run Global loop (time increment)
for cnt = 1:MaxIter
               
    Un = Un1;
    Vn = Vn1;
    An = An1;

    %Newton Raphson loop
    NRc = 0;
    ConvFlag = 0;
    
    while (NRc<25 && ConvFlag == 0)
        NRc = NRc+1;
        Un2 = Un1;
        Vn2 = Un1;
        An2 = Un1;
        Lam2 = Lam;
        [ConvFlag, Exflag, Un1, Lam, Vn1, An1] = NRIteration(Un1, Un, Vn, An, Exflag, Lam);

        if ConvFlag == 2
            ConvFlag = 0;
            Un1 = Un2;
            Vn1 = Vn2;
            An1 = An2;
        end
    end

    %Save Results
    Ures(:,cnt) = Un1;
    Vres(:,cnt) = Vn1;
    Lres(cnt) = Lam(3);
    
    fprintf('time: %2.4f s | Iteration: %d/%d | Total NR steps %d | Flag %d\n',dt*cnt,cnt,MaxIter,NRc, Exflag );

end

save('Results.mat', 'Ures', 'Vres', 'Lres');