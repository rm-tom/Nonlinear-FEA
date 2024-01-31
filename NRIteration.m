function [Rflag, Exflag, Un1, Lam, Vn1, An1] = NRIteration(Un1, Un, Vn, An, Exflag, Lam)

    % This is a function to do one iteration of Newton Raphson Method
    % Output flag are
    %   Rflag = This values returns 1 when the residual is below the
    %   required tolerance
    %   Exflag = This value return 1 when the next iteration should impose
    %   the contrainst

    % Debugging outputs
    %Op1 = sprintf('\t -> %d,  %e',Exflag,Lam(1));
        
    %Defining persistent matrices
    persistent GlobalM
    persistent CV
        
    if isempty(GlobalM)
        %GlobalM = getMassM();
        GlobalM = Assemble('Mass');
    end
    if isempty(CV)
        CV = importdata('CVs.mat');
    end
    
    %Initializing values
    Rflag = 0;
    if Exflag 
        Vn(CV.CN) = 0;
    end
    An1 = 1/CV.bet/CV.dt^2*(Un1 - Un - Vn*CV.dt) - An*(1-2*CV.bet)/(2*CV.bet);
    
    % Forming Residual
    if Exflag
        Res1 = GlobalM*An1 + Assemble('StrDiv', Un1) + CV.duk*Lam;
        Res2 = CV.Gf*(Un1(CV.CN) - CV.umax);
        Res = [Res1;Res2];
    else
        Res = GlobalM*An1 + Assemble('StrDiv', Un1);
    end
    
    %fprintf('\t\t\t\t Residual %e \n',norm(Res));
        
    %Solving Normal System or Extended System
    if (norm(Res)<CV.MaxTol)
        Rflag = 1;
    else
        GlobalK = Assemble('Stiff',Un1, CV.bet*CV.dt^2);
        if Exflag
            K = [GlobalK, CV.duk;CV.duk', zeros(length(CV.CN),length(CV.CN))];
            Dels = -K\Res;
            Un1 = Un1 + Dels(1:length(Un1));
            Lam = Lam + Dels(length(Un1)+1:end);
            
        else
            du = GlobalK\(-Res);
            Un1 = Un1 + du;
        end
    end
    
    %Output Values
    An1 = 1/CV.bet/CV.dt^2*(Un1 - Un - Vn*CV.dt) - An*(1-2*CV.bet)/(2*CV.bet);
    Vn1 = Vn + CV.dt/2*(An + An1);
    
    %Update Exflag (Extended system yes/no flag)
    
    %Positive Pressure yes/no (Prflag
    Prflag = 0;
    if max(Lam) >= 0
        Prflag = 1;
    end
    
    %Penetration yes/no (Peflag)
    g = Un1(CV.CN) - CV.umax;
    Peflag = 0;
    if max(g) >=0
        Peflag = 1;
    end

     g =-g;
% 
%     if max(g)>0
%         Peflag = 2;
%     end
    
%     if Peflag == 2
%         ExflagNew = 1;
%     elseif Prflag == 1 || Peflag == 1
%         ExflagNew = 1;
%     elseif Prflag == 0
%         ExflagNew = 0;
%     end

    if min(g)<-1e-10 || (max(Lam) > 0 && min(g)<1e-10)
        ExflagNew = 1;
    else
        ExflagNew =0;
    end

%     if Prflag == 1 && Peflag == 1
%         ExflagNew = 1;
%     else
%         ExflagNew = 0;
%     end

    
    
    if abs(ExflagNew - Exflag) < 1e-10
        Exflag = ExflagNew;
    else
        Rflag  = 2;
        Exflag = ExflagNew;
    end

    % Debugging outputs
    %Op2 = sprintf('\t %d,  %e, %e',Exflag,max(Lam));
    %disp(strcat(Op1, Op2));
    
end