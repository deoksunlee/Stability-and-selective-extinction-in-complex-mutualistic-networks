function Xout = LVstationary6(Amatrix,c,m,Tmax,ExtinctTh)

% Amatrix NP x NA matrix representing mutualistic interactions
% c : competition strength
% m : mutualism strength
% Tmax : time length for numerical integration
% ExtinctTh (set to 10^-5) : If X(i)<ExtinctTh, species i is considered extinct

% LV equation : dx_i/dt = x_i (1 + B_ij x_j) (i=1,...,S = NP+NA)
% Interaction matrix : B = - I - c (J^(0) - I) + m A
% Competition: J^(0) = J^PP '+' J^AA (direct sum)
% Mutualism: A = Amatrix '+' Amatrix' (for PA and AP, respectively)

% Xout(1,1:S) : Xinf from solving numerically LV with Amatrix
% Xout(2,1:S) : Xstar from -B^{-1} with no species extinct 
% Xout(3,1:S) : Xstareff from -(Beff)^-1 with selected species extinct
% Xout(4,1:S) : Xtinf from solving LV with Atmatrix (annealed)
% Xout(5,1:S) : Xtstar from -tildeB^{-1} with no species extinct 
% Xout(6,1:S) : Xtstareff from -(tildeBeff)^-1 with selected species extinct


[NP NA] = size(Amatrix);
S = NP + NA;

% Growth factor %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
GF = ones(S,1);

% B0matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B0matrix = -(1-c)*eye(S,S)- c*ones(S,S);
B0matrix(1:NP, NP+(1:NA))=0;
B0matrix(NP+(1:NA),1:NP)=0;

% Bmatrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Bmatrix = B0matrix;
Bmatrix(1:NP,NP+(1:NA))=m*Amatrix;
Bmatrix(NP+(1:NA),1:NP)=m*Amatrix';

% Atmatrix (Annealed approximation) %%%%%%%%%%%%
kP = sum(Amatrix,2);
kA = sum(Amatrix,1);
L = sum(kP);
Atmatrix = kP*kA/L;

% Btmatrix (in terms of annealed A)%%%%%%%%%%%%%
Btmatrix = B0matrix;
Btmatrix(1:NP,NP+(1:NA))=m*Atmatrix;
Btmatrix(NP+(1:NA),1:NP)=m*Atmatrix';

Bcell = {Bmatrix, Btmatrix};
Acell = {Amatrix,Atmatrix};

for i=1:2  % original and annealed
   
    B = Bcell{i};
    A = Acell{i};
    % Computation: Xinf/Xtinf %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    X0 = ones(1,S);
    X=X0;
    odefun(X) = @(t,X) X.*(GF+ B*X);
    tspan = [0 Tmax];
    %[t,X]=ode45(odefun,tspan,X0);
    [t,X]=ode15s(odefun,tspan,X0);
    row = 1+3*(i-1);
    Xout(row,:) = X(end,:);
    
    % Computation: Xstar /Xtstar %%%%%%%%%%%%%%%%%%%%%%%%%%
    Binv = inv(B);
    row = 2+3*(i-1);
    Xout(row,:) = -(Binv*GF)';

    % Computation: Xstareff / Xtstareff %%%%%%%%%%%%%%%%%%%
    Xeff = Xout(row,:);
    Surv = 1:S;
    NewExtinct = [];
    
    while true
        % Xstareff updated %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Xeff(NewExtinct)=0;
        Surv=setdiff(Surv,NewExtinct);
        % Surv is always non-empty.
        
        Beff = B(Surv,Surv);
        GFeff =GF(Surv);
        Beffinv = inv(Beff);
        Xeff(Surv)=-(Beffinv*GFeff)';
        
        % NewExtinct updated %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        SurvP = Surv(Surv<=NP);
        SurvA = Surv(Surv>NP);
        Aeff = A(SurvP,SurvA - NP);
        
        kPeff = sum(Aeff,2);
        kAeff = sum(Aeff,1);
        [~,Pkmax] = max(kPeff); 
        [~,Akmax] = max(kAeff);
        
        Xref = [repmat(Xeff(SurvP(Pkmax)),[1 NP]) ...
            repmat(Xeff(SurvA(Akmax)),[1 NA])];
        
        NewExtinct = find(Xeff.*Xref<0);
        
        if isempty(NewExtinct)
            break;
        end
        
    end
    row = 3+3*(i-1);
    Xout(row,:) = Xeff;

end
end
