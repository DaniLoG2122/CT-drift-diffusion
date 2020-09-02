function solstruct=DriftanalyseCT(sol)

v2struct(sol.params);
x=sol.x;
t=sol.t;
xpoints=length(x);%%%%%%%%%%%%%%%%%%%%careful
xnm=1e7*x;
xmax=max(x);
for ii=1:1:layers
    if(ii==1)
        if(symm==1)
            BM=ones(length(t), xpoints)*diag(x <= layer{1}.tp)+ones(length(t), xpoints)*diag(x >= xmax-layer{1}.tp);
        else
            BM=ones(length(t), xpoints)*diag(x <= layer{1}.tp);
        end
        nstatM = (-layer{1}.NA)*BM  + (layer{1}.ND)*BM   ;
        niM = layer{1}.ni*BM ;
        EAM = layer{1}.EA*BM ;
        IPM = layer{1}.IP*BM ;
        EiM = layer{1}.Ei*BM  ;
        mu_e=layer{1}.mue*BM ;
        mu_h=layer{1}.mup*BM  ;
        N0CM = layer{1}.N0C*BM ;
        N0VM = layer{1}.N0V*BM ;
        Krad=layer{1}.krad*BM;
        Ni = layer{1}.ni*BM ;
        Ni2 = layer{1}.ni*layer{1}.ni*BM ;
        Nt = layer{1}.nt*BM ;
        Pt = layer{1}.pt*BM ;
        Taun=layer{1}.taun*BM ;
        Taup=layer{1}.taup*BM ;
        Cnp=layer{1}.Cn*BM;
        enp=layer{1}.en*BM;
    else
        if(symm==1)
            BM=ones(length(t), xpoints)*diag(layer{ii}.XL<=x & x <= layer{ii}.XR)+ones(length(t), xpoints)*diag(xmax-layer{ii}.XL>=x & x > xmax-layer{ii}.XR);
        else
            BM=ones(length(t), xpoints)*diag(layer{ii}.XL<=x & x <= layer{ii}.XR);
        end
        nstatM = nstatM+(-layer{ii}.NA)*BM  + (layer{ii}.ND)*BM   ;
        niM = niM+layer{ii}.ni*BM ;
        EAM = EAM+layer{ii}.EA*BM ;
        IPM = IPM+layer{ii}.IP*BM ;
        EiM = EiM+layer{ii}.Ei*BM  ;
        mu_e= mu_e+layer{ii}.mue*BM ;
        mu_h= mu_h+layer{ii}.mup*BM  ;
        N0CM = N0CM+layer{ii}.N0C*BM ;
        N0VM = N0VM+layer{ii}.N0V*BM ;
        Krad=Krad+layer{ii}.krad*BM;
        Ni = Ni+layer{ii}.ni*BM ;
        Nt = Nt+layer{ii}.nt*BM ;
        Pt = Pt+layer{ii}.pt*BM ;
        Taun=Taun+layer{ii}.taun*BM ;
        Taup=Taup+layer{ii}.taup*BM ;
        Cnp=Cnp+layer{ii}.Cn*BM;
        enp=enp+layer{ii}.en*BM;
    end
end
%%%%% Analysis %%%%%
% split the solution into its component parts (e.g. electrons, holes and efield)
u1 = sol.sol(:,:,1); %%% electron density
u2 = sol.sol(:,:,2); %%% hole density
u3 = sol.sol(:,:,3); %%% CT states
u4 = sol.sol(:,:,4); %%% electric potential
u5 = sol.sol(:,:,5); %%% exciton formatiom
% plot the output
rhoc = (u1 + u2  );     % Net charge density calculated from adding individual charge densities
V=u4-EAM;
Ecb = EAM-u4;                                 % Conduction band potential
Evb = IPM-u4;                                 % Valence band potential
Efn = real(Ecb+(kB*T/q)*log((u1)./N0CM));        % Electron quasi-Fermi level
Efp = real(Evb-(kB*T/q)*log(u2./N0VM));        % Hole quasi-Fermi level
Phin = real(EiM+(kB*T/q)*log(u1./niM)-EAM);     % Chemical Potential electrons
Phip = real(EiM-(kB*T/q)*log(u2./niM)-EAM);     % Chemical Potential holes
Phi = Phin - Phip;
% Rec=Krad*((u1.*u2)-(Ni^2))+(((u1.*u2)-Ni^2)/((Taun.*(u2+Pt)) + (Taup.*(u1+Nt))));
%Rec=max(0,Krad.*((u1.*u2)-(Ni2))+(((u1.*u2)-Ni2)./((Taun.*(u2+Pt)) + (Taup.*(u1+Nt)))));
if discretetrap==1
    Recn= Krad.*((u1.*u2)-niM.^2)-u1.*Cnp.*(layer{2}.NTD-u3)+enp.*u3;
    Recp= Krad.*((u1.*u2)-niM.^2)-u2.*Cnp.*u3+enp.*(layer{2}.NTD-u3);
else
    Recn=(Krad.*((u1.*u2)-(niM.^2))+(((u1.*u2)-niM.^2)./((Taun.*(u2+Pt)) + (Taup.*(u1+Nt)))));%./u1;%
    Recp=(Krad.*((u1.*u2)-(niM.^2))+(((u1.*u2)-niM.^2)./((Taun.*(u2+Pt)) + (Taup.*(u1+Nt)))));%./u2;
end
% Recomination Rate - NEEDS SORTING 24/03/2016
% NTS! Electrodes do not have same dielectric constant as pervoskite so
% this is incorrect
Vint = Vbi - (V(:, round(xpoints)) - V(:, 1));           % Internal voltage
% ATM => Total charge density calculated by adding densities is correct.
% Pot from solution u4 and field from u4 should be correct - should be
% verified...
Potc = 0;                          % Potential calculated from Field
Fp = -gradient(u4(end, :), x);                      % Electric field calculated from u4
Potp = u4(end, :);
rhoctot = (trapz(x, rhoc, 2)-trapz(x, rhoc(1,:), 2))/xmax;   % Net charge
% Recomination Rate - NEEDS SORTING 24/03/2016
Urec = 0;% (krad((u1*u2)-ni^2) + sn*(u1-htln0)).*nBM + (krad((u1*u2)-ni^2)+ sp*(u2-etlp0)).*pBM;
if symm==1
    Voc = Efn(:, round(xpoints/2)) - Efp(:, 1) ;                   % Open Circuit Voltage
    if pulseon == 1                                         % AC coupled mode
        
        Voc = Voc - Voc(1, :);
        Voc = Voc*1000;
        t = (t-(tstart+pulselen))*1e6;
    end
else
    Voc = Efp(:, 1 ) - Efn(:, xpoints) ;                  % Open Circuit Voltage
end
Voc0=Voc(1, :);

ntot = trapz(x, u1, 2);     % Total
ptot = trapz(x, u2, 2);
nCTtot = trapz(x, u3, 2);     % Total
Extot = trapz(x, u5, 2);
%     rhoctot=  trapz(x, rhoc, 2)/xmax;
% Calculates current at every point and all times
if calcJ == 1
    
    % find the internal current density in the device
    Jndiff = zeros(length(t), length(x));
    Jndrift = zeros(length(t), length(x));
    Jpdiff = zeros(length(t), length(x));
    Jpdrift= zeros(length(t), length(x));
    Jpart = zeros(length(t), length(x));
    Jtot = zeros(length(t));
    
    for jj=1:length(t)
        
        tj = t(jj);
        
        [nloc,dnlocdx] = pdeval(0,x,u1(jj,:),x);
        [ploc,dplocdx] = pdeval(0,x,u2(jj,:),x);
        [Vloc, Floc] = pdeval(0,x,u4(jj,:),x);
        
        % Particle currents
        Jndiff(jj,:) = (mu_e(1,:).*dnlocdx*kB*T)*(1000*e);
        Jndrift(jj,:) = (-mu_e(1,:).*nloc.*Floc)*(1000*e);
        
        Jpdiff(jj,:) = (-mu_h(1,:).*dplocdx*kB*T)*(1000*e);
        Jpdrift(jj,:) = (-mu_h(1,:).*ploc.*Floc)*(1000*e);
        
        
        % Particle current
        Jpart(jj,:) = Jndiff(jj,:) + Jndrift(jj,:) + Jpdiff(jj,:) + Jpdrift(jj,:);
        
        % Electric Field
        Floct(jj,:) = Floc;
        
    end
    
    
    % Currents at the boundaries (should be the same)
    Jpartl = -Jpart(:,1);
    if symm==1
        Jpartr = -Jpart(:,round(xpoints/2));
        
        % Displacement Current at right hand side
        Fend = (Floct(:, round(xpoints/2)));
        Jdispr =0;%; -(e*1000)*eppn*gradient(Floct(:, round(xpoints/2)), t);
        
        Jtotr = Jpartr + Jdispr;
    else
        Jpartr = -Jpart(:,end);
        
        % Displacement Current at right hand side
        Fend = (Floct(:, end));
        Jdispr = 0;%-(e*1000)*eppn*gradient(Floct(:, end), t);
        
        Jtotr = Jpartr + Jdispr;
    end
    
end

% Calculates currents only for right hand x points at all times
if calcJ == 2
    
    % find the internal current density in the device
    Jndiff = zeros(length(t), 1);
    Jndrift = zeros(length(t), 1);
    Jpdiff = zeros(length(t), 1);
    Jpdrift= zeros(length(t), 1);
    Jpart = zeros(length(t), 1);
    
    for jj=1:length(t)
        
        [nloc,dnlocdx] = pdeval(0,x,u1(jj,:),x(end));
        [ploc,dplocdx] = pdeval(0,x,u2(jj,:),x(end));
        [Vloc, Floc] = pdeval(0,x,u4(jj,:),x(end));
        
        % Particle currents
        Jndiff(jj) = (mu_e(1,end )*kB*T*dnlocdx)*(1000*e);
        Jndrift(jj) = (-mu_e(1,end)*nloc.*Floc)*(1000*e);
        
        Jpdiff(jj) = (-mu_h(1,end)*kB*T*dplocdx)*(1000*e);
        Jpdrift(jj) = (-mu_h(1,end)*ploc.*Floc)*(1000*e);
        
        
        
        % Particle current
        Jpart(jj) = Jndiff(jj) + Jndrift(jj) + Jpdiff(jj) + Jpdrift(jj) ;
        
        % Electric Field
        Floct(jj) = Floc;
        
    end
    
    % Currents at the boundaries
    Jpartr = Jpart';
    
    
    Jdispr = 0;%(e*1000)*eppn*gradient(Floct, t);
    
    Jtotr = Jpartr + Jdispr;
    
end
if sol.params.BC ==5;
    x0=[0.1,5];
    fun = @(x,xdata) x(1)*sin(2*pi*sol.params.wAC*xdata+x(2));
    xVoc = lsqcurvefit(fun,x0,t',(Voc-Voc(1)));
    
    x0=[1,5];
    errorfitI=100;
    while(errorfitI>10)
        xI = lsqcurvefit(fun,x0,t',(Jtotr-(max(Jtotr)+min(Jtotr))/2));
        errorfitI=sum(((fun(xI,t')-(Jtotr-(max(Jtotr)+min(Jtotr))/2))./max(Jtotr)).^2);
        x0(1)=x0(1)*100;
    end
    Zi=xVoc(1)/xI(1)*(exp(j*(xVoc(2)-xI(2))));
    Capacitance=(imag(1./Zi)/2/pi/sol.params.wAC);
    %
    %             figure(1)
    %             plot(t',(Jtotr-(max(Jtotr)+min(Jtotr))/2),t',fun(xI,t'),'*')
    
end
%--------------------------------------------------------------------------------------

% Readout solutions to structure
solstruct=sol;
solstruct.n = u1(end, :)'; solstruct.p = u2(end, :)';  solstruct.V = V(end, :)'; solstruct.nct = u3(end, :)';
solstruct.nexc = u5(end, :)';		solstruct.Urec = Urec;
solstruct.Ecb = Ecb(:, :)'; solstruct.Evb = Evb(:, :)'; solstruct.Efn = Efn(:, :)'; solstruct.Efp = Efp(:, :)'; ...
    solstruct.rho = rhoc; solstruct.Field = Fp; solstruct.Voc = Voc; solstruct.xnm = xnm';solstruct.nt = u3(end, :)';
solstruct.Recn = Recn;solstruct.rhoctot=rhoctot;solstruct.nCTtot=nCTtot;solstruct.Extot=Extot;
solstruct.t=t;solstruct.Recp = Recp;solstruct.Voc0=Voc0;
if calcJ ~= 0
    solstruct.Jtotr = Jtotr; solstruct.Jpartr = Jpartr'; solstruct.Jdispr = Jdispr';
    solstruct.Jndiff=Jndiff;solstruct.Jndrift=Jndrift;solstruct.Jpdiff=Jpdiff;solstruct.Jpdrift=Jpdrift;
    solstruct.Jpart = Jpart;
end
if sol.params.BC ==5;
    solstruct.Zi=Zi;
    solstruct.Capacitance=Capacitance;
    solstruct.errorfitI=errorfitI;
    
end
% Store params
end