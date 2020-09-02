function solstruct = pndriftHCT(varargin)
localtime = cputime;
% Look up pdepe solver
% Requires v2struct toolbox for unpacking parameters structure
% IMPORTANT! Currently uses parameters from pnParamsH - all variables must be
% declared in pnDrift (line 52)


%
% Piers Barnes last modified (09/01/2016)
% Phil Calado last modified (07/07/2016)

% This version allows a previous solution to be used as the input
% conditions. If there is no input argument asssume default flat background
% condtions. If there is one argument, assume it is the previous solution
% to be used as the initial conditions. If there are two input arguments,
% assume that first are the x points from the previous solution, and the
% second is the previous solution.

set(0,'DefaultLineLinewidth',2);
set(0,'DefaultAxesFontSize',24);
set(0,'DefaultFigurePosition', [600, 400, 640, 400]);
% Declare Variables - rememeber to add them here if adding to parameters list
% The scopng rules for nested and anonymous functions require that all variables used within the function be present in the text of the code.
% [BC, Int,T,Vapp,Vbi,figson,htlp0, kB,m, OM,pii, pn, pp, pulseint, pulselen,pulseon,q,ti, tp, tn, t0,tmax,...
%     tmesh_type,tpoints, xmax,tstart, symm,layer,layers,gridx]  = deal(0);
[BL,BC, Bn ,calcJ, e, EA,Eg,Ei,Etetl, Ethtl, IP,IC,ilt, Int,N0,NA,ND,NI,PhiA,...
    PhiC,T,Tn,Vapp,Vbi,cn,cp,Efnside, Efpside, edge,ep,epoints,epp0,eppp,eppi,eppn,...
    et,etln0,etlp0,fastrec,figson,htln0,htlp0, kB,kext, krad,kradetl, kradhtl,m,mobset,...
    mobseti, mue_p,muh_p, mue_i, muh_i, mue_n, muh_n, mui, ni, ntetl, nthtl,OM, pedge,...
    pii, pn, pp, ptetl, pthtl,pulseint, pulselen,pulseon,q,se,sn, sp, side, etlsn, etlsp,...
    htlsn, htlsp, taun_etl, taun_htl, taup_etl, taup_htl, ti, tinter, tp, tn, t0,taun,...
    taup,tmax, tmesh_type,tpoints, v, varlist, varstr, wn, wp, wscr, x0,xmax,xmesh_type,...
    xpoints, tedge, thosrh ,tstart,p0,taup_il,taun_il,EAn,EAp,Egn,Egp,Eip,Ein,nin,nip,Xi,DEApi,DEAin,...
    DIPpi, DIPin, DN0pi,DN0in,Xip,IPn,IPp,N0i,N0n,N0p,PhiCV,PhiAV,symm,mu,kBT0,layer,layers,gridx,equilibrium,discretetrap,Vtransient,wAC,AbsTol,RelTol]  = deal(0);
if isempty(varargin)
    params = pnParamsHCT;                         % Calls Function EParams and stores in sturcture 'params'
    v2struct(params);   
     x = gridx;
%     plot(x)
%     input('toz');
elseif length(varargin) == 1
    % Call input parameters function
    icsol = varargin{1, 1}.sol;
    
    x = varargin{1, 1}.x;
    params = pnParamsHCT;                         % Calls Function EParams and stores in sturcture 'params'   
    v2struct(params);  
    
    
     
elseif length(varargin) == 2 

    if varargin{1, 1}.sol == 0
        params = varargin{2};
        v2struct(params); 
        x = gridx;
%             plot(x)
%     input('toz');
    else
        Voc=varargin{1, 1}.Voc(end);
        
        icsol = varargin{1, 1}.sol;
        x = varargin{1, 1}.x;
        params = varargin{2};
        
        v2struct(params);   
    end

end
                      


% define solution mesh either logarithmically or linearly spaced points



xnm = x*1e7;     % x in nm for plotting 

% Define the x points to give the initial

icx = x;

% Time mesh
if tmesh_type == 1
    t =  [linspace(0,tmax,tpoints)];
    t= abs(t);
    options = odeset('MaxOrder',5,'MaxStep',t0*100,'InitialStep',t0*10,'RelTol',1e-3,'AbsTol',1e-6,'NonNegative',[1 1 1 0 1],'Stats','on') ;

elseif tmesh_type == 2
    if pulseon == 1 && tstart>t0 && pulselen<tmax
        % Tweak solver options - limit maximum time step size during integration.

       options = odeset('MaxOrder',5,'MaxStep',t0*1000,'InitialStep',t0,'RelTol',RelTol,'AbsTol',AbsTol,'NonNegative',[1 1 1 0 1],'Stats','on') ;
 
        t = [logspace(log10(t0),log10(tstart),tpoints/20) - t0,logspace(log10(tstart+t0/100),log10(2*tstart+t0/100),tpoints*9.5/20)-t0,logspace(log10(2*tstart+2*t0/100),log10(tmax),tpoints*9.5/20)-t0];

    else 
            t = logspace(log10(t0),log10(tmax),tpoints) - t0;
            % Tweak solver options - limit maximum time step size during integration.

        options = odeset('MaxOrder',5,'MaxStep',t0*10,'InitialStep',t0,'RelTol',RelTol,'AbsTol',AbsTol,'NonNegative',[1 1 1 0 1],'Stats','on') ;
    end
end


% Generation
genspace = linspace(0,tn+tp,pii);

if OM == 1 && Int ~= 0; %OM = Optical Model %%% Int == ligth intensity
    % Beer-Lambert - Currently requires solution in the workspace
    Gx1S = evalin('base', 'BL1Sun');                    % 1 Sun generation profile
    Gx1S = Gx1S';
    GxLas = evalin('base', 'BL638');
    GxLas = GxLas';
   
elseif OM == 2 && Int ~= 0;
    % Call Transfer Matrix code: [Gx1, Gx2] = TMPC1(layers, thicknesses, activeLayer1, activeLayer2)
    [Gx1S, GxLas] = TMPC1({'SiO2' 'P3HT' 'SiO2'}, [pp pii pn], 2, 2);
    Gx1S = Gx1S';
    GxLas = GxLas';
  
end


xpoints=length(x);
xmax=max(x);

impottime=cputime-localtime 
localtime = cputime;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%555
% Call solver - inputs with '@' are function handles to the subfunctions
% below for the: equation, initial conditions, boundary conditions
sol = pdepe(m,@pdex4pde,@pdex4ic,@pdex4bc,x,t,options); %%%% m is the symmetry of the problem, first imput is the function itself, then initial conditions, boundary conditions

assignin('base', 'sol', sol);
SizesSol=size(sol);
solstruct.params = params;  solstruct.tspan=t ;solstruct.x = x; solstruct.t = t(1:SizesSol(1)); solstruct.sol = sol;  
solstruct = DriftanalyseCT(solstruct); %% change the info into a matrix version from the pdepe function that gives you a tensor
timespent=cputime-localtime
solstruct.timespent=timespent;
assignin('base', 'sol', solstruct);
if figson==1
    driftplotCT(solstruct);
end 
% --------------------------------------------------------------------------
% Set up partial differential equation (pdepe) (see MATLAB pdepe help for details of c, f, and s)
function [c,f,s] = pdex4pde(x,t,u,DuDx)
sim=0;

if symm==1
    if x>=xmax/2
        
        x=xmax -x;
        sim=1;
        if(x<0)
            x=0;
        end
    end 
end 

for kk=1:1:layers
         
        if(x<=layer{kk}.XR && x>=layer{kk}.XL )
            break;
        end 
end 



%if side == 1
     
% Beer Lambert or Transfer Matrix 1 Sun
if Int ~= 0 && OM ==1 || Int ~= 0 && OM == 2
     
      if x > tp && x < (tp+ti) 
          g = Int*interp1(genspace, Gx1S, (x-tp));
      else
          g = 0;
      end
 
%     Add pulse
    if pulseon == 1
        if  t >= 10e-6 && t < pulselen + 10e-6
           if x > tp && x < (tp+ti)
                lasg = pulseint*interp1(genspace, GxLas, (x-tp));
                g = g + lasg;
           end
        end
   end
  
% Uniform Generation
elseif OM == 0

          if Int ~= 0      
               g = Int*2.5e21;
          else
               g = 0;
          end

            % Add pulse
        if pulseon == 1 
            if  t >= tstart && t < pulselen + tstart %&& x<=layer{kk}.XL+4*layer{kk}.tinterR
                g = g+pulseint*2.5e21;
            end
        end
elseif OM == 3

          if Int ~= 0      
               g = Int*2.5e21;
            
          else
               g = 0;
          end

            % Add pulse
        if pulseon == 1 
            if  t >= tstart && t < pulselen + tstart && x<=layer{kk}.XR-layer{kk}.tinterR
                g = g+pulseint*2.5e21*(1-exp((tstart-t)*1e12));
            end
        end


else
        g = 0;
        
end

if layer{kk}.int==0
    g=0;
    
end 



% Prefactors set to 1 for time dependent components - can add other
% functions if you want to include the multiple trappng model
c = [1
     1
     1
     0
     1];

f = [  (layer{kk}.mue*((u(1))*-DuDx(4)+kB*T*DuDx(1)));
       (layer{kk}.mup*((u(2))*DuDx(4)+kB*T*DuDx(2)));
       0;
       DuDx(4);
       0;];     
   

 if(x<layer{kk}.XL+layer{kk}.XiL && kk>1 && x>layer{kk}.XL) %%equations 21 and 22 -> transport across layers
     N0C=layer{kk-1}.N0C+(layer{kk}.N0C-layer{kk-1}.N0C)*(x-layer{kk}.XL+layer{kk}.XiL)/2/layer{kk}.XiL;
     N0V=layer{kk-1}.N0V+(layer{kk}.N0V-layer{kk-1}.N0V)*(x-layer{kk}.XL+layer{kk}.XiL)/2/layer{kk}.XiL;
     layer{kk}.DN0CL=(log(N0C/layer{kk-1}.N0C))/2/layer{kk}.XiL;
     layer{kk}.DN0VL=(log(N0V/layer{kk-1}.N0V))/2/layer{kk}.XiL;
     
     f = [(layer{kk}.mue*((u(1))*(-DuDx(4)+((-1)^sim)*layer{kk}.DEAL-((-1)^sim)*layer{kk}.DN0CL*kB*T)+kB*T*DuDx(1)));
          (layer{kk}.mup*((u(2))*(DuDx(4)-((-1)^sim)*layer{kk}.DIPL-((-1)^sim)*layer{kk}.DN0VL*kB*T)+kB*T*DuDx(2)));     
           0;
           DuDx(4);
           0;];     
     
 end 
 if(x>layer{kk}.XR-layer{kk}.XiR && kk<layers && x<layer{kk}.XR ) %%equations 21 and 22 -> transport accross layers
     N0C=layer{kk}.N0C+(layer{kk+1}.N0C-layer{kk}.N0C)*(x-layer{kk}.XR+layer{kk}.XiR)/2/layer{kk}.XiR;
     N0V=layer{kk}.N0V+(layer{kk+1}.N0V-layer{kk}.N0V)*(x-layer{kk}.XR+layer{kk}.XiR)/2/layer{kk}.XiR;
     layer{kk}.DN0CL=(log(N0C/layer{kk}.N0C))/2/layer{kk}.XiL;
     layer{kk}.DN0VL=(log(N0V/layer{kk}.N0V))/2/layer{kk}.XiL;
     f = [(layer{kk}.mue*((u(1))*(-DuDx(4)+((-1)^sim)*layer{kk}.DEAR-((-1)^sim)*layer{kk}.DN0CR*kB*T)+kB*T*DuDx(1)));
          (layer{kk}.mup*((u(2))*(DuDx(4)-((-1)^sim)*layer{kk}.DIPR-((-1)^sim)*layer{kk}.DN0VR*kB*T)+kB*T*DuDx(2)));     
           0;
           DuDx(4);
           0;];  
     
 end
if discretetrap==1 %% Uniform generation and srh recombination
   s = [g- layer{kk}.krad*((u(1)*u(2))-(layer{kk}.ni^2))-u(1)*layer{kk}.Cn*(layer{kk}.NTD-u(3))+layer{kk}.en*u(3);%(((u(1)*u(2))-layer{kk}.ni^2)/((layer{kk}.d*(u(2)+layer{kk}.pt)) + (layer{kk}.taup*(u(1)+layer{kk}.nt)))))%-max((u(1)-htln0)/thosrh,(u(2)-htlp0)/thosrh)%-max((u(1)-htln0)/taun_htl,(u(2)-htlp0)/taup_htl);%g - kradhtl*((u(1)*u(2))-(ni^2)) - (((u(1)*u(2))-ni^2)/((taun_htl*(u(2)+ptetl)) + (taup_htl*(u(1)+ntetl))));
        g- layer{kk}.krad*((u(1)*u(2))-(layer{kk}.ni^2))-u(2)*layer{kk}.Cp*u(3)+layer{kk}.ep*(layer{kk}.NTD-u(3)); %+(((u(1)*u(2))-layer{kk}.ni^2)/((layer{kk}.taun*(u(2)+layer{kk}.pt)) + (layer{kk}.taup*(u(1)+layer{kk}.nt)))))%-max((u(1)-htln0)/thosrh,(u(2)-htlp0)/thosrh)%-max((u(1)-htln0)/taun_htl,(u(2)-htlp0)/taup_htl) ;%g - kradhtl*((u(1)*u(2))-(ni^2)) - (((u(1)*u(2))-ni^2)/((taun_htl*(u(2)+ptetl)) + (taup_htl*(u(1)+ntetl))));
        u(1)*layer{kk}.Cn*(layer{kk}.NTD-u(3))-layer{kk}.en*u(3)-u(2)*layer{kk}.Cp*u(3)+layer{kk}.ep*(layer{kk}.NTD-u(3));
        (q/layer{kk}.epp)*(-u(1)+u(2)-layer{kk}.NA+layer{kk}.ND-u(3));];
else
%     if (u(1)<1e8 || u(2)<1e8)
% 
%         if (u(1)<1e8)
%             f = [  0;
%                 (layer{kk}.mup*((u(2))*DuDx(4)+kB*T*DuDx(2)));
%                 0;
%                 DuDx(4);
%                 0;];
%             s = [layer{kk}.kdis*u(3)+ layer{kk}.kfor*((layer{kk}.ni^2));
%                 layer{kk}.kdis*u(3)+ layer{kk}.kfor*((layer{kk}.ni^2));
%                 layer{kk}.kdisexc*u(5)-layer{kk}.kfor*((layer{kk}.ni^2))-layer{kk}.kdis*u(3)-layer{kk}.krec*u(3);
%                 (q/layer{kk}.epp)*(+u(2)-layer{kk}.NA+layer{kk}.ND);
%                 g-layer{kk}.kdisexc*u(5);];
%         end
%         if (u(2)<1e8)
%             f = [  (layer{kk}.mue*((u(1))*-DuDx(4)+kB*T*DuDx(1)));
%                 0;
%                 0;
%                 DuDx(4);
%                 0;];
%             s = [layer{kk}.kdis*u(3)+ layer{kk}.kfor*((layer{kk}.ni^2));
%                 layer{kk}.kdis*u(3)+ layer{kk}.kfor*((layer{kk}.ni^2));
%                 layer{kk}.kdisexc*u(5)-layer{kk}.kfor*((layer{kk}.ni^2))-layer{kk}.kdis*u(3)-layer{kk}.krec*u(3);
%                 (q/layer{kk}.epp)*(-u(1)-layer{kk}.NA+layer{kk}.ND);
%                 g-layer{kk}.kdisexc*u(5);];
%         end
%     else 
    if u(1)<0 %% Approximation ni^2 = p*n
         u(1)=(layer{kk}.ni^2)/u(2);
    end
    if u(2)<0 
         u(2)=(layer{kk}.ni^2)/u(1);
    end
    if   u(3)<0 
%         u(3)=layer{kk}.kfor*((u(1)*u(2))-(layer{kk}.ni^2))/(layer{kk}.kdis+layer{kk}.krec);
% pause (1)
    end
      s = [layer{kk}.kdis*(u(3))- layer{kk}.kfor*((u(1)*u(2))-(layer{kk}.ni^2));
        layer{kk}.kdis*(u(3))- layer{kk}.kfor*((u(1)*u(2))-(layer{kk}.ni^2));
         layer{kk}.kdisexc*(u(5))+layer{kk}.kfor*((u(1)*u(2))-(layer{kk}.ni^2))-(layer{kk}.kdis+layer{kk}.krec)*(u(3));%-layer{kk}.kfor*layer{kk}.ni^2/(layer{kk}.kdis+layer{kk}.krec));%%%%layer{kk}.kdisexc*(u(5))+layer{kk}.kfor*((u(1)*u(2))-(layer{kk}.ni^2))-(layer{kk}.kdis+layer{kk}.krec)*(u(3)-layer{kk}.kfor*layer{kk}.ni^2/(layer{kk}.kdis+layer{kk}.krec));
        (q/layer{kk}.epp)*(-u(1)+u(2)-layer{kk}.NA+layer{kk}.ND);
        g-layer{kk}.kdisexc*(u(5))-layer{kk}.krecexc*(u(5));];%abs

%     end
end 
% % x
% if x<8.4184e-07
%     
%     f
% else 
% %     input('tot')
% end 
    if equilibrium==1  %solution for the potential
%         c = [0
%              0
%              0
%              0];
         f = [0;
              0;     
              0;
              DuDx(4);
              0;];  
          s = [0; %%equation 16 
                  0;
                  0;
                  (q/layer{kk}.epp)*(-u(1)+u(2)-layer{kk}.NA+layer{kk}.ND);
                  0;];
    end 
 
% 
% x
% input('toz')
end           %pdex4pde

% --------------------------------------------------------------------------

% Define initial conditions.
function u0 = pdex4ic(x)
    for ii=1:1:layers
         
        if(x<layer{ii}.XR)
            break;
        end 
    end 
    if length(varargin) == 0 | varargin{1, 1}.sol == 0

        if symm==1
            if x>=xmax/2
                x=xmax -x;

                if(x<0)
                    x=0;
                end

            end 
        end 
    %ii=max(find(icx==x))
        toz=layer{1}.ni*layer{1}.ni/layer{1}.p0;
        u0  = [layer{1}.ni*layer{1}.ni/layer{1}.p0/exp(-(x/xmax)*log(layer{1}.p0/toz));
               layer{1}.p0*exp(-(x/xmax)*log(layer{1}.p0/toz));
               0;
               x/xmax*Vbi;
               0;]; 
        if layers~=1
            if ii==1
                XL=0;
                XR=layer{ii}.tp;
                previouSno=0;
                previouSpo=0;
                previouSNA=0;
                previouSND=0;
                previouSu0=0 ;
                Nextno=layer{ii+1}.n0;
                Nextpo=layer{ii+1}.p0;
            elseif (ii==layers)
                XL=layer{ii}.XL;
                XR=layer{ii}.XR;
                previouSno=layer{ii-1}.n0;
                previouSpo=layer{ii-1}.p0;
                previouSNA=layer{ii-1}.NA;
                previouSND=layer{ii-1}.ND;
                previouSu0=((layer{ii-1}.NA-layer{ii-1}.ND)*layer{ii-1}.wl^2)/2/(2*layer{ii-1}.epp);
                Nextno=layer{ii}.n0;
                Nextpo=layer{ii}.p0;
            else
                XL=layer{ii}.XL;
                XR=layer{ii}.XR;
                previouSno=layer{ii-1}.n0;
                previouSpo=layer{ii-1}.p0;
                previouSNA=layer{ii-1}.NA;
                previouSND=layer{ii-1}.ND;
                previouSu0=((layer{ii-1}.NA-layer{ii-1}.ND)*layer{ii-1}.wl^2)/2/(2*layer{ii-1}.epp);
                Nextno=layer{ii+1}.n0;
                Nextpo=layer{ii+1}.p0;
            end 
            if x < XL +layer{ii}.wl
                   u0 = [layer{ii}.n0*(1+(x-XL)^2/(layer{ii}.wl^2))/2+previouSno*(1-((x-XL)^2)/(layer{ii}.wl^2))/2;                          %ni*exp((Efnside - (-q*((((q*NA)/(2*eppp))*(x-tp+wp)^2))))/(kB*T));
                         layer{ii}.p0*(1+(x-XL)^2/layer{ii}.wl^2)/2+previouSpo*(1-(x-XL)^2/layer{ii}.wl^2)/2;
                         0;
                         x/xmax*Vbi;
                         0;];

            elseif  x <= XL +layer{ii}.tp-layer{ii}.wr
               u0 = [layer{ii}.n0;
                    layer{ii}.p0;
                    0;
                    x/xmax*Vbi;
                    0;]; 
            elseif x <= XL +layer{ii}.tp
               u0 = [layer{ii}.n0*(1+(x-XR)^2/layer{ii}.wr^2)/2+Nextno*(1-(x-XR)^2/layer{ii}.wr^2)/2;                          %ni*exp((Efnside - (-q*((((q*NA)/(2*eppp))*(x-tp+wp)^2))))/(kB*T));
                     layer{ii}.p0*(1+(x-XR)^2/layer{ii}.wr^2)/2+Nextpo*(1-((x-XR)^2/layer{ii}.wr^2))/2;
                     0;
                     x/xmax*Vbi;
                     0;];
    %         end 
            end 
                u0 = [u0(1);                          %ni*exp((Efnside - (-q*((((q*NA)/(2*eppp))*(x-tp+wp)^2))))/(kB*T));
                     u0(2);
                      0;%layer{ii}.NTD*u0(1)/(u0(1)+u0(2));
                     x/xmax*Vbi;
                     0;];
        end 
    
    
    %}

    elseif length(varargin) == 1 
        % insert previous solution and interpolate the x points
        Vapp0=varargin{1, 1}.params.Vapp;
        Vappx=(Vapp-Vapp0)*(x/xmax);

        u0 = [abs(interp1(icx,icsol(end,:,1),x));
              abs(interp1(icx,icsol(end,:,2),x));
              abs(interp1(icx,icsol(end,:,3),x));
              interp1(icx,icsol(end,:,4),x);
              abs(interp1(icx,icsol(end,:,5),x));];

    elseif   max(max(max(varargin{1, 1}.sol))) ~= 0
    % insert previous solution and interpolate the x points
        u0 = [  abs(interp1(icx,icsol(end,:,1),x));
                abs(interp1(icx,icsol(end,:,2),x));
               abs( interp1(icx,icsol(end,:,3),x));
                interp1(icx,icsol(end,:,4),x);
               abs( interp1(icx,icsol(end,:,5),x));];
    end
%     figure(1)
%     hold on
% 
%     semilogy(x,u0(1),'*')

end

% --------------------------------------------------------------------------

% --------------------------------------------------------------------------

% Define boundary condtions, refer pdepe help for the precise meaning of p
% and you l and r refer to left and right.
% in this example I am controlling the flux through the boundaries using
% the difference in concentration from equilibrium and the extraction
% coefficient.
function [pl,ql,pr,qr] = pdex4bc(xl,ul,xr,ur,t)
% if pulseon ==1   %for CE measurment
%     if  t >= Biasoff 
%                Vapp = 0;
%     end
% end 
% Zero current 
% t
% tho=1e-9;
% Vapp0=varargin{1, 1}.params.Vapp;
% -Vapp0+(Vapp0-Vapp)*(1-exp(-t/tho));
% t
if BC == 0
    
pl = [0;
      0;
      0;
       -ul(4);
       0;];
      
ql = [1; 
      1; 
      1;
         0;
         1;];
        
pr = [0;
      0;
      0;
       -ur(4) + Vbi - Vapp;
       0;];  
    
qr = [1; 
      1; 
      1;
      0;
      1;];


% Fixed charge at the boundaries- contact in equilibrium with etl and htl
elseif BC == 1
%     -1e-0*(ul(1) - layer{1}.n0)
pl = [0;%-1e-0*(ul(1) - layer{1}.n0);
      (ul(2)- layer{1}.p0);
      0;
         -ul(4);
         0;];
      
ql = [1; 
      0; 
      1;
      0;
      1;];
%         -1e-0*(ur(2) - layer{layers}.p0)
        
pr = [(ur(1)- layer{layers}.n0);
      0;%-1e-0*(ur(2) - layer{layers}.p0);
      0; 
     -ur(4)+Vbi+(Voc-Vapp)*(1-exp(-t/t0/100));
     0;];  
    
qr = [0;
      1; 
      1;
      0;
      1;];
  
% Non- selective contacts - equivalent to infinite surface recombination
% velocity for minority carriers
elseif BC == 2
    
pl = [ul(1) - layer{1}.n0;
      ul(2) - layer{1}.p0;
     0;
      -ul(4);
      0;];
      
ql = [0;  
      0;
     1;
      0;
     1;];
        
pr = [ur(1) - layer{layers}.n0;
      ur(2) - layer{layers}.p0;
      0;
      -ur(4)+Vbi+Voc+(-Voc-Vapp)*(1-exp(-t/t0/100));%(1-exp(-t/t0/100));%
      0;];  %+(Vapp0-Vapp)*(1-exp(-t/tho))

qr = [0;
      0;
      1;
      0;
      1;];
% For CE
elseif BC == 3
    
pl = [ul(1) - layer{1}.n0;
      ul(2) - layer{1}.p0;
      ul(3);
      -ul(4);
       ul(5);];
      
ql = [0;  
      0;
      0;
      0;
      0;];
        
pr = [ur(1) - layer{layers}.n0;
      ur(2) - layer{layers}.p0;
       ur(3);
      -ur(4)+Vbi-(Voc+Vapp)*(1-exp(-t/t0/10))+Voc;
       ur(5);] ; %+(Vapp0-Vapp)*(1-exp(-t/tho))
  
qr = [0;
      0;
      0;
      0;
      0;];
% Vbi+(Voc-Vapp)*(exp(-t/t0))

% Open Circuit. Problem - doesn't give zero field at the rh boundary even
% though it should 
elseif BC == 4
pl = [ul(1) - layer{1}.n0;
      ul(2) - layer{1}.p0;
      0;
      ul(4);
      0;];
      
ql = [0; 
      0;
      1;
      0;
      1;];
        
pr = [ur(1) - layer{1}.n0;
      ur(2) - layer{1}.p0;
     0;
      ur(4);
      0;];  
    
qr = [0; 
      0; 
      1;
      0;
      1;];

%impedance Boundary conditions  
elseif BC == 5
    
pl = [ul(1) - layer{1}.n0;
      ul(2) - layer{1}.p0;
      0;
      -ul(4);
      0];
      
ql = [0;  
      0;
      1;
      0];
        
pr = [ur(1) - layer{layers}.n0;
      ur(2) - layer{layers}.p0;
     0;
      -ur(4)+Vbi-Vapp+params.Vtransient*sin(2*pi*params.wAC*t);];  %+(Vapp0-Vapp)*(1-exp(-t/tho))
    
qr = [0;
      0;
      1;
      0];
elseif BC == 6
    
pl = [ul(1) - layer{1}.n0;
      ul(2) - layer{1}.p0;
      ul(3);
      ul(4);
       ul(5);];
      
ql = [0;  
      0;
      0;
      0;
      0;];
        
pr = [ur(1) - layer{layers}.n0;
      ur(2) - layer{layers}.p0;
      ur(3);
      -ur(4)+Vbi-Vapp;
      ur(5);];  %+(Vapp0-Vapp)*(1-exp(-t/tho))

    
qr = [0;
      0;
     0;
      0;
      0;];
end

end

end


