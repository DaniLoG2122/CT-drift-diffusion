%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Int=[0.079969/1.0227, 0.24997/1.0227, 0.3977/1.0227, 0.63635/1.0227, 1.0227/1.0227, 1.6023/1.0227, 2.9319/1.0227, 4.34/1.0227]*4.95;
Tempmat=137:20:297;
mobep=[linspace(0.1,1.1,5)*1e-4;linspace(0.1,1.1,5)*1e-4;linspace(0.1,1.1,5)*1e-4;linspace(0.9,1.9,5)*1e-4;linspace(1.9,2.9,5)*1e-4;linspace(2.9,3.9,5)*1e-4;linspace(3.9,4.9,5)*1e-4;linspace(4.9,5.9,5)*1e-4];%logspace(-5,-3,5);
kdis=logspace(9,11,5);
kdisexc=logspace(10,12,5);
kfor=logspace(-10,-7,5);
krec=logspace(9,11,5);
krecexc=logspace(10,12,5);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ints=1;
tt=1;
mm=1;
kd=1;
kde=1;
kf=1;
kr=1;
kre=1;
name=1;
% maxNumCompThreads(1);
%%%%%%%%%%%%Equilibrium sol%%%%%%%%%%%
p=hstruct(pnParamsHCT).data;
p.T=Tempmat(tt);
p=hstruct(p).data;
sol_eq=EquilibratePNHCT(0,p);
results.name=name;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p.Int=0;
p.BC=2;
p.Vapp =0;
p.figson =0;
p.Vtransient=0;
p.wAC=0;
p.symm=0;
p.pulseon=0;
p.tmax=1e-2;
p.tmesh_type = 2;
p.layer{2}.muee=mobep(tt,mm);
p.layer{2}.mupp=mobep(tt,mm);
p.layer{2}.kdis=kdis(kd); %exp(((p.T-297)*kdis(kd))/(297*1000*p.T*p.kB))*1e11;
p.layer{2}.kdisexc=kdisexc(kde); %exp(((p.T-297)*kdisexc(kde))/(297*1000*p.T*p.kB))*1e12;     
p.layer{2}.kfor=kfor(kf);
p.layer{2}.krec=krec(kr);
p.layer{2}.krecexc=krecexc(kre);
p=hstruct(p).data;
%%%Build a structure in which the result are going to be stored%%%%%
sol_eq=pndriftHCT(sol_eq,p);
%%%%%%%With this solution we can calculate the solution for Jsc%%%%%%%%
p.Int=Int(ints); %4.95;
%%%%%Note the inetnesity was 0 so now we have generation%%%%%%
    p.BC=2;
    p.Vapp =0;
    p.figson =0;
    p.Vtransient=0;
    p.wAC=0;
    p.tmax=1e-2;
    p=hstruct(p).data;
    p.tmesh_type = 2;
    %%% Logarithmic mess
    p.tpoints = 1000;
    sol_Jsc=pndriftHCT(sol_eq,p);
    p.tmax=5e-2;
    %%% If ptmax is set to values below 5e-2 there are spikes popping up in the
    %%% plots
    %%%%%%Note that the time reduces to calculate a quicker solution%%%%%%
    results.J=1;
    while (results.J<0)<1 %Initially it will be 0 because J=1 (set)
        p.BC=2;
        p.symm=0;
        p.Vapp =1.5; %%% Bias voltage is applied now
        p.figson =0;
        p.Vtransient=0;
        p.wAC=0;
        p.tmax= p.tmax/2; %%% Reduces the time in every iteration
        p.tmesh_type = 2;
        p.tpoints = 1000;
        p=hstruct(p).data;
        sol_JV=pndriftHCT(sol_Jsc,p);
        results.V=-sol_JV.Voc;
        results.J=sol_JV.Jtotr;
    end
%%% Once we get the results, we can store the values of the voltage and the
%%% current, and the PCE and FF
results.p=p;
results.VV=results.V(1:find(results.J<0,1));
results.JJ=results.J(1:find(results.J<0,1));
results.Jsc=sol_JV.Jtotr(1,1);
results.Voc=interp1(results.JJ,results.VV,0,'spline');
results.PCE=max(results.VV.*results.JJ);
results.FF=max(results.VV.*results.JJ)./results.Voc/results.Jsc;
save(['results' num2str(name) '.mat'],'results')