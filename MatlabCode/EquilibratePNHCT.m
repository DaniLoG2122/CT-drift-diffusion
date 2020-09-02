function y=EquilibratePNHCT(Vapp,p)

% Swicth off mobility and set time to short time step 
%p=p1;

% p.discretetrap=0;
% EAN0=p.EAn;
%p.NI = a;
%%% Sets diferent parameters to calculate the initial solution from sol=0
%%% as initial guess
p.tpoints = 50;
p.equilibrium=1;
p.symm=0;
p.pulseon=0;
p.mobset=1;
p.BC=0;
p.Vapp=Vapp;
p.Int=0;
p.thosrh=5e13;
p.figson = 0;
for ii=1:1:p.layers
    p.layer{ii}.krec  = 1; % [cm3 s-1] Bulk Radiative Recombination coefficient [nominally 1e-10]
    p.layer{ii}.kfor = 1e-15;   % [s] SRH time constant for electrons
    p.layer{ii}.kdis = 1;  
end 
sol.sol=0;
toz=hstruct(p);
p=toz.data;
p.tmax = 1e-2; %%% To calculate really quickly
p.t0 = p.tmax/1e3;
% Run with initial solution

sol = pndriftHCT (sol,p);
%%% Now we have the initial solution, so we can calculate another one with
%%% more precession
p.equilibrium=0;
p.mobset=1;
p.BC=2;
p.figson = 0;
p.tmax = 1e-3; %%% More time than before to calculate an extended solution
p.t0 = p.tmax/1e6;
toz=hstruct(p);
p=toz.data;
sol = pndriftHCT(sol, p);

p.mobset=1;
p.tmax = 1e-6;
p.t0 = p.tmax/1e3;
p.Int=0;


sol_eq = pndriftHCT(sol, p);

for ii=2:1:2%p.layers
    p.layer{ii}.krec  = 1e10; % [cm3 s-1] Bulk Radiative Recombination coefficient [nominally 1e-10]
    p.layer{ii}.kfor = 1e-13;   % [s] SRH time constant for electrons
    p.layer{ii}.kdis = 1;  
end 
p.BC=2;
p.Int=0;
% 
p.calcJ = 1; %%% now J can be calculated with a more stable solution
p.tmax = 1e-3;
p.t0 = p.tmax/1e3;
toz=hstruct(p);
p=toz.data;

sol_eq = pndriftHCT(sol_eq, p);


y=sol_eq;


% assignin('base', 'sol_eq', sol_eq);

end
