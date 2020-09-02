function solution=changep(value,parameter)
p=hstruct(pnParamsHCT).data;
if (strcmp(parameter,'T'))
    ;
else sol_eq=EquilibratePNHCT(0,p);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%Soluciton for as many parameters we include to change%%%%%
for j=1:1:length(value)
    %%% Equilibrium solution
    p=hstruct(pnParamsHCT).data;
    if (strcmp(parameter,'T'))
        p.T=value(j);
        p=hstruct(p).data;
        sol_eq=EquilibratePNHCT(0,p);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Now we need to re run the solution with the parameters we want and
    %%% taking the initial guess as sol_eq
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
    %%%%%%%%%%%%Re-run hstruct to save the changes into p%%%%%%%%%%%%%%%
    if (strcmp(parameter,'mobilitye'))
        p.layer{2}.muee=value(j);
        p.layer{2}.mupp=value(j);
    elseif (strcmp(parameter,'mobilityp'))
        p.layer{2}.mupp=value(j);
    elseif (strcmp(parameter,'kdisexc'))
        p.layer{2}.kdisexc=value(j);
    elseif (strcmp(parameter,'kdis'))
        p.layer{2}.kdis=value(j);
    elseif (strcmp(parameter,'kfor'))
        p.layer{2}.kfor=value(j);
    elseif (strcmp(parameter,'krec'))
        p.layer{2}.krec=value(j);
    elseif (strcmp(parameter,'krecexc'))
        p.layer{2}.krecexc=value(j);
%     elseif (strcmp(parameter,'T'))
%         p.T=value(j);
    end
    p=hstruct(p).data;
    %%%Build a structure in which the result are going to be stored%%%%%
    sol_eq=pndriftHCT(sol_eq,p);
    %%%%%%%With this solution we can calculate the solution for Jsc%%%%%%%%
    p.Int=4.95;
    if (strcmp(parameter,'Int'))
        p.Int = value(j);
    end
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
    solution{j}.results.J=1;
    while (max(solution{j}.results.J<0)<1) %Initially it will be 0 because J=1 (set)
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
%         if (strcmp(parameter,'T'))
%             p.tmax=p.tmax/10;
%             p.t(1,1)=p.t(1,2);
%             p.AbsTol=1e-5;
%             p.RelTol=1e-2;
%         end
        solution{j}.sol_JV=pndriftHCT(sol_Jsc,p);
        solution{j}.results.V=-solution{j}.sol_JV.Voc;
        solution{j}.results.J=solution{j}.sol_JV.Jtotr;
    end
    %%% Once we get the results, we can store the values of the voltage and the
    %%% current, and the PCE and FF
    solution{j}.results.VV=solution{j}.results.V(1:find(solution{j}.results.J<0,1));
    solution{j}.results.JJ=solution{j}.results.J(1:find(solution{j}.results.J<0,1));
    solution{j}.results.Jsc=solution{j}.sol_JV.Jtotr(1,1);
    solution{j}.results.Voc=interp1(solution{j}.results.JJ,solution{j}.results.VV,0,'spline');
    solution{j}.results.PCE=max(solution{j}.results.VV.*solution{j}.results.JJ);
    solution{j}.results.FF=max(solution{j}.results.VV.*solution{j}.results.JJ)./solution{j}.results.Voc/solution{j}.results.Jsc;
%     %% Finally the solution for the Voc can be calculated
%     p.Vapp =solution{j}.results.Voc;
%     %%% We use the Voc calculated befores as an input for this step
%     p.Int=3;
%     p.figson =0;
%     p.Vtransient=0;
%     p.wAC=0;
%     p.tmax=1e-5;
%     p.tmesh_type = 2;
%     p.tpoints = 1000;
%     p=hstruct(p).data;
%     solution{j}.sol_Voc=pndriftHCT(sol_Jsc,p); %solution{j}.sol_Voc=pndriftHCT(solution{j}.sol_Jsc,p);
end
end