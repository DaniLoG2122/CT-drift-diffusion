function changepplot(solution)
figure('Name','JV');
for j=1:1:length(solution)
    subplot(1,length(solution),j)
    plot (-solution{j}.sol_JV.Voc, -solution{j}.sol_JV.Jtotr);
    xlabel('Voltage [V]');
    ylabel('Current [mA/cm-2]');
    ylim([-(solution{j}.results.Jsc + 0.1*solution{j}.results.Jsc),0])
    xlim([0,solution{j}.results.Voc])
    title(strcat('T = ', num2str(solution{j}.sol_JV.params.T,'%10.0f\n'),'K'));% & Int = ', num2str(ans{i}.sol_JV.params.Int,'%10.2f\n')))
    %legend(strcat('krec = ',num2str(solution{j}.sol_JV.params.layer{2}.krec,'%10.0e\n')),'Location','northwest')
end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
figure('Name','Band Diagram')
for j=1:1:length(solution)
    subplot(1,length(solution),j)
    tiii=length(solution{j}.sol_JV.t);
    xpoints=length(solution{j}.sol_JV.x);%%%%%%%%%%%%%%%%%%%%careful
    xnm=1e7*solution{j}.sol_JV.x;
    %set(FigHandle, 'units','normalized','position',[.1 .1 .4 .4]);
    if solution{j}.sol_JV.params.pulseon==1
    bb=find(solution{j}.sol_JV.t>0, 1 );
%         plot (xnm, sol.Efn(:,tiii), '--', xnm, sol.Efp(:,tiii), '--', xnm, sol.Ecb(:,tiii), xnm, sol.Evb(:,tiii));
    plot (xnm, solution{j}.sol_JV.Efn(:,1), '--', xnm, solution{j}.sol_JV.Efp(:,1), '--', xnm, solution{j}.sol_JV.Ecb(:,1), xnm, solution{j}.sol_JV.Evb(: ,1));
    else
    plot (xnm, solution{j}.sol_JV.Efn(:,tiii), '--', xnm, solution{j}.sol_JV.Efp(:,tiii), '--', xnm, solution{j}.sol_JV.Ecb(:,tiii), xnm, solution{j}.sol_JV.Evb(:,tiii));
    end
    legend('E_{fn}', 'E_{fp}', 'CB', 'VB');
    set(legend,'FontSize',16);
    xlabel('Position [nm]');
    ylabel('Energy [eV]');
    if solution{j}.sol_JV.params.symm==1
    xlim([0, xnm(round(xpoints/2))]);
    else
    xlim([0, xnm(end)]);
    end
    %ylim([-2.5, 0.5]);
    set(legend,'FontSize',14);
    set(legend,'EdgeColor',[1 1 1]);
    grid on;
%     x1=[0 300 300 0]; 
%     y1=[-1.35 -1.35 0.43 0.43];
%     x2=[300 400 400 300]; 
%     y2=[-1.35 -1.35 0.43 0.43];
%     x3=[400 700 700 400];
%     y3=[-1.35 -1.35 0.43 0.43];
%     ylim([-1.35, 0.43]);
%     patch(x1,y1,'cyan','FaceAlpha',.1)
%     patch(x2,y2,'yellow','FaceAlpha',.1)
%     patch(x3,y3,'black','FaceAlpha',.1)
%     legend('E_{fn}', 'E_{fp}', 'CB', 'VB','Electron transport layer', 'Active layer', 'Hole transfer layer');
    drawnow;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Name','Detailed current');
for j=1:1:length(solution)
    subplot(1,length(solution),j)
    if solution{j}.sol_JV.params.calcJ == 1
    plot(xnm,solution{j}.sol_JV.Jndiff(end, :),xnm,solution{j}.sol_JV.Jndrift(end, :),xnm,solution{j}.sol_JV.Jpdiff(end, :),xnm,solution{j}.sol_JV.Jpdrift(end, :),xnm,solution{j}.sol_JV.Jpart(end, :));
    legend('Jn diff','Jn drift','Jp diff','Jp drift','Total J');
    xlabel('Position [nm]');
    ylabel('Current Density [mA cm^-2]');
    set(legend,'FontSize',12);
    set(legend,'EdgeColor',[1 1 1]);
    if solution{j}.sol_JV.params.symm==1
    xlim([0, xnm(round(xpoints/2))]);
    else
    xlim([0, xnm(end)]);
    end
    grid on;
    drawnow;
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%