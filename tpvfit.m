function y=tpvfit(sol)
t=sol.t/1e6;
Voc=sol.Voc;
nce=sol.rhoctot-sol.rhoctot(1);
Rstart=find(t<-sol.params.pulselen , 1, 'last' );
Rend=find(t<0, 1, 'last' );
last=length(t);
[maxvoc,Rend]=max(Voc);
Rend=Rend+20;
options = optimset('Display','off');
% last=find(Voc(Rend:last)<(maxvoc/3), 1 )+Rend;
F1 = @(x,xdata)x(1)*(1-exp(-x(2)*xdata));
x0=[30e6,1e-9];
lb = [];
ub = [];
[xr,~,~,~,~] = lsqcurvefit(F1,x0,t(Rstart:Rend)'-t(Rstart),(Voc(Rstart:Rend))*1e6,lb,ub,options);
F2 = @(x,xdata)x(1)*(exp(-x(2)*xdata));
x0=[3e7,1e-9];
[xNce,resnorm,~,~,~] = lsqcurvefit(F2,x0,t(Rend:last)',(nce(Rend:last))/1e8,lb,ub,options);
x0=[3e7,1e-9];
[xf,resnorm,~,~,~] = lsqcurvefit(F2,x0,t(Rend:last)',(Voc(Rend:last))*1e6,lb,ub,options);

last1=last;
    resnorm1=1;

while resnorm1>0.01
    last1=int32(last1-length(t)/300);
    [xf,resnorm,~,~,~] = lsqcurvefit(F2,x0,t(Rend:last1)',(Voc(Rend:last1))*1e6,lb,ub,options);
    resnorm1=sum(abs(F2(xf,t(Rend:last1)')/1e6-abs(Voc(Rend:last1)))./(Voc(Rend:last1)))/(last1-Rend) ;

end 
    [xf2,resnorm,~,~,~] = lsqcurvefit(F2,x0,t(last1:last)',(Voc(last1:last))*1e6,lb,ub,options);

while resnorm1>0.001
    
    [xf2,resnorm,~,~,~] = lsqcurvefit(F2,x0,t(last1:last)',(Voc(last1:last))*1e6,lb,ub,options);
    last1=int32(last1+length(t)/100);
        resnorm1=sum(abs(F2(xf2,t(last1:last)')/1e6-abs(Voc(last1:last)))./(Voc(last1:last)))/(last-last1) ;
    
end 
% 
x0=[3e7,1e-9];
    [xNce2,resnorm,~,~,~] = lsqcurvefit(F2,x0,t(Rend:last)',(nce(Rend:last))/1e8,lb,ub,options);
    
figure(1)
plot(t(Rstart:Rend)'-t(Rstart),F1(xr,t(Rstart:Rend)'-t(Rstart))/1e6,t(Rstart:Rend)'-t(Rstart),(Voc(Rstart:Rend)),'*')
 xlabel('time [s]');
    ylabel('Voltage [mV]');
    legend({'fit' 'data'})
figure(2)
semilogy(t(Rend:last)',F2(xf,t(Rend:last)')/1e6,t(Rend:last)',(Voc(Rend:last)),'*')
 xlabel('time [s]');
    ylabel('Voltage [mV]');
    legend({'fit' 'data'})
    figure(3)
semilogy(t(Rend:last)',F2(xNce,t(Rend:last)'),t(Rend:last)',(nce(Rend:last))/1e8,'*')
 xlabel('time [s]');
    ylabel('excess charge density [cm-3]');
    legend({'fit' 'data'})
y.risetaoVoc=1/xr(2);
y.falltaoVoc=1/xf(2);
y.falltao2Voc=1/xf2(2);
y.falltaonce=1/xNce(2);
y.falltaonce2=1/xNce2(2);
y.resnorm=sum(abs(F2(xf,t(Rend:last)')/1e6-(Voc(Rend:last)))./(Voc(Rend:last)))/(last-Rend);
end 