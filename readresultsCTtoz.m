% clear all;
name=0;
name0=0;
for kk=1:1:5
    for ii=1:1:5
        for tt=1:1:50
            
            name0=name0+1;
            if name0>0
                if (exist(['results' num2str(name0) '.mat']))
                    load(['results' num2str(name0) '.mat'])
                    %                  load(['results' num2str(name) '.mat'],'sol_JV');
                    
                    %                  load(['results' num2str(name) '.mat'],'Voc');
                    if (isfield(results,'Voc'))
                        if isempty(results.Voc)
                        else
                            name=name+1;
                            table(name,8)=name;
                            table(name,1)=results.p.layer{1, 2}.krec;
                            table(name,2)=results.p.layer{1, 2}.kfor;
                            table(name,3)=results.p.layer{1, 2}.kdis;
                            table(name,4)=results.p.layer{1, 2}.kdisexc;
                            table(name,5)=results.p.Int;
                            table(name,7)=results.p.T;
                            table(name,6)=results.p.layer{1, 2}.muee;
                            datanum=8;
                            table(name,datanum+1)=results.Jsc;
                            table(name,datanum+2)=results.Voc;
                            try
                                V=results.V;
                                J=results.J;
                                JV{name}.V=V;
                                JV{name}.J=J;
                                vmax=find(V>table(name,datanum+2),1);
                                table(name,datanum+3)=max(V(1:vmax).*J(1:vmax));
                                table(name,datanum+4)=max(V(1:vmax).*J(1:vmax))./table(name,datanum+1)/table(name,datanum+2);
                                CE{name}=results.CE;
                                table(name,datanum+6)=results.nce2;
                                TPV{name}=results.TPV;
                                maxnce=find(TPV{name}.nce==max(TPV{name}.nce),1);
                                f=fit(TPV{name}.t(maxnce:end)',TPV{name}.nce(maxnce:end)./max(TPV{name}.nce),'exp2');
                                if f.a>f.c
                                    table(name,datanum+7)=-1/f.b*1e-6;
                                    table(name,datanum+8)=-1/f.d*1e-6;
                                else
                                    table(name,datanum+7)=-1/f.d*1e-6;
                                    table(name,datanum+8)=-1/f.b*1e-6;
                                end
                                table(name,datanum+9)=f.a;
                                table(name,datanum+10)=f.c;
                                maxTPV=find(TPV{name}.Voc==max(TPV{name}.Voc),1);
%                                 f=fit(TPV{name}.t(maxTPV:end)',TPV{name}.nce(maxTPV:end)./max(TPV{name}.nce),'exp1');
%                                 table(name,datanum+11)=-1/f.b*1e-6;
                                maxTPV=maxTPV+50;
                                f=fit(TPV{name}.t(maxTPV:end)',TPV{name}.Voc(maxTPV:end)./max(TPV{name}.Voc),'exp1');
                                
                                table(name,datanum+11)=-1/f.b*1e-6;
                                table(name,datanum+12)=max(results.TPV.Voc);
                                table(name,datanum+13)=results.TPV.nCTtot(1);
                                table(name,datanum+14)=results.TPV.nCTtot(1)*results.p.layer{1, 2}.krec*1.61917e-19*1000*1e-5;
                                table(name,datanum+15)=1.1+0.026*log(table(name,datanum+1)/(results.p.layer{1, 2}.krec*results.p.layer{1, 2}.kfor/results.p.layer{1, 2}.kdis))+0.026*log(6.25e-15);
                                table(name,datanum+16)=(results.p.layer{1, 2}.krec*results.p.layer{1, 2}.kfor/results.p.layer{1, 2}.kdis);
                                
                            catch ME
                                fprintf(ME.message);
                                continue;  % Jump to next iteration of: for i
                            end
                        end
                        %
                    end
                end
            end
        end
    end
end

