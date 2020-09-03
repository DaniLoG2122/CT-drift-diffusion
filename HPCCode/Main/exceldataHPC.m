n=140625;
clear bad
bad=[];
% bad will store the results that either don't exist (1) or are empty (0)
for i=80526
    name=strcat('results',num2str(i),'.mat');
    A=exist(name);
    s=dir(name);
    the_size=s.bytes;
    if A==2 || the_size~=0
        load(name)
        if numel(fieldnames(results)) == 10
            ;
        else bad=[bad;i,0];
            continue
        end
        intensity(i,1)=results.p.Int*(1.0227/4.95);
        temperature(i,1)=results.p.T;
        mobility(i,1)=results.p.layer{2}.mue;
        kfor(i,1)=results.p.layer{2}.kfor;
        kdis(i,1)=results.p.layer{2}.kdis;
        kdisexc(i,1)=results.p.layer{2}.kdisexc;
        krec(i,1)=results.p.layer{2}.krec;
        krecexc(i,1)=results.p.layer{2}.krecexc;
        Jsc(i,1)=results.Jsc;
        Voc(i,1)=results.Voc;
        FiF(i,1)=results.FF;
        Index(i,1)=i;
    else bad=[bad;i,1];
    end
end
header = {'Intensity','Temperature', 'Mobility' 'Bfor' , 'kdis', 'kdisexc', 'krec', 'krecexc', 'Jsc', 'Voc', 'FF', 'Index'};
data=[intensity,temperature,mobility,kfor,kdis,kdisexc,krec,krecexc,Jsc,Voc,FiF,Index];
header2 = {'Index', 'Existence'};
if size(bad) ~=0
    data2=[bad(:,1),bad(:,2)];
    info2=array2table(data2,'VariableNames',header2);
    writetable(info2,strcat('check',num2str(temperature(1,1)),'.xlsx'))
end
% Check if you have created an Excel file previously or not 
info=array2table(data,'VariableNames',header);
writetable(info,'results.xlsx')