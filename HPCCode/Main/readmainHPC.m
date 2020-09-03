% Read txt into cell A
fid = fopen('mainHPC.m','r');
i = 1;
tline = fgetl(fid);
A{i} = tline;
while ischar(tline)
    i = i+1;
    tline = fgetl(fid);
    A{i} = tline;
end
fclose(fid);
name=0;
for ints=1:1:8
    for tt=1:1:9
        for mm=1:1:5
            for kd=1:1:5
                for kde=1:1:5
                    for kf=1:1:5
                        for kr=1:1:5
                            for kre=1:1:5
                                name=name+1;
                                B=A;
                                % Change cell A
                                B{11} =  strrep(A{11},'=1',['=',num2str(ints)]);
                                B{12} =  strrep(A{12},'=1',['=',num2str(tt)]);
                                B{13} =  strrep(A{13},'=1',['=',num2str(mm)]);
                                B{14} =  strrep(A{14},'=1',['=',num2str(kd)]);
                                B{15} =  strrep(A{15},'=1',['=',num2str(kde)]);
                                B{16} =  strrep(A{16},'=1',['=',num2str(kf)]);
                                B{17} =  strrep(A{17},'=1',['=',num2str(kr)]);
                                B{18} =  strrep(A{18},'=1',['=',num2str(kre)]);
                                B{19} =  strrep(A{19},'=1',['=',num2str(name)]);
                                % Write cell A into txt
                                fid = fopen(['main' num2str(name) '.m'], 'w');
                                for i = 1:numel(B)
                                    if B{i+1} == -1
                                        fprintf(fid,'%s', B{i});
                                        break
                                    else
                                        fprintf(fid,'%s\n', B{i});
                                    end
                                end
                                fclose(fid);
                            end
                        end
                    end
                end
            end
        end
    end
end


