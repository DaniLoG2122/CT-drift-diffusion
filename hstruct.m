classdef hstruct < handle
	properties
		data
	end
	
	methods
		function obj = hstruct(data)
			obj.data = data;
            obj.data.t0 = obj.data.tmax/1000;
            if obj.data.BC==5
                obj.data.t0 = obj.data.tmax/1000;
            end 
            if obj.data.pulseon==1
                obj.data.t0 = obj.data.pulselen/100;
            end 
            for i=1:1:obj.data.layers
                %obj.data.layer{i}.PhiC = obj.data.layer{i}.EA -obj.data.layer{i}.PhiCV;       %n doping
                %obj.data.layer{i}.PhiA = obj.data.layer{i}.IP + obj.data.layer{i}.PhiAV;     % p doping
                obj.data.layer{i}.Eg = obj.data.layer{i}.EA-obj.data.layer{i}.IP;                      % Band Gap
                obj.data.layer{i}.ND = obj.data.layer{i}.N0C*exp((-obj.data.layer{i}.PhiCV)/(obj.data.kB*300));                       
                obj.data.layer{i}.NA = obj.data.layer{i}.N0V*exp((-obj.data.layer{i}.PhiAV)/(obj.data.kB*300));
                  if  obj.data.mobset == 0
                         obj.data.layer{i}.mue = 0;     % electron mobility
                         obj.data.layer{i}.mup = 0;     % hole mobility
                  else 
                         obj.data.layer{i}.mue = obj.data.layer{i}.muee;     % electron mobility
                         obj.data.layer{i}.mup = obj.data.layer{i}.mupp;     % hole mobility
                  end 
                % Doping concentration and band bending
                    obj.data.layer{i}.ni = sqrt(obj.data.layer{i}.N0C*obj.data.layer{i}.N0V)*(exp(-obj.data.layer{i}.Eg/(2*obj.data.kB*obj.data.T)));          % Intrinsic carrier density
                    obj.data.layer{i}.n0 = obj.data.layer{i}.ni;   % Background density electrons in etl/n-type
                    obj.data.layer{i}.p0 = obj.data.layer{i}.ni;     % Background density holes in etl/n-type
                    obj.data.layer{i}.c0 =obj.data.layer{i}.n0*1e-1 ;
                    obj.data.layer{i}.Phi=0;
                    %obj.data.layer{i}.ND = 0;%obj.data.layer{i}.N0C*exp((obj.data.layer{i}.PhiC-obj.data.layer{i}.EA)/(obj.data.kB*obj.data.T));                       
                    %obj.data.layer{i}.NA = 0;%obj.data.layer{i}.N0V*exp((obj.data.layer{i}.IP-obj.data.layer{i}.PhiA)/(obj.data.kB*obj.data.T));      
                if (obj.data.layer{i}.PhiCV>0)
                    obj.data.layer{i}.n0 = obj.data.layer{i}.ND;     % Background density electrons in etl/n-type
                    obj.data.layer{i}.PhiC = obj.data.layer{i}.EA - obj.data.kB*obj.data.T*log(obj.data.layer{i}.N0C/obj.data.layer{i}.n0);
                    obj.data.layer{i}.p0 = obj.data.layer{i}.N0V*exp((obj.data.layer{i}.IP-obj.data.layer{i}.PhiC)/(obj.data.kB*obj.data.T));     % Background density holes in etl/n-type
                    obj.data.layer{i}.c0 =obj.data.layer{i}.p0*1e-1  ;
                    obj.data.layer{i}.Phi=obj.data.layer{i}.PhiC;
                    %obj.data.layer{i}.ND = obj.data.layer{i}.N0C*exp((obj.data.layer{i}.PhiC-obj.data.layer{i}.EA)/(obj.data.kB*obj.data.T));                       
                    obj.data.layer{i}.NA = 0;%obj.data.layer{i}.N0V*exp((obj.data.layer{i}.IP-obj.data.layer{i}.PhiA)/(obj.data.kB*obj.data.T));      
                elseif (obj.data.layer{i}.PhiAV>0)
                    obj.data.layer{i}.p0 = obj.data.layer{i}.NA; %obj.data.layer{i}.N0V*exp((obj.data.layer{i}.IP-obj.data.layer{i}.PhiA)/(obj.data.kB*obj.data.T));     % Background density holes in etl/n-type
                    obj.data.layer{i}.PhiA = obj.data.layer{i}.IP + obj.data.kB*obj.data.T*log(obj.data.layer{i}.N0V/obj.data.layer{i}.p0);
                    obj.data.layer{i}.n0 = obj.data.layer{i}.N0C*exp((obj.data.layer{i}.PhiA-obj.data.layer{i}.EA)/(obj.data.kB*obj.data.T));     % Background density electrons in etl/n-type
                    obj.data.layer{i}.c0 =obj.data.layer{i}.n0*1e-1 ;
                    obj.data.layer{i}.Phi=obj.data.layer{i}.PhiA;
                    obj.data.layer{i}.ND = 0;%obj.data.layer{i}.N0C*exp((obj.data.layer{i}.PhiC-obj.data.layer{i}.EA)/(obj.data.kB*obj.data.T));                       
                    %obj.data.layer{i}.NA = obj.data.layer{i}.N0V*exp((obj.data.layer{i}.IP-obj.data.layer{i}.PhiA)/(obj.data.kB*obj.data.T));      
                else  
                    %obj.data.layer{i}.n0 = obj.data.layer{i}.ni;   % Background density electrons in etl/n-type
                    %obj.data.layer{i}.p0 = obj.data.layer{i}.ni;     % Background density holes in etl/n-type
                    %obj.data.layer{i}.c0 =obj.data.layer{i}.n0*1e-1 ;
                    %obj.data.layer{i}.Phi=0;
                    obj.data.layer{i}.ND = 0;%obj.data.layer{i}.N0C*exp((obj.data.layer{i}.PhiC-obj.data.layer{i}.EA)/(obj.data.kB*obj.data.T));                       
                    obj.data.layer{i}.NA = 0;%obj.data.layer{i}.N0V*exp((obj.data.layer{i}.IP-obj.data.layer{i}.PhiA)/(obj.data.kB*obj.data.T));      
                end 
                % Intrinsic Fermi Energy
                obj.data.layer{i}.Ei = 0.5*((obj.data.layer{i}.EA+obj.data.layer{i}.IP));
                % Trap properties
                obj.data.layer{i}.nt =  obj.data.layer{i}.ni*exp(( obj.data.layer{i}.Ete- obj.data.layer{i}.Ei)/( obj.data.kB* obj.data.T));     % Density of CB electrons when Fermi level at trap state energy
                obj.data.layer{i}.pt =  obj.data.layer{i}.ni*exp(( obj.data.layer{i}.Ei- obj.data.layer{i}.Eth)/( obj.data.kB* obj.data.T));     % Density of VB holes when Fermi level at trap state energy
                obj.data.layer{i}.en=1/obj.data.layer{i}.NTA/obj.data.layer{i}.taun*obj.data.layer{i}.nt;
                obj.data.layer{i}.Cn=1/obj.data.layer{i}.NTA/obj.data.layer{i}.taun;
                obj.data.layer{i}.ep=1/obj.data.layer{i}.NTD/obj.data.layer{i}.taup*obj.data.layer{i}.pt;
                obj.data.layer{i}.Cp=1/obj.data.layer{i}.NTD/obj.data.layer{i}.taup;
                %%%%%%%%%%%%heterojunction section 
                if(i>1)
                obj.data.layer{i}.DEAL=(obj.data.layer{i}.EA-obj.data.layer{i-1}.EA)/2/obj.data.layer{i}.XiL;
                obj.data.layer{i}.DIPL=(obj.data.layer{i}.IP-obj.data.layer{i-1}.IP)/2/obj.data.layer{i}.XiL;
                obj.data.layer{i}.DN0CL=(log(obj.data.layer{i}.N0C/obj.data.layer{i-1}.N0C))/2/obj.data.layer{i}.XiL;
                obj.data.layer{i}.DN0VL=(log(obj.data.layer{i}.N0V/obj.data.layer{i-1}.N0V))/2/obj.data.layer{i}.XiL;

                end 
                if(i<obj.data.layers)
                obj.data.layer{i}.DEAR=(obj.data.layer{i+1}.EA-obj.data.layer{i}.EA)/2/obj.data.layer{i}.XiR;
                obj.data.layer{i}.DIPR=(obj.data.layer{i+1}.IP-obj.data.layer{i}.IP)/2/obj.data.layer{i}.XiR;
                obj.data.layer{i}.DN0CR=(log(obj.data.layer{i+1}.N0C/obj.data.layer{i}.N0C))/2/obj.data.layer{i}.XiR;
                obj.data.layer{i}.DN0VR=(log(obj.data.layer{i+1}.N0V/obj.data.layer{i}.N0V))/2/obj.data.layer{i}.XiR;
                end 
            end 
            
        %%%%%%%%%%%%%%Grid      
        x=0;
        for ii=1:1:obj.data.layers
          if(ii==1)
               x=[linspace(0,obj.data.layer{ii}.XiL,obj.data.layer{ii}.XipL),...
                    linspace((obj.data.layer{ii}.tinterL/max(1,obj.data.layer{ii}.epointsL))+obj.data.layer{ii}.XiL,obj.data.layer{ii}.XiL+obj.data.layer{ii}.tinterL,obj.data.layer{ii}.epointsL),...
                    linspace(obj.data.layer{ii}.XiL+obj.data.layer{ii}.tinterL+(obj.data.layer{ii}.tinterL/max(1,obj.data.layer{ii}.epointsL)),obj.data.layer{1}.tp-obj.data.layer{1}.tinterR,obj.data.layer{1}.pp),...
                    linspace(obj.data.layer{1}.tp-obj.data.layer{1}.tinterR+(obj.data.layer{1}.tinterR/max(1,obj.data.layer{1}.epointsR)),obj.data.layer{1}.tp-obj.data.layer{1}.XiR,obj.data.layer{1}.epointsR),...
                    linspace(obj.data.layer{1}.tp-obj.data.layer{1}.XiR+(obj.data.layer{1}.XiR/max(1,obj.data.layer{1}.XipR)),obj.data.layer{1}.tp,obj.data.layer{1}.XipR)];
                obj.data.layer{1}.XL=0;
                obj.data.layer{1}.XR=max(x);

            else
                obj.data.layer{ii}.XL=max(x)+(obj.data.layer{ii}.XiL/obj.data.layer{ii}.XipL);
                x=[x,linspace(max(x)+(obj.data.layer{ii}.XiL/max(1,obj.data.layer{ii}.XipL)),max(x)+obj.data.layer{ii}.XiL,obj.data.layer{ii}.XipL),...
                    linspace(max(x)+(obj.data.layer{ii}.tinterL/max(1,obj.data.layer{ii}.epointsL))+obj.data.layer{ii}.XiL,max(x)+obj.data.layer{ii}.tinterL,obj.data.layer{ii}.epointsL),...
                    linspace(max(x)+obj.data.layer{ii}.XiL+obj.data.layer{ii}.tinterL+(obj.data.layer{ii}.tinterL/max(1,obj.data.layer{ii}.epointsL)),max(x)+obj.data.layer{ii}.tp-obj.data.layer{ii}.tinterR,obj.data.layer{ii}.pp),....
                    linspace(max(x)+obj.data.layer{ii}.tp-obj.data.layer{ii}.tinterR+(obj.data.layer{ii}.tinterR/max(1,obj.data.layer{ii}.epointsR)),max(x)+obj.data.layer{ii}.tp-obj.data.layer{ii}.XiR,obj.data.layer{ii}.epointsR),...
                    linspace(max(x)+obj.data.layer{ii}.tp-obj.data.layer{ii}.XiR+(obj.data.layer{ii}.XiR/max(1,obj.data.layer{ii}.XipR)),max(x)+obj.data.layer{ii}.tp,obj.data.layer{ii}.XipR)];
                obj.data.layer{ii}.XR=max(x);
            end 
        end 
           
        obj.data.gridx=x;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if(obj.data.layers==1)
                obj.data.Vbi=2;
                 obj.data.layer{1}.n0 = obj.data.layer{1}.N0C*exp(-0.1/(obj.data.kB*obj.data.T));     % Background density electrons in etl/n-type
                    obj.data.layer{1}.p0 = obj.data.layer{1}.N0V*exp(-0.1/(obj.data.kB*obj.data.T));     % Background density holes in etl/n-type
               
            else
            obj.data.Vbi=obj.data.layer{obj.data.layers}.Phi-obj.data.layer{1}.Phi;
            end 
            end

	end
end