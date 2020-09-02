function [layer,layers]=readlayers(filename)

%this function reads the excel table and stores the information in a layer
%structure.

delimiterIn = ',';
headerlinesIn = 1;
data = importdata(filename,delimiterIn,headerlinesIn);
layers=length(data.data(:,1));
epp0 = 552434;        % e^2 eV^-1 cm^-1 -Checked (02-11-15)

for ii=1:1:layers
    layer{ii}.epp=data.data(ii,1)*epp0; % Dielectric constant
    layer{ii}.EA = data.data(ii,2);        % Conduction band energy
    layer{ii}.IP = data.data(ii,3);     % Valence band energy
    layer{ii}.PhiCV = data.data(ii,4);   % n doping
    layer{ii}.PhiAV = data.data(ii,5);     % p doping
    layer{ii}.N0C=data.data(ii,6);      % effective Density Of States at the conduction band
    layer{ii}.N0V=data.data(ii,7);      % effective Density Of States at the Valence band
    layer{ii}.muee=data.data(ii,8);        %mobility
    layer{ii}.mupp=data.data(ii,9);        %mobility
    layer{ii}.krad  = data.data(ii,10); % [cm3 s-1] Bulk Radiative Recombination coefficient [nominally 1e-10]
    layer{ii}.taun = data.data(ii,11);   % [s] SRH time constant for electrons
    layer{ii}.taup = data.data(ii,12);  
    layer{ii}.Ete = data.data(ii,13);  % ((EA-IP)/2+IP)-0.2;      % Deep trap energy- currently a guess!
    layer{ii}.Eth = data.data(ii,14);
    layer{ii}.NTA = data.data(ii,15);
    layer{ii}.NTD = data.data(ii,16);
    layer{ii}.tp =data.data(ii,17)*1e-7;   %  layer thickness
    layer{ii}.pp =data.data(ii,18);       % number layer points
    layer{ii}.tinterL = data.data(ii,19)*1e-7;      % Interfacial region thickness
    layer{ii}.epointsL = data.data(ii,20);        % No. of points for electrode [for log mesh]
    layer{ii}.XiL =data.data(ii,21)*1e-7;            % Interfacial region thickness Heterojunction
    layer{ii}.XipL =data.data(ii,22);              % No. of points for electrode [for l     
    layer{ii}.tinterR = data.data(ii,23)*1e-7;      % Interfacial region thickness
    layer{ii}.epointsR = data.data(ii,24);        % No. of points for electrode [for log mesh]
    layer{ii}.XiR =data.data(ii,25)*1e-7;            % Interfacial region thickness Heterojunction
    layer{ii}.XipR =data.data(ii,26);    
    layer{ii}.wr =data.data(ii,27)*1e-7;%ND/(ND+NA)*sqrt(2*eppp/q*Vbi*(1/ND));%50e-7;  %((-ti*NA*q) + ((NA^0.5)*(q^0.5)*(((ti^2)*NA*q) + (4*eppi*Vbi))^0.5))/(2*NA*q);
    layer{ii}.wl =data.data(ii,28)*1e-7;%ND/(ND+NA)*sqrt(2*eppp/q*Vbi*(1/ND));%50e-7;  %((-ti*NA*q) + ((NA^0.5)*(q^0.5)*(((ti^2)*NA*q) + (4*eppi*Vbi))^0.5))/(2*NA*q);
    layer{ii}.int=data.data(ii,29);
    layer{ii}.kdisexc=data.data(ii,30);
    layer{ii}.kdis=data.data(ii,31);
    layer{ii}.kfor=data.data(ii,32);
    layer{ii}.krec=data.data(ii,33);
    layer{ii}.krecexc=data.data(ii,34);
end

end
    
