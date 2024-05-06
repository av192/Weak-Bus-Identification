% Calling the excel sheets and getting the line and bus datas 
clc;
% Information about the bus matrix
  % nd   V     Ang.    Pg        Qg       PL        QL      Gs       jBs   Type
  % (1) (2)    (3)     (4)      ( 5)      (6)       (7)     (8)      (9)   (10)
  % Colum 11: if the bus has shunt element =1, if it hasnt shunt element =0
  % bt; 1->slack, 2->PV 3->PQ bus
bus=[1  1.060  0.000   0.2324    -0.169    0.000    0.000    0.000    0.000  1   0.0;
     2  1.045  0.000   0.4000    0.4240    0.217    0.127    0.000    0.000  2   0.0;
     3  1.01  0.000   0.0000    0.2340    0.942    0.190    0.000    0.000  2   0.0;
     4  1.00  0.000   0.0000    0.0000    0.478    0.039    0.000    0.000  3   0.0;
     5  1.00  0.000   0.0000    0.0000    0.076    0.016    0.000    0.000  3   0.0;
     6  1.07  0.000   0.0000    0.1220    0.112    0.075    0.000    0.000  2   0.0;
     7  1.00  0.000   0.0000    0.0000    0.000    0.000    0.000    0.000  3   0.0;
     8  1.09  0.000   0.0000    0.1740    0.000    0.000    0.000    0.000  2   0.0;
     9  1.00  0.000   0.0000    0.0000    0.295    0.166    0.000    0.190  3   1.0;
     10 1.00  0.000   0.0000    0.0000    0.090    0.058    0.000    0.000  3   0.0;
     11 1.00  0.000   0.0000    0.0000    0.035    0.018    0.000    0.000  3   0.0;
     12 1.00  0.000   0.0000    0.0000    0.061    0.016    0.000    0.000  3   0.0;
     13 1.00  0.000   0.0000    0.0000    0.135    0.058    0.000    0.000  3   0.0;
     14 1.00  0.000   0.0000    0.0000    0.149    0.050    0.000    0.000  3   0.0];
 
 %Information about the line matrix
%COL 1.-  From bus
%COL 2.-  to bus
%COL 3.-  R P.U.
%COL 4.-  Xl P.U.
%COL 5.-  Xc (parallel) P.U.
%COL 6.-  Type of line: 0==Line ; tap value if its a transformer
%COL 7.- phase shifter angle
line=[1  2   0.01938   0.05917   0.02640    0.000   0.00;
      2  3   0.04699   0.19797   0.02190    0.000   0.00;
      2  4   0.05811   0.17632   0.01870    0.000   0.00;
      1  5   0.05403   0.22304   0.02460    0.000   0.00;
      2  5   0.05695   0.17388   0.01700    0.000   0.00;
      3  4   0.06701   0.17103   0.01730    0.000   0.00;
      4  5   0.01335   0.04211   0.00640    0.000   0.00; 
      5  6   0.00000   0.25202   0.00000    0.932   0.00;
      4  7   0.00000   0.20912   0.00000    0.978   0.00;
      7  8   0.00000   0.17615   0.00000    0.000   0.00;
      4  9   0.00000   0.55618   0.00000    0.969   0.00;
      7  9   0.00000   0.11001   0.00000    0.000   0.00;
      9  10  0.03181   0.08450   0.00000    0.000   0.00;
      6  11  0.09498   0.19890   0.00000    0.000   0.00;
      6  12  0.12291   0.25581   0.00000    0.000   0.00;
      6  13  0.06615   0.13027   0.00000    0.000   0.00;
      9  14  0.12711   0.27038   0.00000    0.000   0.00;
      10 11  0.82050   0.19207   0.00000    0.000   0.00;
      12 13  0.22092   0.19988   0.00000    0.000   0.00;
      13 14  0.17093   0.34802   0.00000    0.000   0.00];


j=sqrt(-1);  i = sqrt(-1);

fb = line(:,1);      % From bus number
tb = line(:,2);      % To bus number
R = line(:,3);       % Resistance, R
X = line(:,4);       % Reactance, X...
b = j*line(:,5);     % Shunt Admittance, B/2 

Z = R + j*X;             % Z matrix
y= 1./Z;   %branch admittance

nbranch=length(line(:,1)); % no. of branches
nbus = max(max(fb), max(tb));  % no. of buses

%  Forming the Y Bus Matrix

for n = 1:nbranch
Ybus=zeros(nbus,nbus);     % initialize Ybus to zero
               
% Formation of the off diagonal elements
for k=1:nbranch;
       Ybus(fb(k),tb(k))=Ybus(fb(k),tb(k))-y(k);
       Ybus(tb(k),fb(k))=Ybus(fb(k),tb(k));
    end
end

% Formation of the diagonal elements
for  n=1:nbus
     for k=1:nbranch
         if fb(k)==n | tb(k)==n
         Ybus(n,n) = Ybus(n,n)+y(k) + b(k);
            else, end
     end
end
fprintf('\n\n\t\t\t --------[POWER FLOW SOLUTIONS]---------');
fprintf('\n\n\t\t\t             <--METHODS-->');
fprintf('\n\n\t *1) GAUSS SIDEL METHOD  \t   *2) NEWTON RAPHSON METHOD \t');
fprintf('\n\n\t          [Tolerance: 0.001]     [Acceleration Factor:1.6]')


casi=input('\n\n   ENTER THE METHOD YOU PREFER :');

switch(casi)
    
    case 1 
% Gauss-Seidel method

nbuses=length(bus(:,1)); % number of buses of the electric power system
V=bus(:,2); Vprev=V; % Initial bus voltages
The=bus(:,3); % Initial bus angles
% Net power (Generation - Load)
P=bus(:,4)-bus(:,6);     
Q=bus(:,5)-bus(:,7); 
% ++++++++++++++ First, compute the addmitance matrix Ybus ++++++++++++++++
[Y] = Y_admi(line,bus,nbuses); % function to get the admittance matrix 

% +++++++++++++++++++++++ Start iterative process +++++++++++++++++++++++++
tolerance=1; 
iteration=0;
st=clock; % start the iteration time clock
while (tolerance > 1e-8)
    for k=2:nbuses
        PYV=0;
        for i=1:nbuses
            if k ~= i
                PYV = PYV + Y(k,i)* V(i);  % Vk * Yik
            end
        end
        if bus(k,10)==2 % PV bus
            % Estimate Qi at each iteration for the PV buses
            Q(k)=-imag(conj(V(k))*(PYV + Y(k,k)*V(k)));
        end
        V(k) = (1/Y(k,k))*((P(k)-j*Q(k))/conj(V(k))-PYV); % Compute bus voltages
        if bus(k,10) == 2 % For PV buses, the voltage magnitude remains same, but the angle changes
            V(k)=abs(Vprev(k))*(cos(angle(V(k)))+j*sin(angle(V(k))));
        end
    end
    iteration=iteration+1; % Increment iteration count
    tolerance = max(abs(abs(V) - abs(Vprev))); % Tolerance at the current iteration
    Vprev=V;
end
ste=clock; % end the iteration time clock

% Now, we have found V and Theta, then, we can compute the power flows

% +++++++++++++++++++++++++++++ Power flow ++++++++++++++++++++++++++++++++
% currents at each node
I=Y*V;
% Power at each node
S=V.*conj(I); % Complex power
for k=1:nbuses
    if bus(k,10)==1
        % Real and reactive generation at the Slack bus
        Pgen(k)=real(S(k));
        Qgen(k)=imag(S(k));
    end
    if bus(k,10)==2
        % Real and reactive generation at the PV buses
        Pgen(k)=real(S(k))+bus(k,6);
        Qgen(k)=imag(S(k))+bus(k,7);
    end
    if bus(k,10)==3
        Pgen(k)=0;
        Qgen(k)=0;
    end
end
% calculate the line flows and power losses
FromNode=line(:,1);
ToNode=line(:,2);
nbranch = length(line(:,1)); % number of branches
for k=1:nbranch
    a=line(k,6);   % Find out if is a line or a transformer, a=0 -> line, a=1 -> transformer, 0<a<1 -> Transformer
    switch a    % for both cases use the pi model
        case 0  %if its a line a=0
            b=1i*line(k,5);
            suceptancia(k,1)=b/2; 
            suceptancia(k,2)=b/2;
        otherwise %if its a transformer
            Zpq=line(k,3)+1i*line(k,4);
            Ypq=Zpq^-1;
            suceptancia(k,1)=(Ypq/a)*((1/a)-1); 
            suceptancia(k,2)=Ypq*(1-(1/a));
    end
end
% Define admmitance of lines
r = line(:,3);
rx = line(:,4);
z = r + j*rx;
y = ones(nbranch,1)./z;
% Define complex power flows
Ss = V(FromNode).*conj((V(FromNode) - V(ToNode)).*y ...
   + V(FromNode).*suceptancia(:,1)); % complex flow of the sending buses
Sr = V(ToNode).*conj((V(ToNode) - V(FromNode)).*y ...
   + V(ToNode).*suceptancia(:,2)); % complex low of the receiving buses

% Define active and reactive power flows
Pij=real(Ss);
Qij=imag(Ss);
Pji=real(Sr);
Qji=imag(Sr);

% Active power lossess
P_loss=sum(Pij+Pji);

% Reactive power lossess
Q_loss=sum(Qij+Qji);

fprintf('\n\n\t\t\t\t\t GAUSS SIDEL SOLUTION');

fprintf('\n\n\t 1) Y BUS');
fprintf('\n\t 2) LINE FLOW SOLUTION');
fprintf('\n\t 3) LINE LOSSES SOLUTION');
fprintf('\n\t 4) EXIT');

opt=input('\n\n Choose your option : ');

if(opt==1)
    
%  DISPLAYING Y BUS
fprintf('                               Y BUS  \n\n')

display(Ybus);

end
if (opt==2)
    
%  DISPLAYING POWER FLOW SOLUTIONS UPTO A VALUE OF 3 DECIMAL PLACES

% disp(tech)
fprintf('                               %g Iterations  \n\n', iteration)
head =['    Bus   Voltage   Angle     ------Load------     ---Generation--- '
       '    No.   Mag.      Degree      MW       Mvar        MW       Mvar  '
       '                                                                    '];
disp(head)
ywz=[  bus(:,1)    abs(V)  (180/pi)*angle(V)  bus(:,6)*100 bus(:,7)*100 Pgen'*100  Qgen'*100];
disp(ywz)

end


if(opt==3)

% CALCULATING LINE FLOW LOSSES 

SLT = 0;
fprintf('\n')
fprintf('                           Line Flow and Losses \n\n')
fprintf('     --Line--            Power at bus & line flow     \n')
fprintf('     from      to          MW      Mvar     MVA            \n')

l=1:1:length(line(:,1));
xy=[l' FromNode ToNode Pij Qij];
yx=[l' ToNode  FromNode Pji Qji];
disp(xy)
disp(yx)

if(opt==4)
fprintf('\n\n\t    Have a good day! ');
end
end

    case 2

%  Newton-Raphson method

basemva = 100;        %Base MVA 
tolerance = 0.0001;   %Tolerance 
mi = 80;              %Maximum Iterations

% Keys for check purposes
ns=0; ng=0; Vm=0; delta=0; yload=0; deltad=0;

nbus = length(bus(:,1));
for k=1:nbus
n=bus(k,1);             % Bus Number
bt(n)=bus(k,2);         % Bus type
Vm(n)=bus(k,3);         % Magnitude of bus voltage
delta(n)=bus(k, 4);     % Bus voltage Angle
Pd(n)=bus(k,5);         % Power Demand
Qd(n)=bus(k,6);         % Reactive Power Demand  
Pg(n)=bus(k,7);         % Power generated
Qg(n) = bus(k,8);       % Reactive Power Generated
end
svc_loc=14;     % SVC location
svc_value=0.0;  % SVC Size
svc=[svc_loc,svc_value];
% svc=[];       % For Base case
lfData{3}=100;  % Load Percentage 

fprintf('\n\n\t\t\t\t\t NEWTON RAPHSON SOLUTION');

fprintf('\n\n\t 1) Y BUS');
fprintf('\n\t 2) LINE FLOW SOLUTION');
fprintf('\n\t 3) LINE LOSSES SOLUTION');
fprintf('\n\t 4) EXIT');

opt=input('\n\n Choose your option : ');

if(opt==1)

%  DISPLAYING Y BUS
fprintf('Y BUS  \n\n')
display(Ybus);

elseif(opt==4)
    fprintf('\n\n\tHave a good day! \n');
else
    % for 14-bus system
    % pf=run_pf(nbus,svc,lfData,opt);  

    % for 30 bus system
    pf=run_pf_30(30,svc,lfData,opt);
    
end

end