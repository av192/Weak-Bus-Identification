function flow_result=run_pf(nbus, svc, lfData, display)
    finalres=load_flow (nbus,svc, lfData);
    result=power_flow(nbus,finalres,svc);

       flow_result.bus=[finalres{1} finalres{2} result{3} result{4} result{1} result{2} result{5} result{6}];
       flow_result.totalPower=[sum(result{3}) sum(result{4}) sum(result{1}) sum(result{2}) sum(result{5}) sum(result{6})];
      line_loss=line_flow(nbus,finalres);
       flow_result.line=[line_loss{1} line_loss{2}];
        flow_result.lineLoss=line_loss{3};
         flow_result.baseMVA=finalres{5};
         data=case14();                          % Load IEEE14 Data
    nl = length(data.line);
    V = finalres{1};
    z = data.line(:,3) + sqrt(-1)*data.line(:,4);
    Q = result{4}/data.baseMVA;
    fb = data.line(:,1);
    tb = data.line(:,2);
    

         %%  Fast Voltage stability index
for ii=1:nl
      m=fb(ii); 
      n=tb(ii);                                   
      Delta(ii)=finalres{2}(m)-finalres{2}(n);
      Theta=angle(z)*(180/pi);
     
      LSIrij(m,n)=(4*imag(z(ii))*Q(n))/((V(m)*sin(Theta(ii)-Delta(ii)))^2); % Line Stability Index
      FVSIrij(m,n)=(4*((abs(z(ii)))^2)*Q(n))/(((V(m))^2)*imag(z(ii))); % Fast voltage Stability Index
end

    if display==2
        head=[
        '======================================================================================='
        '|     Bus Data                                                                        |'
        '======================================================================================='
        '    Bus  Voltage  Angle    ------Load------    ---Generation---      -----Power-----   '
        '    No.  Mag.     Degree     MW       Mvar       MW       Mvar        MW        Mvar   '
        '                                                                                       '];
        disp(head)
        for n=1:nbus
             fprintf(' %5g', n), fprintf(' %7.3f', finalres{1}(n)),
             fprintf(' %8.3f', finalres{2}(n)), fprintf(' %9.3f', result{3}(n)),
             fprintf(' %9.3f', result{4}(n)),  fprintf(' %9.3f', result{1}(n)),
             fprintf(' %9.3f ', result{2}(n)), fprintf(' %9.3f ', result{5}(n)),fprintf(' %9.3f\n', result{6}(n))
        end
        fprintf('   -----------------------------------------------------------------------------------\n')
        fprintf('    Total              ')
        fprintf(' %9.3f', sum(result{3})), fprintf(' %9.3f', sum(result{4})),
        fprintf(' %9.3f', sum(result{1})), fprintf(' %9.3f', sum(result{2})), fprintf('  %9.3f', sum(result{5})),
        fprintf('  %9.3f\n\n', sum(result{6}))
        
        FVSIrij = sparse(FVSIrij);
        LSIrij = sparse(LSIrij);
        fprintf('   -----------------------------------------------------------------------------------\n')
        fprintf('    FVSI            \n ')
       disp(FVSIrij);
       fprintf('   -----------------------------------------------------------------------------------\n')
        fprintf('    LSI            \n ')
       disp(LSIrij);
        % Store Result 
        flow_result.bus=[finalres{1} finalres{2} result{3} result{4} result{1} result{2} result{5} result{6}];
        flow_result.totalPower=[sum(result{3}) sum(result{4}) sum(result{1}) sum(result{2}) sum(result{5}) sum(result{6})];
        line_loss=line_flow(nbus,finalres);
    

    elseif  display==3

        head2=[
        '============================================================================== '
        '|     Line Data                                                               |'
        '============================================================================== '];
        disp(head2);
        
        disp("   ----Line---   ---From>>To---       ---To>>From---       ---Line Loss---"); 
        
        disp("    from  to     MW       Mvar        MW        Mvar        MW        Mvar");
        disp("                                                                          ");
        nl=line_loss{1}(:,1); nr=line_loss{1}(:,2);
        Pij=line_loss{2}(:,1);
        Qij=line_loss{2}(:,2);
        Pji=line_loss{2}(:,3);
        Qji=line_loss{2}(:,4);
        Ploss=line_loss{2}(:,5);
        Qloss=line_loss{2}(:,6);
        total_ploss=line_loss{3}(1);
        total_qloss=line_loss{3}(2);

        for n=1:length(nl)
            fprintf(' %5.0f', nl(n)), fprintf(' %4.0f', nr(n)),
            fprintf(' %9.3f', Pij(n)), fprintf(' %9.3f', Qij(n)), fprintf('  %9.3f', Pji(n)),
            fprintf('  %9.3f', Qji(n)), fprintf('  %9.3f', Ploss(n)), fprintf('  %9.3f\n', Qloss(n)),
        end
        fprintf('   ------------------------------------------------------------------------\n')
        fprintf('    Total                                             '), fprintf(' %10.3f', total_ploss),
        fprintf(' %9.3f\n', total_qloss)
        flow_result.line=[line_loss{1} line_loss{2}];
        flow_result.lineLoss=line_loss{3};
        flow_result.baseMVA=finalres{5};
    
   

     end
end

function finalout=load_flow(nbus, svc, lfData)
    
    Ybus=admittance(nbus, svc);             % Admittance Matrix
    data=case14();                          % Load IEEE14 Data
    busdata=data.bus;                       % Bus Data
    multiple = 1.0;
    load_percentage=lfData{3};               % Load Percentage
    busdata (:,7) =busdata(:,7) + (busdata (:, 7) * ((load_percentage-100)/100));   % Real power update
    
    basemva=100;
    type_bus=busdata (:,2);
    voltage=busdata(:,3);
    volt_angle=zeros (length (voltage),1);
    Pg=busdata (:,5)/basemva;               % Generator real power
    Qg=busdata(:,6)/basemva;                % Generator reactive power
    Pl=busdata(:,7)/basemva;                % Load real power
    Ql=multiple*busdata(:,8)/basemva;                % Load reactive power
    Qmin=busdata(:,9)/basemva;              % Minimum reactive power
    Qmax=busdata(:,10)/basemva;             % Maximum reactive power
    Pf=Pg - Pl;                             % final Real Power
    Qf=Qg - Ql;                             % final Reactive Power
    Psp=Pf;
    Qsp=Qf;
    
    conductance=real (Ybus);
    susceptance=imag (Ybus);
    pq_loc=find (type_bus==3);                      % PQ bus location
    npq=length (pq_loc);                            % Number of PQ bus 
    
    tolerance = 1;
    iter = 1;                                       % iteration starting
    
    while (tolerance > 1e-5)
        Pf = zeros (nbus, 1);
        Qf = zeros (nbus, 1); 
        % Calculate P and Q
        for i = 1:nbus
            for k = 1:nbus
               Pf (i) = Pf(i) + voltage (i) * ...
                   voltage (k) * (conductance (i, k) *cos (volt_angle (i)-volt_angle(k))+...
                   susceptance(i,k)*sin (volt_angle (i)-volt_angle(k)));
               Qf (i) = Qf(i) + voltage (i) * ...
                   voltage (k) * (conductance (i, k) *sin (volt_angle (i)-volt_angle(k))-...
                   susceptance(i,k)*cos (volt_angle (i)-volt_angle(k)));
            end
        end
        
        %Checking Q-limit violations..
        if iter <=7 && iter> 2                      % Only checked up to 7th iterations..
           for n = 2: nbus
                if type_bus (n)== 2
                    QG = Qf (n) +Ql(n);
                    if QG < Qmin (n)
                        voltage(n)= voltage(n)+ 0.01;
                    elseif QG > Qmax (n)
                        voltage(n)= voltage(n)- 0.01;
                    end
                end
           end
        end
        % Calculate change from specified value
        dPa = Psp-Pf;
        dQa = Qsp-Qf;
        k = 1;
        dQ=zeros (npq, 1);
        for i = 1:nbus
            if type_bus (i)== 3
                dQ (k, 1)= dQa (i);
                k = k+1;
            end
        end
        dP = dPa (2:nbus);
        mismatch_vec= [dP; dQ];
        jacob1=zeros (nbus-1, nbus-1);
        for i = 1: (nbus-1)
            m =i+1;
            for k=1: (nbus-1)
                n = k+1;
                if n == m
                    for n = 1:nbus
                        jacob1(i, k)= jacob1(i, k) + voltage(m)*voltage(n)*...
                            (-conductance(m, n) *sin (volt_angle(m) -volt_angle(n))+...
                              susceptance(m, n)*cos (volt_angle(m)- volt_angle(n)));
                    end
                    jacob1(i, k)=jacob1(i, k)-voltage(m) ^2*susceptance(m, m);
                else
                    jacob1(i, k) = voltage(m)*voltage(n)*(conductance(m, n)*...
                        sin(volt_angle(m)-volt_angle(n))-susceptance (m, n)*...
                        cos(volt_angle(m)-volt_angle(n)));
                end
            end
        end
        jacob2=zeros (nbus-1, npq);
        for i = 1: (nbus-1)
            m = i+1;
            for k = 1:npq
                n = pq_loc(k);
                if n == m
                    for n = 1:nbus
                        jacob2 (i, k) = jacob2(i, k) + voltage(n)*...
                        (conductance(m, n)*cos(volt_angle(m) -volt_angle(n))+ ...
                        susceptance(m, n)*sin (volt_angle(m) -volt_angle(n)));
                    end
                    jacob2(i, k) = jacob2(i, k) + voltage(m)*conductance(m, m);
                else
                    jacob2(i, k) = voltage(m)*(conductance(m, n)*... 
                        cos(volt_angle(m) -volt_angle(n)) + susceptance(m,n)*...
                        sin(volt_angle(m)-volt_angle(n)));
                end
            end
        end
        
        jacob3=zeros (npq, nbus-1);
        for i = 1:npq
            m = pq_loc(i);
            for k = 1: (nbus-1)
                n = k+1;
                if n == m
                    for n=1:nbus
                        jacob3(i, k)= jacob3(i, k) + voltage(m)*...
                            voltage(n) * (conductance(m, n)*...
                            cos(volt_angle(m) - volt_angle(n)) + susceptance(m, n)*...
                            sin(volt_angle(m) - volt_angle(n)));
                    end
                    jacob3 (i, k) = jacob3 (i, k) - voltage(m)^2*conductance(m, m);
                else
                     jacob3 (i, k) = voltage(m) * voltage(n)*...
                        (-conductance(m, n) * cos(volt_angle(m) - volt_angle(n))- ...
                        susceptance(m, n) * sin (volt_angle(m) -volt_angle(n)));
                end
            end
        end
        jacob4=zeros(npq,npq);
        for i = 1:npq
            m = pq_loc(i);
            for k = 1:npq
                n= pq_loc (k);
                if n == m
                    for n = 1:nbus
                        jacob4(i, k) = jacob4(i, k) + voltage(n) *...
                            (conductance(m, n) * sin (volt_angle(m) -volt_angle(n))...
                            - susceptance(m, n) * cos(volt_angle(m) -volt_angle(n)));
                    end 
                    jacob4(i, k) = jacob4(i, k)- voltage(m) * susceptance(m,m);
                else
                    jacob4(i, k)= voltage(m) * (conductance(m, n)*...
                        sin(volt_angle (m) -volt_angle(n)) - susceptance(m, n)*...
                        cos(volt_angle (m) -volt_angle(n)));
                end
            end
        end
        final_jacob= [jacob1 jacob2; jacob3 jacob4];
        corr_vec=final_jacob\mismatch_vec;
        volt_ang_chang=corr_vec(1:nbus-1);
        change_volt_mag=corr_vec(nbus:end);
        volt_angle(2:nbus)= volt_ang_chang + volt_angle(2:nbus);
        k = 1;
        
        for i = 2:nbus
            if type_bus(i) == 3
                voltage(i) = change_volt_mag(k) + voltage (i); 
                k = k+1;
            end
        end
            iter=iter + 1;
            tolerance=max(abs(mismatch_vec));
    end
    finalout{1}=real(voltage);          % Real voltage magnitude
    finalout{2}=180/pi * volt_angle;    % Angle in degree
    finalout{3}=voltage;                % Complex Voltage
    finalout{4}=volt_angle;             % Angle in Radian
    finalout{5}=basemva;
end

function lossRes=line_flow(nbus,lfresult)

    % calculate the line flows and power losses
    if nbus==14
        data=case14();              % Load IEEE14 bus system data
    elseif nbus==30
        data=case30();              % Load IEEE30 bus system data
    end
    linedata=data.line;         % Load Line Data
    nl=linedata(:,1);           % From Bus 
    nr=linedata(:,2);           % To bus 
    nbr=length(linedata(:,1));  % Number of branches
    voltage = lfresult{3};      % Voltage in p.u
    angle = lfresult{4};
    basemva=lfresult{5};
    complex_volt = voltage.*cos(angle) + 1i*voltage.*sin(angle);    % Complex voltage in PU
    suceptanc=zeros(nbr,2);
    for k=1:nbr
        a=linedata(k,6);    % Find out if is a line or a transformer, a=0 -> line, a=1 -> transformer, 0<a<1 -> Transformer
        switch a            % for both cases use the pi model
            case 0          % if its a line a=0
                b=1i*linedata(k,5);
                suceptanc(k,1)=b/2; 
                suceptanc(k,2)=b/2;
            otherwise       %if its a transformer
                Zpq=linedata(k,3)+1i*linedata(k,4);
                Ypq=Zpq^-1;
                suceptanc(k,1)=(Ypq/a)*((1/a)-1); 
                suceptanc(k,2)=Ypq*(1-(1/a));
        end
    end
    % Define admmitance of lines
    r = linedata(:,3);
    x = linedata(:,4);
    z = r + 1i*x;
    y = ones(nbr,1)./z;
    % Define complex power flows
    Ss = complex_volt(nl).*conj((complex_volt(nl) - complex_volt(nr)).*y ...
       + complex_volt(nl).*suceptanc(:,1));     % complex flow of the sending buses
    Sr = complex_volt(nr).*conj((complex_volt(nr) - complex_volt(nl)).*y ...
       + complex_volt(nr).*suceptanc(:,2));     % complex low of the receiving buses
    
    % Define active and reactive power flows
    Pij=real(Ss);
    Qij=imag(Ss);
    Pji=real(Sr);
    Qji=imag(Sr);
    
    % Active power lossess
    Ploss=Pij+Pji;
    total_Ploss=sum(Pij+Pji);
    
    % Reactive power lossess
    Qloss=Qij+Qji;
    total_Qloss=sum(Qij+Qji);
    
    % Store Result
    lossRes{1}=[nl nr];
    lossRes{2}=basemva*[Pij Qij Pji Qji Ploss Qloss];
    lossRes{3}=basemva*[total_Ploss total_Qloss];
end

function result = power_flow(nbus, lfresult, svc)
data=case14();
    basemva=data.baseMVA;           % Base MVA
    busdata=data.bus;              % Bus Data

    % Get admittance matrix
    Ybus = admittance(nbus, svc);

    % Initial Generation power
    Pgen=zeros(nbus,1);
    Qgen=zeros(nbus,1);
    Pl = busdata(:, 7) / basemva;   % Load real power
    Ql = busdata(:, 8) / basemva;   % Load reactive power
    % Extract voltage and angle from load flow results
    voltage = lfresult{3};
    angle = lfresult{4};

    % Calculate bus power
    complex_volt = voltage.*cos(angle) + 1i*voltage.*sin(angle);
    % currents at each node
    I=Ybus*complex_volt;
    % Power at each node
    S=complex_volt.*conj(I); % Complex power
    Pbus=real(S);
    Qbus=imag(S);
    for k=1:nbus
        if busdata(k,2)==1
            % Real and reactive generation at the Slack bus
            Pgen(k)=real(S(k));
            Qgen(k)=imag(S(k));
        end
        if busdata(k,2)==2
            % Real and reactive generation at the PV buses
            Pgen(k)=real(S(k))+Pl(k);
            Qgen(k)=imag(S(k))+Ql(k);
        end
        if busdata(k,2)==3
            Pgen(k)=0;
            Qgen(k)=0;
        end
    end
    result{1}=Pgen*basemva;
    result{2}=Qgen*basemva;
    result{3}=Pl*basemva;
    result{4}=Ql*basemva;
    result{5}=Pbus*basemva;
    result{6}=Qbus*basemva;
end

function Ybus = admittance(nbus, svc)
    if nbus==14
        data=case14();              % Load IEEE14 bus system data
    elseif nbus==30
        data=case30();              % Load IEEE30 bus system data
    end
    linedata = data.line;
    nl = linedata (:,1);       % From Bus
    nr = linedata (:,2);       % To Bus
    r = linedata (:,3);        % Line Resistance
    x = linedata (:,4);        % Line Reactance 
    b = linedata (:,5);        % Ground Admittance or Susceptance
    
    if (~isempty(svc))
        loc_value = round(svc (1:1));
        svc_value = svc(1+1: end);
        for km=1:length(loc_value)
            QG=svc_value(km);
            b(loc_value(km)) = b(loc_value(km) ) + (QG);
        end
    end
    a = linedata(:,6);             % Tranformer Tap(a)
    z = complex(r, x);              % Line Impedance in complex form
    y = 1./z;                       % Line Admittance in complex form
    b = complex(0, b);              % Line Susceptance in complex form
    
    nbus = max(max(nl), max(nr));   % Number of Bus
    nbr = length(nl);                % Number of Line

    % Initial Admittance Matrix
    Ybus = zeros(nbus, nbus);
    
    for k=1:nbr
        Ybus (nl(k), nr(k)) = Ybus(nl(k), nr(k)) - y(k)/a(k);
        Ybus (nr(k), nl(k)) = Ybus(nl(k), nr(k));
    end
    for m = 1:nbus
        for n = 1:nbr
            if nl (n) == m
            Ybus(m, m) = Ybus(m, m) + y(n)/(a(n)^2) + b(n);
            elseif nr(n) == m
                Ybus(m, m) = Ybus(m,m) + y(n) + b(n);
            end
        end
    end
end

