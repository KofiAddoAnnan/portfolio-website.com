clc
clear all
% In this program, the Gauss-Seidel Method is used to obtain the voltage
% values of a three(3) bus system.Also, the slack bus power, line flows and
% line losses are found.

fprintf('<strong>Name:</strong> Kofi Addo Annan\n')
fprintf('<strong>Index Number:</strong> PG%d%d%d%d%d%d%d\n',7,3,4,6,5,2,3)
fprintf('<strong>Student ID:</strong> %d%d%d%d%d%d%d%d\n',2,0,9,8,7,8,5,2)
fprintf('<strong>Course:</strong> EE%d%d%d - Energy System Analysis\n',5,5,3)
fprintf('<strong>Programme:</strong> MPhil. Power Systems Engineering\n\n')

%Objective of this program
fprintf('<strong>Objective of this program</strong>\n')
fprintf('<strong>Using the Guass-Seidel Method, write a MATLAB program to obtain:</strong>\n')
fprintf('1.) The voltages at bus 2 and 3 of the three(3)-bus network.\n')
fprintf('2.) Calculate the slack bus power, line flows and line losses.\n')
fprintf('<strong>A tolerance of 0.001 is to be used for both the real and imaginary parts of the voltages.</strong>\n')

%Define the data for the bus system
Busbar=['1';'2';'3'];
Voltage_pu = [0.99+0i; 1+0i; 1+0i]; %Assume a value of V2=V3=1+0i
Generation_pu = [0; 0.25-0.1i; 0.4+0.1i]; %Assume a 0 value for slack bus
Load_pu = [1+2i; 0.2+0.1i; 0];

data = table(Busbar,Voltage_pu,Generation_pu,Load_pu)

fprintf('<strong>Obtain the admittance matrix from the data given:</strong>\n')
%Create the admittance matrix from data
Y12=2-4i;
Y13=1-2i;
Y23=1-2i;

Ybus=zeros(3);
%row one
Ybus(1,1)= Y12+Y13;
Ybus(1,2)=-(Y12);
Ybus(1,3)=-(Y13);
%row 2
Ybus(2,1)=-(Y12);
Ybus(2,2)= Y12+Y23;
Ybus(2,3)=-(Y23);
%row 3
Ybus(3,1)=-(Y13);
Ybus(3,2)=-(Y23);
Ybus(3,3)=Y13+Y23;

%Display Ybus
Ybus

%Determine the bus 2 specified apparent power(s)
fprintf('<strong>Determine the bus 2 specified injected apparent power(S2):</strong>\n')
S2 = Generation_pu(2) - Load_pu(2)

%Real Power P2
P2 = real(S2)
%Reactive Power Q2
Q2 = imag(S2)

%Determine the bus 3 specified apparent power(s)
fprintf('<strong>Determine the bus 3 specified injected apparent power(S3):</strong>\n')
S3 = Generation_pu(3) - Load_pu(3)

%Real Power P3
P3 = real(S3)
%Reactive Power Q3
Q3 = imag(S3)

%Define the tolerance value
tolerance = 0.001;

%Gauss-Seidel Iteration
for k = 1:50


V2_new=1/Ybus(2,2)*(conj(S2)/conj(Voltage_pu(2))-Ybus(2,1)*Voltage_pu(1)-Ybus(2,3)*Voltage_pu(3));

V3_new=1/Ybus(3,3)*(conj(S3)/conj(Voltage_pu(3))-Ybus(3,1)*Voltage_pu(1)-Ybus(3,2)*V2_new);


fprintf('After %d Iteration',k)
V2_new
V3_new


error = max(abs([real(V2_new)-real(Voltage_pu(2)) imag(V2_new)-imag(Voltage_pu(2)) real(V3_new)-real(Voltage_pu(3)) imag(V3_new)-imag(Voltage_pu(3))]))
if error < tolerance
    break;
end
Voltage_pu(2) = V2_new;
Voltage_pu(3) = V3_new;
end

%Display the iterated voltage values
fprintf('<strong>Display the iterated voltage values for bus 1 and 2 respectively:</strong>\n')
disp('Final Voltages After Iteration:')
disp(['V2: ', num2str(abs(V2_new)), ' <', num2str(angle(V2_new)*180/pi),' ', 'degrees',' ', 'or',' ',num2str(V2_new)])
disp(['V3: ', num2str(abs(V3_new)), ' <', num2str(angle(V3_new)*180/pi),' ', 'degrees',' ', 'or',' ',num2str(V3_new)])

%Slack bus power calculaton
fprintf('<strong>Slack bus power calculation:</strong>\n')
disp('S(P+jQ) =V * conj(I) or S(P-jQ) = conj(V) * I')
I1=(Ybus(1,1)*Voltage_pu(1)) + (Ybus(1,2)*V2_new) + (Ybus(1,3)*V3_new);
S1 = Voltage_pu(1)*conj(I1);
S1_conj = conj(Voltage_pu(1))*I1;
disp('Slackbus Power:')
disp(['S1_(P+jQ): ', num2str(S1),' ', 'or',' ', 'S1_(P-jQ): ',num2str(S1_conj)])

%Determine the line flows of the network
fprintf('<strong>Determine the line flows of the network:</strong>\n')
I12=Ybus(1,2)*(Voltage_pu(1)-V2_new);
I21=-I12;
I13=Ybus(1,3)*(Voltage_pu(1)-V3_new);
I31=-I13;
I23=Ybus(2,3)*(V2_new-V3_new);
I32=-I23;

S12= Voltage_pu(1)*conj(I12);
S21= V2_new*conj(I21);
S13= Voltage_pu(1)*conj(I13);
S31= V3_new*conj(I31);
S23= V2_new*conj(I23);
S32= V3_new*conj(I32);
disp('Line Flows:')
disp(['S12: ', num2str(S12)])
disp(['S21: ', num2str(S21)])
disp(['S13: ', num2str(S13)])
disp(['S31: ', num2str(S31)])
disp(['S23: ', num2str(S23)])
disp(['S32: ', num2str(S32)])

%Determine the line losses of the network
fprintf('<strong>Determine the line losses of the network:</strong>\n')
SLoss12= S12+S21;
SLoss13= S13+S31;
SLoss23= S23+S32;
disp('Line Losses:')
disp(['SLoss12: ', num2str(SLoss12)])
disp(['SLoss13: ', num2str(SLoss13)])
disp(['SLoss23: ', num2str(SLoss23)])
