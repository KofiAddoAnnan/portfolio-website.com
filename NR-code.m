clc
clear
% In this program, the Newton-Raphson method of iteration is used on a two
% bus system.

fprintf('<strong>Name:</strong> Kofi Addo Annan\n')
fprintf('<strong>Index Number:</strong> PG%d%d%d%d%d%d%d\n',7,3,4,6,5,2,3)
fprintf('<strong>Student ID:</strong> %d%d%d%d%d%d%d%d\n',2,0,9,8,7,8,5,2)
fprintf('<strong>Course:</strong> EE%d%d%d - Energy System Analysis\n',5,5,3)
fprintf('<strong>Programme:</strong> MPhil. Power Systems Engineering\n\n')

%Objective of this program
fprintf('<strong>Objective of this program</strong>\n')
fprintf('<strong>Using the Newton-Raphson Method for two interconnected generating stations, Determine:</strong>\n')
fprintf('1.) The voltage angle and the generated reactive power Qg at station 2.\n')
fprintf('2.) The generated active and reactive power Pg and Qg at station 1.\n')
fprintf('<strong>A tolerance of 0.001 is to be used for power mismatch.</strong>\n')

%Define the data for the bus system
Bus=['1';'2'];
Voltage_pu = [1+0i; 1+0i]; %Assume a value of V2=1+0i where angle of V2 is zero(0)
Generation_pu = [0+0i; 20+0i]; %Assume a 0 for all unspecified values
Load_pu = [25+15i; 15+5i];

data = table(Bus,Voltage_pu,Generation_pu,Load_pu)

%Create the admittance matrix from data
fprintf('<strong>Obtain the admittance matrix from the data given:</strong>\n')
Z12 = 0.005+0.05i;
Y12 = inv(Z12);
fprintf('The Ybus admittance matrix in polar form:\n')
fprintf('Magnitude and angle values are placed in alternative columns:\n')
Ybus = zeros(2,4);
Ybus(1,1)= abs(Y12) ;
Ybus(1,2)= angle(Y12)*180/pi;
Ybus(1,3)= abs(-Y12);
Ybus(1,4)= angle(-Y12)*180/pi;
Ybus(2,1)= abs(-Y12);
Ybus(2,2)= angle(-Y12)*180/pi;
Ybus(2,3)= abs(Y12); 
Ybus(2,4)= angle(Y12)*180/pi;

%Display the Ybus matrix
Ybus

fprintf('<strong>First Objective: Determine the voltage angle(<) and the generated reactive power(QG2) of Station 2.</strong>\n')
%Determine the injected P2 at Station 2 in per unit
P2 = real(Generation_pu(2))-real(Load_pu(2));
Mag_V2 = abs(Voltage_pu(2));
Angle_V2 = angle(Voltage_pu(2))*180/pi;
Mag_V1 = abs(Voltage_pu(1))
Angle_V1 = angle(Voltage_pu(1))*180/pi;
%Determine P2 calculated
fprintf('Find the calculated P2 to obtain the power mismatch matrix:\n')
%Convert all degrees values into radians before the cos function
radYbus_2_2 = Ybus(2,2)*pi/180;
radYbus_2_4 = Ybus(2,4)*pi/180;
radYbus_1_2 = Ybus(1,2)*pi/180;
radYbus_1_4 = Ybus(1,4)*pi/180;
radAngle_V2 = Angle_V2 *pi/180;
radAngle_V1 = Angle_V1*pi/180;

tolerance = 0.001

for k = 1:5
P2_calculated = Mag_V2*Ybus(2,1)*abs(Voltage_pu(1))*cos(radYbus_2_2+radAngle_V1-radAngle_V2) + Mag_V2^2*(Ybus(2,3))*cos(radYbus_2_4);
Power_Mismatch_P2 = P2 - P2_calculated;

%Differentiate to obtain the Jacobian matrix
fprintf('Differentiate P2_calculated to find the Jacobian matrix:\n')
J11 = Mag_V2*Ybus(2,1)*abs(Voltage_pu(1))*sin(radYbus_2_2+radAngle_V1-radAngle_V2);

fprintf('Find the value the correction matrix/change in angle which is the form:\n')
fprintf('<strong>[Power Mismatch] = [Jacobian Matrix] * [Correction Matrix]</strong>\n')
fprintf('<strong>[Correction Matrix] = [Jacobian Matrix]^-1 * [Power Mismatch]</strong>\n\n')

fprintf('<strong>After %d Iteration</strong>\n',k)
Correction_matrix = pinv(J11) * Power_Mismatch_P2;

Iterated_AngleV2_radians = radAngle_V2 + Correction_matrix

Iterated_AngleV2_degrees = Iterated_AngleV2_radians * 180/pi

error = abs(Power_Mismatch_P2)
if error < tolerance
    break;
end

radAngle_V2 = Iterated_AngleV2_radians;
end

%Determine the generated reactive power for Station 2
Injected_Q2 = -Mag_V2*Ybus(2,1)*abs(Voltage_pu(1))*sin(radYbus_2_2+radAngle_V1-radAngle_V2) - Mag_V2^2*(Ybus(2,3))*sin(radYbus_2_4)
Generated_Q2 = Injected_Q2 + imag(Load_pu(2))

%Determine the generated active power Pg at Station 1
Injected_P1 = Mag_V1*Ybus(1,1)*Mag_V1*cos(radYbus_1_2) + Mag_V1*Ybus(1,3)*Mag_V2*cos(radYbus_1_4+radAngle_V2-radAngle_V1)
Generated_P1 = Injected_P1 + real(Load_pu(1))

%Determine the generated reactive power Qg at Station 1
Injected_Q1 = -Mag_V1*Ybus(1,1)*Mag_V1*sin(radYbus_1_2) - Mag_V1*Ybus(1,3)*Mag_V2*sin(radYbus_1_4+radAngle_V2-radAngle_V1)
Generated_Q1 = Injected_Q1 + imag(Load_pu(1))

%Display the results
%Display the iterated voltage angle value
fprintf('<strong>Display the iterated voltage angle for Station 2:</strong>\n')
disp('Final Voltage Angle After Iteration:')
disp(['V2: ', num2str(abs(Voltage_pu(2))), ' <', num2str(Iterated_AngleV2_degrees),' ', 'degrees',' ', 'or',' ', 'V2: ', num2str(abs(Voltage_pu(2))), ' <', num2str(Iterated_AngleV2_radians),' ', 'radians' ])
%Display the generated reactive power for station 2
fprintf('<strong>Display the generated reactive power Qg for Station 2:</strong>\n')
disp('Generated Reactive Power For Station 2:')
disp(['Qg2: ', num2str(Generated_Q2), 'p.u'])
%Display the generated active power Pg at Station 1
fprintf('<strong>Display the generated active power Pg for Station 1:</strong>\n')
disp('Generated Active Power For Station 1:')
disp(['Pg1: ', num2str(Generated_P1), 'p.u'])
%Display the generated reactive power Qg at Station 1
fprintf('<strong>Display the generated reactive power Qg for Station 1:</strong>\n')
disp('Generated Reactive Power For Station 1:')
disp(['Qg1: ', num2str(Generated_Q1), 'p.u'])