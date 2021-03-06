% EE 4753/EE 5243 Analysis of Power Systems - N. Gatsis
% Gauss-Seidel on 2-bus network
% Bus 1 is slack; bus 2 is PQ

% Bus admittance matrix: y12 has magnitude 5, angle -80 degrees
Y=[5*exp(-1i*deg2rad(80)), 5*exp(1i*deg2rad(100)); ...
   5*exp(1i*deg2rad(100)), 5*exp(-1i*deg2rad(80))];

% Specified quantities
V1=1; % slack bus
PL2=2; % Load at bus 2
QL2=1;

% Set up P and Q injections
P2=-PL2;
Q2=-QL2;

% initialization
V2=1;

% max number of iterations and initial iteration index
maxIter=40;
ind=1;

% vector of iterates
V2iter=zeros(1,maxIter+1);
V2iter(1)=V2;

V2temp=0;
tol=1e-4;
while (ind<=maxIter)&&(abs(V2-V2temp)>=tol),
    V2temp=V2;
    V2=(1/Y(2,2))*((P2-1i*Q2)/V2'-Y(2,1)*V1);
    ind=ind+1;
    V2iter(ind)=V2;
end

S1=V1*(Y(1,1)*V1+Y(1,2)*V2)';
P1=real(S1);
Q1=imag(S1);


figure(3)
plot(sqrt(conj(V2iter(:,1:ind)).*V2iter(:,1:ind))') % Plot the magnitude 
title('Magnitude of V_2(i) (pu)') 
xlabel('Iteration index i')
figure(4)
plot(atan2(imag(V2iter(:,1:ind))',real(V2iter(:,1:ind))')) % Plot the phase
title('Angle of V_2(i) (rad)')
xlabel('Iteration index i')