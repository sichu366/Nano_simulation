% EE 4753/EE 5243 Analysis of Power Systems - N. Gatsis
% Gauss-Seidel on an N-bus network; general code applied to N=4 example
% Bus 1 is slack bus; buses 2, 3, ..., N are PQ buses
clear;
N=4;

P=[0; 0.5; -1; 0.3]; % P(1) will be found at the end
Q=[0; -0.2; 0.5; -0.1]; % Q(1) will be found at the end
V=zeros(N,1); % Vector of voltages for N buses
              % V(1) will be constant; 
              % V(2:N) will be computed in each iteration
V1=1.04+0i;
V(1)=V1;

% Define the bus admittance matrix
Y=[3-9*1i       -2+6*1i         -1+3*1i       0;
    -2+6*1i     3.666-11*1i     -0.666+2*1i    -1+3*1i;
    -1+3*1i     -0.666+2*1i      3.666-11*1i   -2+6*1i;           
    0           -1+3*1i          -2+6*1i       3-9*1i];


V(2:N)=1+0i; % initilization: flat start

maxIter=150; % Maximum number of iterations allowed
vIter=zeros(N,maxIter+1); % Matrix to save my iterates
vIter(:,1)=V; % Save the initialization
tol=1e-4; % Desired relative error between two successive iterates
Vtemp=zeros(N,1); % Vector to save the voltages before they are changed
                  % It is needed for checking the convergence criterion
ind=1; % iteration index

while (ind<=maxIter)&&(max(abs(V-Vtemp))>=tol)
    %ind
    Vtemp=V; % Save previous voltages and start the iteration
    V(2)=1/Y(2,2)*((P(2)-(Q(2)*1i))/V(2)'-Y(2,1)*V(1)-Y(2,3)*V(3)-Y(2,4)*V(4));
    V(3)=1/Y(3,3)*((P(3)-(Q(3)*1i))/V(3)'-Y(3,1)*V(1)-Y(3,2)*V(2)-Y(3,4)*V(4));
    V(4)=1/Y(4,4)*((P(4)-(Q(4)*1i))/V(4)'-Y(4,1)*V(1)-Y(4,2)*V(2)-Y(4,3)*V(3));
    ind=ind+1; % advance iteration index
    vIter(:,ind) = V; % save the latest iterate
end

S1=V(1)*(Y(1,:)*V)';
P1=real(S1);
Q1=imag(S1);
vIter(:,1:ind)
figure(1)
plot(real(vIter(:,1:ind))') % Plot the real part of the voltage iterates
title('Re[V(i)]')
xlabel('Iteration index i')
figure(2)
plot(imag(vIter(:,1:ind))') % Plot the imaginary part of the voltage iterates
title('Im[V(i)]')
xlabel('Iteration index i')
figure(3)
plot(sqrt(conj(vIter(:,1:ind)).*vIter(:,1:ind))') % Plot the magnitude 
title('Magnitude of V(i)') 
xlabel('Iteration index i')
figure(4)
plot(atan2(imag(vIter(:,1:ind))',real(vIter(:,1:ind))')) % Plot the phase
title('Angle of V(i)')
xlabel('Iteration index i')

