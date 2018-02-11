% EE 4753/EE 5243 Analysis of Power Systems - N. Gatsis
% Gauss-Seidel on an N-bus network; general code applied to N=4 example
% Bus 1 is slack bus; buses 2, 3, ..., N are PQ buses

N=4;

P=[0; -0.8; -0.6; -0.5]; % P(1) will be found at the end
Q=[0; -0.1; -0.1; -0.2]; % Q(1) will be found at the end
V=zeros(N,1); % Vector of voltages for N buses
              % V(1) will be constant; 
              % V(2:N) will be computed in each iteration
V1=1+0i;
V(1)=V1;

% Series admittances (1 is connected to 2, 2 to 3, 3 to 4, 4 to 1)
y12=0.04-5i; % series admittance of line (1,2)
y23=0.05-6i; 
y34=0.056-8.66i; % series admittance of line (3,4)
y14=0.06-7i;

% Define the bus admittance matrix
Y=[y12+y14       -y12         0        -y14;
    -y12        y12+y23     -y23         0;
    0            -y23       y23+y34     -y34;           
    -y14           0          -y34       y14+y34];


V(2:N)=1+0i; % initilization: flat start

maxIter=150; % Maximum number of iterations allowed
vIter=zeros(N,maxIter+1); % Matrix to save my iterates
vIter(:,1)=V; % Save the initialization
tol=1e-4; % Desired relative error between two successive iterates
Vtemp=zeros(N,1); % Vector to save the voltages before they are changed
                  % It is needed for checking the convergence criterion
ind=1; % iteration index
while (ind<=maxIter)&&(max(abs(V-Vtemp))>=tol),
    %ind
    Vtemp=V; % Save previous voltages and start the iteration
    for n=2:N,
        V(n)=(1/Y(n,n))*((P(n)-(Q(n)*1i))/(V(n)')-...
            Y(n,[1:(n-1),(n+1):N])*V([1:(n-1),(n+1):N]));
        fprintf('Iteration %i, Updated Bus %i\n', ind, n);
%        keyboard  % This will pause execution and return to keyboard 
                   % so you can examine the value of different variables
                   % (see help for 'keyboard' command)
                   % Type 'dbcont' to continue execution 
                   % Type 'dbquit' to exit execution
    end
    ind=ind+1; % advance iteration index
    vIter(:,ind) = V; % save the latest iterate
end

S1=V(1)*(Y(1,:)*V)';
P1=real(S1);
Q1=imag(S1);

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

