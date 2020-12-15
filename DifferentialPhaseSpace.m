function Y_predicted = DifferentialPhaseSpace(x,m,tau)

N = length(x);

%Enlarge the time series xn, 1 <= n <= N to -m <= n <= N + m + 1 
N_new = 2 * (m + 1) + N;
 
xn = zeros(1,N_new);
xn(1:m+1) = (x(2-[-m:0]) - x(1));
xn(m+2:N_new-(m+1)) = x;
xn(N_new-m:N_new) = x(2*N-[N+1:N_new-(m+1)]) - x(N);

%xn = x;

%N_new = length(xn);

%N = N_new - 2 * (m+1);

%x = xn((m+1)+1:end-m-1);


%The polynomial filter to estimate the numerical derivatives 

%The polynomial in the form of matrix
T = zeros((m+1)*2+1, (m+1)+1);
T(:,1) = 1;
T(:,2:end) = ((meshgrid(-(m+1):m+1,1:(m+1))'*tau).^meshgrid(1:(m+1),-(m+1):m+1))./factorial(meshgrid(1:(m+1),-(m+1):m+1));


U = inv(T'*T)*T'; %U is an ((m + 2) x (2m + 3)) matrix

B = xn(meshgrid(1:N,0:2*m+2)+meshgrid(0:2*m+2,1:N)');
Z = U*B(:,1:N); %Estimation of derivation by the least square method
Z(1,:) = x';


 %Filter to correct the estimated results Z

Y = Z(1:m,:); %The state matrix (the derivative matrix)

Tp =  triu(ones(m,m)); %The downward transmitting matrix

Tp = Tp.*tau.^(meshgrid(0:m-1,0:m-1) + (-1)*meshgrid(0:m-1,0:m-1)')./factorial(abs(meshgrid(0:m-1,0:m-1) + (-1)*meshgrid(0:m-1,0:m-1)'));

I = eye(m); %The identity matrix

Tp2 = Tp*Tp;

Kb = inv(I + Tp'*Tp + Tp2'*Tp2);

W1b = Kb;
W2b = Kb*Tp';
W3b = Kb*Tp2';

Tn =  triu(ones(m,m)); %The upward transmitting matrix

Tn = Tn.*(-tau).^(meshgrid(0:m-1,0:m-1) + (-1)*meshgrid(0:m-1,0:m-1)')./factorial(abs(meshgrid(0:m-1,0:m-1) + (-1)*meshgrid(0:m-1,0:m-1)'));

Tn2 = Tn*Tn;

Kf = inv(I + Tn'*Tn + Tn2'*Tn2);

W1f = Kf;
W2f = Kf * Tn';
W3f = Kf * Tn2';

if N >= 4

    Y_predicted = zeros(m, N); %Reconstructed differential phase space
      
	%Downward filtering
	Y_predicted(:,1:2) = (W1b*Y(:,meshgrid(0,1:2)+meshgrid(1:2,0)')) + (W2b*Y(:,meshgrid(1,1:2)+meshgrid(1:2,1)')) + (W3b*Y(:,meshgrid(2,1:2)+meshgrid(1:2,2)'));
	
	%Forward filtering
    Y_predicted(:,3:N) = (W1f*Y(:,meshgrid(0,3:N)+meshgrid(3:N,0)')) + (W2f*Y(:,meshgrid(-1,3:N)+meshgrid(3:N,-1)')) + (W3f*Y(:,meshgrid(-2,3:N)+meshgrid(3:N,-2)'));
    
    Y_predicted(1,:) = Y(1,:);
    

else
    
    Y_predicted = Y;
    
end

end
