function Phi = PHI(X, Y, kernel,  c, alpha, d, q, Q, beta, a)


    N = size(X,1);
    
    M = N+1;
    
    switch kernel
        
        case 'Linear'
    
            Phi = X*Y' + c;
    
        case 'Polynomial'
    
            Phi = (alpha*X*Y' + c).^d;
    
        case 'Gauss'
    
            NormSqrd = sum(X.^2,2)*ones(1,N) + ones(N,1)*sum(Y.^2,2)' - 2*(X*Y');
            Phi = exp(-(NormSqrd/(2*(q^2))));
    
        case 'Exponential'

            NormSqrd = sum(X.^2,2)*ones(1,N) + ones(N,1)*sum(Y.^2,2)' - 2*(X*Y');
            Phi = exp(-(sqrt(NormSqrd)/(2*(q^2))));

        case 'Laplacian'
    
            NormSqrd = sum(X.^2,2)*ones(1,N) + ones(N,1)*sum(Y.^2,2)' - 2*(X*Y');
            Phi = exp(-(sqrt(NormSqrd)/q));
            
        case 'Sigmoid'
         
             Phi = tan(alpha*X*Y' + c); 
             
        case 'Rational_Quadratic'
     
             NormSqrd = sum(X.^2,2)*ones(1,N) + ones(N,1)*sum(Y.^2,2)' - 2*(X*Y');
             Phi = 1 - (NormSqrd ./(NormSqrd + c));
             
        case 'Multiquadric'
    
            NormSqrd = sum(X.^2,2)*ones(1,N) + ones(N,1)*sum(Y.^2,2)' - 2*(X*Y');
            Phi = sqrt(NormSqrd + (c^2));
        
        case 'Inverse_Multiquadric'
            
            NormSqrd = sum(X.^2,2)*ones(1,N) + ones(N,1)*sum(Y.^2,2)' - 2*(X*Y');
            Phi = 1 ./ sqrt(NormSqrd + (c^2));
        
        case 'Circular'

            NormSqrd = sum(X.^2,2)*ones(1,N) + ones(N,1)*sum(Y.^2,2)' - 2*(X*Y');
            ind = find((sqrt(NormSqrd)<q)==1);
            Phi = zeros(N*N,1);
            Phi(ind) = ((2/pi)*acos(-sqrt(NormSqrd(sqrt(NormSqrd)<q))/q)-(2/pi)*(sqrt(NormSqrd(sqrt(NormSqrd)<q))/q).*sqrt(1 - ((sqrt(NormSqrd(sqrt(NormSqrd)<q))/q).^2)));
            Phi = reshape(Phi,[N,N]);
        
        case 'Spherical'
    
           NormSqrd = sum(X.^2,2)*ones(1,N) + ones(N,1)*sum(Y.^2,2)' - 2*(X*Y');
           ind = find((sqrt(NormSqrd)<q)==1);
           Phi = zeros(N*N,1);
           Phi(ind) = 1 - (3/2)*sqrt(NormSqrd(sqrt(NormSqrd)<q))/q + (1/2)*((sqrt(NormSqrd(sqrt(NormSqrd)<q))/q).^3);
           Phi = reshape(Phi,[N,N]);

        case 'Wave'
    
           NormSqrd = sum(X.^2,2)*ones(1,N) + ones(N,1)*sum(Y.^2,2)' - 2*(X*Y');
           ind = find((sqrt(NormSqrd)~=0)==1);
           Phi = zeros(N*N,1);
           Phi(ind) = (Q./sqrt(NormSqrd(sqrt(NormSqrd)~=0))).*sin(sqrt(NormSqrd(sqrt(NormSqrd)~=0))./Q);  
           Phi = reshape(Phi,[N,N]);
    
        case 'Power'
    
           NormSqrd = sum(X.^2,2)*ones(1,N) + ones(N,1)*sum(Y.^2,2)' - 2*(X*Y');
           Phi = -(sqrt(NormSqrd).^d);     
    
        case 'Log'
    
           NormSqrd = sum(X.^2,2)*ones(1,N) + ones(N,1)*sum(Y.^2,2)' - 2*(X*Y');
           Phi = -log((sqrt(NormSqrd).^d) + 1);   
           
        case 'Cauchy'
     
           NormSqrd = sum(X.^2,2)*ones(1,N) + ones(N,1)*sum(Y.^2,2)' - 2*(X*Y');
           Phi = 1 ./ (1 + ((sqrt(NormSqrd).^2)./(q^2)));  
           
           
        case 'T-Student'
    
           NormSqrd = sum(X.^2,2)*ones(1,N) + ones(N,1)*sum(Y.^2,2)' - 2*(X*Y');
           Phi = 1 ./ (1 + (sqrt(NormSqrd).^d));     
 
        end
%     elseif (strcmp(kernel,'Spline')==true)
%     
%            Phi = 1;
%            
%            for i=1:N
%                
%                 Phi = Phi * (1 + X(i)*Y(i) + X(i)*Y(i)*min([X(i) Y(i)]) - ((X(i) + Y(i)) / 2)*(min([X(i) Y(i)])^2) + (min([X(i) Y(i)])^3) / 3);
%            
%            end
%            
%     elseif (strcmp(kernel,'Cauchy')==true)
%     
%            Phi = 1 / (1 + ((norm(X-Y)^2)/(q^2)));     
%     
%     elseif (strcmp(kernel,'Chi-Square')==true)
%     
%            Phi = 0;
%            
%            for i=1:N
%            
%                 Phi = Phi + ((X(i)-Y(i))^2) / ((1/2)*(X(i)+Y(i)));
%            
%            end
%            
%            Phi = 1 - Phi;
%            
%    elseif (strcmp(kernel,'Histogram_Intersection')==true)
%     
%            Phi = 0;
%            
%            for i=1:N
%                
%                Phi = Phi + min([X(i) Y(i)]);
%                
%            end
%            
%     elseif (strcmp(kernel,'Generalized_Histogram')==true)
%     
%            Phi = 0;
%            
%            for i=1:N
%                
%                Phi = Phi + min([(abs(X(i))^alpha) (abs(Y(i))^beta)]);
%                
%            end
%     
%     elseif (strcmp(kernel,'T-Student')==true)
%     
%            Phi = 1 / (1 + (norm(X-Y)^d));     
%     
%     elseif (strcmp(kernel, 'Wavelet')==true)
%           
%           Phi = 1;
%           
%           for i = 1:N
%             
%                 Phi = Phi * (cos(1.75*((X(i)-c)/a))*exp(-(((X(i)-c)/a)^2)/2))*(cos(1.75*((X(i)-c)/a))*exp(-(((X(i)-c)/a)^2)/2));
%               
%           end
%           
%    end
%     
    
   Phi = [ones(N,1) Phi];
            

end