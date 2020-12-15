function [Mi, Sel, Beta] = WRVM(Phi, t)

N = size(Phi, 1);

M = size(Phi, 2);

Max_Real_Step = 1000;

Q2 = var(t) * 0.1;

Sel = zeros(M,1);

Beta(1:N,1) = 1/Q2;

B = diag(Beta);
  
Alpha(1:M,1) = exp(1000);

Arg = (Phi'*t)./diag(sqrt(Phi'*Phi));

[val, i] = max(Arg);

Sel(i) = 1;

Alpha(i) = (norm(Phi(:,i))^2) / ((norm(Phi(:,i)'*t)^2) / (norm(Phi(:,i))^2) - Q2);

A = diag(Alpha(Sel==1));

Mi = (Phi(:,Sel==1)'*B*Phi(:,Sel==1)+A)\Phi(:,Sel==1)'*B*t; 

S = diag(Phi'*B*Phi - Phi'*B*Phi(:,Sel==1)*((Phi(:,Sel==1)'*B*Phi(:,Sel==1)+A)\Phi(:,Sel==1)')*B*Phi);

Q = Phi'*B*t - Phi'*B*Phi(:,Sel==1)*((Phi(:,Sel==1)'*B*Phi(:,Sel==1)+A)\Phi(:,Sel==1)')*B*t;

s = S;
s(Sel==1) = Alpha(Sel==1).*S(Sel==1) ./ (Alpha(Sel==1)-S(Sel==1));

q = Q;
q(Sel==1) = Alpha(Sel==1).*Q(Sel==1) ./ (Alpha(Sel==1)-S(Sel==1));

Real_Step = 0;

R = 5;

i = 1;

for step=1:Max_Real_Step
    
    Qi = (q(i)^2) - s(i);
    
    Did_Nothing = 'false';
	
	Sel_Old = Sel;
    
    if (Qi > 0)
        
        if Alpha(i) < exp(1000)
           Did_Nothing = 'reestimate';
 
        elseif Alpha(i) == exp(1000)
            
            Sel(i) = 1;
            Did_Nothing = 'add';

        end
    else
        
        if Alpha(i) < exp(1000)
            
            Sel(i) = 0;
            Did_Nothing = 'remove';
        
        else
            
            Did_Nothing = 'true';
            
        end
        
    end
    
    if strcmp(Did_Nothing,'true') == false
        
        Real_Step = Real_Step + 1;
        
        Alpha_Old = Alpha;
    
        if strcmp(Did_Nothing,'reestimate') == true
            
              Alpha(i) = (s(i)^2)/((q(i)^2) - s(i));
    
         elseif strcmp(Did_Nothing,'add') == true
        
              Alpha(i) = (s(i)^2)/((q(i)^2) - s(i));
    
        elseif strcmp(Did_Nothing,'remove') == true
        
              Alpha(i) = exp(1000);
        
        end
        
        if mod(Real_Step, R) == 0 %Reestimate noise variance every 5 step
		
		   y = Phi(:,Sel_Old==1)*Mi;
		
		   r = t-y;
    
           l = 3.7773*median(abs(r));
    
           Beta(1:N,1) = Ksi(r,l) / l;

		   B = diag(Beta);
        
        end
    
        A = diag(Alpha(Sel==1));

        Mi = (Phi(:,Sel==1)'*B*Phi(:,Sel==1)+A)\Phi(:,Sel==1)'*B*t;
        
        S = diag(Phi'*B*Phi - Phi'*B*Phi(:,Sel==1)*((Phi(:,Sel==1)'*B*Phi(:,Sel==1)+A)\Phi(:,Sel==1)')*B*Phi);

        Q = Phi'*B*t - Phi'*B*Phi(:,Sel==1)*((Phi(:,Sel==1)'*B*Phi(:,Sel==1)+A)\Phi(:,Sel==1)')*B*t;

        s = S;
        s(Sel==1) = Alpha(Sel==1).*S(Sel==1) ./ (Alpha(Sel==1)-S(Sel==1));

        q = Q;
        q(Sel==1) = Alpha(Sel==1).*Q(Sel==1) ./ (Alpha(Sel==1)-S(Sel==1));
        
        if ~isempty(find(Sel==1, 1)) && ~isempty(find(Sel==0, 1))
            
             logAlpha = abs(log(Alpha(Sel==1)) - log(Alpha_Old(Sel==1)));
             Qi = (q(Sel==0).^2) - s(Sel==0);
             
             if (max(logAlpha) < 10^(-6)) && (max(Qi) <= 0)
          
                break;
           
             end

        elseif ~isempty(find(Sel==1, 1)) && isempty(find(Sel==0, 1))
            
            logAlpha = abs(log(Alpha) - log(Alpha_Old));
            
            if (max(logAlpha) < 10^(-6))
          
                break;
           
            end
            
        end
        
    end %Did Nothing: No

    i = i+1;
    
    if i > M
        
        i = 1;
        
    end

end

end