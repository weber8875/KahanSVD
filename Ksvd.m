function[S]=Ksvd(A)

n=size(A,1);

U = eye(n);
H = A;
V = eye(n);

for k=1:n-1
    x = H(k:n, k);
    [u, sigma] = housevec(x);
    beta = -2 / (u' * u);
    H(k,k) = sigma;
    H(k+1:n,k) = 0;

    v = U(k:n,k:n);
    U(k:n,k:n) = v + beta*(v*u)*u';
    for j = k+1:n
        v = H(k:n, j);
        H(k:n, j) = v + beta * (u' * v) * u;
    end

    %============================

    if  k < n
        y = H(k, k+1:n)'; 
        [u,sigma] = housevec(y);
        beta = -2 /(u'*u);
        H(k,k+1)=sigma;
        H(k,k+2:n)=0;
        for i=k+1:n
        v=H(i,k+1:n);
        H(i,k+1:n)= v + beta*(v*u)*u';
        end

        v = V(k:n,k+1:n);
        V(k:n,k+1:n)= v + beta*(v*u)*u';
       

    end
end
     for i =1:50
        u = [H(1,1)^2,H(1,1)*H(1,2)];
        [c, s] = givens(u);
    
    
        for j = 1:n
           w = c*H(j,1)-s*H(j,2);
           x = s*H(j,1)+c*H(j,2);
    
            H(j,1) = w;
            H(j,2) = x;
        end  
    
        %============================
    
        u=[H(1,1),H(2,1)];
        [c, s] = givens(u);
    
        for j = 1:n
            y = c*H(1,j)-s*H(2,j);
            z = s*H(1,j)+c*H(2,j);
    
            H(1,j) = y;
            H(2,j) = z;
        end
    
        %============================
    
        for i=2:n-1
    
            u=[H(i-1,i),H(i-1,i+1)];
            [c, s] = givens(u);
    
            for j = 1:n
                w = c*H(j,i)-s*H(j,i+1);
                x = s*H(j,i)+c*H(j,i+1);
    
                H(j,i) = w;
                H(j,i+1) = x;
    
            end
    
            u=[H(i,i),H(i+1,i)];
            [c, s] = givens(u);
    
            for j = 1:n
                y = c*H(i,j)-s*H(i+1,j);
                z = s*H(i,j)+c*H(i+1,j);
    
                H(i,j) = y;
                H(i+1,j) = z;
            end
    
        end
        
    end
    H;
    S=abs(diag(H));
    

end


