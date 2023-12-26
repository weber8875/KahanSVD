function[U,V]=Ksvd2(A)

n=size(A,1);

U = eye(n);
H = A;
V = eye(n);

for k=1:n
    x = H(k:n, k);
    [u, sigma] = housevec(x);
    beta = -2 / (u' * u);
    H(k,k) = sigma;
    H(k+1:n,k) = 0;

    R = U(1:n,k:n);
    U(1:n,k:n) = R + beta*(R*u)*u';

    for j = k+1:n
        v = H(k:n, j);
        H(k:n, j) = v + beta * (u' * v) * u;
    end

    

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

        T = V(k+1:n,2:n);
        V(k+1:n,2:n)= T + beta * u * u' * T;
       

    end
end
for s =1:20
    
    u = [H(1,1)^2,H(1,1)*H(1,2)];
    [c, s] = givens(u);
    
    for j = 1:n
       w = c*H(j,1)-s*H(j,2);
       x = s*H(j,1)+c*H(j,2);
    
       H(j,1) = w;
       H(j,2) = x;          

       w2 = c*V(1,j)-s*V(2,j);
       x2 = s*V(1,j)+c*V(2,j);

       V(1,j) = w2;   
       V(2,j) = x2;
    end    
    
    %============================
    
    u=[H(1,1),H(2,1)];
    [c, s] = givens(u);
    
    for j = 1:n
        y = c*H(1,j)-s*H(2,j);
        z = s*H(1,j)+c*H(2,j);
    
        H(1,j) = y;
        H(2,j) = z;

        y2 = c*U(j,1)-s*U(j,2);
        z2 = s*U(j,1)+c*U(j,2);

        U(j,1) = y2;
        U(j,2) = z2;

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
                
                w2 = c*V(i,j)-s*V(i+1,j);
                x2 = s*V(i,j)+c*V(i+1,j);

                V(i,j) = w2;
                V(i+1,j) = x2;
            end

        u=[H(i,i),H(i+1,i)];
        [c, s] = givens(u);
    
            for j = 1:n
                y = c*H(i,j)-s*H(i+1,j);
                z = s*H(i,j)+c*H(i+1,j);
    
                H(i,j) = y;
                H(i+1,j) = z;

                y2 = c*U(j,i)-s*U(j,i+1);
                z2 = s*U(j,i)+c*U(j,i+1);

                U(j,i) = y2;
                U(j,i+1) = z2;
            end      
    
    end
end


V=V'




