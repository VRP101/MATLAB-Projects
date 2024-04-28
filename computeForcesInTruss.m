function f = computeForcesInTruss(loads)

         a = 1/sqrt(2);

         A = zeros(13,13); % Initialize A matrix

         b = zeros(13,1); % Initialize b vector

         % from Eqn 1

         A(1,2)=-1;

         A(1,6)=1;

         % from Eqn 2

         A(2,3)=1;

         % from Eqn 3

         A(3,1)=-a;

         A(3,4)=1;

         A(3,5)=a;

         % from Eqn 4 

         A(4,1)=-a;

         A(4,3)=-1;

         A(4,5)=-a;

         % from Eqn 5

         A(5,4)=-1;

         A(5,8)=1;

         % from Eqn 6

         A(6,7)=1;

         % from Eqn 7

         A(7,5)=-a;

         A(7,6)=-1;

         A(7,9)=a;

         A(7,10)=1;

         % from Eqn 8

         A(8,5)=a;

         A(8,7)=1;

         A(8,9)=a;

         % from Eqn 9

         A(9,10)=-1;

         A(9,13)=1;

         % from Eqn 10

         A(10,11)=1;

         % from Eqn 11

         A(11,8)=-1;

         A(11,9)=-a;

         A(11,12)=a;

         % from Eqn 12

         A(12,9)=-a;

         A(12,11)=-1;

         A(12,12)=-a;

         % from Eqn 13

         A(13,13)=-1;

         A(13,12)=-a;

        

         b(2)=loads(1);

         b(8)=loads(2);

         b(10)=loads(3);

        

         f = solveLinearSystem(A,b);

end

 



function [L, U, P] = luDecomposition(A)

    [n_row, n_column] = size(A); 

    L=eye(n_column); 

    P=eye(n_column); 

    U=A;

    for k=1:n_row-1

        pivot=max(abs(U(k:n_row,k)));

        for j=k:n_row

            if(abs(U(j,k))==pivot)

                ind=j;

                break;

            end

        end



        % swap rows in U

        temp = U(k,k:n_row);

        U(k,k:n_row)=U(ind,k:n_row);

        U(ind,k:n_row)=temp;



        % swap rows in L

        temp = L(k,1:k-1);

        L(k,1:k-1)=L(ind,1:k-1);

        L(ind,1:k-1)=temp;



        % swap rows in P

        temp = P(k,:);

        P(k,:)=P(ind,:);

        P(ind,:)=temp;



        for j=k+1:n_row

            scale=U(j,k)/U(k,k);

            U(j,k:n_row)=U(j,k:n_row)-scale*U(k,k:n_row);

            L(j,k)=scale;

        end

    end

end



 %% step2: Forward substitution



function y = forwardSubstitution(L,b,P)

    N = size(L,2);

    y = zeros(N,1);

    b = P*b; % Rearrange rows in b matrix based on P matrix

    if L(1,1)~=0

       y(1)=b(1)/L(1,1);

    else

       y(1)=0;

    end

    for i =2:N

        if L(i,i)~=0

           y(i)=(b(i)-L(i,1:i-1)*y(1:i-1))/L(i,i);

        else

           y(i)=0;

        end

    end

end



%% Step3: Back substitution



function x = backSubstitution(U,y)

    N = size(U,2);

    x = zeros(N,1);

    if U(N,N)~=0

       x(N)=y(N)/U(N,N);

    else

       x(N)=0;

    end

    for i=1:N-1

       if U(i,i)~=0

           x(N-i)=(y(N-i)-U(N-i,N-i+1:N)*x(N-i+1:N))/U(N-i,N-i);

       else

           x(N-i)=0;

       end

    end

end



%% Linear solve



function x = solveLinearSystem(A,b)

   [L,U,P] = luDecomposition(A);

   y = forwardSubstitution(L,b,P);

   x = backSubstitution(U,y);



end