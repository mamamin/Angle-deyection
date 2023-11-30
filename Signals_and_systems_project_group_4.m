clear 'all';clc;
load("400109768.mat");
N=40;
M=2;
L=4;
my_angles=zeros(M,1000);
for j=1:1000
    data=recieve(:,j);
    data=reshape(data,1,40);
    X0=zeros(N-L,L);
    for i=1:N-L
        X0(i,:)=data(i:i+L-1);
    end
    X1=zeros(N-L,L);
    for i=1:N-L
        X1(i,:)=data(i+1:i+L);
    end
    X0H=X0';
    X1H=X1';
    my_inv=inv(X0H*X0);
    z_m=eig(my_inv*(X0H*X1));
    k=0;
    for o=1:L
        if(abs(z_m(o))>0.9 && abs(z_m(o))<1.1)
           k=k+1;
           my_angles(k,j)=asin(-atan(imag(z_m(o))/real(z_m(o)))/pi)/pi*180;
        end
    end
    if(my_angles(2,j)>my_angles(1,j))
        a=my_angles(2,j);
        my_angles(2,j)=my_angles(1,j);
        my_angles(1,j)=a;
    end
end
subplot(1,2,1);
histogram(my_angles(1,:));
title("angle of first source");
xlabel("degree");
ylabel("number of snapshots");
subplot(1,2,2);
histogram(my_angles(2,:));
title("angle of second source");
xlabel("degree");
ylabel("number of snapshots");