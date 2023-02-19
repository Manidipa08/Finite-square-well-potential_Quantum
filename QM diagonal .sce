//27/7/22
clc
clear
n=1001 //input("enter number of datapoints : ")
h=6.63D-34//plank constant
m=9.1D-31//mass of electrion
e=1.6D-19
xt_1=-0.2//input("enter xt_1 : ")
//xmin = input("enter the value of xmin : ")
xmin=2*xt_1
//xmax = input("enter the value of xmax : ")
xt_2=0.2 //input("enter xt_2 : ")
xmax=2*xt_2
vt=xt_2
//v0=-100
v0 =-100//input("enter the potential : ")
x=linspace(xmin,xmax,n)//range
new = x(2:n-1)*1D-9
xnew = new'
dx=(x(2)-x(1))*1D-9
C=-((h/(2*%pi))^2)/(2*m*(dx^2)*e)
V=zeros(1,n)
for i=1:n
    if abs(x(i))-vt<=1D-6
        V(i)=v0
    end
end
//C=2
//V=1
A=eye(n-2,n-2)
v=diag(V(2:n-1))
//disp(A)
//----------------------By inbuilt command-------------------------
tic()
D=(-2*C)*ones(n-2,1)
A1=diag(D)
//disp(A1)
C1=C*ones(n-3,1)
A2=diag(C1,1)
//disp(A2)
A3=diag(C1,-1)
//disp(A3)
//Hermitian matrix (tridiagonal)
H=A1+A2+A3+v
//disp("By inbuilt : ",H)
[a1 a2]=spec(H)
Z=spec(H)
a2=diag(a2)
//disp("Corresponding eigen vectors : ",a1)
//disp("Eigen values : ",a2)
counter = 0
for i=1:n-2
    if a2(i)<0
        counter=counter+1
        no_EV(i)=a2(i)
    end
end
disp("No. of Bound States : "+string(counter))
disp("No. of Eigen Values : ",no_EV)

//--------------------------Normalizations-------------------------
vnew = V(2:n-1)
for i = 1:counter
    eigenvec = a1(:,i).*a1(:,i)
    Nrm(i)=inttrap(xnew,eigenvec)
    wf(:,i)=(1/sqrt(Nrm(i))).*a1(:,i)
    neigen(:,i) = wf(:,i).*wf(:,i)
    neigen1(:,i) = xnew.*neigen(:,i)
    neigen2(:,i) = xnew.*neigen1(:,i)
    neigen3(:,i) = vnew'.*neigen(:,i)
    ex_X(i)=inttrap(xnew,neigen1(:,i))
    ex_X2(i)= inttrap(xnew,neigen2(:,i))
    ex_V(i) = inttrap(xnew,neigen3(:,i))
end
for i = 1:n-3
    for j = 1:counter
        x2(i) = (xnew(i)+xnew(i+1))/2
        mid_psi(i,j) = (wf(i,j)+wf(i+1,j))/2
        diff_psi(i,j) = (wf(i+1,j)-wf(i,j))/dx
    end
end

for i = 1:n-4
    for j = 1:counter
        x3(i) = (x2(i) + x2(i+1))/2
        mid2_psi(i,j) = (wf(i+2,j) + 2*wf(i+1,j)+ wf(i,j))/4
        diff2_psi(i,j) = (wf(i+2,j) - 2*wf(i+1,j) + wf(i,j))/(dx*dx)
    end
end
//----------------Uncertainty check--------------------------------------
Un = (4*%pi)/h

for i=1:counter
    y2(:,i) = mid_psi(:,i).*diff_psi(:,i)
    y3(:,i) = mid2_psi(:,i).*diff2_psi(:,i)
    ex_p(i) = -1*%i*(h/(2*%pi))*inttrap(x2,y2(:,i))
    ex_p2(i) = -1*(h/(2*%pi))**2*inttrap(x3,y3(:,i))
    sig_x(i) = sqrt(ex_X2(i) - (ex_X(i)*ex_X(i)))
    sig_p(i) = sqrt(ex_p2(i) - (ex_p(i)*ex_p(i)))
    ex_K(i) = ex_p2(i)/(2*m*e)
    un(i) = Un*(sig_x(i) * sig_p(i))
end
//-----------------total energy------------------------
E = ex_V + ex_K
disp("Expectation value <x> =",ex_X)
disp("Expectation value <x2> =",ex_X2)
disp("Expectation value <p> =",ex_p)
disp("Expectation value <p2> =",ex_p2)
disp("Expectation value <V> =",ex_V)
disp("Expectation value <KE> =",ex_K)
disp("Total Energy for all bound state values <V>+<KE> = ",E)
disp("Standard deviation of x =",sig_x)
disp("Standard deviation of p =",sig_p)
disp("Uncertanity Product (hbar/2) = ",un)
for i=1:counter
    if abs(E(i)-Z(i))<=0.01
    disp("Energy Eigenvalues Correct for boundstate"+string(i))
    else
    disp("Energy Eigenvalues Incorrect for boundstate"+string(i))
    end
end
disp("time required : ",toc())
//------Plotting wavefunction vs xnew-----------------
//subplot(221)
//for i=1:counter
    //plot(xnew,wf(:,i))
//end
//subplot(222)
//for i=1:counter
    //plot(xnew,neigen)
//end
//N_sq = inttrap(xnew,eigenvec)
//Plotting of wavefunction------
figure
for i=1:counter
    subplot(2,2,i)
    plot(xnew,wf(:,i))
    m=gca()
    m.x_ticks=tlist(["ticks","locations","labels"],[-2.5000D-10,-1.000D-10,0.000D-10,1.000D-10,2.5000D-10],["-2.5e-10","-1e-10","0e-10","1e-10","2.5e-10"])
    title("Wavefunction plotting for bound state"+string(i),"color","black","Fontsize","5","Fontname","2")
    xlabel("X----->","color","brown","Fontsize","3")
    ylabel("Wavefunction----->","color","brown","Fontsize","3")
    xgrid()
end
figure
for i=1:counter
    subplot(2,2,i)
    plot(xnew,neigen(:,i))
    m=gca()
    m.x_ticks=tlist(["ticks","locations","labels"],[-2.5000D-10,-1.000D-10,0.000D-10,1.000D-10,2.5000D-10],["-2.5e-10","-1e-10","0e-10","1e-10","2.5e-10"])
    title("|U(x)|^2 plotting for bound state"+string(i),"color","black","Fontsize","5","Fontname","2")
    xlabel("X----->","color","brown","Fontsize","3")
    ylabel("|U(x)|^2----->","color","brown","Fontsize","3")
    xgrid()
end
//Plotting Potential as function of grid points : 
Id=ones(1,n-2)
x_plot = x(2:n-1)
figure
m=gca();
m.data_bounds = [xmin*1D-9,-150;xmax*1D-9,50]
plot(x_plot,vnew)
for i=1:counter
    plot(x_plot',Z(i)*Id,'r')
end
title("Finite Square Well Potential Plot ","color","black","Fontsize","5","Fontname","2")
xlabel("X----->","color","brown","Fontsize","3")
ylabel("Potential----->","color","brown","Fontsize","3")
legend("Potential well","Eigen values")


