clc;
clear all;

N = 64;
M = 10;
B = 10;

Totally=20;

for Count_Number = 1:Totally
    
G = gsp_sensor(N);
G = gsp_compute_fourier_basis(G);

Lumbda = G.e;
U = double(full(G.U));

% Sampling on Graphs of Bandlimited graph signal with Vertex Domain Sampling %
[UnSample_Vertex_Index,Sample_Vertex_Index]=crossvalind('LeaveMOut',N,M);
Sample_Vertex_Location = find(Sample_Vertex_Index~=0);
UnSample_Vertex_Location = find(UnSample_Vertex_Index~=0);

% Random Gaussian Values %
d=1+1*randn(M,1);  

% Define GBL Signal--Bandwidth is B %
U_VB = U(:,1:B);
F_GBL = U_VB*d; 

% Define the Non-GBL Signal %
Alpha=zeros(N,1);

for i=1:N
    Alpha(i,1) = 1-2*Lumbda(i,1)/G.lmax;
end

A_NonGBL=zeros(N,M);
for i=1:M
    A_NonGBL(:,i) = U*Alpha;  
end

F_NonGBL = A_NonGBL*d;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Vertex Sampling and Recovery of GBL Signal%
F_Recovery_GBL_Vertex = U(:,1:B)*pinv(U(Sample_Vertex_Location,1:B))*F_GBL(Sample_Vertex_Location);

% Spectral Sampling and Recovery of GBL Signal% 
F_Recovery_GBL_Spectral= U*diag([ones(B,1);zeros(N-B,1)])*U'*F_GBL;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Vertex Sampling and Recovery of Non-GBL Signal%
Sample_Matrix = zeros(N,M);
for i=1:M
    Sample_Matrix(Sample_Vertex_Location(i,1),i)=1;
end

F_Recovery_NonGBL_Vertex = A_NonGBL*pinv(Sample_Matrix'*A_NonGBL)*(Sample_Matrix'*F_NonGBL);

% Spectral Sampling and Recovery of Non-GBL Signal%
G_function=zeros(N,1);
for i=1:N
    G_function(i,1)=exp(-0.5*Lumbda(i,1));
end
S_T = Sample_Matrix'*diag(G_function)*U';

F_Recovery_NonGBL_Spectral = A_NonGBL*pinv(S_T*A_NonGBL)*S_T*F_NonGBL;


% Compute the MSE Number of Different Method of Graph Sampling %

MSE_GBL_Vertex(Count_Number,1) = immse(F_Recovery_GBL_Vertex,F_GBL);
MSE_GBL_Spectral(Count_Number,1) = immse(F_Recovery_GBL_Spectral,F_GBL);

MSE_NonGBL_Vertex(Count_Number,1) = immse( F_Recovery_NonGBL_Vertex,F_NonGBL);
MSE_NonGBL_Spectral(Count_Number,1) = immse( F_Recovery_NonGBL_Spectral,F_NonGBL);

end

Totally_Result(1,1) = mean(MSE_GBL_Vertex);
Totally_Result(1,2) = mean(MSE_GBL_Spectral);
Totally_Result(1,3) = mean(MSE_NonGBL_Vertex);
Totally_Result(1,4) = mean(MSE_NonGBL_Spectral);

disp(Totally_Result);