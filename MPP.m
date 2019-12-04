clear all;
clc;

load('one_year_data.mat');

mpc=case123_new;
[N,R,X,A,a0,r,x]=readMat(mpc);
N=N-1;
clearvars x
load data.mat;

load_base=(mpc.bus(:,3));
load_base_pc=repmat(load_base,1,size(pc_scaled,1));
load_base_pg=repmat(load_base,1,size(pg_scaled,1));

pc=pc_scaled(:,1:N+1).*load_base_pc';
pg=pg_scaled(:,1:N+1).*load_base_pg';
pg_max_scaled=max(pg);
pg_max_scaled=repmat(pg_max_scaled',1,size(pg,1));

alpha=0.7; %penetration
scale=1;

sg_max=1.1*alpha*pg_max_scaled(1:N+1,:);
qg_max=(realsqrt((sg_max.^2)-((alpha*pg').^2)));

pg=scale*(pg(:,1:N))';
pc=scale*(pc(:,1:N))';
qc=scale*(qc(:,1:N))';
qg_max=scale*(qg_max(1:N,:));
ind_zero=find(sum(qg_max')==0);
K=[1:1:N];
ind_nonzero=setdiff(K,ind_zero);
ind_zero1=ind_zero+240;
ind_zero2=ind_zero+360;

lambda=1000;
mu=100;
%epsilon=0.03;
Z=zeros(N,N);
I=eye(N);
Z2=zeros(2*N,2*N);
I2=eye(2*N);
B1=eye(N);
S=B1;
S(ind_zero,:)=[];
B1(ind_nonzero,:)=[];

eye_mod = eye(N);
eye_mod(1,:)=[];
eye_mod(:,1)=-1;

B=[B1 zeros(length(ind_zero),2*N) zeros(size(B1,1),1)];

v_lower = 0.97;
v_upper = 1.03;

td=0.5;

H=2*[(td*R+(1-td)*X'*X) Z Z (1-td)*X*ones(N,1);
    Z 0.0003*I Z zeros(N,1);
    Z Z 0.0003*I zeros(N,1);
    (1-td)*ones(1,N)*X' zeros(1,2*N) (1-td)*N];

c=[zeros(N,1);lambda*ones(N,1);lambda*ones(N,1) ;-(1-td)*2*N];

A=[-X -I Z -ones(N,1);
    X Z -I ones(N,1);
    -S zeros(length(ind_nonzero),2*N) zeros(size(S,1),1);
    S zeros(length(ind_nonzero),2*N) zeros(size(S,1),1);
    Z -I Z zeros(N,1);
    Z Z -I zeros(N,1);
    zeros(1,3*N) 1;
    -zeros(1,3*N) -1];

C= (1-td)*[2*X'*R Z Z;
    Z Z Z;
    Z Z Z;
    2*ones(1,N)*R zeros(1,2*N)];

E=[R Z Z;
    -R Z Z;
    zeros(length(ind_nonzero),N) -S zeros(length(ind_nonzero),N);
    zeros(length(ind_nonzero),N) zeros(length(ind_nonzero),N) S;
    Z Z Z;
    Z Z Z;
    zeros(1,3*N);
    zeros(1,3*N)];

b=[-v_lower*ones(N,1);
    v_upper*ones(N,1);
    zeros(2*N,1);
    zeros(2*length(ind_nonzero),1);
    1.05;
    -0.95];

f=zeros(length(ind_zero),1);

Hi = inv(H);

T=24*30*12;
V_nc=R*(alpha*pg(:,1:T)-pc(:,1:T))-X*(qc(:,1:T));
j=0;

indic_degen = [];
tic
for i=1:1:T
    %yalmip clear
    flag(i) = 0;
    pgi=pg(:,i);
    pci=pc(:,i);
    pi=alpha*pgi-pci;
    qci=qc(:,i);
    qgi_max=qg_max(1:N,i);
    teta=[pi;-qgi_max-qci;qgi_max-qci];
    if (i==1 || j==0)
        flag(i)=0;
    else
        for k=1:size(poly.p_coeff,2)
            indic.p{k}=poly.p_coeff{k}*teta - poly.p_bias{k} ;
            indic.p{k}=indic.p{k} <= 1e-4;
            sum(indic.p{k})==size(temp.A_bar{k},1)
            
            indic.d{k}=temp.G1{k}*teta + temp.w1{k};
            indic.d{k}=indic.d{k} >= -1e-4;
            sum(indic.d{k})==size(temp.G1{k},1)
            if (sum(indic.p{k})>=size(temp.A_bar{k},1) && sum(indic.d{k})>=size(temp.G1{k},1))
                x_sol(:,i)=temp.M{k}*teta+temp.r{k};
                flag(i)=k;
                break
            else
                flag(i)=0;
                continue
            end
        end
    end
    
    if (flag(i)==0)
        j=j+1;
        flag(i) = j;
        
        %tic
        [x(:,j),res(j),D(:,j),D_eq(:,j),C_eq(:,i),constraints(:,j)]=MPP_LM_CVX(teta,A,H,E,b,N,c,B,f,C);
         %time (j) = toc;
        
        %         [cons_sort,idx] = sort(abs(constraints(:,j)));
        %         D_sort = D(idx,j);
        %
        %         close all
        %         plot(cons_sort,'o')
        %         hold
        %         plot(D_sort,'p')
        %         set(gca, 'YScale', 'log')
        %         grid on
        
        
        %close all
        ind_all = (1:size(A,1))';
        
        ind_active{j} = find(abs(constraints(:,j))<=1e-7);
        ind_non_slack = ind_all;
        ind_non_slack(ind_active{j})=[];
        ind_non{j} = ind_non_slack;
        
        temp.A_tilde{j}=A(ind_active{j},:);
        temp.A_bar{j}=A(ind_non{j},:);
        temp.b_tilde{j}=b(ind_active{j});
        temp.b_bar{j}=b(ind_non{j});
        temp.E_tilde{j}=E(ind_active{j},:);
        temp.E_bar{j}=E(ind_non{j},:);
        
        temp.J{j}=[temp.A_tilde{j};B];
        
        if(rank(temp.J{j}*Hi*temp.J{j}') < size(temp.J{j}*Hi*temp.J{j}',1))
            indic_degen(end+1) = i;
            x_sol(:,i)=x(:,j);
            ind_active{j} = [];
            ind_non_slack = [];
            ind_non_slack(ind_active{j})=[];
            ind_non{j} = [];
            
            temp.A_tilde{j}=[];
            temp.A_bar{j}=[];
            temp.b_tilde{j}=[];
            temp.b_bar{j}=[];
            temp.E_tilde{j}=[];
            temp.E_bar{j}=[];
            
            temp.J{j}=[];
            j = j-1;
            continue
        end
        buffer_G=-inv(temp.J{j}*Hi*temp.J{j}')*(temp.J{j}*Hi*C+[temp.E_tilde{j};zeros(size(B,1),size(teta,1))]);
        temp.G1{j}=buffer_G(1:length(ind_active{j}),:);
        temp.G2{j}=buffer_G(length(ind_active{j})+1:end,:);
        temp.M{j}=-Hi*(C + (temp.A_tilde{j})'*temp.G1{j} + B'*temp.G2{j} );
        buffer_w=-inv(temp.J{j}*Hi*temp.J{j}')*(temp.J{j}*Hi*c + [temp.b_tilde{j};f]);
        temp.w1{j}=buffer_w(1:length(ind_active{j}));
        temp.w2{j}=buffer_w(length(ind_active{j})+1:end);
        temp.r{j}=-Hi*(c+(temp.A_tilde{j})'*temp.w1{j} + B'*temp.w2{j} );
        
        poly.p_coeff{j}=temp.A_bar{j}*temp.M{j}-temp.E_bar{j};
        poly.p_bias{j}=temp.b_bar{j}-temp.A_bar{j}*temp.r{j};
        
        
        x_sol(:,i)=x(:,j);
    end
    bound(:,i)=E*teta+b;
    sol(:,i)=A*x_sol(:,i);
    q(:,i)=x_sol(1:N,i);
end
toc

s_l=x_sol(N+1:2*N,:);
s_u=x_sol(2*N+1:3*N,:);
DV=R*(alpha*pg(:,1:T)-pc(:,1:T))+X*q;

qg=q-qc(:,1:size(q,2));



%%
% H=[R Z Z;Z 0.0003*I Z;Z Z 0.0003*I];
% c=[zeros(N,1);lambda*ones(N,1);lambda*ones(N,1)];
% A=[-X -I Z;X Z -I;
%     -S zeros(length(ind_nonzero),2*N);S zeros(length(ind_nonzero),2*N);
%     Z -I Z;Z Z -I];
%
% E=[R Z Z;-R Z Z; zeros(length(ind_nonzero),N) -S zeros(length(ind_nonzero),N)
%     zeros(length(ind_nonzero),N) zeros(length(ind_nonzero),N) S;Z Z Z;Z Z Z];
% b=[(epsilon)*ones(2*N,1);zeros(2*N,1);zeros(2*length(ind_nonzero),1)];
% f=zeros(length(ind_zero),1);

% clear all;
% clc;
%
% load('one_year_data.mat');
%
% mpc=case123_new;
% [N,R,X,A,a0,r,x]=readMat(mpc);
% N=N-1;
% clearvars x
% load_base=(mpc.bus(:,3));
% load_base_pc=repmat(load_base,1,size(pc_scaled,1));
% load_base_pg=repmat(load_base,1,size(pg_scaled,1));
% pc=pc_scaled(:,1:N+1).*load_base_pc';
% pg=pg_scaled(:,1:N+1).*load_base_pg';
% pg_max_scaled=max(pg);
% pg_max_scaled=repmat(pg_max_scaled',1,size(pg,1));
%
% pf = 0.9 + (0.95-0.9).*rand(size(pc,2),1);
% pf=repmat(pf',size(pc,1),1);
% qc=pf.*pc;
% save('data.mat','pg','pc','qc')

%[x(:,j),res(j),D(:,j),D_eq(:,j),C_eq(:,i),constraints(:,j)]=MPP_LM(teta,A,H,E,b,N,c,B,f);
% H=2*[(td*R+(1-td)*X'*X) Z Z (1-td)*X'*ones(N,1);
%      Z 0.0003*I Z zeros(N,1);
%      Z Z 0.0003*I zeros(N,1);
%      (1-td)*ones(1,N)*X zeros(1,2*N) (1-td)*1];