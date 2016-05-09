%%% "Volatility and Correlation-Based Systemic Risk Measures in the US Market"
%%% do file for Granger-causality tests

%Index Cohesion Force

w=30:15:300; %windows used

ICF_base=zeros(4735,19); %initiating the output vectors
F_ra_ICF=zeros(19,2);
F_rd_ICF=zeros(19,2);
F_va_ICF=zeros(19,2);
F_vd_ICF=zeros(19,2);
c_v=zeros(19,3);


DV=SP500(max(w)+1:5035,1); %selecting the dependent variables
DV(:,2)=VIX(max(w)+1:5035,1);

for i = 1:length(w)
    
    [t,m]=size(SP500stocks);
    temp=zeros(t-w(i),1);

    for j=1:t-w(i)
    
        correlmat=corr(SP500stocks(j:j+w(i),:)); %build correlation matrix
        pcorrmat = partialcorr(SP500stocks(j:j+w(i),:),SP500(j:j+w(i),1)); %build p.correlation matrix
        temp(j,1)=sum(sum(correlmat-eye(m)))/sum(sum(pcorrmat-eye(m))); %create the ICF
    
    end
    
    ICF_base(:,i)=temp(max(w)-w(i)+1:length(temp),1); %setting the indicators to the same time period
    
    [F_ra_ICF(i,1),c_v(i,1),c_v(i,2),c_v(i,3)] = granger_cause(DV(1:1963,1),ICF_base(1:1963,i),20);
    F_ra_ICF(i,2)=(F_ra_ICF(i,1)>c_v(i,1))+(F_ra_ICF(i,1)>c_v(i,2))+(F_ra_ICF(i,1)>c_v(i,3));
    
    [F_rd_ICF(i,1),c_v(i,1),c_v(i,2),c_v(i,3)] = granger_cause(DV(1964:4735,1),ICF_base(1964:4735,i),20);
    F_rd_ICF(i,2)=(F_rd_ICF(i,1)>c_v(i,1))+(F_rd_ICF(i,1)>c_v(i,2))+(F_rd_ICF(i,1)>c_v(i,3));

    [F_va_ICF(i,1),c_v(i,1),c_v(i,2),c_v(i,3)] = granger_cause(DV(1:1963,2),ICF_base(1:1963,i),20);
    F_va_ICF(i,2)=(F_va_ICF(i,1)>c_v(i,1))+(F_va_ICF(i,1)>c_v(i,2))+(F_va_ICF(i,1)>c_v(i,3));
    
    [F_vd_ICF(i,1),c_v(i,1),c_v(i,2),c_v(i,3)] = granger_cause(DV(1964:4735,2),ICF_base(1964:4735,i),20);
    F_vd_ICF(i,2)=(F_vd_ICF(i,1)>c_v(i,1))+(F_vd_ICF(i,1)>c_v(i,2))+(F_vd_ICF(i,1)>c_v(i,3));
    
end
clear temp t pcorrmat correlmat c_v i j m



%%%% Absorption Ratio

w=500;

AR_base=zeros(4535,1);

F_ra_AR=zeros(1,2);
F_rd_AR=zeros(1,2);
F_va_AR=zeros(1,2);
F_vd_AR=zeros(1,2);

F_ra_DR=zeros(1,2);
F_rd_DR=zeros(1,2);
F_va_DR=zeros(1,2);
F_vd_DR=zeros(1,2);

c_v=zeros(1,3);


DV=SP500(max(w)+1:5035,1);
DV(:,2)=VIX(max(w)+1:5035,1);

for i = 1:length(w)
    
    [t,m]=size(SP500stocks);
    temp=zeros(t-w(i),1);

    for j=1:t-w(i)
        eigenmat=sort(eig(cov(SP500stocks(j:j+w(i),:))));
        eigenmat=eigenmat./sum(eigenmat);
        temp(j,1)=sum(eigenmat(round(m*4/5):m,1));
    end
    
    AR_base(:,1)=temp(max(w)-w(i)+1:length(temp),1);
    Delta_AR=absorption_ti(temp);
    
    [F_ra_AR(i,1),c_v(i,1),c_v(i,2),c_v(i,3)] = granger_cause(DV(1:1963,1),AR_base(1:1963,i),20);
    F_ra_AR(i,2)=(F_ra_AR(i,1)>c_v(i,1))+(F_ra_AR(i,1)>c_v(i,2))+(F_ra_AR(i,1)>c_v(i,3));
    
    [F_ra_DR(i,1),c_v(i,1),c_v(i,2),c_v(i,3)] = granger_cause(DV(253:2215,1),Delta_AR(1:1963,1),20);
    F_ra_DR(i,2)=(F_ra_DR(i,1)>c_v(i,1))+(F_ra_DR(i,1)>c_v(i,2))+(F_ra_DR(i,1)>c_v(i,3));
    
    [F_rd_AR(i,1),c_v(i,1),c_v(i,2),c_v(i,3)] = granger_cause(DV(1964:4535,1),AR_base(1964:4535,i),20);
    F_rd_AR(i,2)=(F_rd_AR(i,1)>c_v(i,1))+(F_rd_AR(i,1)>c_v(i,2))+(F_rd_AR(i,1)>c_v(i,3));

    [F_rd_DR(i,1),c_v(i,1),c_v(i,2),c_v(i,3)] = granger_cause(DV(2216:4535,1),Delta_AR(1964:4283,1),20);
    F_rd_DR(i,2)=(F_rd_DR(i,1)>c_v(i,1))+(F_rd_DR(i,1)>c_v(i,2))+(F_rd_DR(i,1)>c_v(i,3));
    
    [F_va_AR(i,1),c_v(i,1),c_v(i,2),c_v(i,3)] = granger_cause(DV(1:1963,2),AR_base(1:1963,i),20);
    F_va_AR(i,2)=(F_va_AR(i,1)>c_v(i,1))+(F_va_AR(i,1)>c_v(i,2))+(F_va_AR(i,1)>c_v(i,3));

    [F_va_DR(i,1),c_v(i,1),c_v(i,2),c_v(i,3)] = granger_cause(DV(1:1963,2),Delta_AR(1:1963,1),20);
    F_va_DR(i,2)=(F_va_DR(i,1)>c_v(i,1))+(F_va_DR(i,1)>c_v(i,2))+(F_va_DR(i,1)>c_v(i,3));
    
    [F_vd_AR(i,1),c_v(i,1),c_v(i,2),c_v(i,3)] = granger_cause(DV(1964:4535,2),AR_base(1964:4535,i),20);
    F_vd_AR(i,2)=(F_vd_AR(i,1)>c_v(i,1))+(F_vd_AR(i,1)>c_v(i,2))+(F_vd_AR(i,1)>c_v(i,3));

    [F_vd_DR(i,1),c_v(i,1),c_v(i,2),c_v(i,3)] = granger_cause(DV(2216:4535,2),Delta_AR(1964:4283,1),20);
    F_vd_DR(i,2)=(F_vd_DR(i,1)>c_v(i,1))+(F_vd_DR(i,1)>c_v(i,2))+(F_vd_DR(i,1)>c_v(i,3));
    
end
clear temp t eigenmat c_v i j m AR_base

%Eigenvalue Entropy

w=30:15:300;

EE_base=zeros(4735,19);
DE_base=zeros(4734,19);

F_ra_EE=zeros(19,2);
F_rd_EE=zeros(19,2);
F_va_EE=zeros(19,2);
F_vd_EE=zeros(19,2);

F_ra_DE=zeros(19,2);
F_rd_DE=zeros(19,2);
F_va_DE=zeros(19,2);
F_vd_DE=zeros(19,2);

c_v=zeros(19,3);


DV=SP500(max(w)+1:5035,1);
DV(:,2)=VIX(max(w)+1:5035,1);

for i = 1:length(w)
    
    [t,~]=size(SP500stocks);
    temp=zeros(t-w(i),1);

    for j=1:t-w(i)
    
        eigenmat=eig(corr(SP500stocks(j:j+w(i),:)));
        eigenmat=eigenmat/sum(eigenmat);
        temp(j,1)=-sum(eigenmat.*log(eigenmat));
        temp2=temp(2:length(temp))-temp(1:length(temp)-1);
    
    end
    
    EE_base(:,i)=temp(max(w)-w(i)+1:length(temp),1);
    DE_base(:,i)=temp2(max(w)-w(i)+1:length(temp2),1);
    
    [F_ra_EE(i,1),c_v(i,1),c_v(i,2),c_v(i,3)] = granger_cause(DV(1:1963,1),EE_base(1:1963,i),20);
    F_ra_EE(i,2)=(F_ra_EE(i,1)>c_v(i,1))+(F_ra_EE(i,1)>c_v(i,2))+(F_ra_EE(i,1)>c_v(i,3));
    [F_ra_DE(i,1),c_v(i,1),c_v(i,2),c_v(i,3)] = granger_cause(DV(1:1963,1),DE_base(1:1963,i),20);
    F_ra_DE(i,2)=(F_ra_DE(i,1)>c_v(i,1))+(F_ra_DE(i,1)>c_v(i,2))+(F_ra_DE(i,1)>c_v(i,3));
    
    [F_rd_EE(i,1),c_v(i,1),c_v(i,2),c_v(i,3)] = granger_cause(DV(1964:4735,1),EE_base(1964:4735,i),20);
    F_rd_EE(i,2)=(F_rd_EE(i,1)>c_v(i,1))+(F_rd_EE(i,1)>c_v(i,2))+(F_rd_EE(i,1)>c_v(i,3));
    [F_rd_DE(i,1),c_v(i,1),c_v(i,2),c_v(i,3)] = granger_cause(DV(1964:4734,1),DE_base(1964:4734,i),20);
    F_rd_DE(i,2)=(F_rd_DE(i,1)>c_v(i,1))+(F_rd_DE(i,1)>c_v(i,2))+(F_rd_DE(i,1)>c_v(i,3));

    [F_va_EE(i,1),c_v(i,1),c_v(i,2),c_v(i,3)] = granger_cause(DV(1:1963,2),EE_base(1:1963,i),20);
    F_va_EE(i,2)=(F_va_EE(i,1)>c_v(i,1))+(F_va_EE(i,1)>c_v(i,2))+(F_va_EE(i,1)>c_v(i,3));
    [F_va_DE(i,1),c_v(i,1),c_v(i,2),c_v(i,3)] = granger_cause(DV(1:1963,2),DE_base(1:1963,i),20);
    F_va_DE(i,2)=(F_va_DE(i,1)>c_v(i,1))+(F_va_DE(i,1)>c_v(i,2))+(F_va_DE(i,1)>c_v(i,3));
    
    [F_vd_EE(i,1),c_v(i,1),c_v(i,2),c_v(i,3)] = granger_cause(DV(1964:4735,2),EE_base(1964:4735,i),20);
    F_vd_EE(i,2)=(F_vd_EE(i,1)>c_v(i,1))+(F_vd_EE(i,1)>c_v(i,2))+(F_vd_EE(i,1)>c_v(i,3));
    [F_vd_DE(i,1),c_v(i,1),c_v(i,2),c_v(i,3)] = granger_cause(DV(1964:4734,2),DE_base(1964:4734,i),20);
    F_vd_DE(i,2)=(F_vd_DE(i,1)>c_v(i,1))+(F_vd_DE(i,1)>c_v(i,2))+(F_vd_DE(i,1)>c_v(i,3));

end
clear temp temp2 t pcorrmat correlmat c_v i j m

corr_base=zeros(3,3,3);
corr_base(1:3,1:3,1)=corr([ICF_base(453:4735,3),Delta_AR(1:4283,1),real(EE_base(453:4735,3))]);
corr_base(1:3,1:3,2)=corr([ICF_base(453:4735,9),Delta_AR(1:4283,1),real(EE_base(453:4735,9))]);
corr_base(1:3,1:3,3)=corr([ICF_base(453:4735,19),Delta_AR(1:4283,1),real(EE_base(453:4735,19))]);

corr_past=zeros(3,3,3);
corr_past(1:3,1:3,1)=corr([ICF_base(2403:4735,3),Delta_AR(1951:4283,1),real(EE_base(2403:4735,3))]);
corr_past(1:3,1:3,2)=corr([ICF_base(2403:4735,9),Delta_AR(1951:4283,1),real(EE_base(2403:4735,9))]);
corr_past(1:3,1:3,3)=corr([ICF_base(2403:4735,19),Delta_AR(1951:4283,1),real(EE_base(2403:4735,19))]);