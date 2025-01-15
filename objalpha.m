function ff=objalpha(x)
global K  w  P_c B zeta Htill
NN=x(1:K)';
P=x(K+1:2*K)';
ff=(B*sum(w.*NN.*log2(1+Htill.*P./NN)))/(P_c+zeta*sum(P));
ff=-ff;
end