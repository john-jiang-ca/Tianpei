function [ symOut ] = Partial_MIC_recursive( y_prev, H_prev, list_prev, W_prev, sym_prev, sym_total, SNRd, M, pav, stage, maxStage, symConstell, tol)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
%Multistage Interference Cancellation (MIC) recursive algorithm
%   the stage is calculated from 0(initial stage) to maxStage-1
%initial stage y_prev=y, H_prev=H, list_prev=list, W_prev=zeros(Nt)
%% split reliable and unreliable symbols as well as the corresponding channels from the previous stage
sym_reliable=[];
sym_unreliable=[];
H_reliable=[];
H_unreliable=[];
for count=1:length(W_prev)
    if(W_prev(count)==1)
        sym_reliable=[sym_reliable; sym_prev(count)];
        sym_total(list_prev(count))=sym_prev(count);
        H_reliable=[H_reliable, H_prev(:,count)];
    else
        sym_unreliable=[sym_unreliable, sym_prev(count)];
        H_unreliable=[H_unreliable, H_prev(:,count)];
    end
end

%% if all the estimation of previous stage is reliable stop
if(isempty(sym_unreliable))
    symOut=sym_total;
    return
end
%% update list
list_current=[];
for count=1:length(W_prev)
    if(W_prev(count)==0)
        list_current=[list_current, list_prev(count)];
    end
end
% list_current=list_prev;
%% partial parellel interference cancellation

    if(isempty(H_reliable))
        y_current=y_prev;
    else
y_current=y_prev-H_reliable*sym_reliable;
    end

%% Detection of unreliable symbols
Nt=length(H_unreliable(1,:));
I=eye(Nt);
sym_tmp_u=(I/(H_unreliable'*H_unreliable+SNRd^(-1)*I))*H_unreliable'*y_current;
[sym_tmp]=Rectangular_QAM_slicer(sym_tmp_u,M, pav);
%% check reliability
[W_current]=weightCal(sym_tmp_u, sym_tmp, tol, pav, M);
sym_current=sym_tmp;
if(stage==maxStage)
    for count=1:length(sym_current)
    sym_total(list_current(count))=sym_current(count);
    end
    symOut=sym_total;
    return
end

%% recursive to next stage
[symOut]=Partial_MIC_recursive(y_current, H_unreliable, list_current, W_current, sym_current, sym_total, SNRd, M, pav, stage+1, maxStage, symConstell, tol);
end


%% adaptive  weight function 
function [W]=weightCal(s, s_Q, tol, pav, M)
 d=sqrt(3*pav/(2*(M-1)));
 L=length(s);
 W=zeros(L,1);
 for count=1:L
     if(abs(real(s(count))-real(s_Q(count)))>tol*d||abs(imag(s(count))-imag(s_Q(count)))>tol*d)
         W(count)=0;
     else
         W(count)=1;
     end
 end

end

function [X_hat] = Rectangular_QAM_slicer(X, M, pav)    %slicer
%need to be modified

% sq10=sqrt(10);
if(M==2)
    d=sqrt(pav);
else
d=sqrt(3*pav/(2*(M-1)));
end

N=length(X);
for i=1:N
if(M==2)
    
    if(real(X(i))<=0)
        REAL=-d;
    else
        REAL=d;
    end
    X(i)=complex(REAL,0);
elseif(M==4)
        if(real(X(i))<=0)
        REAL=-d;
    else
        REAL=d;
        end
        if(imag(X(i))<=0)
        IMAG=-d;
    else
        IMAG=d;
        end
        X(i)=complex(REAL, IMAG);
elseif(M==16)
    
        if(real(X(i))<=-2*d)
        REAL=-3*d;
        elseif (real(X(i))>-2*d&&real(X(i))<=0)
        REAL=-d;
        elseif(real(X(i))>0&&real(X(i))<=2*d)
            REAL=d;
        elseif(real(X(i)>2*d))
            REAL=3*d;
        end
        if(imag(X(i))<=-2*d)
        IMAG=-3*d;
        elseif (imag(X(i))>-2*d&&imag(X(i))<=0)
        IMAG=-d;
        elseif(imag(X(i))>0&&imag(X(i))<=2*d)
            IMAG=d;
        elseif(imag(X(i)>2*d))
            IMAG=3*d;
        end
        X(i)=complex(REAL, IMAG);
end

end
X_hat=X;
end


