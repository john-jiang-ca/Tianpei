
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
        symR=-d;
    else
        symR=d;
    end
    X(i)=complex(symR,0);
elseif(M==4)
        if(real(X(i))<=0)
        symR=-d;
    else
        symR=d;
        end
        if(imag(X(i))<=0)
        symI=-d;
    else
        symI=d;
        end
        X(i)=complex(symR, symI);
elseif(M==16)
    
        if(real(X(i))<=-2*d)
        symR=-3*d;
        elseif (real(X(i))>-2*d&&real(X(i))<=0)
        symR=-d;
        elseif(real(X(i))>0&&real(X(i))<=2*d)
            symR=d;
        else(real(X(i)>2*d));
            symR=3*d;
        end
        
        if(imag(X(i))<=-2*d)
        symI=-3*d;
        elseif (imag(X(i))>-2*d&&imag(X(i))<=0)
        symI=-d;
        elseif(imag(X(i))>0&&imag(X(i))<=2*d)
            symI=d;
        else(imag(X(i)>2*d));
            symI=3*d;
        end
        X(i)=complex(symR, symI);
end

end
X_hat=X;
end