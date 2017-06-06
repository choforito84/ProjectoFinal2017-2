#funciones auxiliares
module auxiliares
export girar,girarones,girarzeros,rotacional,circulo
function girar(A, n::Int64; axis::Int64=1)
    if n<=0
        n=abs(n)
        if axis==1
            B=copy(A[(n+1):end,:])
            B=vcat(B,A[1:n,:])
            return B
        elseif axis==2
            B=copy(A[:,(n+1):end])
            B=hcat(B,A[:,1:n])
            return B
        end
    elseif n>=0
        if axis==1
            B=copy(A[1:(end-n),:])
            B=vcat(A[(end-n+1):end,:],B)
            return B
        elseif axis==2
            B=copy(A[:,1:(end-n)])
            B=hcat(A[:,(end-n+1):end],B)
            return B
        end
    end
end
function girarzeros(A, n::Int64; axis::Int64=1)
    if n<=0
        n=abs(n)
        if axis==1
            B=zeros(size(A[(n+1):end,:]))
            B=vcat(B,A[1:n,:])
            return B
        elseif axis==2
            B=zeros(size(A[:,(n+1):end]))
            B=hcat(B,A[:,1:n])
            return B
        end
    elseif n>=0
        if axis==1
            B=zeros(size(A[1:(end-n),:]))
            B=vcat(A[(end-n+1):end,:],B)
            return B
        elseif axis==2
            B=zeros(size(A[:,1:(end-n)]))
            B=hcat(A[:,(end-n+1):end],B)
            return B
        end
    end
end
function girarones(A, n::Int64; axis::Int64=1)
    if n<=0
        n=abs(n)
        if axis==1
            B=ones(size(A[(n+1):end,:]))
            B=vcat(B,A[1:n,:])
            return B
        elseif axis==2
            B=ones(size(A[:,(n+1):end]))
            B=hcat(B,A[:,1:n])
            return B
        end
    elseif n>=0
        if axis==1
            B=ones(size(A[1:(end-n),:]))
            B=vcat(A[(end-n+1):end,:],B)
            return B
        elseif axis==2
            B=ones(size(A[:,1:(end-n)]))
            B=hcat(A[:,(end-n+1):end],B)
            return B
        end
    end
end
function rotacional(ux,uy;h=1)
    if size(ux)!=size(uy)
        error("El tamaño de las matrices debe coincidir")
    end
    resultado=girar(uy,1,axis=1)-girar(uy,-1,axis=1) -(girar(ux,1,axis=2)-girar(ux,-1,axis=2))
    resultado=resultado[2:(end-1),2:(end-1)]
    return resultado
end

function circulo(A;r::Float64=1.0,C=[0.0,0.0])
    if abs(norm(A-C))<=r
        return true
    else
        return false
    end
end 
end