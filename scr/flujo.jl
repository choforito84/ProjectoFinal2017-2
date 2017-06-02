#funciones auxiliares
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
#todas las matrices tendrán la siguiente orientación: la primera entrada (renglón) representa la coordenada x y la segunda entrada (columna) la coordenada en y. 
#Aunque este es el orden usual, en una matriz este orden está transpuesto, por lo que es importante la aclaración
#Para mayor claridad, le daremos a nuestro sistema coordenado una orientación repecto a los puntos cardinales de manera usual, 
#es decir, el plano y>0 es el norte, el plano y<0 es el sur, el plano x>0 el este y x<0 el oeste


dir=Dict(""=>1,"N"=>2,"S"=>3,"E"=>4,"O"=>5,"NE"=>6,"NO"=>7,"SE"=>8,"SO"=>9)
#parámetros del flujo
dims=(30,10)
iter=5
#c=k/h
c=0.5
omega=1
U=[1,0] 
E=[[0,0],
[0,1],
[0,-1],
[1,0],
[-1,0],
[1,1],
[-1,1],
[1,-1],
[-1,-1]]
W=[4/9,
1/9,1/9,1/9,1/9,
1/36,1/36,1/36,1/36]
N=[zeros(dims) for i in 1:9]
#inicializar densidades numéricas
for i in 1:9
    N[i]=W[i]*(ones(dims) + 3/c*dot(E[i],U) + 9/2*(dot(E[i],U)/c)^2 -3/(2*c)*norm(U)^2)
end
#inicializar barreras
Barr=falses(dims)
#Barrera de prueba lateral
Barr[convert(Int64,floor(dims[1]/4)), convert(Int64,floor(dims[2]/3)) : convert(Int64,floor(dims[2]*2/3)) ]=true
#Localizamos puntos vecinos a las barreras para cada dirección
B=[falses(dims) for i in 1:9]
for i in 2:9
    B[i]=girar(Barr,E[i][1],axis=1)
    B[i]=girar(B[i],E[i][2],axis=2)
end
function mover()
    #movemos las partículas en cada dirección hacia su dirección dada
    for i in 1:9
        N[i]=girar(N[i],E[i][1],axis=1)
        N[i]=girar(N[i],E[i][2],axis=2)
    end
    #Añadimos a densidad de los vecinos a la barrera las partículas que colisionaron con esta barrera
    N[dir["N"]][B[dir["N"]]]=N[dir["S"]][Barr]
    N[dir["S"]][B[dir["S"]]]=N[dir["N"]][Barr]
    N[dir["E"]][B[dir["E"]]]=N[dir["O"]][Barr]
    N[dir["O"]][B[dir["O"]]]=N[dir["E"]][Barr]
    
    N[dir["NE"]][B[dir["NE"]]]=N[dir["SO"]][Barr]
    N[dir["NO"]][B[dir["NO"]]]=N[dir["SE"]][Barr]
    N[dir["SE"]][B[dir["SE"]]]=N[dir["NO"]][Barr]
    N[dir["SO"]][B[dir["SO"]]]=N[dir["NE"]][Barr]
end

function colisionar()
    global rho, ux, uy
    NEQ=[zeros(dims) for i in 1:9]
    rho=sum([N[i] for i in 1:9])
    ux=[E[j][1]*N[j] for j in 1:9]
    ux=sum(ux)./rho
    uy=[E[j][2]*N[j] for j in 1:9]
    uy=sum(uy)./rho
    #println(uy)
    for i in 1:9
        NEQ[i]=rho.*W[i].*(1+3*(E[i][1]*ux+E[i][2]*uy)+9/2*(E[i][1]*ux+E[i][2]*uy).^2 -3/2*(ux.^2+uy.^2))
        N[i]=N[i]+omega*(NEQ[i]-N[i])
    end
end
colisionar()
for t in 1:iter
    mover()
    colisionar()
end
using Plots
gr()
heatmap(rho)