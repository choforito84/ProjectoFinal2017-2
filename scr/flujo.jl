#funciones auxiliares
#module fluir
using Plots
s=Dates.time()
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

#todas las matrices tendrán la siguiente orientación: la primera entrada (renglón) representa la coordenada x y la segunda entrada (columna) la coordenada en y. 
#Aunque este es el orden usual, en una matriz este orden está transpuesto, por lo que es importante la aclaración
#Para mayor claridad, le daremos a nuestro sistema coordenado una orientación repecto a los puntos cardinales de manera usual, 
#es decir, el plano y>0 es el norte, el plano y<0 es el sur, el plano x>0 el este y x<0 el oeste


dir=Dict(""=>1,"N"=>2,"S"=>3,"E"=>4,"O"=>5,"NE"=>6,"NO"=>7,"SE"=>8,"SO"=>9)
#parámetros del flujo
dims=(300,100)
#c=k/h
#obtención de parámetros del flujo mediante el argument parser. Restringimos los valores para obtener un flujo estable
if length(ARGS)>=1
    omega=ARGS[1]
elseif length(ARGS)>=2
    vel=ARGS[2]
elseif length(ARGS)>=3
    barrera=ARGS[3]
else
    omega=0.04
    vel=0.12
    barrera="o"
end
if omega>0.3 || omega<0.005
    error("Flujo inestable para viscosidad cinemática ν = $(omega)")    
end
if omega>0.2 || vel<0.005
    error("Flujo inestable para velocidad de tunel v= $(vel)")    
end
c=1
omega=1/(3*omega+0.5)
U=vel*[1,0] 
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
#Escoger geometría de la barerra
if barrera=="|"
    Barr[convert(Int64,floor(dims[1]/4)), convert(Int64,floor(dims[2]/3)) : convert(Int64,floor(dims[2]*2/3)) ]=true
elseif barrera=="o"
     Barr=[circulo([x,y],r=5.0,C=[0.0,0.0]) for x in linspace(-20,40,dims[1]),y in linspace(-10,10,dims[2])]
end

#Localizamos puntos vecinos a las barreras para cada dirección
B=[falses(dims) for i in 1:9]
for i in 2:9
    B[i]=girar(Barr,E[i][1],axis=1)
    B[i]=girar(B[i],E[i][2],axis=2)
end


function mover(N)
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
    return N
end
function colisionar(N)
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
    #Forzamiento debido a condiciones de frontera
 	N[dir["N"]][1,:] = 1/9 * (1 + 3*norm(U) + 4.5*norm(U)^2 - 1.5*norm(U)^2)
	N[dir["S"]][1,:] = 1/9 * (1 - 3*norm(U) + 4.5*norm(U)^2 - 1.5*norm(U)^2)
	N[dir["NE"]][1,:] = 1/36 * (1 + 3*norm(U) + 4.5*norm(U)^2 - 1.5*norm(U)^2)
	N[dir["SE"]][1,:] = 1/36 * (1 + 3*norm(U) + 4.5*norm(U)^2 - 1.5*norm(U)^2)
	N[dir["NO"]][1,:] = 1/36 * (1 - 3*norm(U) + 4.5*norm(U)^2 - 1.5*norm(U)^2)
	N[dir["SO"]][1,:] = 1/36 * (1 - 3*norm(U) + 4.5*norm(U)^2 - 1.5*norm(U)^2)

    return (N,rho,ux,uy)
end
#construir directorio para la salida y crear log
dest=Dates.format(Dates.now(),"HH;MM;SS dd-mm-Y")

mkdir("../gifs/$(dest)")

write("../gifs/$(dest)/log.txt")
println(dest)
open("../gifs/$(dest)/log.txt","w") do f
        write(f,"Fluid simulation with Lattice-Boltzmann. Date: $(Dates.format(Dates.now(),"HH:MM:SS dd-mm-Y")) \n \n")
        write(f,"Viscosity factor = $omega \n")
        write(f,"Fluid speed = $U \n")
    end

iter=1200
println("Making simulation \n")

"""
anim= @animate for t in 1:iter 
    N=mover(N)
    A=colisionar(N)
    global N=A[1]
    global rho=A[2]
    global ux=A[3]
    global uy=A[4]
    heatmap(rho',zlim=(0.8,1.2))
end
"""
F=[]
push!(F,(N,zeros(dims),zeros(dims),zeros(dims)))
for i in 1:iter
    N=mover(F[i][1])
    push!(F,colisionar(N))
end

println("Building animations \n")
anim= @animate for t in 2:10:iter 
    heatmap((F[t][2])',clim=(0.8,1.2),color=:curl,aspect_ratio=:equal)
end
gif(anim,"../gifs/$(dest)/densidad.gif",fps=30)

anim= @animate for t in 2:10:iter 
    heatmap((rotacional(F[t][3],F[t][4]))',clim=(-0.25,0.25),color=:curl,aspect_ratio=:equal)
end
gif(anim,"../gifs/$(dest)/rot.gif",fps=30)

anim= @animate for t in 2:10:iter
    heatmap(sqrt(F[t][3].^2 + F[t][4].^2)',clim=(0.0,vel*1.5),color=:curl,aspect_ratio=:equal)
end
gif(anim,"../gifs/$(dest)/vel.gif",fps=30)
anim= @animate for t in 2:10:iter  
    heatmap((F[t][3])',clim=(-vel,vel),color=:curl,aspect_ratio=:equal)
end
gif(anim,"../gifs/$(dest)/velx.gif",fps=30)
anim= @animate for t in 2:10:iter 
    heatmap((F[t][4])',clim=(-vel,vel),color=:curl,aspect_ratio=:equal)
end
gif(anim,"../gifs/$(dest)/vely.gif",fps=30)

println("Time spent creating animations: $(round(Dates.time()-s,2)) segs")


#export girar, mover, colisionar, dir, dims, c, omega, U, E, W, N, Barr, B
#end

