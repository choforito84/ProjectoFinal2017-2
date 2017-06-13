include("./aux.jl")
using auxiliares
s=Dates.time()

#todas las matrices tendrán la siguiente orientación: la primera entrada (renglón) representa la coordenada x y la segunda entrada (columna) la coordenada en y. 
#Aunque este es el orden usual, en una matriz este orden está transpuesto, por lo que es importante la aclaración
#Para mayor claridad, le daremos a nuestro sistema coordenado una orientación repecto a los puntos cardinales de manera usual, 
#es decir, el plano y>0 es el norte, el plano y<0 es el sur, el plano x>0 el este y x<0 el oeste


dir=Dict(""=>1,"N"=>2,"S"=>3,"E"=>4,"O"=>5,"NE"=>6,"NO"=>7,"SE"=>8,"SO"=>9)
#parámetros del flujo
dims=(300,100)
#obtención de parámetros del flujo mediante el argument parser. Restringimos los valores para obtener un flujo estable
if length(ARGS)>=1
    omega=parse(Float64,ARGS[1])
else 
    omega=0.02
end
if length(ARGS)>=2
    vel=parse(Float64,ARGS[2])
else
    vel=0.1
end
if length(ARGS)>=3
    barrera=ARGS[3]
else
    barrera="|"
end
if length(ARGS)>=4
    geometria=ARGS[4]
else
    geometria="normal"
end
#arrojar errores para valores fuera de parámetros normales
if omega>0.3 || omega<0.005
    error("Flujo inestable para viscosidad cinemática ν = $(omega)")    
end
if vel>0.2 || vel<0.005
    error("Flujo inestable para velocidad de tunel v= $(vel)")    
end
if barrera!="|" && barrera!="o"
    error("$(barrera) no es una barrera aceptable")
end

if geometria!="normal" && geometria!="anillo" && geometria!="toro"
    error("$(geometria) no es una geometria aceptable")
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

#inicializamos particulas virtuales para rastrear su movimiento
particulasx=linspace(0,dims[1],12)[2:(end-1)]
particulasy=linspace(0,dims[2],12)[2:(end-1)]
particulas=[ [px,py] for px in particulasx, py in particulasy]

#distintas definiciones de la funcion de movimiento, dependiendo de la geometría que le queramos dar al problema


if geometria=="normal"
    function mover(N)
        #movemos las partículas en cada dirección hacia su dirección dada
        for i in 1:9
            N[i]=girar(N[i],E[i][1],axis=1)
            N[i]=girar(N[i],E[i][2],axis=2)
        end
        #insertando particulas por el lado izquierdo
        N[dir["E"]][1,:]=1/9*(1+3*norm(U)+3*norm(U)^2)
        N[dir["SE"]][1,:]=1/36*(1+3*norm(U)+3*norm(U)^2)
        N[dir["NE"]][1,:]=1/36*(1+3*norm(U)+3*norm(U)^2)
        #insertando particulas por el lado derecho
        N[dir["O"]][end,:]=1/9*(1+3*norm(U)+3*norm(U)^2)
        N[dir["SO"]][end,:]=1/36*(1+3*norm(U)+3*norm(U)^2)
        N[dir["NO"]][end,:]=1/36*(1+3*norm(U)+3*norm(U)^2)

        #Asignando valores fijos para la frontera superior e inferior como condiciones de frontera
        N[dir["E"]][:,1] = 1/9 * (1 + 3*norm(U) + 4.5*norm(U)^2 - 1.5*norm(U)^2)
        N[dir["O"]][:,1] = 1/9 * (1 - 3*norm(U) + 4.5*norm(U)^2 - 1.5*norm(U)^2)
        N[dir["NE"]][:,1] = 1/36 * (1 + 3*norm(U) + 4.5*norm(U)^2 - 1.5*norm(U)^2)
        N[dir["SE"]][:,1] = 1/36 * (1 + 3*norm(U) + 4.5*norm(U)^2 - 1.5*norm(U)^2)
        N[dir["NO"]][:,1] = 1/36 * (1 - 3*norm(U) + 4.5*norm(U)^2 - 1.5*norm(U)^2)
        N[dir["SO"]][:,1] = 1/36 * (1 - 3*norm(U) + 4.5*norm(U)^2 - 1.5*norm(U)^2)
        N[dir["E"]][:,end] = 1/9 * (1 + 3*norm(U) + 4.5*norm(U)^2 - 1.5*norm(U)^2)
        N[dir["O"]][:,end] = 1/9 * (1 - 3*norm(U) + 4.5*norm(U)^2 - 1.5*norm(U)^2)
        N[dir["NE"]][:,end] = 1/36 * (1 + 3*norm(U) + 4.5*norm(U)^2 - 1.5*norm(U)^2)
        N[dir["SE"]][:,end] = 1/36 * (1 + 3*norm(U) + 4.5*norm(U)^2 - 1.5*norm(U)^2)
        N[dir["NO"]][:,end] = 1/36 * (1 - 3*norm(U) + 4.5*norm(U)^2 - 1.5*norm(U)^2)
        N[dir["SO"]][:,end] = 1/36 * (1 - 3*norm(U) + 4.5*norm(U)^2 - 1.5*norm(U)^2)

        #Añadimos a densidad de los vecinos a la barrera las partículas que colisionaron con esta barrera
        N[dir["N"]][B[dir["N"]]]=N[dir["S"]][Barr]
        N[dir["S"]][B[dir["S"]]]=N[dir["N"]][Barr]
        N[dir["E"]][B[dir["E"]]]=N[dir["O"]][Barr]
        N[dir["O"]][B[dir["O"]]]=N[dir["E"]][Barr]
        
        N[dir["NE"]][B[dir["NE"]]]=N[dir["SO"]][Barr]
        N[dir["NO"]][B[dir["NO"]]]=N[dir["SE"]][Barr]
        N[dir["SE"]][B[dir["SE"]]]=N[dir["NO"]][Barr]
        N[dir["SO"]][B[dir["SO"]]]=N[dir["NE"]][Barr]
        
        #ajustando valores en la frontera izquierda y derecha
        N[dir["O"]][end,:]= N[dir["O"]][end-1,:]
        N[dir["SO"]][end,:]=N[dir["SO"]][end-1,:]
        N[dir["NO"]][end,:]=N[dir["NO"]][end-1,:]

        return N
    end


elseif geometria=="anillo"
    function mover(N)
        #movemos las partículas en cada dirección hacia su dirección dada
        for i in 1:9
            N[i]=girar(N[i],E[i][1],axis=1)
            N[i]=girar(N[i],E[i][2],axis=2)
        end

        #Asignando valores fijos para la frontera superior e inferior como condiciones de frontera
        N[dir["E"]][:,1] = 1/9 * (1 + 3*norm(U) + 4.5*norm(U)^2 - 1.5*norm(U)^2)
        N[dir["O"]][:,1] = 1/9 * (1 - 3*norm(U) + 4.5*norm(U)^2 - 1.5*norm(U)^2)
        N[dir["NE"]][:,1] = 1/36 * (1 + 3*norm(U) + 4.5*norm(U)^2 - 1.5*norm(U)^2)
        N[dir["SE"]][:,1] = 1/36 * (1 + 3*norm(U) + 4.5*norm(U)^2 - 1.5*norm(U)^2)
        N[dir["NO"]][:,1] = 1/36 * (1 - 3*norm(U) + 4.5*norm(U)^2 - 1.5*norm(U)^2)
        N[dir["SO"]][:,1] = 1/36 * (1 - 3*norm(U) + 4.5*norm(U)^2 - 1.5*norm(U)^2)
        N[dir["E"]][:,end] = 1/9 * (1 + 3*norm(U) + 4.5*norm(U)^2 - 1.5*norm(U)^2)
        N[dir["O"]][:,end] = 1/9 * (1 - 3*norm(U) + 4.5*norm(U)^2 - 1.5*norm(U)^2)
        N[dir["NE"]][:,end] = 1/36 * (1 + 3*norm(U) + 4.5*norm(U)^2 - 1.5*norm(U)^2)
        N[dir["SE"]][:,end] = 1/36 * (1 + 3*norm(U) + 4.5*norm(U)^2 - 1.5*norm(U)^2)
        N[dir["NO"]][:,end] = 1/36 * (1 - 3*norm(U) + 4.5*norm(U)^2 - 1.5*norm(U)^2)
        N[dir["SO"]][:,end] = 1/36 * (1 - 3*norm(U) + 4.5*norm(U)^2 - 1.5*norm(U)^2)

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


elseif geometria=="toro"
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
    return (N,rho,ux,uy)
end
#construir directorio para la salida y crear log
dest=Dates.format(Dates.now(),"HH;MM;SS dd-mm-Y")
mkdir("../outputs/$(dest)")
write("../outputs/$(dest)/log.txt")
open("../outputs/$(dest)/log.txt","w") do f
        write(f,"Simulación de flujo unidireccional con algoritmo de Lattice-Boltzmann. Fecha: $(Dates.format(Dates.now(),"HH:MM:SS dd-mm-Y")) \n \n")
        write(f,"Viscosidad cinemática = $(omega) \n")
        write(f,"Velocidad del fluido = $U \n")
        write(f,"Barrera: $(barrera)\n")
        write(f,"Geometría de la superficie = $geometria \n")
    end




#iteración importante para la simulación
iter=7200
println("Realizando simulación \n")


F=[]
P=[]
push!(F,(N,zeros(dims),zeros(dims),zeros(dims)))
push!(P,particulas)
for k in 1:iter
    for i in 1:10
        for j in 1:10
            x=convert(Int64,floor(particulas[i,j][1]))
            y=convert(Int64,floor(particulas[i,j][2]))
            particulas[i,j]+= [F[k][3][x,y],F[k][4][x,y]]
            if particulas[i,j][1]>=299
                particulas[i,j][1]=2
            end
            if particulas[i,j][1]<=2
                particulas[i,j][1]=2
            end
            if particulas[i,j][2]>=99
                particulas[i,j][2]=99
            end
            if particulas[i,j][2]<=2
                particulas[i,j][2]=2
            end
        end
    end
    push!(P,copy(particulas))
    N=mover(F[k][1])
    push!(F,colisionar(N))
end




println("Construyendo animaciones \n")

using Plots

anim= @animate for t in 2:30:iter
    heatmap((F[t][2])',clim=(0.8,1.2),color=:curl,aspect_ratio=:equal,xlabel="x",ylabel="y",title="densidad")
end
gif(anim,"../outputs/$(dest)/densidad.gif",fps=30)

anim= @animate for t in 2:30:iter
 
    heatmap((rotacional(F[t][3],F[t][4]))',clim=(-0.25,0.25),color=:curl,aspect_ratio=:equal,xlabel="x",ylabel="y",title="vorticidad")
end
gif(anim,"../outputs/$(dest)/rot.gif",fps=30)

anim= @animate for t in 2:30:iter

    heatmap(sqrt(F[t][3].^2 + F[t][4].^2)',clim=(0.0,vel*1.5),color=:curl,aspect_ratio=:equal,xlabel="x",ylabel="y",title="rapidez")
end
gif(anim,"../outputs/$(dest)/vel.gif",fps=30)
anim= @animate for t in 2:30:iter
  
    heatmap((F[t][3])',clim=(-vel,vel),color=:curl,aspect_ratio=:equal,xlabel="x",ylabel="y",title="velocidad en x")
end
gif(anim,"../outputs/$(dest)/velx.gif",fps=30)
anim= @animate for t in 2:30:iter
    heatmap((F[t][4])',clim=(-vel,vel),color=:curl,aspect_ratio=:equal,xlabel="x",ylabel="y",title="velocidad en y")
end
gif(anim,"../outputs/$(dest)/vely.gif",fps=30)


anim= @animate for t in 2:30:iter
    X=[P[t][i,j][1] for i in 1:10 for j in 1:10]
    Y=[P[t][i,j][2] for i in 1:10 for j in 1:10]
    scatter(X,Y,xlim=(0.0,dims[1]),ylim=(0,dims[2]))
end
gif(anim,"../outputs/$(dest)/part.gif",fps=30)


println("Tiempo para realizar animaciones: $(round(Dates.time()-s,2)) segs")



