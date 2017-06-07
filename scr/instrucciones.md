Para correr la simulación de manera Hay que ejecutar el siguiente comando desde una terminal:

`julia simulacion.jl visc vel bar geo`
Este comando corre la simulación con los siguientes parámetros
* `visc`: Velocidad inicial macroscópica del fluido en la dirección x positiva, debe de estar entre 0.005 y 0.2 
* `vel`: Velocidad inicial macroscópica del fluido en la dirección x positiva, debe de estar entre 0.005 y 0.2 
* `bar`: geometría de la barrera que impacta el fluido. Puede elegirse entre `o` para un circulo y `|` para un segmento de recta vertical.
* `geo`: geometría de las condiciones de frontera. Puede elegirse entre condiciones normales de un túnel de viento cerrado (`normal`), condiciones periodicas en el extremo x y -x (`anillo`) y condiciones periodicas en los cuatro lados (`toro`)

Para correr la simulación con condiciones por default, simplemente ejecutar `julia simulacion.jl`