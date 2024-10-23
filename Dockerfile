# Amigooo te dejo unos comentarios para que sepas que hace cada cosa


# La imagen se contruye en dos etapas: builder y runner,
# la idea es que el contexto en el que compilas y corres el programa es
# distinto, ya que habitualmente tenés muchos más requisitos para compilar.

# Uso una imagen de C
FROM gcc:12.4-bookworm as builder
# Elijo el dir de trabajo y copio archivos
WORKDIR /usr/src/coconut
COPY . .
# Instalo dependencias y compilo
RUN apt install -y libcurl-dev
RUN gcc -o coconut coconut.c -lcurl

# Uso ahora una imagen de linux mucho más chiquita
FROM debian:11-slim as runner
WORKDIR /usr/src/coconut
# Copio las dependencias y libs necesarias para que ande el programa
# desde builder
COPY --from=builder /lib/x86_64-linux-gnu/* /lib/x86_64-linux-gnu/
COPY --from=builder /lib64/* /lib64/
# Copio el ejecutable desde builder
COPY --from=builder /usr/src/coconut/coconut /usr/src/coconut/coconut
# Agrego un usuario para que quien use esta imagen no tenga privilegios de root
# sin saberlo
RUN useradd -ms /bin/bash cocouser
USER cocouser
# El entrypoint es el comando que se corre implicitamente al correr la imagen
# el resto de cosas que le agregues al comando `docker run IMAGEN` son sus argumentos.
ENTRYPOINT /usr/src/coconut/coconut 
# Para probarlo podes correr 
# `docker run egonik/cocodocker ARGUMENTOS` 
# La primera vez tiene que descargar todo asi que va a tardar un toque.
# Como todo se corre en una imagen virtualizada los archivos de salida van a quedar
# dentro del contenedor. Para poder acceder podes usar un volumen, que monta una carpeta de la
# imagen dentro de tu compu. Para eso tenes que agregar -v /carpeta/de/tu/pc:/carpeta/del/contenedor
# el path que tenes que usar es el completo.

# Te tendría que quedar algo asi:
# `docker run -v /tu/carpeta/lalalal/:/usr/src/coconut ARGUMENTOS`
