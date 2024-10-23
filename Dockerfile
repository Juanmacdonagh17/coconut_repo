# The image is built in two stages: builder and runner.
# The idea is that the context in which you compile and run the program is
# different since you usually have many more requirements for compiling.

# use a C image
FROM gcc:12.4-bookworm as builder
# choose the working directory and copy files
WORKDIR /usr/src/coconut
COPY . .
#  install dependencies and compile
RUN apt install -y libcurl-dev
RUN gcc -o coconut coconut.c -lcurl

#  much smaller Linux image

FROM debian:11-slim as runner
WORKDIR /usr/src/coconut

#  copy the necessary dependencies and libs for the program to run from the builder

COPY --from=builder /lib/x86_64-linux-gnu/* /lib/x86_64-linux-gnu/
COPY --from=builder /lib64/* /lib64/
#  copy the executable from the builder
COPY --from=builder /usr/src/coconut/coconut /usr/src/coconut/coconut
#  add a user so that whoever uses this image doesn't have root privileges without knowing it

RUN useradd -ms /bin/bash cocouser
USER cocouser

# The entrypoint is the command that runs implicitly when running the image
# the rest of the things you add to the `docker run IMAGE` command are its arguments.

ENTRYPOINT /usr/src/coconut/coconut 

# To test it, you can run 
# `docker run user/cocodocker ARGUMENTS` 
# The first time, it has to download everything, so it will take a bit.
# Since everything runs in a virtualized image, the output files will remain
# inside the container. To access them, you can use a volume, which mounts a folder from the
# image inside your computer. For that, you have to add -v /your/pc/folder:/container/folder
# you need to use the full path.

# You should end up with something like this:
# `docker run -v /your/folder/lalalal/:/usr/src/coconut ARGUMENTS`