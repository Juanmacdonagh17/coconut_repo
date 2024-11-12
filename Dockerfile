# builder
FROM gcc:12.4-bookworm as builder
WORKDIR /usr/src/coconut
COPY . .
# curl and gcc
RUN apt-get update && apt-get install -y libcurl-dev && rm -rf /var/lib/apt/lists/*
RUN gcc -o coconut coconut.c -lcurl

# Runner
FROM debian:11-slim as runner

# add a user so that whoever uses this image doesn't have root privileges without knowing it
RUN useradd -ms /bin/bash cocouser

# switch to the user's home directory (cocouser)
WORKDIR /home/cocouser

# install required runtime libraries
RUN apt-get update && apt-get install -y libcurl4 && rm -rf /var/lib/apt/lists/*

# copy the executable from the builder
COPY --from=builder /usr/src/coconut/coconut /home/cocouser/coconut

# change ownership to cocouser
RUN chown -R cocouser:cocouser /home/cocouser

USER cocouser

ENTRYPOINT ["/home/cocouser/coconut"]
