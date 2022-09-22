FROM dealii/dealii:latest

USER root

RUN apt-get update \
 && apt-get install -y vim ack colordiff gdb less mlocate source-highlight ripgrep \
 && apt-get clean autoclean \
 && rm -rf /var/lib/apt/lists/*

RUN updatedb
