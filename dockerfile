# Use the official Ubuntu 22.04 LTS image as the base
FROM ubuntu:22.04

# Prevent interactive prompts during package installation
ENV DEBIAN_FRONTEND=noninteractive

# Set the working directory inside the container
WORKDIR /app

# Copy the current directory (where the Dockerfile is located) to the /app directory in the container
COPY . /app

# Update the package list and install dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    build-essential \
    g++-12 \
    libeigen3-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# Add the PATH modification to .bashrc
RUN LINE='export PATH="/usr/include/eigen3:$PATH"' && \
    grep -qxF "$LINE" ~/.bashrc || echo "$LINE" >> ~/.bashrc

# Set the default command to run bash (optional, you can modify as needed)
CMD ["/bin/bash"]