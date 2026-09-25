FROM ubuntu
ENV DEBIAN_FRONTEND=noninteractive

# Install any needed packages specified in requirements.txt
RUN apt-get update && apt-get install -y git cmake g++ libcgal-dev libgmp3-dev

# Set the working directory to /app
WORKDIR /app

# Download and compile TetWild
RUN git clone https://github.com/wildmeshing/fTetWild.git --recursive
WORKDIR /app/fTetWild/build
RUN cmake .. && make

WORKDIR /data

ENTRYPOINT ["/app/fTetWild/build/FloatTetwild_bin"]
