
FROM ubuntu:latest

WORKDIR /build

COPY . .

RUN apt-get update -qq && \
    apt-get install -y zlib1g-dev \
                    git \
                    cmake \
                    build-essential \
                    python3

RUN rm -rf build; mkdir build; cd build; cmake ..; make; make install

WORKDIR /bin

CMD ['bash']