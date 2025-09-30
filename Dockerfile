# if building on M1 uncomment the line 2, and comment out line 3
#FROM --platform=linux/amd64 ubuntu as dynadjust-build
FROM ubuntu as dynadjust-build

RUN mkdir -p /app/DynAdjust/ && mkdir -p /opt/dynadjust/geoid_file
COPY . /app/DynAdjust/
WORKDIR /app

RUN apt-get update && \
    DEBIAN_FRONTEND="noninteractive" apt-get install -y --no-install-recommends \
    libxerces-c-dev \
    xsdcxx \
    libboost-program-options-dev \
    libopenblas-dev \
    liblapacke-dev \
    cmake \
    make \
    g++ \
    wget && \
    rm -rf /var/lib/apt/lists/*

RUN cd ./DynAdjust && \
    cmake -S dynadjust -B build -DBUILD_TESTING=ON && \
    cmake --build build --parallel && \
    ctest --test-dir build

RUN wget --no-check-certificate "https://s3-ap-southeast-2.amazonaws.com/geoid/AUSGeoid/AUSGeoid2020_20180201.gsb" -P "/opt/dynadjust/geoid_file"
