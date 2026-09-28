# if building on M1 uncomment the line 2, and comment out line 3
#FROM --platform=linux/amd64 ubuntu  as dynadjust-build
FROM ubuntu as dynadjust-build
RUN mkdir -p /app/DynAdjust/ && mkdir -p /opt/dynadjust/geoid_file
COPY . /app/DynAdjust/
WORKDIR /app
RUN apt-get update &&\
 apt-get install -y --no-install-recommends apt-utils &&\
 apt-get -y --no-install-recommends install sudo &&\
 DEBIAN_FRONTEND="noninteractive" &&\
 apt update &&\
 apt-get install -y tzdata &&\
 apt-get install -y wget &&\
 apt-get install -y software-properties-common &&\
 apt-get install -y --no-install-recommends gnupg2 &&\
 apt-get install -y --no-install-recommends p7zip &&\
 apt-get install -y --no-install-recommends libboost-system-dev libboost-filesystem-dev libboost-timer-dev libboost-thread-dev libboost-program-options-dev &&\
 apt-get install -y --no-install-recommends libxerces-c-dev &&\
 apt-get install -y --no-install-recommends cmake &&\
 apt-get install -y --no-install-recommends make &&\
 apt-get install -y --no-install-recommends g++ &&\
 apt-get install -y --no-install-recommends xsdcxx
RUN wget -qO /tmp/intel-oneapi-key.pub https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS.PUB &&\
 gpg --dearmor --output /usr/share/keyrings/oneapi-archive-keyring.gpg /tmp/intel-oneapi-key.pub &&\
 rm /tmp/intel-oneapi-key.pub &&\
 echo "deb [signed-by=/usr/share/keyrings/oneapi-archive-keyring.gpg] https://apt.repos.intel.com/oneapi all main" | tee /etc/apt/sources.list.d/oneAPI.list > /dev/null
RUN apt-get update &&\
 apt-get install -y intel-oneapi-compiler-dpcpp-cpp intel-oneapi-mkl-devel intel-oneapi-tbb-devel
RUN apt-get -y install cpio &&\
 cd ./DynAdjust &&\
 chmod +x ./resources/make_dynadjust_gcc.sh &&\
 ./resources/make_dynadjust_gcc.sh --auto
RUN wget "https://s3-ap-southeast-2.amazonaws.com/geoid/AUSGeoid/AUSGeoid2020_20180201.gsb" -P "/opt/dynadjust/geoid_file"
