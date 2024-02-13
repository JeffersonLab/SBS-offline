FROM almalinux:9.2


ARG APP_VERSION
ARG REPO_NAME
ARG ANALYZER_REPO_NAME=analyzer
ARG ANALYZER_VERSION

SHELL ["/bin/bash", "-c"]
ADD http://pki.jlab.org/JLabCA.crt /etc/pki/ca-trust/source/anchors/JLabCA.crt
RUN update-ca-trust
RUN dnf update -y
RUN dnf -y install 'dnf-command(config-manager)'
RUN dnf -y install  epel-release  dnf-plugins-core
RUN dnf config-manager --set-enabled crb


RUN dnf -y install epel-release  git g++ cmake gcc-c++ make root root-mathcore root-montecarlo-eg \
                               root-mathmore root-gui root-hist root-physics root-genvector

ADD https://github.com/JeffersonLab/${ANALYZER_REPO_NAME}/archive/refs/tags/${ANALYZER_VERSION}.tar.gz .
RUN tar -xvf ${ANALYZER_VERSION}.tar.gz && rm ${ANALYZER_VERSION}.tar.gz
WORKDIR "/${ANALYZER_REPO_NAME}-${ANALYZER_VERSION}"
RUN cmake -DCMAKE_INSTALL_PREFIX=/usr/local/analyzer -B BUILD -S /${ANALYZER_REPO_NAME}-${ANALYZER_VERSION}/
RUN cmake --build BUILD -j8
RUN cmake --install BUILD
ENV CMAKE_INSTALL_PREFIX="/usr/local/analyzer"
ENV ANALYZER="/usr/local/analyzer"
ENV PATH="/usr/local/analyzer/bin:/usr/bin/root:$PATH"
ENV LD_LIBRARY_PATH="/usr/local/analyzer/lib64:$LD_LIBRARY_PATH"
ENV ROOTSYS="/usr/"
RUN git clone https://github.com/JeffersonLab/SBS-offline.git --branch ${APP_VERSION} /SBS-offline
RUN mkdir /build
WORKDIR "/build"
RUN cmake -DCMAKE_INSTALL_PREFIX=/usr/local/sbs-offline  -S /SBS-offline 
RUN make install
#script to set up the environment for SBS-offline
ENV SBS="/usr/local/sbs-offline"
ENV SBSOFFLINE="/usr/local/sbs-offline"
ENV PATH="/usr/local/sbs-offline/bin:$PATH"
ENV LD_LIBRARY_PATH="/usr/local/sbs-offline/lib64:$LD_LIBRARY_PATH"
