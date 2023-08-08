########################################################
#        Renku install section - do not edit           #

# For finding latest versions of the base image see
# https://github.com/SwissDataScienceCenter/renkulab-docker
FROM renku/renkulab-r:4.2.0-0.18.1 as builder

# RENKU_VERSION determines the version of the renku CLI
# that will be used in this image. To find the latest version,
# visit https://pypi.org/project/renku/#history.
ARG RENKU_VERSION=2.6.1

# Install renku from pypi or from github if a dev version
RUN if [ -n "$RENKU_VERSION" ] ; then \
        source .renku/venv/bin/activate ; \
        currentversion=$(renku --version) ; \
        if [ "$RENKU_VERSION" != "$currentversion" ] ; then \
            pip uninstall renku -y ; \
            gitversion=$(echo "$RENKU_VERSION" | sed -n "s/^[[:digit:]]\+\.[[:digit:]]\+\.[[:digit:]]\+\(rc[[:digit:]]\+\)*\(\.dev[[:digit:]]\+\)*\(+g\([a-f0-9]\+\)\)*\(+dirty\)*$/\4/p") ; \
            if [ -n "$gitversion" ] ; then \
                pip install --force "git+https://github.com/SwissDataScienceCenter/renku-python.git@$gitversion" ;\
            else \
                pip install --force renku==${RENKU_VERSION} ;\
            fi \
        fi \
    fi
#             End Renku install section                #
########################################################

FROM renku/renkulab-r:4.2.0-0.18.1


LABEL maintainer="Swiss Data Science Center <info@datascience.ch>"

USER root
SHELL [ "/bin/bash", "-c", "-o", "pipefail" ]

# Add VSCode support for a more pleasant coding and debugging experience for .py files. More details about issues, comments, and automatically installing extensions on https://renku.discourse.group/t/using-visual-studio-code-in-renkulab-interactive-sessions/249.
RUN curl -s https://raw.githubusercontent.com/SwissDataScienceCenter/renkulab-docker/master/scripts/install-vscode.sh | bash

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update -y \
    && apt-get install -yq --no-install-recommends \
        dbus-x11 \
        firefox \
        net-tools \
        less \
        xfce4 \
        xfce4-panel \
        xfce4-session \
        xfce4-settings \
        xorg \
        xubuntu-icon-theme \
        gnome-terminal \
        fonts-dejavu \
        git-gui \
        gitk \
        emacs \
    && apt-get autoremove --purge \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/* \
    && rm -rf /tmp/* \
    && find /var/log -type f -exec cp /dev/null \{\} \;

#################################################################
# Install noVNC

ENV novnc_version=1.1.0

RUN curl -sSfL https://github.com/novnc/noVNC/archive/v${novnc_version}.tar.gz | tar xvz -C /opt/ && \
    chmod a+rX -R /opt/noVNC-${novnc_version}

COPY --chown=root:root vnc_renku.html /opt/noVNC-${novnc_version}
COPY --chown=root:root renku-48x48.png /opt/noVNC-${novnc_version}/app/images/icons

COPY --chown=root:root ui.js /opt/noVNC-${novnc_version}/app/ui.js
COPY --chown=root:root base.css /opt/noVNC-${novnc_version}/app/styles/base.css
COPY --chown=root:root fonts /opt/noVNC-${novnc_version}/app/styles/fonts

#################################################################
# Install TigerVNC

ENV tigervnc_version=1.9.0

RUN curl -sSfL https://sourceforge.net/projects/tigervnc/files/stable/${tigervnc_version}/tigervnc-${tigervnc_version}.x86_64.tar.gz/download | tar -zxf - -C /usr/local --strip=2

#################################################################
# Add desktop icons
# Add wallpaper and force the default to point to the renku one

COPY Git-GUI.desktop /home/jovyan/Desktop/
COPY gitk.desktop /home/jovyan/Desktop/
COPY iTol.desktop /home/jovyan/Desktop/
COPY .git_icon.png /home/jovyan/Desktop/
COPY renku_background_dots.png /usr/share/backgrounds/renku/
RUN chmod +x /home/jovyan/Desktop/gitk.desktop && \
    chmod +x /home/jovyan/Desktop/Git-GUI.desktop && \
    mkdir -p /usr/share/backgrounds/renku/ && \
    ln -s -f /usr/share/backgrounds/renku/renku_background_dots.png /usr/share/backgrounds/xfce/xfce-stripes.png

#################################################################
# Install the jupyter extensions


RUN mamba install -y jupyter-server-proxy numpy websockify -c conda-forge \
    && jupyter labextension install @jupyterlab/server-proxy \
    && mamba clean -y --all

COPY jupyter_notebook_config.py /home/jovyan/.jupyter/jupyter_notebook_config.py

COPY post-init.sh /post-init.sh

# install the R dependencies
## make the renv install script and renv.lock file
## available in the working dir and run the install
COPY .renv_install.sh .
COPY renv.lock .
COPY install.R /tmp

RUN bash .renv_install.sh
## ensure renv lock is in the project directory
COPY renv.lock /home/rstudio/renv.lock
RUN R -f /tmp/install.R

# To apply a custom RStudio config uncomment the line below
# ENV RSTUDIO_CONFIG_HOME=/home/rstudio/work/renv/.rstudio_config_dir

## Clean up the /home/rstudio directory to avoid confusion in nested R projects
RUN rm /home/rstudio/.Rprofile; rm /home/rstudio/renv.lock


#################################################################
# Add from here your custom installations
#################################################################

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get -y update && \
  apt-get clean && \
  apt-get install -y --no-install-recommends \
    apt-utils \
    bash-completion \
    vim \
    jq \
    sudo \
    wget \
    curl \
    unzip \
    bzip2 \
    g++ \
    libncurses5-dev \
    libncursesw5-dev \
    parallel \
    libgit2-dev \
    tk-dev \
    jq \
    curl \
    librsvg2-2 \
    librsvg2-dev \
    librsvg2-common \
    zlib1g \
    python3-opencv && \
    apt-get -y clean && \
    apt-get -y autoclean  && \
    apt-get -y autoremove

## add root user
ARG DOCKER_USER=polyploid
ARG DOCKER_PASSWORD=polyploid

RUN useradd ${DOCKER_USER}
RUN usermod -aG sudo ${DOCKER_USER}
RUN echo ${DOCKER_USER}:${DOCKER_PASSWORD} | chpasswd

## if necessary
# RUN su - ${DOCKER_USER}


#################################################################
#################################################################

### Install Octopus
# Get dependencies
RUN apt-get -y update \
    && apt-get -y install \
        build-essential \
        libboost-all-dev \
        libgmp-dev \
        cmake \
        libhts-dev \
        python3-pip \
        git \
    && pip3 install distro


# Install Octopus
ARG THREADS=4
ARG CPU=haswell
WORKDIR /tmp
RUN git clone https://github.com/luntergroup/octopus.git && \
    octopus/scripts/install.py \
    --threads $THREADS \
    --architecture $CPU

# COPY octopus/bin/* /opt/conda/bin/.
# ENV PATH=octopus/bin:${PATH}
RUN mv /tmp/octopus/bin/octopus /usr/local/bin/


# EAGLE-RC
# https://github.com/tony-kuo/eagle


WORKDIR /tmp
RUN git clone https://github.com/tony-kuo/eagle.git
WORKDIR /tmp/eagle
RUN git clone https://github.com/samtools/htslib.git
WORKDIR /tmp/eagle/htslib
RUN git submodule update --init --recursive
WORKDIR /tmp/eagle
RUN make

RUN mv /tmp/eagle/eagle /usr/local/bin/
RUN mv /tmp/eagle/eagle-rc /usr/local/bin/

# install hifiasm
WORKDIR /tmp
RUN git clone https://github.com/chhylp123/hifiasm && \
    cd hifiasm && make && \
    mv /tmp/hifiasm/hifiasm /usr/local/bin/

# Add user to the sudoers
RUN adduser ${NB_USER} sudo && \
    echo "${NB_USER} ALL=(ALL:ALL) ALL" >> /etc/sudoers
ARG NB_PASSWORD=rstudio
RUN echo ${NB_USER}:${NB_PASSWORD} | chpasswd




COPY requirements.txt env/assembly.yml env/environment.yml env/wgd.yml env/ksrates.yml env/REC_env.yml env/snp_calling.yml env/updog_install.R /tmp/
# update conda base environment
RUN mamba env update -q -f /tmp/environment.yml && \
    /opt/conda/bin/pip install -r /tmp/requirements.txt --no-cache-dir

RUN chown -R  ${NB_USER}:${NB_USER} /home/rstudio
RUN chown -R ${NB_USER}:${NB_USER} /opt/conda/

USER ${NB_USER}


# ### installing WGDv1
RUN conda create -n wgd python=3.8 && \
    mamba env update -n wgd -q -f /tmp/wgd.yml 
WORKDIR /tmp
RUN git clone --depth 1 https://github.com/rimjhimroy/wgd.git && \
    cd wgd && \
    /opt/conda/envs/wgd/bin/pip install numpy && \
    /opt/conda/envs/wgd/bin/pip install .



### installing ksrates
RUN conda create -n ksrates python=3.8 && \
    mamba env update -n ksrates -q -f /tmp/ksrates.yml 
WORKDIR /tmp
RUN git clone --depth 1 https://github.com/VIB-PSB/ksrates && \
    cd ksrates && \
    /opt/conda/envs/ksrates/bin/pip install numpy==1.23.0 && \
    /opt/conda/envs/ksrates/bin/pip install .


# install the snp_calling conda environment which has GATK, Freebayes, Updog
RUN conda create -n snp_calling && \
    mamba env update -n snp_calling -f /tmp/snp_calling.yml && \
    /opt/conda/envs/snp_calling/bin/R -f /tmp/updog_install.R

# install the assembly conda environment 
RUN conda create -n assembly && \
    mamba env update -n assembly -f /tmp/assembly.yml

    
    

### installing Reconciliation environment REC_env
RUN conda create -n REC_env python=3.8 && \
    mamba env update -n REC_env -q -f /tmp/REC_env.yml

RUN mamba clean -y --all && \
    mamba env export -n "root"

# install the python dependencies
COPY requirements.txt /tmp/
RUN pip3 install -r /tmp/requirements.txt --no-cache-dir && \
    rm -rf ${HOME}/.renku/venv

WORKDIR /home/rstudio
COPY .Rprofile .   
COPY --from=builder ${HOME}/.renku/venv ${HOME}/.renku/venv
