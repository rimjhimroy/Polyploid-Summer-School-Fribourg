########################################################
#        Renku install section - do not edit           #

FROM renku/renkulab-bioc:RELEASE_3_17-0.17.0 as builder

# RENKU_VERSION determines the version of the renku CLI
# that will be used in this image. To find the latest version,
# visit https://pypi.org/project/renku/#history.
ARG RENKU_VERSION=2.5.0

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

FROM renku/renkulab-bioc:RELEASE_3_17-0.17.0

# Add VSCode support for a more pleasant coding and debugging experience for .py files. More details about issues, comments, and automatically installing extensions on https://renku.discourse.group/t/using-visual-studio-code-in-renkulab-interactive-sessions/249.
RUN curl -s https://raw.githubusercontent.com/SwissDataScienceCenter/renkulab-docker/master/scripts/install-vscode.sh | bash

# Uncomment and adapt if code is to be included in the image
COPY src /code/src

# Uncomment and adapt if your R or python packages require extra linux (ubuntu) software
# e.g. the following installs apt-utils and vim; each pkg on its own line, all lines
# except for the last end with backslash '\' to continue the RUN line
#
USER root
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    apt-utils \
    bash-completion \
    vim \
    jq \
    less \
    wget \
    curl \
    unzip \
    bzip2 \
    python3-opencv && \
    apt-get -y clean && \
    apt-get -y autoclean  && \
    apt-get -y autoremove

USER ${NB_USER}

# install the R dependencies
COPY install.R /tmp/
RUN R -f /tmp/install.R

# install the python dependencies
COPY requirements.txt assembly.yml environment.yml wgd.yml ksrates.yml  /tmp/
RUN conda create -n assembly
RUN mamba env update -n assembly -q -f /tmp/assembly.yml
RUN mamba env update -q -f /tmp/environment.yml && \
    /opt/conda/bin/pip install -r /tmp/requirements.txt --no-cache-dir && \
    mamba clean -y --all && \
    mamba env export -n "root"

COPY --from=builder ${HOME}/.renku/venv ${HOME}/.renku/venv
