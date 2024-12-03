FROM mambaorg/micromamba:1.5.8-alpine3.19
# FROM continuumio/miniconda3:24.1.2-0

USER root
# RUN addgroup -S molgpka && adduser -S -G molgpka molgpka

ENV LANG=C.UTF-8 LC_ALL=C.UTF-8
ENV PATH /opt/conda/bin:/opt/conda/envs/env/bin:/opt/micromamba/bin:/opt/micromamba/envs/env/bin:$PATH
ENV PYTHONPATH /src/src

RUN apk add --upgrade apk-tools
RUN apk upgrade --available
RUN apk add git

ARG CONDA_ENV="MolGpka"

# USER molgpka

COPY environment.yml /src/environment.yml
RUN micromamba create -n $CONDA_ENV -c conda-forge rdkit python=3.11
RUN micromamba install -n $CONDA_ENV -f /src/environment.yml

RUN micromamba install -n $CONDA_ENV -c pytorch -c conda-forge torch-scatter

RUN micromamba run -n $CONDA_ENV pip install numpy==2.0.1

# RUN micromamba remove -n $CONDA_ENV numpy && \
# 	micromamba install -n $CONDA_ENV -c conda-forge numpy

# RUN micromamba clean -p -t -l --trash -y

COPY . /src/

# USER root

COPY docker/flask/uwsgi.ini /etc/uwsgi/uwsgi.ini
# COPY docker/flask/start_flask.sh /src/start_flask.sh

RUN chmod 755 /src/docker/flask/start_flask.sh

WORKDIR /src
# EXPOSE 8080
# RUN chown -R molgpka:molgpka /src

# USER molgpka

# #ENTRYPOINT ["tail", "-f", "/dev/null"]

CMD ["sh", "/src/docker/flask/start_flask.sh"]

# CMD ["micromamba", "run", "-n", "MolGpka", "python", "/src/cts_molgpka_flask.py"]
