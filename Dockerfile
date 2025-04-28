# FROM mambaorg/micromamba@sha256:7358b481a67025cf0356acf9999cca0625b8c6883977d947b164c173471dfcda
FROM micromamba:0.24.0
COPY --chown=$MAMBA_USER:$MAMBA_GROUP env.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes
# add base environment to path to ensure this works in Singularity
# without this, singularity exec and shell will not behave correctly
ENV PATH "$MAMBA_ROOT_PREFIX/bin:$PATH"
COPY transcriptome.py /usr/local/bin
COPY annofilter_junctions.py /usr/local/bin
