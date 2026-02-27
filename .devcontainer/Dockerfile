# last update: 2025-06-03
FROM quay.io/jupyter/datascience-notebook:julia-1.11.4

WORKDIR /workdir
EXPOSE 8888

# install JuliaFormatter for VSCode extension
RUN julia -e 'using Pkg; Pkg.add("JuliaFormatter")'
# install project-specific packages
COPY Project.toml Manifest.toml ./
RUN julia --project=. -e 'using Pkg; Pkg.instantiate()'
# set default to the current project
ENV JULIA_PROJECT=@.
# enable multi-threading
ENV JULIA_NUM_THREADS=12
