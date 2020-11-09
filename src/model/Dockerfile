FROM continuumio/miniconda3:4.8.2
LABEL description="Image for pytorch model"

ENV PYTHONDONTWRITEBYTECODE=true 

# update conda environment with packages and clean up conda installation by removing 
# conda cache/package tarbarlls and python bytecode
COPY environment.yml .
RUN conda env update -f environment.yml \
    && conda clean -afy \
    && find /opt/conda/ -follow -type f -name '*.pyc' -delete 
CMD [ "/bin/bash" ] 