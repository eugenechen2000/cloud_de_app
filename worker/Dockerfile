FROM rocker/r-ver:4.3.3

RUN apt-get update && apt-get install -y python3 python3-pip && rm -rf /var/lib/apt/lists/*
WORKDIR /app
COPY requirements.txt .
RUN pip3 install --no-cache-dir -r requirements.txt
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org')" && \
    R -e "BiocManager::install(c('DESeq2','org.Hs.eg.db','AnnotationDbi'), ask=FALSE, update=FALSE)" && \
    R -e "install.packages(c('ggplot2','ggrepel','msigdbr','fgsea'), repos='https://cloud.r-project.org')"
COPY . .
CMD ["python3", "-m", "rq", "worker", "de-jobs", "--url", "redis://redis:6379/0"]
