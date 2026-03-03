FROM continuumio/miniconda3:latest
ENV DEBIAN_FRONTEND=noninteractive

# Install system dependencies
RUN apt-get update && apt-get install -y \
	wget \
	bzip2 \
	time \
	&& apt-get clean && rm -rf /var/lib/apt/lists/*

# Copy environment file and create conda environment
COPY pseudoDB_env.yaml /tmp/pseudoDB_env.yaml
RUN conda env create -n gatk3 -f /tmp/pseudoDB_env.yaml && \
	conda clean -afy
ENV PATH="/opt/conda/envs/gatk3/bin:$PATH"

# Install GATK 3.8 into the environment
SHELL ["conda", "run", "-n", "gatk3", "/bin/bash", "-c"]
RUN wget https://storage.googleapis.com/gatk-software/package-archive/gatk/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2 && \
	tar -jxvf GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2 && \
	mv GenomeAnalysisTK-3.8-1-0-gf15c1c3ef gatk3.8 && \
	gatk-register gatk3.8/GenomeAnalysisTK.jar && \
	rm -rf GenomeAnalysisTK-3.8-1-0-gf15c1c3ef.tar.bz2 gatk3.8

# Copy pipeline script
COPY script/pseudoDB.py /opt/pseudoDB/pseudoDB.py

# Run
ENTRYPOINT ["/usr/bin/time", "--verbose", "/opt/conda/envs/gatk3/bin/python3", "/opt/pseudoDB/pseudoDB.py"]