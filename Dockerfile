FROM ubuntu:24.04
WORKDIR /home

# Add the gold standard TB reference - H37Rv
# md5sum of ref.fa is 0e8f3acd2c92f20c42105e014adbfe97
# Taken from clockwork ref prep, but this is only H37Rv stuff -- no decontamination reference here!
RUN mkdir ref
COPY ./Ref.H37Rv.tar ./ref/
RUN cd ./ref/ && tar -xvf Ref.H37Rv.tar

# Add the "known lineage SRA tuberculosis tree", which should ONLY be used for debugging, folks!
RUN mkdir example_tree
COPY ./data/for_debugging_only__tb_7K_noQC_diffs_mask2ref.L.fixed.pb ./example_tree/

# Pre-install setup
SHELL ["/bin/bash", "-c"]
ENV APT_KEY_DONT_WARN_ON_DANGEROUS_USAGE=DontWarn
ENV DEBIAN_FRONTEND=noninteractive
USER root
RUN apt-get update && apt-get install -yq --no-install-recommends \
    build-essential \
	ca-certificates \ 
	curl \
	git \
	pigz \
	sudo \
	tree \
	wget \
	vim \
	zip \
	&& rm -rf /var/lib/apt/lists/*

# Conda stuff
ENV CONDA_DIR=/opt/conda
ENV PATH=${CONDA_DIR}/bin:${PATH}
RUN wget https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O /tmp/miniforge.sh && \
    /bin/bash /tmp/miniforge.sh -b -p ${CONDA_DIR} && \
    rm /tmp/miniforge.sh
RUN conda init bash

# Install BTE
USER root
WORKDIR /home
RUN git clone --recurse-submodules https://github.com/aofarrel/BTE.git
WORKDIR /home/BTE
RUN mamba env update -n base --file bte.yml && \
    mamba clean --all -y
RUN export CONDA_PREFIX=/opt/conda && \
    python3 setup.py build_ext && \
    python3 setup.py install

# faSomeRecords and faSize are needed for the UShER WDL workflow 
WORKDIR /home/kentsource
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faSomeRecords
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faSize

# Install UShER
# Previous versions of this image just pulled master, but now I'm pinning to latest tagged commit. UShER doesn't
# have a lot of tagged commits though so I might revert this: https://github.com/yatisht/usher/compare/v0.6.6...master
# Astute readers will know that BTE comes with a packaged version of UShER's source. Too bad! We're not using it!
RUN chmod 775 *
WORKDIR /home
RUN git clone https://github.com/yatisht/usher.git
RUN git checkout $(git describe --tags $(git rev-list --tags --max-count=1))
WORKDIR usher
RUN ./install/installUbuntu.sh 
ENV PATH="/home/usher/build:/home/kentsource:${PATH}"

# Bonus content! Huzzah!
# NOTE: six must be installed before ete3 or else ete3 won't work properly
RUN pip install six numpy ete3 pandas polars requests

# Install PhyloDM
ENV PYO3_USE_ABI3_FORWARD_COMPATIBILITY=1
RUN pip install phylodm

# Add my scripts
RUN mkdir /scripts/
COPY ./find_clusters.py /scripts/
COPY ./process_clusters.py /scripts/
COPY ./summarize_changes.py /scripts/
RUN wget -P /scripts https://raw.githubusercontent.com/aofarrel/diffdiff/main/diffdiff.py
RUN wget -P /scripts -o marcs_incredible_script.pl https://gist.githubusercontent.com/aofarrel/a638f2ff05f579193632f7921832a957/raw/baa77b4f6afefd78ae8b6a833121a413bd359a5e/marcs_incredible_script

