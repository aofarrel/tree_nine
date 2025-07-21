# EXPERIMENTAL - In usher-plus:0.6.4_5 we swapped the base image from ubunutu:22.04 to miniconda3:22.11.1 which may
# cause unexpected issues. If you are having issues with usher-plus:0.6.4_5, try ashedpotatoes/usher-plus:0.6.4_4
# or, for pure UShER, ashedpotatoes/usher-mirror:0.6.4
# Why did I do this? I needed to include BTE, which pretty much requires conda. I looked into installing miniconda
# on the UShER base image, but the overall experience didn't exactly spark joy. So, we're going to handle BTE's
# requirements first, and then install UShER on top of that.
FROM continuumio/miniconda3:22.11.1


# Install BTE
USER root
WORKDIR /home
RUN git clone --recurse-submodules https://github.com/aofarrel/BTE.git
WORKDIR BTE
RUN conda env create -f bte.yml
SHELL ["conda", "run", "-n", "bte", "/bin/bash", "-c"]
RUN python3 setup.py build_ext && \
 python3 setup.py install && \
 ln -s /home/BTE/build/lib.linux-x86_64-cpython-310/bte.cpython-310-x86_64-linux-gnu.so /opt/conda/lib/python3.10/site-packages

# UShER pre-install setup
SHELL ["/bin/bash", "-c"]
ENV APT_KEY_DONT_WARN_ON_DANGEROUS_USAGE=DontWarn
ENV DEBIAN_FRONTEND=noninteractive
USER root
RUN apt-get update && apt-get install -yq --no-install-recommends git wget ca-certificates sudo

# faSomeRecords and faSize are needed for the UShER WDL workflow 
WORKDIR /HOME/kentsource
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faSomeRecords
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faSize

# Install UShER
# Astute readers will know that BTE comes with a packaged version of UShER's source. Too bad! We're not using it!
RUN chmod 775 *
WORKDIR /HOME
RUN git clone https://github.com/yatisht/usher.git
######RUN git checkout $(git describe --tags `git rev-list --tags --max-count=1`) --> no.
WORKDIR usher
RUN ./install/installUbuntu.sh 
ENV PATH="/HOME/usher/build:/HOME/kentsource:${PATH}"

# Add the gold standard TB reference - H37Rv
# md5sum of ref.fa is 0e8f3acd2c92f20c42105e014adbfe97
# Taken from clockwork ref prep, but this is only H37Rv stuff -- no decontamination reference here!
RUN mkdir ref
COPY ./Ref.H37Rv.tar ./ref/
RUN cd ./ref/ && tar -xvf Ref.H37Rv.tar

# Add the "known lineage SRA tree", which should ONLY be used for debugging, folks!
RUN mkdir example_tree
COPY ./data/for_debugging_only__tb_7K_noQC_diffs_mask2ref.L.fixed.pb ./example_tree/

# Bonus content! Huzzah!
RUN apt-get install -y tree vim zip pigz
RUN pip install six numpy ete3 pandas polars requests # six is a prereq for ete3 that doesn't get installed when installing ete3
RUN mkdir /scripts/
COPY ./find_clusters.py /scripts/
COPY ./process_clusters.py /scripts/
COPY ./summarize_changes.py /scripts/
RUN wget -P /scripts https://raw.githubusercontent.com/aofarrel/diffdiff/main/diffdiff.py
RUN wget -P /scripts -o marcs_incredible_script.pl https://gist.githubusercontent.com/aofarrel/a638f2ff05f579193632f7921832a957/raw/baa77b4f6afefd78ae8b6a833121a413bd359a5e/marcs_incredible_script

