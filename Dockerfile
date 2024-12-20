FROM ubuntu:22.04
ENV APT_KEY_DONT_WARN_ON_DANGEROUS_USAGE=DontWarn
ENV DEBIAN_FRONTEND=noninteractive
USER root
RUN apt-get update && apt-get install -yq --no-install-recommends \
    git wget \
    ca-certificates \
    sudo python3 python3-pip
# faSomeRecords and faSize are needed for the UShER WDL workflow 
WORKDIR /HOME/kentsource
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faSomeRecords
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faSize
RUN chmod 775 *
WORKDIR /HOME
RUN git clone https://github.com/yatisht/usher.git 
WORKDIR usher
## Checkout latest release
#RUN git checkout $(git describe --tags `git rev-list --tags --max-count=1`)
RUN ./install/installUbuntu.sh 
## set the path
ENV PATH="/HOME/usher/build:/HOME/kentsource:${PATH}"

# add the TB reference used by clockwork-plus and others (not committed to this repo since it's a bit big)
# see gs://ucsc-pathogen-genomics-public/tb/ref/clockwork, md5sum should be fca996be5de559f5f9f789c715f1098b
RUN mkdir ref
COPY ./Ref.H37Rv.tar ./ref/
RUN cd ./ref/ && tar -xvf Ref.H37Rv.tar

# add the "known lineage" SRA tree
RUN mkdir example_tree
COPY ./data/tb_alldiffs_mask2ref.L.fixed.pb ./example_tree/