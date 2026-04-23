# version: ashedpotatoes/usher-plus:0.6.6_rev6

# Hardcoded-for-reproducibility stuff you may eventually want to update:
# * UShER v0.6.6
# * BTE commit a9631c63e7bdd246ddd61ca5ecfddc93dad038cf
# * A bunch of my cluster-related scripts (kind of; see bottom of Dockerfile)
# * Marc Perry's persistent cluster IDs perl script

# Regarding ARM (as of mid-April 2026):
# * Running on Terra (GCP) hardware is top priority; other hardwares are secondary considerations
# * This image is built on x86 hardware
# * BTE and UShER rely on tbb and as such kinda sorta don't support building on ARM (UShER has a workaround iirc)
# * Nevertheless, my testing indicates BTE and UShER within x86 Docker images run fine on ARM via Rosetta
# * polars has an ARM version, an x86 version, and a lts-x86 version
# * An x86-compiled Docker image will contain x86 polars and will **NOT** work on an ARM machine because it contains
#   instructions Rosetta doesn't currently support, but the lts-x86 version of polars SHOULD work since it lacks
#   those instructions
# * So, in summary, getting all three to play nicely on ARM seems to require polars-lts-cpu, which forces polars to use
#   older x86 instructions that Rosetta can handle

# Base image: In usher-plus:0.6.4_5 we swapped the base image from ubuntu:22.04 to miniconda3:22.11.1 because I now rely on
# BTE, which essentially requires conda to install correctly. I looked into installing miniconda on the UShER base image
# instead of vice versa, but didn't exactly spark joy, so we're starting with a miniconda base image now.
FROM continuumio/miniconda3:22.11.1

# Essentially required if you don't want conda to take a bazillion years (verbose in case it gets stuck anyway)
ENV CONDA_VERBOSITY=2
RUN conda install -n base conda-libmamba-solver
RUN conda config --set solver libmamba

WORKDIR /HOME

# Install the usual suspects
# git, wget, ca-certificates, and sudo required for UShER
RUN apt-get update 
RUN apt-get install -yq --no-install-recommends tree zip pigz git wget ca-certificates sudo

# Install BTE
# Most recent bugfix currently doesn't have a tagged commit, so we will reset to a9631c6 (most recent commit on main
# as of mid-April 2026) in case a breaking change ever gets pushed
RUN git clone --recurse-submodules https://github.com/jmcbroome/BTE.git
WORKDIR BTE
RUN git fetch origin && git reset --hard a9631c6
RUN conda env create -f bte.yml
SHELL ["conda", "run", "-n", "bte", "/bin/bash", "-c"]
RUN python3 setup.py build_ext && \
 python3 setup.py install && \
 ln -s /HOME/BTE/build/lib.linux-x86_64-cpython-310/bte.cpython-310-x86_64-linux-gnu.so /opt/conda/lib/python3.10/site-packages

# UShER seems to require this be made explicit
SHELL ["/bin/bash", "-c"]
ENV DEBIAN_FRONTEND=noninteractive
USER root

# faSomeRecords and faSize are needed for the UShER WDL workflow 
WORKDIR /HOME/kentsource
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faSomeRecords
RUN wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/faSize

# Install UShER
RUN chmod 775 *
WORKDIR /HOME
RUN git clone https://github.com/yatisht/usher.git
WORKDIR usher
RUN git checkout v0.6.6
RUN ./install/installUbuntu.sh 
ENV PATH="/HOME/usher/build:/HOME/kentsource:${PATH}"

# Python installations
# * six must be installed before ete3 or else ete3 won't work properly
# * If you swap out `polars-lts-cpu` with `polars` and compile this image on x86, the resulting
#   Docker image will not run polars on ARM via Rosetta... so don't do that! Yes, there is an
#   ARM version of polars but as noted above we can't use it without breaking UShER.
RUN pip install six numpy ete3 pandas polars-lts-cpu requests phylodm taxoniumtools

# Dump my stuff (mostly cluster-related) into one folder to maintain the illusion of me being organized
WORKDIR /HOME
RUN mkdir ash
WORKDIR ash

# H37Rv tuberculosis reference genome as processed by clockwork reference_prepare (aka refprep)
# Mirror: gs://ucsc-pathogen-genomics-public/tb/ref/clockwork-v0.11.3/Ref.H37Rv.tar
# Tarball contains an md5 file; ref.fa is 0e8f3acd2c92f20c42105e014adbfe97
RUN mkdir ref
COPY ./Ref.H37Rv.tar ./ref/
RUN cd ./ref/ && tar -xvf Ref.H37Rv.tar
WORKDIR ../ash

# Add a test tree that should ONLY be used for debugging/quick tests
# Mirror: gs://ucsc-pathogen-genomics-public/tb/tree/alldiffs_mask2ref.L.fixed.pb
# MD5: c2c2afee2ad51527376ef2c1b1e00d4f
# Provenance: This is one of the oldest versions of the TB tree Lily and Ash worked on, made up of ~7000
# SRA samples of "known" lineages (we now know some of the lineages were declared by authors/uploaders 
# incorrectly as confirmed by tree placement and TBProfiler). It was generated on an ancient version
# of Ash's myco/Tree Nine workflows and has very little QC filtering.
RUN mkdir example_tree
COPY ./data/for_debugging_only__tb_7K_noQC_diffs_mask2ref.L.fixed.pb ./example_tree/

# Add my scripts -- it is expected you're building this Dockerfile within the Tree Nine git, if
# not, get 'em here: github.com/aofarrel/tree_nine
# UShER already has a scripts folder; for clarity these ones will go in /HOME/ash
RUN mkdir scripts
COPY ./find_clusters.py /HOME/ash/scripts/
COPY ./process_clusters.py /HOME/ash/scripts/
COPY ./summarize_changes.py /HOME/ash/scripts/
COPY ./summarize_changes_alt.py /HOME/ash/scripts/
RUN wget -O scripts/extract_long_rows_and_truncate.sh https://raw.githubusercontent.com/aofarrel/tsvutils/refs/tags/0.0.1/extract_long_rows_and_truncate.sh
RUN wget -O scripts/equalize_tabs.sh https://raw.githubusercontent.com/aofarrel/tsvutils/refs/tags/0.0.1/equalize_tabs.sh
RUN wget -O scripts/diffdiff.py https://raw.githubusercontent.com/aofarrel/diffdiff/0.0.9/diffdiff.py
RUN wget -O scripts/marcs_incredible_script_v1.pl https://gist.githubusercontent.com/aofarrel/a638f2ff05f579193632f7921832a957/raw/e6cfc0bfaf34c0daf2e7c6625ed53228b35c3cbb/marcs_incredible_script_v1.pl
RUN wget -O scripts/marcs_incredible_script_v2.pl https://gist.githubusercontent.com/aofarrel/6a458634abbca4eb16d120cc6694d5aa/raw/cb9643864ea0b3c141de6145ac51f89798e77157/marcs_incredible_script_v2.pl

# Prevent a debug print firing every time container spins up
ENV CONDA_VERBOSITY=1
