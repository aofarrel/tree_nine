FROM ashedpotatoes/usher-mirror:0.6.4

# TB reference used by clockwork-plus and others (not committed to this repo since it's a bit big)
# see gs://ucsc-pathogen-genomics-public/tb/ref/clockwork-v0.12.5/, md5sum should be 07c0e4af1df72e33b3e9882204ba0702
RUN mkdir ref
COPY ./Ref.H37Rv.tar ./ref/
RUN cd ./ref/ && tar -xvf Ref.H37Rv.tar

# "known lineage" SRA tree, used for quick debugging
RUN mkdir example_tree
COPY ./data/tb_alldiffs_mask2ref.L.fixed.pb ./example_tree/

RUN apt-get install -y tree
RUN pip install six numpy ete3 pandas polars requests # six is a prereq for ete3 that doesn't get installed when installing ete3
RUN mkdir /scripts/
RUN wget https://raw.githubusercontent.com/aofarrel/diffdiff/main/diffdiff.py
RUN wget https://gist.githubusercontent.com/aofarrel/a638f2ff05f579193632f7921832a957/raw/baa77b4f6afefd78ae8b6a833121a413bd359a5e/marcs_incredible_script
RUN mv marcs_incredible_script /scripts/marcs_incredible_script.pl
RUN mv diffdiff.py /scripts/diffdiff.py