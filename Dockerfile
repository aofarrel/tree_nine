FROM ashedpotatoes/usher-mirror:0.6.4

# TB reference - H37Rv
# md5sum of ref.fa is 0e8f3acd2c92f20c42105e014adbfe97
# taken from clockwork ref prep, but this is only H37Rv stuff -- no decontamination reference here!
RUN mkdir ref
COPY ./Ref.H37Rv.tar ./ref/
RUN cd ./ref/ && tar -xvf Ref.H37Rv.tar

# add the "known lineage SRA tree", which should ONLY be used for debugging, folks!
RUN mkdir example_tree
COPY ./data/for_debugging_only__tb_7K_noQC_diffs_mask2ref.L.fixed.pb ./example_tree/

RUN apt-get install -y tree
RUN pip install six numpy ete3 pandas polars requests # six is a prereq for ete3 that doesn't get installed when installing ete3
RUN mkdir /scripts/
RUN wget https://raw.githubusercontent.com/aofarrel/diffdiff/main/diffdiff.py
RUN wget https://gist.githubusercontent.com/aofarrel/a638f2ff05f579193632f7921832a957/raw/baa77b4f6afefd78ae8b6a833121a413bd359a5e/marcs_incredible_script
RUN mv marcs_incredible_script /scripts/marcs_incredible_script.pl
RUN mv diffdiff.py /scripts/diffdiff.py