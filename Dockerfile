FROM python:stretch
# Paired end trimmer for QIASeq DNA/RNA reads
#
# This Dockerfile will create an image that will run
# a git pull & install for each invocation of the
# trimming tool. That is to ensure the latest code
# is being used.

ARG repository=qiaseq

# Installation of requirements
RUN pip3 install edlib Cython

# Initial installation of the trimming code
RUN git clone https://github.com/${repository}/read-trimmer && cd read-trimmer && python3 setup.py build_ext --inplace

# install the wrapper script that will run updates before the actual trimming
RUN echo '#!/bin/bash\n\
echo "# [$(date)] Running code updates ..."\n\
cd /read-trimmer\n\
git pull origin master\n\
python3 setup.py build_ext --inplace\n\
echo "# [$(date)] Code updated."\n\
python3 /read-trimmer/trimmer/run.py "$@"' > /trim

RUN chmod +x /trim

ENTRYPOINT [ "/trim" ]
CMD [ "--help" ]

