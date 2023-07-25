FROM python:stretch
# Paired end trimmer for QIASeq DNA/RNA reads
#
# This Dockerfile will create an image that will run
# a git pull & install for each invocation of the
# trimming tool. That is to ensure the latest code
# is being used.

ARG repository=aaron-obrien
ARG branch=master

# Installation of requirements
RUN pip3 install edlib Cython

# Initial installation of the trimming code
RUN git clone -b ${branch} https://github.com/${repository}/read-trimmer

# install the wrapper script that will run updates before the actual trimming
RUN echo '#!/bin/bash\n\
echo "# [$(date)] Running code updates ..."\n\
cd /read-trimmer\n\
git pull origin ${branch}\n\
echo "# [$(date)] Code updated."\n\
python3 /read-trimmer/trimmer/setup.py\n\
python3 /read-trimmer/trimmer/run.py "$@"' > /trim

RUN chmod +x /trim

ENTRYPOINT [ "/trim" ]
CMD [ "--help" ]