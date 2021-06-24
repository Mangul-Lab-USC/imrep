FROM ubuntu:18.04

RUN apt-get update --fix-missing
RUN apt-get install -y build-essential \
    wget \
    bzip2 \
    samtools \
    libbam-dev \
    libhts-dev \
    python-minimal

RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/*

ENV PATH="/root/miniconda3/bin:${PATH}"
ARG PATH="/root/miniconda3/bin:${PATH}"

RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && mkdir /root/.conda \
    && bash Miniconda3-latest-Linux-x86_64.sh -b \
    && rm -f Miniconda3-latest-Linux-x86_64.sh

RUN pip install --upgrade pip

WORKDIR /imrep
COPY . .

RUN cd suffix_tree && python ./setup.py install --user

RUN pip install -r requirements.txt

ENTRYPOINT ["python", "./imrep.py"]