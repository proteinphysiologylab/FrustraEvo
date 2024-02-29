FROM proteinphysiologylab/frustraevo:latest

RUN rm -r /root/FrustraEvo 

WORKDIR /opt

RUN git clone https://github.com/FranceCosta/FrustraEvo.git