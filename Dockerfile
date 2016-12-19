FROM python:2.7-slim
RUN mkdir -p /opt/mdtscripts/ /inputs /outputs
WORKDIR /outputs
ADD ./ /opt/mdtscripts/
RUN cd /opt/mdtscripts && \
   pip install -r requirements.txt
ENTRYPOINT ["python","/opt/mdtscripts/vde.py"]
