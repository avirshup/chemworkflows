FROM python:2.7-slim
RUN mkdir -p /opt/chemworkflows/ /inputs /outputs
WORKDIR /outputs
ADD ./ /opt/chemworkflows/
RUN cd /opt/chemworkflows \
 && pip install -e .
ENTRYPOINT ["chemworkflow"]
