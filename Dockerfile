
FROM python:3.8-slim-buster

RUN apt-get update
RUN apt-get install -y gcc python3-dev python3-geopandas
COPY . /SpacialModulesPrueba/
WORKDIR /SpacialModulesPrueba
ADD ./tables /tables
COPY requirements.txt .
RUN pip3 install -r requirements.txt

RUN pip install scipy
RUN pip install statsmodels
ENV PYTHONPATH /SpacialModulesPrueba

CMD python ./examples/main.py