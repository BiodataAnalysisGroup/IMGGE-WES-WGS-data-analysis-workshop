FROM ubuntu:focal

RUN apt update && apt install -y python3

RUN mkdir /home/test

COPY ./hello.py /home/test

WORKDIR /home/test

CMD ["python3","hello.py"]
