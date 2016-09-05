FROM davidsblom/foam-fsi
COPY . /pipeline/source
RUN cd /pipeline/source && travis/run.sh
