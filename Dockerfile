FROM debian:unstable
RUN apt-get clean && apt-get update && apt-get install -y cppcheck
COPY . /pipeline/source
RUN cd /pipeline/source/src/fsi && cppcheck --error-exitcode=1 --enable=performance,portability .
RUN cd /pipeline/source/src/RBFMeshMotionSolver && cppcheck --error-exitcode=1 --enable=performance,portability .
RUN cd /pipeline/source/src/tests && cppcheck --error-exitcode=1 --enable=performance,portability .
RUN cd /pipeline/source/applications && cppcheck --error-exitcode=1 --enable=performance,portability .
