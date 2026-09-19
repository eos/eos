FROM eoshep/release-base:2026-09-17-453bfd7
ARG EOSVERSION
RUN pip3 install eoshep==${EOSVERSION}
