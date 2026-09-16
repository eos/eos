FROM eoshep/release-base:2026-09-17-453bfd7
ARG EOSVERSION
RUN pip3 install --no-cache-dir "eoshep==${EOSVERSION}"
