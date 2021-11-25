ARG BASE
FROM zivgitlab.wwu.io/ag-ohlberger/dune-community/dune-gdt-super/$BASE

COPY . ${SUPERDIR}/${MY_MODULE}/
RUN ls -l ${SUPERDIR}/${MY_MODULE}/
RUN cd ${SUPERDIR} && \
    echo ${SRC_DCTRL} ${BLD} --only=${MY_MODULE} all && \
    ${SRC_DCTRL} ${BLD} --only=${MY_MODULE} all
