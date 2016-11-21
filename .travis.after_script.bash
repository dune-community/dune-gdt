#!/bin/bash

git config --global hooks.clangformat ${CLANG_FORMAT}
PYTHONPATH=${SUPERDIR}/scripts/python/ python3 -c "import travis_report as tp; tp.clang_format_status(\"${TRAVIS_BUILD_DIR}\")"
${SRC_DCTRL} ${BLD} --only=${MY_MODULE} configure
${SRC_DCTRL} ${BLD} --only=${MY_MODULE} make doc
${SUPERDIR}/.travis/init_sshkey.sh ${encrypted_95fb78800815_key} ${encrypted_95fb78800815_iv} keys/dune-community/dune-community.github.io
${SUPERDIR}/.travis/deploy_docs.sh ${MY_MODULE} ${DUNE_BUILD_DIR}
