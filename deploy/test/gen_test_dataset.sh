#!/usr/bin/env bash
#
# @file gen_test_dataset.sh
# @author Mike S Wang
# @brief Automate test dataset generation.
#

TESTPROG="python deploy/test/gen_test_dataset.py"

# Test parameters
KRANGE="0.005 0.105"
RRANGE="50. 150."

# 2-point statistics
MULTIPOLES="0 2"
for degree in ${MULTIPOLES}
do
    ${TESTPROG} --output-tag _gpp \
        --case sim --statistic 'powspec' \
        --multipole $degree --range ${KRANGE}
    ${TESTPROG} --output-tag _gpp \
        --case sim --statistic '2pcf' \
        --multipole $degree --range ${RRANGE}
    ${TESTPROG} --output-tag _lpp \
        --case survey --statistic 'powspec' \
        --multipole $degree --range ${KRANGE}
    ${TESTPROG} --output-tag _lpp \
        --case survey --statistic '2pcf' \
        --multipole $degree --range ${RRANGE}
    ${TESTPROG} \
        --case random --statistic '2pcf-win' \
        --multipole $degree --range ${RRANGE}
done

# 3-point statistics

## Diagonal form
MULTIPOLES="000 202"
for degrees in ${MULTIPOLES}
do
    ${TESTPROG} --output-tag _gpp \
        --case sim --statistic 'bispec' \
        --multipole $degrees --range ${KRANGE}
    ${TESTPROG} --output-tag _gpp \
        --case sim --statistic '3pcf' \
        --multipole $degrees --range ${RRANGE}
    ${TESTPROG} --output-tag _lpp \
        --case survey --statistic 'bispec' \
        --multipole $degrees --range ${KRANGE}
    ${TESTPROG} --output-tag _lpp \
        --case survey --statistic '3pcf' \
        --multipole $degrees --range ${RRANGE}
    ${TESTPROG} \
        --case random --statistic '3pcf-win' \
        --multipole $degrees --range ${RRANGE}
done

## Row form (0th bin)
FORM='row'

MULTIPOLES="000"
for degrees in ${MULTIPOLES}
do
    ${TESTPROG} --output-tag _gpp \
        --case sim --statistic 'bispec' \
        --multipole $degrees --range ${KRANGE} --form ${FORM}
    ${TESTPROG} --output-tag _gpp \
        --case sim --statistic '3pcf' \
        --multipole $degrees --range ${RRANGE} --form ${FORM}
    ${TESTPROG} --output-tag _lpp \
        --case survey --statistic 'bispec' \
        --multipole $degrees --range ${KRANGE} --form ${FORM}
    ${TESTPROG} --output-tag _lpp \
        --case survey --statistic '3pcf' \
        --multipole $degrees --range ${RRANGE} --form ${FORM}
    ${TESTPROG} \
        --case random --statistic '3pcf' \
        --multipole $degrees --range ${RRANGE} --form ${FORM}
done
