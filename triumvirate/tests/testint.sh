#!/bin/bash
# test_intg.sh
# Perform integration tests over all cases.

# ========================================================================
# Set-up
# ========================================================================

# Set script arguments.
help_parameters () {
   echo ""
   echo "Usage: $0 -r -o"
   echo -e "\t-r Recompile code"
   echo -e "\t-o Enable OpenMP"
   exit 1
}

flag_recompile="false"
flag_useomp="false"
while getopts "ro" opt
do
   case ${opt} in
      "r" ) flag_recompile="true" ;;
      "o" ) flag_useomp="true" ;;
      "?" ) help_parameters ;;
   esac
done

# Enter base directory.
cd /home/mikesw/Documents/Bispec/Triumvirate

testdir_input="triumvirate/tests/test_input"
testdir_output="triumvirate/tests/test_output"

# (Re)compile code.
if ${flag_recompile}; then
  make clean; make install useomp=${flag_useomp}
fi

# Change template scripts.
cp triumvirate/tests/testint.py triumvirate/tests/testint_temp.py
sed -i "s/test_params.yml/test_params_temp.yml/g" triumvirate/tests/testint_temp.py


# ========================================================================
# Run
# ========================================================================

test_case () {
  # Determine test case.
  local catalogue_type=${1}
  local space=${2}
  local npoint=${3}

  echo
  echo "Case: ${catalogue_type} ${space} ${npoint}"
  echo

  # Set I/O.
  local test_paramfile_ini="${testdir_input}/params/test_params.ini"
  local test_paramfile_yml="${testdir_input}/params/test_params.yml"

  local test_paramfile_ini_temp="${testdir_input}/params/test_params_temp.ini"
  local test_paramfile_yml_temp="${testdir_input}/params/test_params_temp.yml"

  if [[ $catalogue_type == 'sim' ]]; then
    local data_catalogue_file="test_catalogue_sim.dat"
  elif [[ $catalogue_type == 'survey' ]]; then
    local data_catalogue_file="test_catalogue_data.dat"
  else
    local data_catalogue_file="none"
  fi

  if [[ $catalogue_type == 'random' ]]; then
    local rand_catalogue_file="test_catalogue_rand.dat"
  elif [[ $catalogue_type == 'survey' ]]; then
    local rand_catalogue_file="test_catalogue_rand.dat"
  else
    local rand_catalogue_file="none"
  fi

  # Set parameters.
  if ${flag_useomp}; then
    local tag="omp"
  else
    local tag="nomp"
  fi

  if [[ $catalogue_type == 'sim' ]]; then
    local catalogue_columns="x,y,z,ws"
  else
    local catalogue_columns="x,y,z,nz,ws"
  fi

  if [[ $space == 'fourier' ]]; then
    local bin_min=0.005
    local bin_max=0.205
    local range="[0.005, 0.205]"
  fi
  if [[ $space == 'config' ]]; then
    local bin_min=20.
    local bin_max=820.
    local range="[20., 820.]"
  fi

  if [[ $catalogue_type == 'random' ]]; then
    if [[ $npoint == '2pt' ]]; then
      local statistic="2pcf-win"
    fi
    if [[ $npoint == '3pt' ]]; then
      local statistic="3pcf-win"
    fi
  else
    if [[ $space == 'fourier' ]]; then
      if [[ $npoint == '2pt' ]]; then
        local statistic="powspec"
      fi
      if [[ $npoint == '3pt' ]]; then
        local statistic="bispec"
      fi
    fi
    if [[ $space == 'config' ]]; then
      if [[ $npoint == '2pt' ]]; then
        local statistic="2pcf"
      fi
      if [[ $npoint == '3pt' ]]; then
        local statistic="3pcf"
      fi
    fi
  fi

  # Run.
  python3 application/gen_param_file.py ${test_paramfile_ini} -i \
    -o ${test_paramfile_ini_temp} \
    -p data_catalogue_file ${data_catalogue_file} \
    -p rand_catalogue_file ${rand_catalogue_file} \
    -p catalogue_columns ${catalogue_columns} \
    -p output_tag _${catalogue_type}_cpp${tag} \
    -p catalogue_type ${catalogue_type} \
    -p statistic_type ${statistic} \
    -p bin_min ${bin_min} \
    -p bin_max ${bin_max} \
    -p measurement_dir ${testdir_output} \
    -p num_bins 20

  build/triumvirate ${test_paramfile_ini_temp}

  rm ${test_paramfile_ini_temp}

  python3 application/gen_param_file.py ${test_paramfile_yml} -y \
    -o ${test_paramfile_yml_temp} \
    -p files.data_catalogue ${data_catalogue_file} \
    -p files.rand_catalogue ${rand_catalogue_file} \
    -p tags.output _${catalogue_type}_py${tag} \
    -p catalogue_type ${catalogue_type} \
    -p statistic_type ${statistic} \
    -p range "${range}" \
    -p directories.measurements ${testdir_output} \
    -p num_bins 20

  python3 triumvirate/tests/testint_temp.py

  rm ${test_paramfile_yml_temp}
}

test_case sim fourier 2pt
test_case sim config 2pt
test_case sim fourier 3pt
test_case sim config 3pt

test_case survey fourier 2pt
test_case survey config 2pt
test_case survey fourier 3pt
test_case survey config 3pt

test_case random config 2pt
test_case random config 3pt


# ========================================================================
# Clean-up
# ========================================================================

rm triumvirate/tests/testint_temp.py
