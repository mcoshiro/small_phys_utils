A place to collect various small scripts and utilities used during physics analysis and other tasks.

## Instructions to generate Higgs to Z gamma photon corrections

For more information on these scripts, run with the `--help` flag.

~~~~bash
source set_env.sh
./compile.sh
python39 scripts/hzg/generate_photonidskim.py -i "/net/cms28/cms28r0/oshiro/tnp_tuples/Run2018*.root" -o "/net/cms28/cms28r0/oshiro/tnp_tuples/photonidskim_data_2018.root"
python39 scripts/hzg/generate_photonidskim.py -i "/net/cms28/cms28r0/oshiro/tnp_tuples/DY_LO2018.root" -o "/net/cms28/cms28r0/oshiro/tnp_tuples/photonidskim_simu_2018.root"
python39 scripts/hzg/generate_photon_preselection_weights.py -m "/net/cms28/cms28r0/oshiro/tnp_tuples/DY_LO2018.root" -d "/net/cms28/cms28r0/oshiro/tnp_tuples/Run2018*.root" -e json/photon_presel_eff2018.json -c json/photon_presel_corr2018.json
#python39 scripts/hzg/generate_photon_bkgsubweights.py -i "/net/cms28/cms28r0/oshiro/tnp_tuples/photonidskim_data_2018.root" -j json/bkg_weights_2018.json
python39 scripts/hzg/apply_photon_bkgsubweights.py -i "/net/cms28/cms28r0/oshiro/tnp_tuples/photonidskim_data_2018.root" -j json/bkg_weights_2018.json
python39 scripts/hzg/apply_photon_preweights.py -s generate -m "/net/cms28/cms28r0/oshiro/tnp_tuples/photonidskim_simu_2018.root" -d "/net/cms28/cms28r0/oshiro/tnp_tuples/photonidskim_data_2018.root" -j 2018
python39 scripts/hzg/apply_photon_preweights.py -s validate -m "/net/cms28/cms28r0/oshiro/tnp_tuples/photonidskim_simu_2018.root" -d "/net/cms28/cms28r0/oshiro/tnp_tuples/photonidskim_data_2018.root" -j 2018
python39 scripts/hzg/train_photonid_dnn.py -m "/net/cms28/cms28r0/oshiro/tnp_tuples/photonidskim_simu_2018.root" -d "/net/cms28/cms28r0/oshiro/tnp_tuples/photonidskim_data_2018.root" -r barrel -t barrel2018
python39 scripts/hzg/train_photonid_dnn.py -m "/net/cms28/cms28r0/oshiro/tnp_tuples/photonidskim_simu_2018.root" -d "/net/cms28/cms28r0/oshiro/tnp_tuples/photonidskim_data_2018.root" -r endcap -t endcap2018
python39 scripts/convert_dnn_to_cpp.py -i json/nn_model_barrel2018.h5 -s json/scaler_file_barrel2018 -o dnn_barrel2018.hpp
python39 scripts/convert_dnn_to_cpp.py -i json/nn_model_endcap2018.h5 -s json/scaler_file_endcap2018 -o dnn_endcap2018.hpp
python39 scripts/hzg/validate_photon_dnn.py -m "/net/cms28/cms28r0/oshiro/tnp_tuples/photonidskim_simu_2018.root" -d "/net/cms28/cms28r0/oshiro/tnp_tuples/photonidskim_data_2018.root" -n dnn_barrel2018.hpp -r barrel
~~~~
