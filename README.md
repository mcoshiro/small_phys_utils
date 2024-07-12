A place to collect various small scripts and utilities used during physics analysis and other tasks.

## Instructions to generate Higgs to Z gamma photon corrections

For more information on these scripts, run with the `--help` flag.

~~~~bash
source set_env.sh
./compile.sh
python39 scripts/hzg/generate_photon_preselection_weights.py -m "/data2/oshiro/ntuples/2018/DY_LO2018.root" -d "/data2/oshiro/ntuples/2018/Run2018*.root" -e json/photon_presel_eff2018.json -c json/photon_presel_corr2018.json
python39 scripts/hzg/plot_photon_preselection_weights.py -j json/photon_presel_corr2018.json
python39 scripts/hzg/generate_photonidskim.py -i "/data2/oshiro/ntuples/2018/Run2018*.root" -o "/data2/oshiro/ntuples/2018/photonidskim_data_2018.root"
python39 scripts/hzg/generate_photonidskim.py -i "/data2/oshiro/ntuples/2018/DY_LO.root" -o "/data2/oshiro/ntuples/2018/photonidskim_simu_2018.root"
#python39 scripts/hzg/generate_photon_bkgsubweights.py -i "/data2/oshiro/ntuples/2018/photonidskim_data_2018.root" -j json/bkg_weights_2018.json
python39 scripts/hzg/apply_photon_bkgsubweights.py -i "/data2/oshiro/ntuples/2018/photonidskim_data_2018.root" -j json/bkg_weights_2018.json
python39 scripts/hzg/apply_photon_preweights.py -s generate -m "/data2/oshiro/ntuples/2018/photonidskim_simu_2018.root" -d "/data2/oshiro/ntuples/2018/photonidskim_data_2018.root" -j 2018
python39 scripts/hzg/apply_photon_preweights.py -s validate -m "/data2/oshiro/ntuples/2018/photonidskim_simu_2018.root" -d "/data2/oshiro/ntuples/2018/photonidskim_data_2018.root" -j 2018
python39 scripts/hzg/train_photonid_dnn.py -m "/data2/oshiro/ntuples/2018/photonidskim_simu_2018.root" -d "/data2/oshiro/ntuples/2018/photonidskim_data_2018.root" -r barrel -t barrel2018
python39 scripts/hzg/train_photonid_dnn.py -m "/data2/oshiro/ntuples/2018/photonidskim_simu_2018.root" -d "/data2/oshiro/ntuples/2018/photonidskim_data_2018.root" -r endcap -t endcap2018
python39 scripts/convert_dnn_to_cpp.py -i json/nn_model_barrel2018.h5 -s json/scaler_file_barrel2018 -o dnn_barrel2018.hpp
python39 scripts/convert_dnn_to_cpp.py -i json/nn_model_endcap2018.h5 -s json/scaler_file_endcap2018 -o dnn_endcap2018.hpp
python39 scripts/hzg/validate_photon_dnn.py -m "/data2/oshiro/ntuples/2018/photonidskim_simu_2018.root" -d "/data2/oshiro/ntuples/2018/photonidskim_data_2018.root" -n dnn_barrel2018.hpp -r barrel
python39 scripts/hzg/validate_photon_dnn.py -m "/data2/oshiro/ntuples/2018/photonidskim_simu_2018.root" -d "/data2/oshiro/ntuples/2018/photonidskim_data_2018.root" -n dnn_endcap2018.hpp -r endcap
python39 scripts/hzg/generate_photonid_dnn_syst.py -m "/net/cms11/cms11r0/pico/NanoAODv9UCSB1/htozgamma_kingscanyon_v1/2018/mc/merged_zgmc_llg/merged_pico_llg_ZGToLLG_01J_5f_lowMLL_lowGPt_TuneCP5_13TeV-amcatnloFXFX-pythia8_zgmc_llg_nfiles_77.root" -d "/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/2018/data/merged_zgdata_llg/merged_raw_pico_llg_DoubleMuon_zgdata_llg_nfiles_206.root" -n dnn_barrel2018.hpp -r barrel -s make_plots
python39 scripts/hzg/generate_photonid_dnn_syst.py -m "/net/cms11/cms11r0/pico/NanoAODv9UCSB1/htozgamma_kingscanyon_v1/2018/mc/merged_zgmc_llg/merged_pico_llg_ZGToLLG_01J_5f_lowMLL_lowGPt_TuneCP5_13TeV-amcatnloFXFX-pythia8_zgmc_llg_nfiles_77.root" -d "/net/cms11/cms11r0/pico/NanoAODv9/htozgamma_kingscanyon_v1/2018/data/merged_zgdata_llg/merged_raw_pico_llg_DoubleMuon_zgdata_llg_nfiles_206.root" -n dnn_endcap2018.hpp -r endcap -s make_plots
~~~~
