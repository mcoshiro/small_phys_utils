/**
 * function to train a BDT
 */
void train_bdt() {
  //abduct old code from ttZ tools
  //load TMVA
  TMVA::Tools::Instance();
	//TFile* sig_input_file = TFile::Open("ntuples/train_boost_sig.root");
	//TFile* bak_input_file = TFile::Open("ntuples/train_boost_bak.root");
	//TFile* output_file = TFile::Open("ntuples/bdt_output.root","RECREATE");
	TFile* sig_input_file = TFile::Open("ntuples/train_kinbdt_idmva_nomasscut_sig.root");
	TFile* bak_input_file = TFile::Open("ntuples/train_kinbdt_idmva_nomasscut_bak.root");
	TFile* output_file = TFile::Open("ntuples/bdt_idmva_ptbias_output.root","RECREATE");

	TMVA::Factory *bdt_factory = new TMVA::Factory("zg_idmva_ptbias_bdt",output_file,"!V:ROC:!Correlations:!Silent:Color:!DrawProgressBar:AnalysisType=Classification");
	TMVA::DataLoader bdt_loader("dataset");
	//bdt_loader.AddVariable("higgspt",'F');
	//bdt_loader.AddVariable("higgsdr",'F');
	//bdt_loader.AddVariable("zpt",'F');
	//bdt_loader.AddVariable("phpt",'F');
	bdt_loader.AddVariable("photon_mva",'F');
	bdt_loader.AddVariable("min_dR",'F');
	bdt_loader.AddVariable("max_dR",'F');
	bdt_loader.AddVariable("pt_mass",'F');
	bdt_loader.AddVariable("cosTheta",'F');
	bdt_loader.AddVariable("costheta",'F');
	bdt_loader.AddVariable("phi",'F');
	bdt_loader.AddVariable("photon_res",'F');
	bdt_loader.AddVariable("photon_rapidity",'F');
	bdt_loader.AddVariable("l1_rapidity",'F');
	bdt_loader.AddVariable("l2_rapidity",'F');
	bdt_loader.AddVariable("photon_ptransverse",'F');
	//bdt_loader.AddVariable("photon_pt_mass",'F');
	//bdt_loader.AddVariable("decorr_photon_pt",'F');
  bdt_loader.SetBackgroundWeightExpression("weight");
  bdt_loader.SetSignalWeightExpression("weight");

	bdt_loader.AddSignalTree(static_cast<TTree*>(sig_input_file->Get("tree")));
	bdt_loader.AddBackgroundTree(static_cast<TTree*>(bak_input_file->Get("tree")));

	TCut cut_s, cut_b;
	bdt_loader.PrepareTrainingAndTestTree(cut_s,cut_b,"nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V");
	bdt_factory->BookMethod(&bdt_loader,TMVA::Types::kBDT,"BDT","!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20");
  //bdt_factory->BookMethod(&bdt_loader,TMVA::Types::kBDT,"BDTG","!H:!V:NTrees=2000::BoostType=Grad:Shrinkage=0.1:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=3:MaxDepth=4" );
	bdt_factory->TrainAllMethods();
	bdt_factory->TestAllMethods();
	bdt_factory->EvaluateAllMethods();
	output_file->Close();
	delete bdt_factory;
	//TMVA::TMVAGui("bdt_output.root");
	return 0;
}