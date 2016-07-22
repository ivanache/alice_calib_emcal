#!/usr/bin/gawk {system("root -b -q " FILENAME); exit;}
// -*- mode: c++; -*-

void runCalibEmcal(const char *run_mode = "full",
				   const int lhc_run = 1)
{
	gROOT->ProcessLine(".include $ROOTSYS/include");
	gROOT->ProcessLine(".include $ALICE_ROOT/include");

	// Load base root libraries
	gSystem->Load("libTree");
	gSystem->Load("libGeom");
	gSystem->Load("libVMC");
	gSystem->Load("libSTEERBase");
	gSystem->Load("libPhysics");
	gSystem->Load("libESD");
	gSystem->Load("libAOD");

	// Load analysis framework libraries
	gSystem->Load("libANALYSIS");
	gSystem->Load("libANALYSISalice");
	gSystem->Load("libEMCALUtils");
	gSystem->Load("libPWGPPEMCAL");

	gROOT->ProcessLine(".L AliAnalysisTaskCalibEmcal.cxx+g");

	AliAnalysisManager *mgr = new AliAnalysisManager();

	gROOT->Macro("$ALICE_ROOT/ANALYSIS/macros/train/"
				 "AddESDHandler.C");

	AliAnalysisAlien *plugin =
		new AliAnalysisAlien("pluginCalibEmcal");

	// plugin->AddIncludePath("-I. -I$ROOTSYS/include -I$ALICE_ROOT "
	// "-I$ALICE_ROOT/include -I$ALICE_PHYSICS/include");

	plugin->SetGridWorkingDir("workdir");
	plugin->SetGridOutputDir("outputdir");
	plugin->SetAliPhysicsVersion("v5-08-14-01-1");
	plugin->SetAdditionalLibs(
		"AliAnalysisTaskCalibEmcal.h "
		"AliAnalysisTaskCalibEmcal.cxx");
	plugin->SetAnalysisSource("AliAnalysisTaskCalibEmcal.cxx");
	plugin->SetRunPrefix("000");

	const int run_number_lhc12b[] = {

		177580, 177592, 177597, 177601, 177624, 177681, 177798,
		177799, 177804, 177805, 177810, 177858, 177860, 177861,
		177864, 177869, 177932, 177938, 177942, 178018, 178024,
		178025, 178026, 178028, 178029, 178030, 178031, 178052,
		178053, 178163, 178167, 178220,

		-1
	};

	const int run_number_lhc15j[] = {

		236967, 236970, 237003, 237030, 237050, 237106, 237109,
		237119, 237178, 237257, 237287, 237511, 237513, 237516,
		237646, 237673, 237686, 237700, 237701, 237703, 237704,
		237709, 237712, 237766, 237773, 237778, 237781, 237788,
		237792, 237794, 237805, 237843, 237846,

		-1
	};

	const int *run_number;

	switch (lhc_run) {
	case 1:
		plugin->SetGridDataDir("/alice/data/2012/LHC12b");
		plugin->SetDataPattern("pass1/*/AliESDs.root");
		run_number = run_number_lhc12b;
		break;
	case 2:
		plugin->SetGridDataDir("/alice/data/2015/LHC15j");
		plugin->SetDataPattern("muon_calo_pass1/*/AliESDs.root");
		run_number = run_number_lhc15j;
		break;
	}

	for (const int *r = run_number; *r != -1; r++) {
		plugin->AddRunNumber(*r);
	}

	const char *alien_close_se = gSystem->Getenv("alien_CLOSE_SE");

	if (alien_close_se != NULL) {
		const char *file = mgr->GetCommonFileName();

		plugin->SetDefaultOutputs(kFALSE);
		plugin->SetOutputFiles(Form(
			"%s@%s", file, alien_close_se));
		plugin->SetOutputArchive(Form(
			"log_archive.zip:stdout,stderr@%s "
			"root_archive.zip:%s,*.stat@%s",
			alien_close_se, file, alien_close_se));
	}

	plugin->SetRunMode(run_mode);
	mgr->SetGridHandler(plugin);
	gROOT->Macro("macros/AddAliAnalysisTaskCalibEmcal.C");
	if (mgr->InitAnalysis()) {
		mgr->StartAnalysis("grid");
	}
}
