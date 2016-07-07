#!/usr/bin/gawk {system("root -b -q " FILENAME); exit;}

void runCalibEmcal(void)
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
	plugin->SetAliPhysicsVersion("v5-08-12-01-1");
	plugin->SetAdditionalLibs("AliAnalysisTaskCalibEmcal.h "
							  "AliAnalysisTaskCalibEmcal.cxx");
	plugin->SetAnalysisSource("AliAnalysisTaskCalibEmcal.cxx");
	plugin->SetGridDataDir("/alice/data/2015/LHC15j");
	plugin->SetDataPattern("muon_calo_pass1/*/AliESDs.root");
	plugin->SetRunPrefix("000");

	int run_number[] = {

236892,
#if 0

236893,
236965,
236967,
236968,
236969,
236970,
236972,
236973,
237001,
237003,
237029,
237030,
237049,
237050,
237051,
237056,
237061,
237104,
237106,
237109,
237111,
237112,
237115,
237119,
237174,
237178,
237180,
237244,
237245,
237251,
237253,
237255,
237257,
237259,
237286,
237287,
237288,
237289,
237291,
237330,
237334,
237335,
237336,
237339,
237340,
237342,
237346,
237350,
237351,
237355,
237356,
237358,
237360,
237362,
237364,
237365,
237366,
237368,
237372,
237388,
237391,
237392,
237393,
237394,
237396,
237397,
237400,
237401,
237406,
237408,
237409,
237502,
237504,
237505,
237507,
237511,
237512,
237513,
237515,
237516,
237645,
237646,
237670,
237671,
237673,
237675,
237676,
237678,
237681,
237684,
237686,
237691,
237698,
237699,
237700,
237701,
237703,
237704,
237705,
237706,
237707,
237708,
237709,
237710,
237711,
237712,
237713,
237765,
237766,
237768,
237773,
237777,
237778,
237779,
237780,
237781,
237782,
237787,
237788,
237789,
237790,
237791,
237792,
237793,
237794,
237795,
237796,
237805,
237806,
237842,
237843,
237844,
237845,
237846,
237847,
237945,
237948,
237969,
237971,
237972,
237974,
237978,
237979,
237981,
237982,
237983,
238073,
238091,
238097,
238129,
238130,
238131,
238132,
238133,
238136,
238139,
238140,
238142,
238143,
238144,
238145,
238146,
238147,
238148,
238154,
238159,
238160,
238161,
238164,
238170,
238176,
238179,
238184,
238185,
238186,
238187,
238188,
238395,
238396,
238397,
238398,
238399,
238400,
238429,
238430,
238431,
238432,
238451,
238454,
238455,
238456,
238457,
238458,
238459,
238460,
238466,
238468,
238469,
238470,
238472,
238474,
238570,
238571,
238576,
238578,
238579,
238582,
238583,
238587,
238594,
238595,
238598,
238603,
238604,
238605,
238606,
238607,
238610,
238611,
238614,
238615,
238621,
238622,

#endif
-1 };

	for (int *r = run_number; *r != -1; r++) {
		plugin->AddRunNumber(*r);
	}

	plugin->SetRunMode("terminate");
	mgr->SetGridHandler(plugin);
	gROOT->Macro("macros/AddAliAnalysisTaskCalibEmcal.C");
	if (mgr->InitAnalysis()) {
		mgr->StartAnalysis("grid");
	}
}
