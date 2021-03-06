#!/usr/bin/gawk {system("root -b -q " FILENAME); exit;}
// -*- mode: c++; -*-

#define RUN_PERIOD(x, y) ((x) << 5) + ((y) - 'a')

void runCalibEmcal(const char *run_mode = "full",
				   const int lhc_run_period = RUN_PERIOD(15, 'o'))
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

	const int run_number_lhc10h[] = {

		137135, 137161, 137162, 137230, 137231, 137232, 137235,
		137236, 137243, 137366, 137430, 137431, 137432, 137434,
		137439, 137440, 137441, 137443, 137530, 137531, 137539,
		137541, 137544, 137546, 137549, 137595, 137608, 137638,
		137639, 137685, 137686, 137691, 137692, 137693, 137704,
		137718, 137722, 137724, 137751, 137752, 137844, 137848,
		138190, 138192, 138197, 138201, 138225, 138275, 138364,
		138396, 138438, 138439, 138442, 138469, 138534, 138578,
		138582, 138583, 138621, 138624, 138638, 138652, 138653,
		138662, 138666, 138730, 138732, 138837, 138870, 138871,
		138872, 139028, 139029, 139036, 139037, 139038, 139105,
		139107, 139173, 139309, 139310, 139314, 139328, 139329,
		139360, 139437, 139438, 139465, 139503, 139505, 139507,
		139510,

		-1
	};

	const int run_number_lhc11h[] = {

		170593, 170572, 170388, 170312, 170311, 170309, 170308,
#if 1
		170270, 170269, 170268, 170230, 170228, 170207, 170204,
		170203, 170193, 170155, 170083, 169965, 169923, 169859,
		169858, 169855, 169846, 169838, 169837, 169835, 169591,
		169588, 169557, 169554, 169550, 169512, 169504, 169498,
		169475, 169419, 169417, 169415, 169411, 169238, 169160,
		169156, 169144, 169138, 169035, 168826, 168512, 168511,
		168467, 168464, 168458, 168361, 167988, 167987, 167920,
		167915,
#endif

		-1
	};

	const int run_number_lhc12b[] = {

		177580, 177592, 177597, 177601, 177624, 177681, 177798,
		177799, 177804, 177805, 177810, 177858, 177860, 177861,
		177864, 177869, 177932, 177938, 177942, 178018, 178024,
		178025, 178026, 178028, 178029, 178030, 178031, 178052,
		178053, 178163, 178167, 178220,

		-1
	};

	const int run_number_lhc15f[] = {

		224895, 224896, 224897, 224898, 224930, 224988, 224997,
		225000, 225011, 225026, 225031, 225035, 225050, 225051,
		225052, 225093, 225105, 225106, 225576, 225578, 225579,
		225580, 225582, 225586, 225587, 225702, 225709, 225716,
		225763, 225768, 225799, 225917, 226114, 226115, 226166,
		226440, 226444, 226445, 226466, 226468, 226530,
		226532,

		-1
	};

	const int run_number_lhc15i[] = {

		235716, 235812, 235840, 235887, 235894, 235899, 236149,
		236282, 236332, 236350, 236355, 236358, 236141, 236156,
		236162, 236225, 236232, 236241, 236245, 236391, 236394,
		236442, 236445, 236447, 236454, 236466, 236555, 236559,
		236819, 236820, 236823, 236849,

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

	const int run_number_lhc15o[] = {

		245683, 245700, 245702, 245705, 245829, 245831, 245833,
		245949, 245952, 245954, 246001, 246003, 246037, 246042,
		246052, 246053, 246087, 246089, 246113, 246115,

		-1
	};

	const int run_number_lhc16h[] = {
#if 1
		// Restricted list
		254604, 254606, 254607, 254621, 254629, 254630, 254632,
		254640, 254644, 254646, 254648, 254649, 254651, 254652,
		254653, 254654,
#else
		254604, 254606, 254607, 254621, 254629, 254630, 254632,
		254640, 254644, 254646, 254648, 254649, 254651, 254652,
		254653, 254654, 255249, 255251, 255252, 255253, 255255,
		255256, 255275, 255276, 255350, 255351, 255352, 255418,
		255419, 255420, 255421, 255440, 255463, 255465, 255466,
		255467,
#endif

		-1
	};

	const int *run_number;

	plugin->SetGridDataDir(Form(
		"/alice/data/%d/LHC%02d%c", 2000 + (lhc_run_period >> 5),
		lhc_run_period >> 5, (lhc_run_period & 0x1f) + 'a'));

	switch (lhc_run_period) {
	case RUN_PERIOD(10, 'h'):
		plugin->SetDataPattern("ESDs/pass2/*/AliESDs.root");
		run_number = run_number_lhc10h;
		break;
	case RUN_PERIOD(11, 'h'):
		plugin->SetDataPattern("ESDs/pass1_calo/*/AliESDs.root");
		run_number = run_number_lhc11h;
		break;
	case RUN_PERIOD(12, 'b'):
		plugin->SetDataPattern("pass1/*/AliESDs.root");
		run_number = run_number_lhc12b;
		break;
	case RUN_PERIOD(15, 'f'):
		plugin->SetDataPattern("pass1/*/AliESDs.root");
		run_number = run_number_lhc15f;
		break;
	case RUN_PERIOD(15, 'i'):
		plugin->SetDataPattern("muon_calo_pass2/*/AliESDs.root");
		run_number = run_number_lhc15i;
		break;
	case RUN_PERIOD(15, 'j'):
		plugin->SetDataPattern("muon_calo_pass2/*/AliESDs.root");
		run_number = run_number_lhc15j;
		break;
	case RUN_PERIOD(15, 'o'):
		plugin->SetDataPattern("muon_calo_pass1/*/AliESDs.root");
		run_number = run_number_lhc15o;
		break;
	case RUN_PERIOD(16, 'h'):
		plugin->SetDataPattern("muon_calo_pass1/*/AliESDs.root");
		run_number = run_number_lhc16h;
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
