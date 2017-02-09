// #include "AliAnalysisTaskCalibEmcal.h"

AliAnalysisTaskCalibEmcal *AddAliAnalysisTaskCalibEmcal(TString name = "AliAnalysisTaskCalibEmcal")
{
	AliAnalysisManager *mgr = AliAnalysisManager::GetAnalysisManager();
	AliAnalysisTaskCalibEmcal *task =
		new AliAnalysisTaskCalibEmcal(name.Data());

	AliEMCALRecoUtils *reco_util = task->GetEMCALRecoUtils();
  
	gROOT->LoadMacro("$ALICE_PHYSICS/PWGPP/EMCAL/macros/ConfigureEMCALRecoUtils.C");
  
	ConfigureEMCALRecoUtils(reco_util, kFALSE, kTRUE, kTRUE,
							kFALSE, kFALSE, kFALSE); 

	reco_util->SetNumberOfCellsFromEMCALBorder(0);
	reco_util->SwitchOnRecalibration();
	reco_util->SwitchOnRunDepCorrection();

	TString filename = mgr->GetCommonFileName();

	filename += ":AliAnalysisTaskCalibEmcal";

	mgr->AddTask(task);
	mgr->ConnectInput(task, 0, mgr->GetCommonInputContainer());
	mgr->ConnectOutput(task, 1, mgr->CreateContainer("histogram", TList::Class(), AliAnalysisManager::kOutputContainer, filename.Data()));

	return task;
}
