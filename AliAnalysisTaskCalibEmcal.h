#ifndef ALIANALYSISTASKCALIBEMCAL_H_
#define ALIANALYSISTASKCALIBEMCAL_H_

#include <vector>
#include <TList.h>
#include <TH2D.h>
#include <AliEMCALGeometry.h>
#include <AliEMCALRecoUtils.h>
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskCalibEmcal : public AliAnalysisTaskSE {
private:
	TString _emcal_geometry_name; //!
	TString _trigger_name; //!
	TList *_list; //!
	TH1D *_histogram_run_number_event_triggered; //!
	TH2D *_histogram_cell_id_timing; //!
	TH2D *_histogram_cell_id_timing_9_gev; //!
	TH2D *_histogram_timing_amplitude; //!
	TH2D *_histogram_cell_id_amplitude; //!
	TH2D *_histogram_cell_id_amplitude_scale; //!
	AliEMCALRecoUtils *_reco_utils; //!
	AliEMCALGeometry *_emcal_geometry; //!
	std::vector<TH2D *> _emcal_pi0_calibration_factor; //!
	size_t _ncell; //!
	Double_t _timing_min; //!
	Double_t _timing_max; //!
public:
	AliAnalysisTaskCalibEmcal(void);
	AliAnalysisTaskCalibEmcal(const char *name);
	~AliAnalysisTaskCalibEmcal(void);
	virtual void UserCreateOutputObjects(void);
	Double_t emcal_scale(Int_t cell_id);
	virtual void UserExec(Option_t *);
	AliEMCALRecoUtils *GetEMCALRecoUtils(void);
	ClassDef(AliAnalysisTaskCalibEmcal, 1);
};

#endif // ALIANALYSISTASKCALIBEMCAL_H_
