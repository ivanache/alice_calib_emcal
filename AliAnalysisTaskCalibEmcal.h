#ifndef ALIANALYSISTASKCALIBEMCAL_H_
#define ALIANALYSISTASKCALIBEMCAL_H_

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
	TH2D *_histogram_cell_id_amplitude; //!
	AliEMCALRecoUtils *_reco_utils; //!
	AliEMCALGeometry *_emcal_geometry; //!
public:
	AliAnalysisTaskCalibEmcal(void);
	AliAnalysisTaskCalibEmcal(const char *name);
	~AliAnalysisTaskCalibEmcal(void);
	virtual void UserCreateOutputObjects(void);
	virtual void UserExec(Option_t *);
	AliEMCALRecoUtils *GetEMCALRecoUtils(void);
	ClassDef(AliAnalysisTaskCalibEmcal, 1);
};

#endif // ALIANALYSISTASKCALIBEMCAL_H_
