// -*- mode: c++; -*-

#ifndef ALIANALYSISTASKCALIBEMCAL_H_
#define ALIANALYSISTASKCALIBEMCAL_H_

#include <vector>
#include <deque>
#include <algorithm>
#include <TList.h>
#include <TH2D.h>
#include <AliEMCALGeometry.h>
#include <AliEMCALRecoUtils.h>
#include "AliAnalysisTaskSE.h"

class AliAnalysisTaskCalibEmcal : public AliAnalysisTaskSE {
private:
	TString _emcal_geometry_name; //!
	TString _trigger_name; //!
	std::vector<std::deque<std::vector<std::pair<TLorentzVector, double> > > > _mixed_event; //!
	TList *_list; //!

	TH1D *_histogram_run_number_event_triggered; //!
	TH1D *_histogram_multiplicity_0_5; //!
	TH2D *_histogram_cell_id_timing; //!
	TH2D *_histogram_cell_id_timing_9_gev; //!
	TH2D *_histogram_timing_amplitude; //!
	TH2D *_histogram_cell_id_amplitude; //!
	TH2D *_histogram_cell_id_amplitude_scale; //!
	TH2D *_histogram_cell_id_amplitude_low_mult; //!
	TH2D *_histogram_cell_id_amplitude_high_mult; //!
	TH2D *_histogram_cell_id_pi0_mass; //!
	TH2D *_histogram_cell_id_pi0_mass_mixed; //!
	TH2D *_histogram_cell_id_pi0_mass_recal; //!
	TH2D *_histogram_sm_ieta_pi0_mass; //!
	TH2D *_histogram_sm_ieta_pi0_mass_mixed; //!
	TH2D *_histogram_sm_ieta_pi0_mass_recal; //!
	TH2D *_histogram_sm_iphi_pi0_mass; //!
	TH2D *_histogram_sm_iphi_pi0_mass_mixed; //!
	TH2D *_histogram_sm_iphi_pi0_mass_recal; //!
	TH2D *_histogram_sm_ieta_pi0_mass_same; //!
	TH2D *_histogram_sm_ieta_pi0_mass_same_recal; //!
	TH2D *_histogram_sm_ieta_pi0_mass_cross; //!
	TH2D *_histogram_sm_ieta_pi0_mass_cross_recal; //!
	TH2D *_histogram_pi0_deta_dphi; //!
	TH2D *_histogram_pi0_deta_dphi_ieta_even; //!
	TH2D *_histogram_pi0_deta_dphi_ieta_odd; //!
	TH2D *_histogram_pi0_e_geom_mean_dr; //!
	TH2D *_histogram_pi0_e_geom_mean_dr_ieta_even; //!
	TH2D *_histogram_pi0_e_geom_mean_dr_ieta_odd; //!
	TH2D *_histogram_pi0_correlated_cell_id_e; //!

	AliEMCALRecoUtils *_reco_util; //!
	AliEMCALGeometry *_emcal_geometry; //!
	std::vector<TH2D *> _emcal_pi0_calibration_factor; //!

	size_t _ncell; //!
	size_t _nsm_ieta; //!
	size_t _nsm_iphi; //!
	Double_t _timing_min; //!
	Double_t _timing_max; //!
	Double_t _pi0_e_min; //!
	Double_t _pi0_e_max; //!
	Int_t _pi0_ncell_min; //!
	Double_t _pi0_m02_min; //!
	Double_t _pi0_m02_max; //!
	Double_t _pi0_dtiming_max; //!
public:
	AliAnalysisTaskCalibEmcal(void);
	AliAnalysisTaskCalibEmcal(const char *name);
	AliAnalysisTaskCalibEmcal(const AliAnalysisTaskCalibEmcal &);
	AliAnalysisTaskCalibEmcal &operator=
		(const AliAnalysisTaskCalibEmcal &);
	~AliAnalysisTaskCalibEmcal(void);
	virtual void UserCreateOutputObjects(void);
	Double_t emcal_scale(const Int_t cell_id);
	void init_temperature_calib(const Int_t run_number);
	virtual void UserExec(Option_t *);
	AliEMCALRecoUtils *GetEMCALRecoUtils(void);
	ClassDef(AliAnalysisTaskCalibEmcal, 1);
};

#endif // ALIANALYSISTASKCALIBEMCAL_H_
