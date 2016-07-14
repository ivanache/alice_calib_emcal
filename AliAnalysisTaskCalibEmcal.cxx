#ifndef __CINT__
#include <boost/format.hpp>
#endif // __CINT__
#include <TGeoManager.h>
#include <TFile.h>
#include <AliVEvent.h>
#include <AliESDEvent.h>
#include <AliAODEvent.h>
#include <AliESDHeader.h>
#include <AliVTrack.h>
#include <AliVCluster.h>
#include <AliVCaloCells.h>
#include <AliOADBContainer.h>
#include "AliAnalysisTaskCalibEmcal.h"

namespace {

	void to_sm_ieta_iphi(unsigned int &sm, unsigned int &ieta,
						 unsigned int &iphi, unsigned int n)
	{
		sm = n < 11520 ? n / 1152 :
			n < 12288 ? 10 + (n - 11520) / 384 :
			n < 16897 ? 12 + (n - 12288) / 768 :
			18 + (n - 16896) / 384;

		const unsigned int n0 =
			sm < 10 ? sm * 1152 :
			sm < 12 ? 11520 + (sm - 10) * 384 :
			sm < 18 ? 12288 + (sm - 12) * 768 :
			16896 + (sm - 18) * 384;
		const unsigned int n1 = n - n0;
		const unsigned int nphi =
			sm < 10 ? 24 : sm < 12 ? 8 : sm < 18 ? 24 : 8;

		ieta = 2 * (n1 / (2 * nphi)) + 1 - (n1 % 2);
		iphi = (n1 / 2) % nphi;
	}

	void set_log_axis(TAxis *axis)
	{
		const double axis_range[2] = {
			log10(axis->GetXmin()), log10(axis->GetXmax())
		};
		const size_t nbin = axis->GetNbins();
		const double bin_width =
			(axis_range[1] - axis_range[0]) / nbin;
		double *bin = new double[nbin + 1];

		for (size_t i = 0; i <= nbin; i++) {
			bin[i] = pow(10, axis_range[0] + i * bin_width);
		}
		axis->Set(nbin, bin);
	}

}

ClassImp(AliAnalysisTaskCalibEmcal);

AliAnalysisTaskCalibEmcal::AliAnalysisTaskCalibEmcal(void)
	: AliAnalysisTaskSE(),
	  _emcal_geometry_name("EMCAL_COMPLETE12SMV1_DCAL_8SM"),
	  _trigger_name("CEMC7"), _list(NULL),
	  _histogram_run_number_event_triggered(NULL),
	  _histogram_cell_id_timing(NULL),
	  _histogram_cell_id_timing_9_gev(NULL),
	  _histogram_timing_amplitude(NULL),
	  _histogram_cell_id_amplitude(NULL),
	  _histogram_cell_id_amplitude_scale(NULL),
	  _reco_utils(new AliEMCALRecoUtils),
	  _emcal_geometry(NULL),
	  _emcal_pi0_calibration_factor(),
	  _ncell(17664),
	  _timing_min(300), _timing_max(800) // 560-610 ns in Run-1
{
}

AliAnalysisTaskCalibEmcal::AliAnalysisTaskCalibEmcal(const char *name)
	: AliAnalysisTaskSE(name),
	  _emcal_geometry_name("EMCAL_COMPLETE12SMV1_DCAL_8SM"),
	  _trigger_name("CEMC7"), _list(NULL),
	  _histogram_run_number_event_triggered(NULL),
	  _histogram_cell_id_timing(NULL),
	  _histogram_cell_id_timing_9_gev(NULL),
	  _histogram_timing_amplitude(NULL),
	  _histogram_cell_id_amplitude(NULL),
	  _histogram_cell_id_amplitude_scale(NULL),
	  _reco_utils(new AliEMCALRecoUtils),
	  _emcal_geometry(NULL),
	  _emcal_pi0_calibration_factor(),
	  _ncell(17664),
	  _timing_min(300), _timing_max(800) // 560-610 ns in Run-1
{
	DefineOutput(1, TList::Class());
}

AliAnalysisTaskCalibEmcal::~AliAnalysisTaskCalibEmcal(void)
{
}

void AliAnalysisTaskCalibEmcal::UserCreateOutputObjects(void)
{
	_list = new TList();
	_list->SetOwner(kTRUE);

	/////////////////////////////////////////////////////////////////

	static const unsigned int run_begin = 236892;
	static const unsigned int run_end = 238622;

	_histogram_run_number_event_triggered =
		new TH1D("_histogram_run_number_event_triggered", "",
				 run_end - run_begin + 1, run_begin - 0.5,
				 run_end + 0.5);

	_list->Add(_histogram_run_number_event_triggered);

	/////////////////////////////////////////////////////////////////

	_histogram_cell_id_timing =
		new TH2D("_histogram_cell_id_timing", "",
				 _ncell, -0.5, _ncell - 0.5,
				 1000, 0, 1000);

	_list->Add(_histogram_cell_id_timing);

	/////////////////////////////////////////////////////////////////

	_histogram_cell_id_timing_9_gev =
		new TH2D("_histogram_cell_id_timing_9_gev", "",
				 _ncell, -0.5, _ncell - 0.5,
				 1000, 0, 1000);

	_list->Add(_histogram_cell_id_timing_9_gev);

	/////////////////////////////////////////////////////////////////

	_histogram_timing_amplitude =
		new TH2D("_histogram_timing_amplitude", "",
				 1000, 0, 1000,
				 1000, 0.1, 100);

	_list->Add(_histogram_timing_amplitude);

	/////////////////////////////////////////////////////////////////

	_histogram_cell_id_amplitude =
		new TH2D("_histogram_cell_id_amplitude", "",
				 _ncell, -0.5, _ncell - 0.5,
				 1000, 0.1, 100);
	set_log_axis(_histogram_cell_id_amplitude->GetYaxis());

	_list->Add(_histogram_cell_id_amplitude);

	/////////////////////////////////////////////////////////////////

	_histogram_cell_id_amplitude_scale =
		new TH2D("_histogram_cell_id_amplitude_scale", "",
				 _ncell, -0.5, _ncell - 0.5,
				 1000, 0.1, 100);
	set_log_axis(_histogram_cell_id_amplitude_scale->GetYaxis());

	_list->Add(_histogram_cell_id_amplitude_scale);

	/////////////////////////////////////////////////////////////////

	PostData(1, _list);
}

Double_t AliAnalysisTaskCalibEmcal::emcal_scale(Int_t cell_id)
{
	unsigned int sm;
	unsigned int ieta;
	unsigned int iphi;

	to_sm_ieta_iphi(sm, ieta, iphi, cell_id);

	if (_emcal_pi0_calibration_factor[sm] != NULL) {
		const Int_t root_bin =
			_emcal_pi0_calibration_factor[sm]->
			GetBin(ieta + 1, iphi + 1);

		return _emcal_pi0_calibration_factor[sm]->
			GetBinContent(root_bin);
	}
	else {
		return 1;
	}
}

void AliAnalysisTaskCalibEmcal::UserExec(Option_t *option)
{
	AliVEvent *event = InputEvent();
    AliESDEvent *esd_event = dynamic_cast<AliESDEvent *>(event);
    AliAODEvent *aod_event = dynamic_cast<AliAODEvent *>(event);

	const Int_t run_number = event->GetRunNumber();

	if (_emcal_geometry == NULL) {
#if 1
		TGeoManager::Import(
			"$ALICE_PHYSICS/OADB/EMCAL/geometry_2015.root");
#endif
		_emcal_geometry =
			AliEMCALGeometry::GetInstance(_emcal_geometry_name);
	}

#if 1
	if (_emcal_geometry != NULL) {
		for (size_t i = 0; i < _ncell; i++) {
			_emcal_geometry->PrintCellIndexes(i, 1, "");
		}
	}
#endif

	if (_emcal_pi0_calibration_factor.empty()) {
		TFile *file = TFile::Open(
			"multiplyPi0CalibrationFactors_TextToHisto_"
			"EMCALcoeffs2012pizOnlyNoEspectra_"
			"DCALandThirdsAllOne.root");

		if (file != NULL) {
			for (size_t i = 0; i < 20; i++) {
				_emcal_pi0_calibration_factor.push_back(
					static_cast<TH2D *>(file->Get(
						(boost::format{"EMCALRecalFactors_SM%d"} % i).
						str().c_str())));
			}
		}
	}

    TString trigger_class =
		esd_event ? esd_event->GetFiredTriggerClasses() :
		aod_event ? aod_event->GetFiredTriggerClasses() :
		"";

	if (!trigger_class.Contains(_trigger_name)) {
		return;
	}

	_histogram_run_number_event_triggered->Fill(run_number);

	TRefArray *calo_cluster = new TRefArray();

	event->GetEMCALClusters(calo_cluster);

	AliVCaloCells *emcal_cell = event->GetEMCALCells();

	for (Int_t i = 0; i < calo_cluster->GetEntriesFast() - 1; i++) {
		AliVCluster *c = (AliVCluster *)calo_cluster->At(i);
#if 1
		if (_emcal_geometry != NULL &&
			_reco_utils->ClusterContainsBadChannel(
				_emcal_geometry, c->GetCellsAbsId(),
				c->GetNCells())) {
			continue;
		}
#endif

		const Double_t timing = c->GetTOF() * 1e+9;

		for (Int_t j = 0; j < c->GetNCells(); j++) {
			const Int_t cell_id = c->GetCellsAbsId()[j];
			const Double_t amplitude =
				emcal_cell->GetCellAmplitude(cell_id);

			_histogram_cell_id_timing->Fill(cell_id, timing);
			if (amplitude >= 9.0) {
				_histogram_cell_id_timing_9_gev->Fill(
					cell_id, timing);
			}
			_histogram_timing_amplitude->Fill(timing, amplitude);
		}

		if (!(timing >= _timing_min && timing < _timing_max)) {
			continue;
		}

		for (Int_t j = 0; j < c->GetNCells(); j++) {
			const Int_t cell_id = c->GetCellsAbsId()[j];
			const Double_t amplitude =
				emcal_cell->GetCellAmplitude(cell_id);
			const Double_t scale = 1;

			_histogram_cell_id_amplitude->Fill(cell_id, amplitude);
			_histogram_cell_id_amplitude_scale->Fill(
				cell_id, amplitude * scale);
		}
	}
}

AliEMCALRecoUtils *AliAnalysisTaskCalibEmcal::GetEMCALRecoUtils(void)
{
	return _reco_utils;
}
