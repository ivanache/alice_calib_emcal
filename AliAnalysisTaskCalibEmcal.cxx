#include <AliVEvent.h>
#include <AliESDEvent.h>
#include <AliAODEvent.h>
#include <AliVTrack.h>
#include <AliVCluster.h>
#include <AliVCaloCells.h>
#include "AliAnalysisTaskCalibEmcal.h"

ClassImp(AliAnalysisTaskCalibEmcal);

AliAnalysisTaskCalibEmcal::AliAnalysisTaskCalibEmcal(void)
	: AliAnalysisTaskSE(),
	  _emcal_geometry_name("EMCAL_COMPLETE12SMV1_DCAL_8SM"),
	  _trigger_name("CEMC7"), _list(NULL),
	  _histogram_cell_id_amplitude(NULL),
	  _reco_utils(new AliEMCALRecoUtils), _emcal_geometry(NULL)
{
	_emcal_geometry =
		AliEMCALGeometry::GetInstance(_emcal_geometry_name);

	fprintf(stderr, "%s:%d: _emcal_geometry = %p\n", __FILE__,
			__LINE__, _emcal_geometry);
}

AliAnalysisTaskCalibEmcal::AliAnalysisTaskCalibEmcal(const char *name)
	: AliAnalysisTaskSE(name),
	  _emcal_geometry_name("EMCAL_COMPLETE12SMV1_DCAL_8SM"),
	  _trigger_name("CEMC7"), _list(NULL),
	  _histogram_cell_id_amplitude(NULL),
	  _reco_utils(new AliEMCALRecoUtils), _emcal_geometry(NULL)
{
	_emcal_geometry =
		AliEMCALGeometry::GetInstance(_emcal_geometry_name);

	fprintf(stderr, "%s:%d: _emcal_geometry = %p\n", __FILE__,
			__LINE__, _emcal_geometry);

	DefineOutput(1, TList::Class());
}

AliAnalysisTaskCalibEmcal::~AliAnalysisTaskCalibEmcal(void)
{
}

void AliAnalysisTaskCalibEmcal::UserCreateOutputObjects(void)
{
	_list = new TList();
	_list->SetOwner(kTRUE);

	const size_t ncell = 65536;
	// const size_t ncell = 17632;

	_histogram_cell_id_amplitude =
		new TH2D("_histogram_cell_id_amplitude", "",
				 ncell, -0.5, ncell - 0.5,
				 1000, 0.1, 100);

	TAxis *axis = _histogram_cell_id_amplitude->GetYaxis();
	const double axis_range[2] = {
		log10(axis->GetXmin()), log10(axis->GetXmax())
	};
    const size_t nbin = axis->GetNbins();
    const double bin_width = (axis_range[1] - axis_range[0]) / nbin;
	double *bin = new double[nbin + 1];

	for (size_t i = 0; i <= nbin; i++) {
		bin[i] = pow(10, axis_range[0] + i * bin_width);
	}
    axis->Set(nbin, bin);

	_list->Add(_histogram_cell_id_amplitude);
	PostData(1, _list);
}

void AliAnalysisTaskCalibEmcal::UserExec(Option_t *option)
{
	AliVEvent *event = InputEvent();
    AliESDEvent *esd_event = dynamic_cast<AliESDEvent*>(event);
    AliAODEvent *aod_event = dynamic_cast<AliAODEvent*>(event);
    
    TString trigger_class =
		esd_event ? esd_event->GetFiredTriggerClasses() :
		aod_event ? aod_event->GetFiredTriggerClasses() :
		"";

	if (!trigger_class.Contains(_trigger_name)) {
		return;
	}

	TRefArray *calo_cluster = new TRefArray();

	event->GetEMCALClusters(calo_cluster);

	AliVCaloCells *emcal_cell = event->GetEMCALCells();

	for(Int_t i = 0; i < calo_cluster->GetEntriesFast() - 1; i++) {
		AliVCluster *c = (AliVCluster *)calo_cluster->At(i);

		if (_emcal_geometry != NULL &&
			_reco_utils->ClusterContainsBadChannel(
				_emcal_geometry, c->GetCellsAbsId(),
				c->GetNCells())) {
			continue;
		}

		for (Int_t j = 0; j < c->GetNCells(); j++) {
			const Int_t cell_id = c->GetCellsAbsId()[j];
			const Float_t amplitude =
				emcal_cell->GetCellAmplitude(cell_id);

			_histogram_cell_id_amplitude->Fill(cell_id, amplitude);
		}
	}
}

AliEMCALRecoUtils *AliAnalysisTaskCalibEmcal::GetEMCALRecoUtils(void)
{
	return _reco_utils;
}
