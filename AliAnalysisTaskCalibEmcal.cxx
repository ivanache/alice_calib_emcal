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

	static const unsigned long lcm_range[47] = {
		1UL, /* 1 */ 1UL, /* 2 */ 2UL, /* 3 */ 6UL, /* 4 */ 12UL,
		/* 5 */ 60UL, /* 6 */ 60UL, /* 7 */ 420UL, /* 8 */ 840UL,
		/* 9 */ 2520UL, /* 10 */ 2520UL, /* 11 */ 27720UL,
		/* 12 */ 27720UL, /* 13 */ 360360UL, /* 14 */ 360360UL,
		/* 15 */ 360360UL, /* 16 */ 720720UL, /* 17 */ 12252240UL,
		/* 18 */ 12252240UL, /* 19 */ 232792560UL,
		/* 20 */ 232792560UL, /* 21 */ 232792560UL,
		/* 22 */ 232792560UL, /* 23 */ 5354228880UL,
		/* 24 */ 5354228880UL, /* 25 */ 26771144400UL,
		/* 26 */ 26771144400UL, /* 27 */ 80313433200UL,
		/* 28 */ 80313433200UL, /* 29 */ 2329089562800UL,
		/* 30 */ 2329089562800UL, /* 31 */ 72201776446800UL,
		/* 32 */ 144403552893600UL, /* 33 */ 144403552893600UL,
		/* 34 */ 144403552893600UL, /* 35 */ 144403552893600UL,
		/* 36 */ 144403552893600UL, /* 37 */ 5342931457063200UL,
		/* 38 */ 5342931457063200UL, /* 39 */ 5342931457063200UL,
		/* 40 */ 5342931457063200UL, /* 41 */ 219060189739591200UL,
		/* 42 */ 219060189739591200UL, /* 43 */ 9419588158802421600UL,
		/* 44 */ 9419588158802421600UL,
		/* 45 */ 9419588158802421600UL,
		/* 46 */ 9419588158802421600UL
	};

	double angular_range_reduce(const double x)
	{
		if(!std::isfinite(x)) {
			return x;
		}

		static const double cody_waite_x_max = 1608.4954386379741381;
		static const double two_pi_0 = 6.2831853071795649157;
		static const double two_pi_1 = 2.1561211432631314669e-14;
		static const double two_pi_2 = 1.1615423895917441336e-27;
		double ret;

		if(x >= -cody_waite_x_max && x <= cody_waite_x_max) {
			static const double inverse_two_pi =
				0.15915494309189534197;
			const double k = rint(x * inverse_two_pi);
			ret = ((x - (k * two_pi_0)) - k * two_pi_1) -
				k * two_pi_2;
		}
		else {
			long double sin_x;
			long double cos_x;

			sincosl(x, &sin_x, &cos_x);
			ret = (double)atan2l(sin_x, cos_x);
		}
		if(ret == -M_PI) {
			ret = M_PI;
		}

		return ret;
	}

	void to_sm_ieta_iphi(unsigned int &sm, unsigned int &ieta,
						 unsigned int &iphi, unsigned int n)
	{
		sm = n < 11520 ? n / 1152 :
			n < 12288 ? 10 + (n - 11520) / 384 :
			n < 16896 ? 12 + (n - 12288) / 768 :
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

	unsigned int flat_sm_ieta(const unsigned int sm,
							  const unsigned int ieta)
	{
		return sm < 12 ? sm * 48 + ieta :
			sm < 18 ? 576 + (sm - 12) * 32 + ieta :
			768 + (sm - 18) * 48 + ieta;
	}

	unsigned int flat_sm_iphi(const unsigned int sm,
							  const unsigned int iphi)
	{
		return sm < 10 ? sm * 24 + iphi :
			sm < 12 ? 240 + (sm - 10) * 8 + iphi :
			sm < 18 ? 256 + (sm - 12) * 24 + iphi :
			400 + (sm - 18) * 8 + iphi;
	}

	bool sm_ieta_warm(const unsigned int sm,
					  const unsigned int ieta)
	{
#if 0
		if (sm < 12 || sm >= 18) {
			if ((sm % 2 == 0 && ieta >= 36) ||
				(sm % 2 == 1 && ieta < 12)) {
				return true;
			}
			else if ((sm % 2 == 0 && ieta >= 38) ||
				(sm % 2 == 1 && ieta < 10)) {
				return false;
			}
		}

		if (sm < 7) {
			return (ieta + sm % 2 + 1) % 2 != 0;
		}
		else if (sm == 7) {
			return (ieta + sm % 2 + (ieta < 32 ? 1 : 0)) % 2 != 0;
		}
		else {
			return (ieta + sm % 2) % 2 != 0;
		}
#else
		return (ieta + sm % 2) % 2 != 0;
#endif
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

#define EMCAL_GEOMETRY_NAME	"EMCAL_COMPLETE12SMV1_DCAL_8SM"
//#define EMCAL_GEOMETRY_NAME	"EMCAL_COMPLETE12SMV1"

//#define TRIGGER_NAME	"CEMC7"
#define TRIGGER_NAME	"CINT7"
//#define TRIGGER_NAME	"CTRUE"
//#define TRIGGER_NAME	"CINT7EG1"

#define CLASS_INITIALIZATION							\
	_emcal_geometry_name(EMCAL_GEOMETRY_NAME),			\
		_trigger_name(TRIGGER_NAME),					\
		_mixed_event(22,								\
			std::deque<std::vector<std::pair<			\
				TLorentzVector, double> > >()),			\
		_list(NULL),									\
		_histogram_run_number_event_triggered(NULL),	\
		_histogram_multiplicity_0_5(NULL),				\
		_histogram_cell_id_timing(NULL),				\
		_histogram_cell_id_timing_9_gev(NULL),			\
		_histogram_timing_amplitude(NULL),				\
		_histogram_cell_id_amplitude(NULL),				\
		_histogram_cell_id_amplitude_scale(NULL),		\
		_histogram_cell_id_amplitude_low_mult(NULL),	\
		_histogram_cell_id_amplitude_high_mult(NULL),	\
		_histogram_cell_id_pi0_mass(NULL),				\
		_histogram_cell_id_pi0_mass_mixed(NULL),		\
		_histogram_cell_id_pi0_mass_recal(NULL),		\
		_histogram_sm_ieta_pi0_mass(NULL),				\
		_histogram_sm_ieta_pi0_mass_mixed(NULL),		\
		_histogram_sm_ieta_pi0_mass_recal(NULL),		\
		_histogram_sm_iphi_pi0_mass(NULL),				\
		_histogram_sm_iphi_pi0_mass_mixed(NULL),		\
		_histogram_sm_iphi_pi0_mass_recal(NULL),		\
		_histogram_sm_ieta_pi0_mass_same(NULL),			\
		_histogram_sm_ieta_pi0_mass_same_recal(NULL),	\
		_histogram_sm_ieta_pi0_mass_cross(NULL),		\
		_histogram_sm_ieta_pi0_mass_cross_recal(NULL),	\
		_histogram_pi0_deta_dphi(NULL),					\
		_histogram_pi0_deta_dphi_ieta_even(NULL),		\
		_histogram_pi0_deta_dphi_ieta_odd(NULL),		\
		_histogram_pi0_e_geom_mean_dr(NULL),			\
		_histogram_pi0_e_geom_mean_dr_ieta_even(NULL),	\
		_histogram_pi0_e_geom_mean_dr_ieta_odd(NULL),	\
		_histogram_pi0_correlated_cell_id_e(NULL),		\
		_reco_util(new AliEMCALRecoUtils),				\
		_emcal_geometry(NULL),							\
		_emcal_pi0_calibration_factor(),				\
		_ncell(17664), _nsm_ieta(864), _nsm_iphi(416),	\
		_timing_min(300), _timing_max(800),				\
		_pi0_e_min(4.0), _pi0_e_max(9.0),				\
		_pi0_ncell_min(2),								\
		_pi0_m02_min(0.1), _pi0_m02_max(0.5),			\
		_pi0_dtiming_max(100.0)

AliAnalysisTaskCalibEmcal::AliAnalysisTaskCalibEmcal(void)
	: AliAnalysisTaskSE(), CLASS_INITIALIZATION
{
}

AliAnalysisTaskCalibEmcal::AliAnalysisTaskCalibEmcal(
	const char *name)
	: AliAnalysisTaskSE(name), CLASS_INITIALIZATION
{
	DefineOutput(1, TList::Class());
}

AliAnalysisTaskCalibEmcal::AliAnalysisTaskCalibEmcal(
	const AliAnalysisTaskCalibEmcal &x)
	: AliAnalysisTaskSE(), CLASS_INITIALIZATION
{
}

AliAnalysisTaskCalibEmcal &AliAnalysisTaskCalibEmcal::operator=(
	const AliAnalysisTaskCalibEmcal &x)
{
	return *this;
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

	_histogram_multiplicity_0_5 =
		new TH1D("_histogram_multiplicity_0_5", "",
				 lcm_range[9], 0, lcm_range[9]);

	_list->Add(_histogram_multiplicity_0_5);

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
				 3 * lcm_range[6], 0.1, 100);

	_list->Add(_histogram_timing_amplitude);

	/////////////////////////////////////////////////////////////////

	_histogram_cell_id_amplitude =
		new TH2D("_histogram_cell_id_amplitude", "",
				 _ncell, -0.5, _ncell - 0.5,
				 3 * lcm_range[6], 0.1, 100);
	set_log_axis(_histogram_cell_id_amplitude->GetYaxis());

	_list->Add(_histogram_cell_id_amplitude);

	/////////////////////////////////////////////////////////////////

	_histogram_cell_id_amplitude_scale =
		new TH2D("_histogram_cell_id_amplitude_scale", "",
				 _ncell, -0.5, _ncell - 0.5,
				 3 * lcm_range[6], 0.1, 100);
	set_log_axis(_histogram_cell_id_amplitude_scale->GetYaxis());

	_list->Add(_histogram_cell_id_amplitude_scale);

	/////////////////////////////////////////////////////////////////

	_histogram_cell_id_amplitude_low_mult =
		new TH2D("_histogram_cell_id_amplitude_low_mult", "",
				 _ncell, -0.5, _ncell - 0.5,
				 3 * lcm_range[6], 0.1, 100);
	set_log_axis(_histogram_cell_id_amplitude_low_mult->GetYaxis());

	_list->Add(_histogram_cell_id_amplitude_low_mult);

	/////////////////////////////////////////////////////////////////

	_histogram_cell_id_amplitude_high_mult =
		new TH2D("_histogram_cell_id_amplitude_high_mult", "",
				 _ncell, -0.5, _ncell - 0.5,
				 3 * lcm_range[6], 0.1, 100);
	set_log_axis(_histogram_cell_id_amplitude_high_mult->GetYaxis());

	_list->Add(_histogram_cell_id_amplitude_high_mult);

	/////////////////////////////////////////////////////////////////

	_histogram_cell_id_pi0_mass =
		new TH2D("_histogram_cell_id_pi0_mass", "",
				 _ncell, -0.5, _ncell - 0.5,
				 lcm_range[7], 0, 1);

	_list->Add(_histogram_cell_id_pi0_mass);

	/////////////////////////////////////////////////////////////////

	_histogram_cell_id_pi0_mass_mixed =
		new TH2D("_histogram_cell_id_pi0_mass_mixed", "",
				 _ncell, -0.5, _ncell - 0.5,
				 lcm_range[7], 0, 1);

	_list->Add(_histogram_cell_id_pi0_mass_mixed);

	/////////////////////////////////////////////////////////////////

	_histogram_cell_id_pi0_mass_recal =
		new TH2D("_histogram_cell_id_pi0_mass_recal", "",
				 _ncell, -0.5, _ncell - 0.5,
				 lcm_range[7], 0, 1);

	_list->Add(_histogram_cell_id_pi0_mass_recal);

	/////////////////////////////////////////////////////////////////

	_histogram_sm_ieta_pi0_mass =
		new TH2D("_histogram_sm_ieta_pi0_mass", "",
				 _nsm_ieta, -0.5, _nsm_ieta - 0.5,
				 lcm_range[7], 0, 1);

	_list->Add(_histogram_sm_ieta_pi0_mass);

	/////////////////////////////////////////////////////////////////

	_histogram_sm_ieta_pi0_mass_mixed =
		new TH2D("_histogram_sm_ieta_pi0_mass_mixed", "",
				 _nsm_ieta, -0.5, _nsm_ieta - 0.5,
				 lcm_range[7], 0, 1);

	_list->Add(_histogram_sm_ieta_pi0_mass_mixed);

	/////////////////////////////////////////////////////////////////

	_histogram_sm_ieta_pi0_mass_recal =
		new TH2D("_histogram_sm_ieta_pi0_mass_recal", "",
				 _nsm_ieta, -0.5, _nsm_ieta - 0.5,
				 lcm_range[7], 0, 1);

	_list->Add(_histogram_sm_ieta_pi0_mass_recal);

	/////////////////////////////////////////////////////////////////

	_histogram_sm_iphi_pi0_mass =
		new TH2D("_histogram_sm_iphi_pi0_mass", "",
				 _nsm_iphi, -0.5, _nsm_iphi - 0.5,
				 lcm_range[7], 0, 1);

	_list->Add(_histogram_sm_iphi_pi0_mass);

	/////////////////////////////////////////////////////////////////

	_histogram_sm_iphi_pi0_mass_mixed =
		new TH2D("_histogram_sm_iphi_pi0_mass_mixed", "",
				 _nsm_iphi, -0.5, _nsm_iphi - 0.5,
				 lcm_range[7], 0, 1);

	_list->Add(_histogram_sm_iphi_pi0_mass_mixed);

	/////////////////////////////////////////////////////////////////

	_histogram_sm_iphi_pi0_mass_recal =
		new TH2D("_histogram_sm_iphi_pi0_mass_recal", "",
				 _nsm_iphi, -0.5, _nsm_iphi - 0.5,
				 lcm_range[7], 0, 1);

	_list->Add(_histogram_sm_iphi_pi0_mass_recal);

	/////////////////////////////////////////////////////////////////

	_histogram_sm_ieta_pi0_mass_same =
		new TH2D("_histogram_sm_ieta_pi0_mass_same", "",
				 _nsm_ieta, -0.5, _nsm_ieta - 0.5,
				 lcm_range[7], 0, 1);

	_list->Add(_histogram_sm_ieta_pi0_mass_same);

	/////////////////////////////////////////////////////////////////

	_histogram_sm_ieta_pi0_mass_same_recal =
		new TH2D("_histogram_sm_ieta_pi0_mass_same_recal", "",
				 _nsm_ieta, -0.5, _nsm_ieta - 0.5,
				 lcm_range[7], 0, 1);

	_list->Add(_histogram_sm_ieta_pi0_mass_same_recal);

	/////////////////////////////////////////////////////////////////

	_histogram_sm_ieta_pi0_mass_cross =
		new TH2D("_histogram_sm_ieta_pi0_mass_cross", "",
				 _nsm_ieta, -0.5, _nsm_ieta - 0.5,
				 lcm_range[7], 0, 1);

	_list->Add(_histogram_sm_ieta_pi0_mass_cross);

	/////////////////////////////////////////////////////////////////

	_histogram_sm_ieta_pi0_mass_cross_recal =
		new TH2D("_histogram_sm_ieta_pi0_mass_cross_recal", "",
				 _nsm_ieta, -0.5, _nsm_ieta - 0.5,
				 lcm_range[7], 0, 1);

	_list->Add(_histogram_sm_ieta_pi0_mass_cross_recal);

	/////////////////////////////////////////////////////////////////

	_histogram_pi0_deta_dphi =
		new TH2D("_histogram_pi0_deta_dphi", "",
				 lcm_range[7], -1, 1,
				 lcm_range[7], -1, 1);

	_list->Add(_histogram_pi0_deta_dphi);

	/////////////////////////////////////////////////////////////////

	_histogram_pi0_deta_dphi_ieta_even =
		new TH2D("_histogram_pi0_deta_dphi_ieta_even", "",
				 lcm_range[7], -1, 1,
				 lcm_range[7], -1, 1);

	_list->Add(_histogram_pi0_deta_dphi_ieta_even);

	/////////////////////////////////////////////////////////////////

	_histogram_pi0_deta_dphi_ieta_odd =
		new TH2D("_histogram_pi0_deta_dphi_ieta_odd", "",
				 lcm_range[7], -1, 1,
				 lcm_range[7], -1, 1);

	_list->Add(_histogram_pi0_deta_dphi_ieta_odd);

	/////////////////////////////////////////////////////////////////

	_histogram_pi0_e_geom_mean_dr =
		new TH2D("_histogram_pi0_e_geom_mean_dr", "",
				 3 * lcm_range[6], 0.1, 100,
				 lcm_range[7], 0, 0.2);
	set_log_axis(_histogram_pi0_e_geom_mean_dr->GetXaxis());

	_list->Add(_histogram_pi0_e_geom_mean_dr);

	/////////////////////////////////////////////////////////////////

	_histogram_pi0_e_geom_mean_dr_ieta_even =
		new TH2D("_histogram_pi0_e_geom_mean_dr_ieta_even", "",
				 3 * lcm_range[6], 0.1, 100,
				 lcm_range[7], 0, 0.2);
	set_log_axis(_histogram_pi0_e_geom_mean_dr_ieta_even->GetXaxis());

	_list->Add(_histogram_pi0_e_geom_mean_dr_ieta_even);

	/////////////////////////////////////////////////////////////////

	_histogram_pi0_e_geom_mean_dr_ieta_odd =
		new TH2D("_histogram_pi0_e_geom_mean_dr_ieta_odd", "",
				 3 * lcm_range[6], 0.1, 100,
				 lcm_range[7], 0, 0.2);
	set_log_axis(_histogram_pi0_e_geom_mean_dr_ieta_odd->GetXaxis());

	_list->Add(_histogram_pi0_e_geom_mean_dr_ieta_odd);

	/////////////////////////////////////////////////////////////////

	_histogram_pi0_correlated_cell_id_e =
		new TH2D("_histogram_pi0_correlated_cell_id_e", "",
				 _ncell, -0.5, _ncell - 0.5,
				 3 * lcm_range[6], 0.1, 100);
	set_log_axis(_histogram_pi0_correlated_cell_id_e->GetYaxis());

	_list->Add(_histogram_pi0_correlated_cell_id_e);

	/////////////////////////////////////////////////////////////////

	PostData(1, _list);
}

Double_t AliAnalysisTaskCalibEmcal::emcal_scale(const Int_t cell_id)
{
#if 0
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
#else
	double a = 1;
	double chi2_after = DBL_MAX;

	switch (cell_id) {

case 3: a = 1.22034643392; chi2_after = 53.5384375605;
case 13: a = 0.94574314978; chi2_after = 1.48334575416;
case 15: a = 0.948081827006; chi2_after = 4.48661214047;
case 58: a = 2.21154493162; chi2_after = 2396.21975189;
case 62: a = 1.54056473405; chi2_after = 1882.94857214;
case 74: a = 4.48037492791; chi2_after = 4419.72871493;
case 103: a = 4.19204000918; chi2_after = 42601.3404094;
case 145: a = 0.832939315382; chi2_after = 1.37632490866;
case 177: a = 3.6843494325; chi2_after = 0;
case 198: a = 2.70310791598; chi2_after = 4599.52418006;
case 288: a = 0.925703716728; chi2_after = 15.875037796;
case 289: a = 0.851187283685; chi2_after = 30.2093213993;
case 290: a = 0.93317575551; chi2_after = 14.9802329615;
case 291: a = 0.884436358089; chi2_after = 29.2023727048;
case 292: a = 0.92602167835; chi2_after = 25.2318334539;
case 294: a = 0.888436954003; chi2_after = 27.4773687506;
case 295: a = 0.955319634408; chi2_after = 13.2692383254;
case 296: a = 0.904397392868; chi2_after = 21.4568969738;
case 297: a = 0.960359499094; chi2_after = 11.2866176283;
case 298: a = 0.887104333228; chi2_after = 36.8547158876;
case 299: a = 0.947554858545; chi2_after = 28.7738155293;
case 300: a = 0.930828446977; chi2_after = 16.2744594429;
case 301: a = 0.919139274617; chi2_after = 21.7668925946;
case 302: a = 0.895706324315; chi2_after = 22.9692633103;
case 303: a = 0.854293919514; chi2_after = 50.2432690829;
case 338: a = 0.944807879732; chi2_after = 15.4724787389;
case 340: a = 0.901393368429; chi2_after = 38.4619693522;
case 341: a = 0.948840598312; chi2_after = 23.249744701;
case 342: a = 0.895597896534; chi2_after = 47.4538208481;
case 344: a = 0.86547796881; chi2_after = 38.2731684087;
case 345: a = 0.927524684157; chi2_after = 26.8130920829;
case 346: a = 0.915830729918; chi2_after = 27.3301030255;
case 347: a = 0.96553083784; chi2_after = 14.7914541667;
case 348: a = 0.925961659327; chi2_after = 22.3080632764;
case 349: a = 0.88123397888; chi2_after = 46.5171389396;
case 350: a = 0.953733317297; chi2_after = 19.9281149385;
case 351: a = 0.905205690264; chi2_after = 44.8489651439;
case 353: a = 2.88826461527; chi2_after = 5578.77694346;
case 554: a = 2.8853001502; chi2_after = 2180.27302122;
case 593: a = 2.60283717688; chi2_after = 1653.35013239;
case 594: a = 1.21066238255; chi2_after = 3.23442069261;
case 595: a = 2.77033994909; chi2_after = 1221.69770076;
case 597: a = 2.69608958969; chi2_after = 3831.80243187;
case 599: a = 2.63422276952; chi2_after = 2022.75396939;
case 601: a = 2.62094058997; chi2_after = 1394.41400099;
case 603: a = 2.6513793973; chi2_after = 844.78735893;
case 605: a = 2.75052821308; chi2_after = 1265.56762412;
case 607: a = 2.74726877043; chi2_after = 1587.73999843;
case 640: a = 2.75087483022; chi2_after = 2323.15453805;
case 641: a = 2.94452352499; chi2_after = 3756.63688588;
case 643: a = 2.62739870724; chi2_after = 2244.29025618;
case 644: a = 2.87300880294; chi2_after = 1070.97583482;
case 645: a = 3.10913009877; chi2_after = 2761.7577102;
case 646: a = 2.56834487813; chi2_after = 2219.00836129;
case 647: a = 2.76075264649; chi2_after = 4416.89345634;
case 648: a = 2.68955062046; chi2_after = 1469.01005339;
case 649: a = 2.81866423754; chi2_after = 1780.5437424;
case 650: a = 2.83288527961; chi2_after = 1637.16915777;
case 651: a = 2.89129685381; chi2_after = 2030.85372007;
case 652: a = 2.59961759544; chi2_after = 2651.35136936;
case 653: a = 2.96859130686; chi2_after = 2141.75511932;
case 654: a = 3.06879407667; chi2_after = 1056.66773679;
case 720: a = 2.49181379366; chi2_after = 530.925639223;
case 909: a = 1.02439586623; chi2_after = 0.984800026527;
case 917: a = 1.45113487632; chi2_after = 70.0875376917;
case 1002: a = 2.72139219718; chi2_after = 1516.55661866;
case 1038: a = 2.73177736449; chi2_after = 7148.1165097;
case 1061: a = 0.602255231675; chi2_after = 96.7935320232;
case 1067: a = 0.568147259231; chi2_after = 157.899047003;
case 1204: a = 2.8187148053; chi2_after = 11153.5204838;
case 1212: a = 2.93584560376; chi2_after = 1374.14460373;
case 1222: a = 2.69380900251; chi2_after = 4539.39938617;
case 1276: a = 0.66303002809; chi2_after = 41.51560853;
case 1288: a = 4.99452919287; chi2_after = 4619.42885028;
case 1366: a = 2.47101597171; chi2_after = 1253.45994499;
case 1376: a = 3.71894396838; chi2_after = 10004.2150329;
case 1380: a = 2.44831441884; chi2_after = 4136.81207441;
case 1384: a = 1.53794215024; chi2_after = 34.278836885;
case 1386: a = 3.14020596736; chi2_after = 11224.9753659;
case 1704: a = 3.8624501957; chi2_after = 4765.15546045;
case 1711: a = 2.00063997294; chi2_after = 206.509265993;
case 1738: a = 3.01062765979; chi2_after = 3915.45922266;
case 1836: a = 2.68895491311; chi2_after = 6833.67871383;
case 1837: a = 3.02172318734; chi2_after = 5378.98168611;
case 1838: a = 2.63992100803; chi2_after = 1007.05966649;
case 1844: a = 2.60207156575; chi2_after = 4308.08058304;
case 1860: a = 2.92635377201; chi2_after = 15281.1847526;
case 1892: a = 2.89324683937; chi2_after = 9948.04506847;
case 1963: a = 2.78881421708; chi2_after = 1504.01194317;
case 2014: a = 1.94490991755; chi2_after = 246.430586673;
case 2022: a = 3.66462826824; chi2_after = 30.0196822126;
case 2026: a = 3.9052582765; chi2_after = 11862.0218649;
case 2047: a = 3.346290396; chi2_after = 2536.37655826;
case 2112: a = 0.937234889419; chi2_after = 1.51989075958;
case 2114: a = 1.05637088116; chi2_after = 3.29097444201;
case 2115: a = 0.456558471624; chi2_after = 298.273077406;
case 2116: a = 0.679952598268; chi2_after = 255.919907814;
case 2118: a = 0.914271589411; chi2_after = 3.88637103901;
case 2120: a = 1.12845326201; chi2_after = 2.61074180699;
case 2123: a = 0.469291680416; chi2_after = 1998.43010555;
case 2124: a = 1.09374940176; chi2_after = 3.7058625356;
case 2125: a = 1.0807591508; chi2_after = 1.07546534696;
case 2161: a = 1.42679633749; chi2_after = 16.4491206543;
case 2196: a = 2.4581540003; chi2_after = 2258.59559031;
case 2245: a = 2.56974099508; chi2_after = 1642.35476335;
case 2298: a = 3.88214267074; chi2_after = 54697.7277754;
case 2325: a = 2.60351949865; chi2_after = 1302.40670852;
case 2341: a = 0.957864857627; chi2_after = 1.8904478088;
case 2371: a = 0.983532928691; chi2_after = 2.52834607481;
case 2395: a = 3.19977264211; chi2_after = 11988.6944396;
case 2424: a = 3.02272908734; chi2_after = 12250.1823415;
case 2487: a = 2.84726651144; chi2_after = 10877.3713682;
case 2505: a = 2.37510854883; chi2_after = 1293.83076552;
case 2533: a = 2.07457538473; chi2_after = 692.132688116;
case 2534: a = 2.84788658849; chi2_after = 3465.94843966;
case 2540: a = 2.37423296051; chi2_after = 1993.82651277;
case 2581: a = 2.78854134538; chi2_after = 2724.77632437;
case 2586: a = 2.83140514903; chi2_after = 2271.67335932;
case 2665: a = 2.59361739614; chi2_after = 635.756742087;
case 2778: a = 0.564052138122; chi2_after = 160.802990646;
case 2793: a = 2.22023035744; chi2_after = 84.9869997231;
case 2797: a = 2.98189458227; chi2_after = 2202.44793049;
case 2805: a = 2.60760829455; chi2_after = 4839.9444908;
case 2823: a = 2.20244256906; chi2_after = 152.99321478;
case 2824: a = 2.65463824653; chi2_after = 1709.89715178;
case 2825: a = 2.8178648501; chi2_after = 2306.97691835;
case 2857: a = 2.26445246599; chi2_after = 295.677696755;
case 2884: a = 3.20800321199; chi2_after = 3344.92792248;
case 2888: a = 3.41233640923; chi2_after = 11164.23981;
case 2891: a = 3.10764988949; chi2_after = 13516.8473137;
case 2915: a = 3.11721271699; chi2_after = 10879.9102903;
case 2921: a = 2.70357671132; chi2_after = 1310.96074413;
case 2985: a = 2.64556459167; chi2_after = 3829.42354824;
case 3135: a = 0.762310103263; chi2_after = 16.598868017;
case 3160: a = 1.13275611731; chi2_after = 2.09358271447;
case 3161: a = 2.1269754513; chi2_after = 117.959174324;
case 3196: a = 1.96554598211; chi2_after = 82.9974847453;
case 3223: a = 2.44145258202; chi2_after = 2160.29024665;
case 3236: a = 2.94794634196; chi2_after = 8343.53550628;
case 3244: a = 2.64276122133; chi2_after = 2115.12856663;
case 3259: a = 2.48086053799; chi2_after = 409.822813431;
case 3339: a = 2.84925346213; chi2_after = 11756.4352762;
case 3503: a = 2.32262706884; chi2_after = 1004.90732521;
case 3732: a = 2.827956057; chi2_after = 2544.98493225;
case 3740: a = 2.39566535245; chi2_after = 1094.79753596;
case 3748: a = 2.587578729; chi2_after = 1509.02733258;
case 3754: a = 1.58688591048; chi2_after = 31.2332081747;
case 3772: a = 1.71863050552; chi2_after = 66.6553019675;
case 3796: a = 1.89160568531; chi2_after = 159.295309913;
case 3826: a = 1.04173378869; chi2_after = 213.220039956;
case 3832: a = 0.990647053144; chi2_after = 62.2851230979;
case 3836: a = 1.11023198532; chi2_after = 135.139166181;
case 3838: a = 1.11362583746; chi2_after = 66.7295789663;
case 3839: a = 1.04917314865; chi2_after = 72.6410902395;
case 3854: a = 2.4311625633; chi2_after = 769.143362064;
case 3906: a = 2.1278509167; chi2_after = 552.409677404;
case 3908: a = 2.53627242465; chi2_after = 3844.96812777;
case 3914: a = 2.84870060162; chi2_after = 4448.71652882;
case 3940: a = 2.52079966082; chi2_after = 3243.10088495;
case 3955: a = 0.790231792015; chi2_after = 7.15898580699;
case 3974: a = 2.65668483837; chi2_after = 8217.61843206;
case 4011: a = 2.48074056024; chi2_after = 1570.11294502;
case 4026: a = 0.54547100499; chi2_after = 184.958859544;
case 4027: a = 2.17293201882; chi2_after = 521.173580741;
case 4058: a = 2.93756177912; chi2_after = 4302.81085224;
case 4129: a = 2.49421093243; chi2_after = 882.982836413;
case 4212: a = 2.48959401703; chi2_after = 990.057140406;
case 4230: a = 2.25586430386; chi2_after = 981.934118465;
case 4282: a = 2.60087552206; chi2_after = 995.644986894;
case 4421: a = 2.02204904553; chi2_after = 221.259211271;
case 4516: a = 3.18171441016; chi2_after = 9757.02187543;
case 4530: a = 3.05631940667; chi2_after = 10198.860135;
case 4531: a = 0.876352439319; chi2_after = 2.45353353974;
case 4532: a = 3.15930720246; chi2_after = 8447.59951071;
case 4538: a = 2.98568296035; chi2_after = 5666.66007591;
case 4543: a = 3.15274451609; chi2_after = 13455.1133157;
case 4596: a = 2.4615164993; chi2_after = 2025.09466139;
case 4605: a = 0.778794087532; chi2_after = 9.72683538363;
case 4613: a = 2.75989071706; chi2_after = 2071.3475296;
case 4617: a = 0.977517887361; chi2_after = 1.34994153396;
case 4621: a = 3.56363394855; chi2_after = 17115.5921802;
case 4627: a = 2.7372377233; chi2_after = 2359.9601076;
case 4637: a = 2.38839757927; chi2_after = 480.20533289;
case 4639: a = 0.974478437422; chi2_after = 1.72898173646;
case 4817: a = 3.22980000411; chi2_after = 1910.5079128;
case 4838: a = 1.64437081416; chi2_after = 688.152950736;
case 4846: a = 1.76076585995; chi2_after = 629.853358193;
case 5217: a = 2.8872130242; chi2_after = 1252.63041937;
case 5219: a = 3.92797234079; chi2_after = 252.381420852;
case 5221: a = 2.97391616729; chi2_after = 3700.23850163;
case 5223: a = 3.06177965576; chi2_after = 1304.83711505;
case 5225: a = 3.08913946045; chi2_after = 3056.29497599;
case 5227: a = 3.07375976016; chi2_after = 3006.14350598;
case 5229: a = 3.12914481704; chi2_after = 756.390599011;
case 5231: a = 2.41057797191; chi2_after = 1100.24831805;
case 5264: a = 3.06733056556; chi2_after = 1671.22088915;
case 5265: a = 3.39650355625; chi2_after = 1097.5196157;
case 5266: a = 3.39106886514; chi2_after = 1397.09149279;
case 5267: a = 3.13523846572; chi2_after = 3043.43761439;
case 5268: a = 2.94617633642; chi2_after = 2005.518748;
case 5269: a = 2.93479139086; chi2_after = 2163.82720201;
case 5270: a = 3.27848377052; chi2_after = 1760.39707824;
case 5271: a = 2.78710984128; chi2_after = 4620.66295382;
case 5272: a = 3.57107122427; chi2_after = 1600.60275165;
case 5273: a = 3.03075311841; chi2_after = 3403.32134893;
case 5274: a = 2.96345835011; chi2_after = 2699.44519563;
case 5275: a = 2.8846466958; chi2_after = 5516.19063817;
case 5276: a = 2.66826475769; chi2_after = 1247.75621586;
case 5277: a = 3.01050971354; chi2_after = 3403.70536407;
case 5278: a = 2.30028528639; chi2_after = 1078.20725144;
case 5279: a = 2.60118823778; chi2_after = 869.00257214;
case 5328: a = 0.524446870189; chi2_after = 0;
case 5330: a = 1.88974676628; chi2_after = 39.3608800857;
case 5332: a = 1.95034770259; chi2_after = 18.5552023879;
case 5334: a = 0.720493641748; chi2_after = 0;
case 5336: a = 0.585278383462; chi2_after = 0;
case 5338: a = 1.1667168024; chi2_after = 528.447372798;
case 5340: a = 0.555733309581; chi2_after = 1907.91335988;
case 5430: a = 0.963864962182; chi2_after = 1.22672912645;
case 5448: a = 1.62241586314; chi2_after = 27.9931163236;
case 5608: a = 0.766137330887; chi2_after = 10.184401132;
case 5648: a = 2.87279197265; chi2_after = 1017.23036465;
case 5650: a = 2.83747458917; chi2_after = 1424.98905251;
case 5652: a = 2.80209264994; chi2_after = 963.429832098;
case 5654: a = 2.80912201796; chi2_after = 1230.25847917;
case 5656: a = 2.85452416319; chi2_after = 1376.06391133;
case 5658: a = 2.88618905867; chi2_after = 713.156735411;
case 5660: a = 2.42842370938; chi2_after = 1381.86926256;
case 5663: a = 3.34297137072; chi2_after = 144.12945445;
case 5714: a = 1.82008704054; chi2_after = 43.8145646848;
case 5767: a = 3.08135783674; chi2_after = 4048.09067654;
case 5826: a = 0.401441228831; chi2_after = 1592.06753482;
case 5832: a = 0.996031131446; chi2_after = 15.3095123856;
case 5833: a = 1.03170348623; chi2_after = 10.9121514171;
case 5850: a = 4.89338570323; chi2_after = 1859.78662749;
case 6080: a = 4.55255640138; chi2_after = 617.839112311;
case 6083: a = 4.25800522241; chi2_after = 685.931628599;
case 6085: a = 4.2784900189; chi2_after = 1003.71620865;
case 6094: a = 3.59010601144; chi2_after = 761.448064408;
case 6104: a = 3.34297425606; chi2_after = 20958.8195603;
case 6128: a = 5.23973725775; chi2_after = 2060.61054712;
case 6129: a = 5.0670768619; chi2_after = 1674.77745842;
case 6136: a = 4.51095008032; chi2_after = 5377.8734382;
case 6137: a = 4.98327233486; chi2_after = 3924.68214113;
case 6139: a = 4.87588922416; chi2_after = 995.038952024;
case 6141: a = 4.56753556812; chi2_after = 871.983109924;
case 6275: a = 1.879532007; chi2_after = 89.558194611;
case 6331: a = 3.60447678551; chi2_after = 473.868252918;
case 6340: a = 0.672262900388; chi2_after = 74.7454917296;
case 6481: a = 2.99048794884; chi2_after = 253.246040117;
case 6624: a = 1.37763216313; chi2_after = 58.1833106657;
case 6625: a = 1.34690807139; chi2_after = 40.4029841856;
case 6626: a = 1.37521297043; chi2_after = 93.368572931;
case 6627: a = 1.33238174498; chi2_after = 72.4298136406;
case 6628: a = 1.02532442277; chi2_after = 1.38609296888;
case 6629: a = 1.29417772636; chi2_after = 59.2831406243;
case 6630: a = 1.31429509692; chi2_after = 95.4677351976;
case 6631: a = 1.28848497946; chi2_after = 54.2200621599;
case 6632: a = 1.3593825329; chi2_after = 109.41556721;
case 6633: a = 0.976954097401; chi2_after = 1.30625795213;
case 6634: a = 1.3493524985; chi2_after = 111.261054125;
case 6636: a = 1.33618457089; chi2_after = 108.840506576;
case 6637: a = 0.980620976452; chi2_after = 1.46489150531;
case 6638: a = 1.35394812833; chi2_after = 92.8306598301;
case 6639: a = 1.29538812298; chi2_after = 68.3001296198;
case 6645: a = 0.43553845244; chi2_after = 5381.3946948;
case 6646: a = 0.954820270482; chi2_after = 2.47026647462;
case 6649: a = 0.473322894541; chi2_after = 1310.76121402;
case 6650: a = 0.897237663015; chi2_after = 23.9759570885;
case 6672: a = 1.29595015539; chi2_after = 18.9215352894;
case 6673: a = 1.29659219135; chi2_after = 20.628609483;
case 6674: a = 1.04671991368; chi2_after = 2.4205982009;
case 6676: a = 1.37418321617; chi2_after = 90.9281555188;
case 6677: a = 1.33386571711; chi2_after = 77.8339279016;
case 6678: a = 1.3831768036; chi2_after = 70.6124546914;
case 6680: a = 1.40639438005; chi2_after = 85.0751792211;
case 6682: a = 1.04606756937; chi2_after = 2.55957164857;
case 6683: a = 1.35614859992; chi2_after = 87.3215184489;
case 6684: a = 1.02401635236; chi2_after = 1.1759105824;
case 6685: a = 1.30393074153; chi2_after = 62.5302168064;
case 6686: a = 1.03101397057; chi2_after = 0.834675195711;
case 6701: a = 0.573784103041; chi2_after = 205.432745876;
case 6735: a = 3.28644175673; chi2_after = 3863.08169415;
case 6811: a = 0.859619691795; chi2_after = 41.7596274236;
case 6886: a = 0.961763266558; chi2_after = 3.96305978137;
case 6945: a = 0.969182739224; chi2_after = 1.00068928924;
case 7371: a = 3.92287308505; chi2_after = 1467.2166644;
case 7417: a = 5.33165983737; chi2_after = 5890.4391414;
case 7425: a = 3.5264334742; chi2_after = 7149.1280175;
case 7572: a = 2.94306746263; chi2_after = 1393.60608773;
case 8047: a = 4.42439776858; chi2_after = 3229.37149442;
case 8244: a = 2.21702852742; chi2_after = 365.812367191;
case 8260: a = 1.87658442421; chi2_after = 125.971083878;
case 8275: a = 2.59762810478; chi2_after = 1797.9511821;
case 8340: a = 2.2777979141; chi2_after = 1173.54817404;
case 8352: a = 1.11945866565; chi2_after = 1.83074727398;
case 8356: a = 1.30809801588; chi2_after = 5.49574828799;
case 8357: a = 1.1119693173; chi2_after = 3.62925468095;
case 8420: a = 2.55964957016; chi2_after = 2449.98047095;
case 8610: a = 3.3645372356; chi2_after = 6121.73116885;
case 8807: a = 0.91854899774; chi2_after = 2.77807889834;
case 8809: a = 0.959366017168; chi2_after = 1.47297347038;
case 8813: a = 0.93149705868; chi2_after = 1.42162242209;
case 8855: a = 0.82764821515; chi2_after = 8.1147927438;
case 8904: a = 0.717087824012; chi2_after = 49.1760510501;
case 8907: a = 0.728995146799; chi2_after = 44.1923687576;
case 8911: a = 0.812024588851; chi2_after = 20.2502023954;
case 8916: a = 2.36742715928; chi2_after = 746.276971307;
case 8938: a = 3.17544460066; chi2_after = 6178.07238046;
case 8944: a = 0.754620071668; chi2_after = 38.529345006;
case 8947: a = 0.853191039636; chi2_after = 7.68411420532;
case 9078: a = 2.09699074714; chi2_after = 144.434372618;
case 9200: a = 0.984943467951; chi2_after = 2.92482470555;
case 9243: a = 0.972382821131; chi2_after = 1.28372602046;
case 9269: a = 3.82629537922; chi2_after = 60317.1898903;
case 9291: a = 2.37000810935; chi2_after = 838.183309019;
case 9338: a = 0.541646754905; chi2_after = 170.589035744;
case 9354: a = 1.66417566793; chi2_after = 65.4298162314;
case 9459: a = 1.03122965835; chi2_after = 1.12637535691;
case 9533: a = 2.25368605399; chi2_after = 381.587223005;
case 9598: a = 2.36295305932; chi2_after = 28.4552837344;
case 9703: a = 2.65826994407; chi2_after = 1350.6840056;
case 9706: a = 1.05950912286; chi2_after = 4.34462668645;
case 9769: a = 4.00243769412; chi2_after = 20986.7951018;
case 9798: a = 0.614278263814; chi2_after = 188.381910452;
case 9801: a = 0.455824373002; chi2_after = 700.880461424;
case 9807: a = 0.592456088918; chi2_after = 206.480935034;
case 9872: a = 0.691976968313; chi2_after = 162.122919333;
case 9878: a = 0.680749275946; chi2_after = 73.6640565728;
case 9927: a = 1.85819567519; chi2_after = 56.5269798489;
case 9940: a = 1.47613743115; chi2_after = 14.4605243901;
case 9942: a = 1.57931292264; chi2_after = 17.8489467833;
case 10091: a = 0.690206755052; chi2_after = 17.4446930105;
case 10134: a = 0.674478922009; chi2_after = 153.132501277;
case 10135: a = 0.599244851522; chi2_after = 151.97525757;
case 10138: a = 0.36049538392; chi2_after = 6482.84950289;
case 10143: a = 0.572446825078; chi2_after = 542.263082632;
case 10164: a = 1.79468363469; chi2_after = 203.600033496;
case 10203: a = 3.22449116476; chi2_after = 3355.02904591;
case 10326: a = 2.79336145567; chi2_after = 5356.82616087;
case 10451: a = 3.04319935705; chi2_after = 15213.3329498;
case 10561: a = 4.37702541482; chi2_after = 716.6080561;
case 10566: a = 4.15686440396; chi2_after = 1395.83583019;
case 10567: a = 4.07510433495; chi2_after = 6520.65056431;
case 10570: a = 4.77395438191; chi2_after = 2345.08777862;
case 10572: a = 4.27619008837; chi2_after = 2022.95062446;
case 10575: a = 4.25369474998; chi2_after = 4233.2296251;
case 10609: a = 4.24255654825; chi2_after = 478.700451663;
case 10612: a = 4.28928712546; chi2_after = 383.877545213;
case 10614: a = 4.27750546447; chi2_after = 568.482371567;
case 10615: a = 4.54166600564; chi2_after = 970.477429419;
case 10619: a = 4.72392690734; chi2_after = 2324.2071662;
case 10621: a = 4.74203192265; chi2_after = 816.392302416;
case 10678: a = 1.03405911359; chi2_after = 0.889120459676;
case 10829: a = 0.597461931937; chi2_after = 131.219389756;
case 10831: a = 0.506588399471; chi2_after = 355.285191436;
case 10921: a = 5.16346078221; chi2_after = 23179.6758731;
case 10981: a = 0.41897978542; chi2_after = 1342.56027538;
case 10982: a = 0.440926121775; chi2_after = 977.016435473;
case 10988: a = 0.465691733239; chi2_after = 348.685445744;
case 11024: a = 0.861954770551; chi2_after = 10.0532717765;
case 11034: a = 0.690548759822; chi2_after = 35.9306146831;
case 11042: a = 1.03650803355; chi2_after = 1.03216931054;
case 11044: a = 1.15568977195; chi2_after = 2.04073226627;
case 11048: a = 0.864879249324; chi2_after = 4.14554447495;
case 11050: a = 1.12947944104; chi2_after = 4.00896519713;
case 11052: a = 1.02241944319; chi2_after = 1.19080263056;
case 11054: a = 1.33754217078; chi2_after = 6.28663301967;
case 11091: a = 3.72235993154; chi2_after = 1571.74202323;
case 11093: a = 0.77845969317; chi2_after = 11.5428011955;
case 11095: a = 0.766657131853; chi2_after = 15.2360663054;
case 11097: a = 0.82455423961; chi2_after = 9.63587141526;
case 11099: a = 0.809273517115; chi2_after = 7.51248850852;
case 11102: a = 0.720929972796; chi2_after = 42.8734771103;
case 11141: a = 0.601345399041; chi2_after = 154.670962094;
case 11148: a = 0.806387446972; chi2_after = 13.0659703043;
case 11150: a = 0.775857773156; chi2_after = 25.0324776686;
case 11197: a = 0.828247989587; chi2_after = 11.6778785245;
case 11198: a = 1.28710272795; chi2_after = 4.83798284378;
case 11411: a = 0.691970490183; chi2_after = 71.2385624092;
case 11525: a = 0.917141809369; chi2_after = 6.33607721841;
case 11534: a = 0.944278538125; chi2_after = 3.19318473066;
case 11551: a = 0.969214376603; chi2_after = 2.62815679958;
case 11630: a = 0.671697069089; chi2_after = 394.403034045;
case 11647: a = 0.829435607021; chi2_after = 14.4807808701;
case 11918: a = 0.961042389673; chi2_after = 1.63350264351;
case 12032: a = 0.845458325; chi2_after = 7.68349122523;
case 12042: a = 0.944627141633; chi2_after = 4.65315970629;
case 12043: a = 0.956725193354; chi2_after = 4.41489516729;
case 12044: a = 0.943094245643; chi2_after = 4.3395099733;
case 12045: a = 0.950101022089; chi2_after = 4.96897857427;
case 12046: a = 0.853309226; chi2_after = 10.9041398401;
case 12047: a = 0.95357155415; chi2_after = 3.19379939721;
case 12117: a = 0.77854311265; chi2_after = 14.7179972637;
case 12126: a = 0.908736587657; chi2_after = 2.16731929839;
case 12127: a = 0.911986561348; chi2_after = 9.34897945496;
case 12160: a = 1.05269412539; chi2_after = 22.0670514504;
case 12174: a = 0.986438618393; chi2_after = 10.2290095397;
case 12271: a = 0.963386480137; chi2_after = 2.41966002347;
case 12276: a = 0.879950457333; chi2_after = 13.8632899977;
case 12278: a = 0.912372575791; chi2_after = 8.75675630376;
case 12280: a = 0.913639908116; chi2_after = 9.18224724128;
case 12282: a = 0.896682172483; chi2_after = 6.67293237812;
case 12320: a = 5.45550502726; chi2_after = 1823.5333633;
case 12322: a = 5.47087106233; chi2_after = 2595.79093043;
case 12325: a = 5.8177782481; chi2_after = 3699.97296897;
case 12327: a = 5.01783762099; chi2_after = 7071.87005672;
case 12331: a = 5.6712390139; chi2_after = 5738.82721278;
case 12333: a = 4.63743412726; chi2_after = 5630.88830181;
case 12371: a = 5.61571674514; chi2_after = 3480.41162487;
case 12373: a = 4.47400438754; chi2_after = 8060.06882806;
case 12374: a = 4.77280336453; chi2_after = 9284.93006118;
case 12375: a = 4.8144732306; chi2_after = 7111.47524423;
case 12377: a = 4.42706507631; chi2_after = 6165.71859662;
case 12379: a = 5.44139338079; chi2_after = 1946.28676429;
case 12381: a = 5.11608576059; chi2_after = 2256.24730997;
case 12383: a = 4.28442736213; chi2_after = 1031.37778589;
case 12604: a = 0.690474686621; chi2_after = 209.484499542;
case 12607: a = 0.32759632402; chi2_after = 7657.36850049;
case 12610: a = 0.875160583834; chi2_after = 11.516284423;
case 12613: a = 0.882878054263; chi2_after = 10.787384002;
case 12614: a = 0.437909064154; chi2_after = 2238.97722235;
case 12618: a = 0.591341801417; chi2_after = 827.308113222;
case 12621: a = 0.985077316921; chi2_after = 1.54651109493;
case 12622: a = 0.939818451591; chi2_after = 0.939366981246;
case 12623: a = 0.369349894407; chi2_after = 690.267525906;
case 12783: a = 0.737171186634; chi2_after = 242.681841792;
case 12801: a = 0.704942283773; chi2_after = 73.5112394629;
case 12802: a = 0.619876149315; chi2_after = 272.08931095;
case 12803: a = 0.40069415227; chi2_after = 8625.07235459;
case 12805: a = 0.358054585716; chi2_after = 11380.6568391;
case 12806: a = 0.370484194379; chi2_after = 8114.93934361;
case 12807: a = 0.958267726898; chi2_after = 1.33681768738;
case 12809: a = 0.506515403506; chi2_after = 1454.58948165;
case 12810: a = 0.380985490517; chi2_after = 7259.89184568;
case 12811: a = 0.84933445538; chi2_after = 8.4746432136;
case 12813: a = 0.929203361781; chi2_after = 1.69235827081;
case 12814: a = 0.545079551983; chi2_after = 77.9117270246;
case 12815: a = 0.26179089561; chi2_after = 11093.6633917;
case 12818: a = 0.745449669448; chi2_after = 27.6283172277;
case 12823: a = 0.951340931788; chi2_after = 6.57044742623;
case 12825: a = 0.939310611807; chi2_after = 17.2745294384;
case 12826: a = 0.728996334233; chi2_after = 143.51498549;
case 12829: a = 0.914022398867; chi2_after = 7.33019728832;
case 12848: a = 0.391682636685; chi2_after = 3138.49504093;
case 12850: a = 0.830040751787; chi2_after = 18.2876107671;
case 12853: a = 0.954033384327; chi2_after = 2.47555208834;
case 12856: a = 0.610180279404; chi2_after = 36.9985961262;
case 12858: a = 0.827843702182; chi2_after = 1.54198289088;
case 12862: a = 0.917585177311; chi2_after = 1.35515309336;
case 12870: a = 0.844169698941; chi2_after = 75.8096071099;
case 12871: a = 0.843726143854; chi2_after = 85.0961296301;
case 12874: a = 0.451907409007; chi2_after = 2303.38568566;
case 12879: a = 0.896098658982; chi2_after = 33.0801186428;
case 12914: a = 0.962098115979; chi2_after = 2.67477811666;
case 12921: a = 0.579202046219; chi2_after = 882.855011943;
case 12923: a = 1.05989820376; chi2_after = 3.62643330717;
case 12927: a = 0.531661873381; chi2_after = 1367.56503519;
case 12988: a = 0.713788806179; chi2_after = 23.4338439052;
case 13049: a = 1.89303566409; chi2_after = 200.374765998;
case 13053: a = 3.9791984324; chi2_after = 2704.03280544;
case 13056: a = 0.633578209048; chi2_after = 37.1104438874;
case 13058: a = 0.974932027511; chi2_after = 2.03000138607;
case 13064: a = 0.679411053787; chi2_after = 142.966305344;
case 13068: a = 2.94095040989; chi2_after = 12391.4526719;
case 13120: a = 0.657329778317; chi2_after = 256.883868981;
case 13122: a = 0.41392545494; chi2_after = 5919.80557368;
case 13124: a = 0.938744092185; chi2_after = 4.14183783109;
case 13126: a = 0.346045800074; chi2_after = 5107.80696795;
case 13128: a = 0.706736419616; chi2_after = 100.070060053;
case 13133: a = 0.931587844614; chi2_after = 6.76073757629;
case 13156: a = 0.541514883111; chi2_after = 221.357440648;
case 13248: a = 0.408943615888; chi2_after = 195.516618103;
case 13249: a = 0.400173633758; chi2_after = 192.258453125;
case 13280: a = 0.543670234265; chi2_after = 314.807777316;
case 13282: a = 0.546849781407; chi2_after = 271.917827567;
case 13284: a = 0.52950059451; chi2_after = 350.204899635;
case 13286: a = 0.67621304705; chi2_after = 140.938303845;
case 13287: a = 3.55428858346; chi2_after = 46.3008569576;
case 13288: a = 0.616195556938; chi2_after = 182.465642341;
case 13290: a = 0.551832839882; chi2_after = 271.634179759;
case 13292: a = 0.498108914174; chi2_after = 362.477201125;
case 13293: a = 0.652844254714; chi2_after = 118.571304428;
case 13294: a = 0.383791680487; chi2_after = 182.369923068;
case 13392: a = 0.879385153173; chi2_after = 8.74766708009;
case 13394: a = 0.370688062034; chi2_after = 9131.91349809;
case 13395: a = 0.96645972114; chi2_after = 1.29699350418;
case 13398: a = 0.980799718803; chi2_after = 2.0528035121;
case 13399: a = 0.934755858443; chi2_after = 1.50497636937;
case 13400: a = 0.732298393896; chi2_after = 347.239550814;
case 13402: a = 0.573406322123; chi2_after = 114.348917133;
case 13403: a = 0.587048337202; chi2_after = 613.618890177;
case 13406: a = 0.384376518565; chi2_after = 11950.4750445;
case 13433: a = 0.433319385932; chi2_after = 749.363519133;
case 13457: a = 0.94723494391; chi2_after = 2.99444580864;
case 13458: a = 0.983253641675; chi2_after = 1.27899173398;
case 13460: a = 0.390382913593; chi2_after = 7486.014664;
case 13462: a = 0.639089985186; chi2_after = 72.0512917018;
case 13463: a = 0.977342558356; chi2_after = 2.57873999704;
case 13464: a = 1.05318274264; chi2_after = 9.71983116742;
case 13465: a = 0.740993028236; chi2_after = 78.9877737435;
case 13467: a = 0.450109472588; chi2_after = 3738.7875112;
case 13469: a = 0.969986741033; chi2_after = 1.39299040477;
case 13471: a = 0.714464926481; chi2_after = 156.790990083;
case 13476: a = 0.607990605937; chi2_after = 84.445421131;
case 13477: a = 0.699242890704; chi2_after = 87.1399926555;
case 13479: a = 0.911387483634; chi2_after = 6.15122385632;
case 13481: a = 0.921316586895; chi2_after = 4.12777307614;
case 13485: a = 0.49618678829; chi2_after = 1258.37041477;
case 13487: a = 0.977836135107; chi2_after = 1.19864721881;
case 13508: a = 0.762400757543; chi2_after = 96.7094351923;
case 13509: a = 0.438229331962; chi2_after = 2166.04610879;
case 13510: a = 0.351040775602; chi2_after = 8091.71084641;
case 13512: a = 0.979448969916; chi2_after = 2.6313214593;
case 13513: a = 0.917798901749; chi2_after = 8.10450234323;
case 13517: a = 0.816997260325; chi2_after = 49.7601851045;
case 13521: a = 0.571001571319; chi2_after = 339.998197;
case 13524: a = 0.958916206117; chi2_after = 2.91470094951;
case 13525: a = 0.617406802806; chi2_after = 359.241817756;
case 13553: a = 0.989650759561; chi2_after = 1.48058984544;
case 13558: a = 0.70809995944; chi2_after = 369.644203403;
case 13559: a = 0.70757086014; chi2_after = 55.6025796845;
case 13574: a = 0.808846614025; chi2_after = 18.1843250411;
case 13578: a = 0.529362974612; chi2_after = 297.212343351;
case 13786: a = 0.930988262823; chi2_after = 4.92812559556;
case 13798: a = 0.9431672221; chi2_after = 2.20975132473;
case 13828: a = 0.424914003975; chi2_after = 1921.59075913;
case 13835: a = 0.347628821689; chi2_after = 7267.62460659;
case 13863: a = 0.914400394268; chi2_after = 5.1363190792;
case 13903: a = 0.572527944481; chi2_after = 585.963101939;
case 13921: a = 0.553953852335; chi2_after = 38.0148930234;
case 13926: a = 0.737626321264; chi2_after = 51.4493798485;
case 13933: a = 0.589375343906; chi2_after = 344.587246328;
case 13953: a = 3.54689483044; chi2_after = 32412.3907986;
case 13968: a = 0.753969998339; chi2_after = 14.8004510051;
case 13972: a = 0.854799716617; chi2_after = 29.6821062386;
case 13973: a = 0.886048112372; chi2_after = 26.3715969193;
case 13975: a = 0.95574119589; chi2_after = 16.2025209188;
case 13977: a = 0.551342635618; chi2_after = 540.363152107;
case 13980: a = 0.982089131219; chi2_after = 1.1371317424;
case 13983: a = 0.784103101299; chi2_after = 43.5522625006;
case 14058: a = 0.623085411688; chi2_after = 102.119786036;
case 14060: a = 0.578968061936; chi2_after = 133.786556846;
case 14116: a = 0.943400203798; chi2_after = 1.33509904571;
case 14117: a = 0.723780163618; chi2_after = 109.12372997;
case 14119: a = 0.492631566989; chi2_after = 1349.41163187;
case 14125: a = 0.812927071772; chi2_after = 53.4545947248;
case 14153: a = 1.51239544621; chi2_after = 19.8754144395;
case 14157: a = 0.393964384381; chi2_after = 1986.03252282;
case 14191: a = 5.44041165709; chi2_after = 2374.56062761;
case 14228: a = 1.73138175211; chi2_after = 48.002329272;
case 14229: a = 2.81716896487; chi2_after = 1798.40473552;
case 14234: a = 1.55289897527; chi2_after = 45.4505649234;
case 14236: a = 1.35540954766; chi2_after = 13.5405218453;
case 14239: a = 3.60356190575; chi2_after = 51044.7452796;
case 14249: a = 4.57587734453; chi2_after = 32600.4250975;
case 14303: a = 1.0249587483; chi2_after = 1.86358271905;
case 14305: a = 0.289979530427; chi2_after = 2200.66580263;
case 14306: a = 0.386176778405; chi2_after = 1402.19075444;
case 14308: a = 0.418957442858; chi2_after = 1601.3153694;
case 14309: a = 0.368297909951; chi2_after = 2004.17701599;
case 14311: a = 0.454530964456; chi2_after = 725.804846552;
case 14313: a = 0.73464803844; chi2_after = 70.864797555;
case 14314: a = 0.621755215828; chi2_after = 537.188803845;
case 14317: a = 0.444961451774; chi2_after = 882.803901089;
case 14319: a = 0.527327412609; chi2_after = 406.492780874;
case 14325: a = 0.359276979649; chi2_after = 3521.95606107;
case 14362: a = 0.617142718846; chi2_after = 109.128407182;
case 14363: a = 0.455367051046; chi2_after = 779.356550368;
case 14373: a = 0.530555365224; chi2_after = 266.709675201;
case 14400: a = 0.957829286376; chi2_after = 0.858974782525;
case 14402: a = 0.866036011784; chi2_after = 10.9935719694;
case 14403: a = 0.466663041122; chi2_after = 828.2930631;
case 14409: a = 0.850950434756; chi2_after = 37.6221933271;
case 14438: a = 1.02460524035; chi2_after = 1.63655290039;
case 14498: a = 0.867408548959; chi2_after = 5.1994943166;
case 14501: a = 0.649423749932; chi2_after = 127.442564431;
case 14502: a = 0.609666914082; chi2_after = 126.130580168;
case 14503: a = 0.551109836215; chi2_after = 256.663778583;
case 14505: a = 0.440816042576; chi2_after = 1234.13784645;
case 14544: a = 0.858738477193; chi2_after = 14.2632115353;
case 14545: a = 0.879726522483; chi2_after = 2.69598642715;
case 14549: a = 0.449316053186; chi2_after = 2756.78291894;
case 14550: a = 0.952825964705; chi2_after = 20.1408826304;
case 14553: a = 1.03498428251; chi2_after = 2.91659055776;
case 14554: a = 0.386559112069; chi2_after = 2729.36218015;
case 14557: a = 3.40816372924; chi2_after = 49.4834399741;
case 14558: a = 1.02447677269; chi2_after = 2.34471086921;
case 14592: a = 0.487529132032; chi2_after = 45.3737710917;
case 14594: a = 0.71428195995; chi2_after = 179.003758169;
case 14596: a = 0.801962653802; chi2_after = 50.7940313579;
case 14597: a = 0.865694843319; chi2_after = 40.7769490696;
case 14598: a = 0.602254343648; chi2_after = 617.351481241;
case 14605: a = 0.882746118232; chi2_after = 52.897482905;
case 14610: a = 0.375604424406; chi2_after = 2322.47505796;
case 14616: a = 0.939181136644; chi2_after = 2.47873378305;
case 14617: a = 0.880320619066; chi2_after = 68.608334329;
case 14618: a = 0.464240757873; chi2_after = 1615.98991312;
case 14620: a = 0.424721258518; chi2_after = 3813.80612349;
case 14622: a = 0.81623146007; chi2_after = 19.2804657214;
case 14623: a = 0.90635463626; chi2_after = 41.549788976;
case 14625: a = 1.74247855527; chi2_after = 49.798432475;
case 14638: a = 3.8429348216; chi2_after = 0;
case 14640: a = 0.282218518001; chi2_after = 981.994535532;
case 14647: a = 1.23901801894; chi2_after = 40.4090238121;
case 14660: a = 0.529125942509; chi2_after = 353.123753604;
case 14661: a = 0.852361108031; chi2_after = 5.27623917769;
case 14662: a = 0.35885513314; chi2_after = 3326.88897769;
case 14663: a = 0.547333780444; chi2_after = 471.777021869;
case 14668: a = 0.320584331663; chi2_after = 8334.65689875;
case 14669: a = 0.45923118475; chi2_after = 1002.60273947;
case 14672: a = 0.445656328716; chi2_after = 814.565636619;
case 14676: a = 0.686257499937; chi2_after = 156.224162456;
case 14678: a = 0.97268759589; chi2_after = 3.06285582339;
case 14679: a = 2.74950745238; chi2_after = 6398.91510535;
case 14688: a = 0.381293495915; chi2_after = 1437.79083764;
case 14689: a = 0.430169975503; chi2_after = 736.076011822;
case 14690: a = 0.653041328952; chi2_after = 220.526002267;
case 14691: a = 0.380108285367; chi2_after = 8138.33171852;
case 14694: a = 0.502390654545; chi2_after = 703.859739127;
case 14696: a = 0.747484270103; chi2_after = 34.0224852325;
case 14697: a = 0.51372263596; chi2_after = 1316.15312268;
case 14702: a = 0.493737771069; chi2_after = 1413.03608744;
case 14703: a = 0.849839113115; chi2_after = 7.73209199867;
case 14704: a = 0.793049132207; chi2_after = 114.440457515;
case 14709: a = 0.485995001006; chi2_after = 2222.01018405;
case 14711: a = 0.306333175076; chi2_after = 12948.0083914;
case 14712: a = 0.873482530587; chi2_after = 34.1246726646;
case 14713: a = 0.859839563747; chi2_after = 53.1778451642;
case 14721: a = 0.912359068133; chi2_after = 15.0719817771;
case 14722: a = 1.03886797529; chi2_after = 3.03433541578;
case 14725: a = 1.02644779232; chi2_after = 3.24295384984;
case 14726: a = 0.698712038509; chi2_after = 217.223913416;
case 14727: a = 0.763407046427; chi2_after = 67.1672547512;
case 14731: a = 0.689521578492; chi2_after = 232.402434135;
case 14733: a = 0.657525083488; chi2_after = 241.649867958;
case 14734: a = 0.292665452865; chi2_after = 6095.35277364;
case 14735: a = 0.478851487065; chi2_after = 512.974311124;
case 14736: a = 0.279775145612; chi2_after = 17565.1870701;
case 14737: a = 0.738885138302; chi2_after = 4.61495017521;
case 14738: a = 0.345269233961; chi2_after = 4790.90506499;
case 14739: a = 0.881210161101; chi2_after = 4.79560540755;
case 14741: a = 0.341372293911; chi2_after = 5149.18260949;
case 14742: a = 0.906993221194; chi2_after = 7.25196503649;
case 14744: a = 0.641339232993; chi2_after = 137.022297626;
case 14746: a = 0.985131812467; chi2_after = 15.3793557013;
case 14750: a = 0.660135661254; chi2_after = 552.092831186;
case 14753: a = 0.586198430893; chi2_after = 605.724008052;
case 14754: a = 0.558130700171; chi2_after = 323.013444669;
case 14756: a = 0.899937163495; chi2_after = 8.53299703834;
case 14757: a = 0.41739302317; chi2_after = 4969.81163832;
case 14758: a = 0.526522566864; chi2_after = 2085.76304855;
case 14760: a = 0.323739597015; chi2_after = 25942.1181841;
case 14761: a = 0.601861430008; chi2_after = 174.953949167;
case 14762: a = 0.780644539488; chi2_after = 40.3368727117;
case 14763: a = 0.510186320598; chi2_after = 967.281348231;
case 14764: a = 0.510631896557; chi2_after = 1162.8832448;
case 14765: a = 0.78974958517; chi2_after = 160.221647598;
case 14766: a = 0.415356289847; chi2_after = 4749.79744349;
case 14767: a = 0.486492312555; chi2_after = 515.435941336;
case 14772: a = 1.1399057068; chi2_after = 100.56046569;
case 14773: a = 0.574426700612; chi2_after = 486.542543893;
case 14774: a = 0.778633232829; chi2_after = 45.0377352002;
case 14775: a = 0.427690489992; chi2_after = 2758.42751262;
case 14779: a = 0.397054830642; chi2_after = 7006.95866434;
case 14839: a = 0.457809197865; chi2_after = 793.411447663;
case 14862: a = 5.21771262476; chi2_after = 1163.39262105;
case 14867: a = 2.33246665492; chi2_after = 916.867214538;
case 14874: a = 0.89903312584; chi2_after = 9.76021644488;
case 14977: a = 0.401233682629; chi2_after = 916.119894378;
case 14986: a = 0.640315062494; chi2_after = 300.426916175;
case 14991: a = 0.533081919349; chi2_after = 539.436868266;
case 15009: a = 0.452405356773; chi2_after = 1458.49459087;
case 15012: a = 2.4517185438; chi2_after = 450.3372368;
case 15015: a = 0.478001062514; chi2_after = 601.912203325;
case 15020: a = 0.382781512541; chi2_after = 2088.59415834;
case 15021: a = 0.887868059464; chi2_after = 10.6282344623;
case 15025: a = 0.441216039685; chi2_after = 310.942036949;
case 15027: a = 0.701728474558; chi2_after = 50.126146114;
case 15033: a = 0.519348595677; chi2_after = 390.999880035;
case 15094: a = 0.455395063836; chi2_after = 1054.8813934;
case 15137: a = 0.447363216196; chi2_after = 785.508754998;
case 15140: a = 0.330753048514; chi2_after = 1081.43138361;
case 15143: a = 0.372188789332; chi2_after = 1731.61605129;
case 15149: a = 0.427876536984; chi2_after = 1516.21568335;
case 15152: a = 0.75818870375; chi2_after = 24.3458992402;
case 15154: a = 0.726077584897; chi2_after = 14.4347674427;
case 15157: a = 0.558755870269; chi2_after = 340.700518563;
case 15161: a = 0.511295486292; chi2_after = 328.127106623;
case 15164: a = 0.335129822478; chi2_after = 1787.16258764;
case 15266: a = 0.483219401341; chi2_after = 1837.58988947;
case 15272: a = 0.599088201422; chi2_after = 562.23290064;
case 15280: a = 0.899044055751; chi2_after = 16.6461159963;
case 15282: a = 0.355513446164; chi2_after = 4892.39822012;
case 15284: a = 0.656819982254; chi2_after = 261.359533881;
case 15287: a = 0.925104138108; chi2_after = 2.76496702061;
case 15288: a = 0.402695836282; chi2_after = 6274.22221107;
case 15289: a = 0.984790012886; chi2_after = 1.85078021667;
case 15290: a = 0.823880906075; chi2_after = 5.91536932708;
case 15291: a = 0.407770629155; chi2_after = 5363.45838445;
case 15292: a = 0.966836416894; chi2_after = 1.22141322623;
case 15293: a = 0.962284231651; chi2_after = 1.04900025829;
case 15294: a = 0.449558958737; chi2_after = 2393.5017712;
case 15295: a = 0.49024841214; chi2_after = 2224.11227587;
case 15298: a = 0.89590587052; chi2_after = 3.12528097453;
case 15299: a = 0.472556645674; chi2_after = 2985.44699258;
case 15318: a = 0.370911745745; chi2_after = 6586.86122446;
case 15326: a = 0.365445239955; chi2_after = 8758.57232061;
case 15328: a = 0.708355208634; chi2_after = 241.942523042;
case 15331: a = 0.426416627292; chi2_after = 4616.33633461;
case 15333: a = 0.450042313953; chi2_after = 1537.16447015;
case 15334: a = 5.61534030449; chi2_after = 0;
case 15338: a = 0.425611604593; chi2_after = 5350.64323311;
case 15340: a = 0.949648111793; chi2_after = 3.4288180697;
case 15341: a = 0.88817432577; chi2_after = 4.4983398422;
case 15342: a = 0.942570358389; chi2_after = 1.15229432226;
case 15345: a = 0.454975034825; chi2_after = 1217.73048944;
case 15348: a = 0.450005506967; chi2_after = 2187.02729473;
case 15350: a = 0.584652683947; chi2_after = 709.061552801;
case 15354: a = 0.345357612889; chi2_after = 6759.01991488;
case 15359: a = 0.366470485625; chi2_after = 663.036418321;
case 15367: a = 0.964246835698; chi2_after = 1.28254780643;
case 15398: a = 0.701107680254; chi2_after = 160.877818879;
case 15460: a = 0.734305190086; chi2_after = 19.3946954537;
case 15462: a = 0.963519994631; chi2_after = 1.90043020117;
case 15469: a = 0.712885904655; chi2_after = 105.86920519;
case 15476: a = 0.656707649694; chi2_after = 54.2662067465;
case 15493: a = 0.631238094415; chi2_after = 142.018809443;
case 15526: a = 0.789830867039; chi2_after = 15.2466181547;
case 15670: a = 0.599350703269; chi2_after = 199.217362044;
case 15718: a = 0.640805387052; chi2_after = 156.65625195;
case 15719: a = 2.04922101633; chi2_after = 168.555157984;
case 15754: a = 5.12282889743; chi2_after = 39910.0893788;
case 15777: a = 0.437122402749; chi2_after = 1891.22922065;
case 15783: a = 0.944590967585; chi2_after = 0.899475986559;
case 15785: a = 2.39705199263; chi2_after = 522.276655011;
case 15821: a = 1.87399817198; chi2_after = 927.119617703;
case 15861: a = 0.548165866377; chi2_after = 442.377039274;
case 15870: a = 0.441229477837; chi2_after = 899.700080313;
case 15877: a = 0.567949980255; chi2_after = 316.148977003;
case 15908: a = 0.343546706585; chi2_after = 3228.3123988;
case 15945: a = 0.960598113434; chi2_after = 1.30454775787;
case 15947: a = 0.423501182477; chi2_after = 1096.40204804;
case 15949: a = 0.454979106547; chi2_after = 744.273768196;
case 15955: a = 0.44351078655; chi2_after = 890.737427392;
case 15965: a = 3.83354318087; chi2_after = 7821.16107617;
case 15967: a = 3.70094710977; chi2_after = 2809.64747568;
case 15972: a = 0.395858916214; chi2_after = 2353.55716225;
case 15979: a = 0.489985097153; chi2_after = 421.191199448;
case 16026: a = 0.61615760122; chi2_after = 114.265788666;
case 16054: a = 0.406633776343; chi2_after = 1699.94774914;
case 16088: a = 0.802729644175; chi2_after = 36.721814494;
case 16089: a = 0.81835873618; chi2_after = 76.8491612472;
case 16090: a = 0.944022648553; chi2_after = 7.00188423672;
case 16091: a = 0.392307871718; chi2_after = 5067.63483222;
case 16094: a = 0.787112242841; chi2_after = 90.3090405693;
case 16215: a = 2.04400024471; chi2_after = 136.741959108;
case 16258: a = 2.03218686653; chi2_after = 353.480812893;
case 16268: a = 2.0729419305; chi2_after = 262.439424783;
case 16304: a = 0.583250815117; chi2_after = 1273.08692931;
case 16306: a = 0.920002225148; chi2_after = 4.02415781261;
case 16307: a = 0.905704026549; chi2_after = 21.8882967791;
case 16308: a = 0.875566291449; chi2_after = 14.4896835632;
case 16317: a = 0.788044897178; chi2_after = 71.637199397;
case 16388: a = 0.352186092611; chi2_after = 9919.96663363;
case 16428: a = 0.760357246133; chi2_after = 200.269996339;
case 16429: a = 0.369189957277; chi2_after = 6827.95966965;
case 16440: a = 0.577659019649; chi2_after = 277.493965407;
case 16443: a = 0.365266003645; chi2_after = 11207.0716588;
case 16446: a = 0.536882558283; chi2_after = 1270.73775964;
case 16465: a = 0.673224247143; chi2_after = 92.1649948428;
case 16468: a = 0.781064113108; chi2_after = 66.1913415387;
case 16469: a = 0.932995395638; chi2_after = 3.19312114949;
case 16472: a = 0.884306979038; chi2_after = 4.11193629821;
case 16473: a = 0.539467490121; chi2_after = 618.781937733;
case 16475: a = 0.796452164901; chi2_after = 102.547424803;
case 16476: a = 0.477104155007; chi2_after = 2438.38428493;
case 16477: a = 0.795954336726; chi2_after = 101.028746285;
case 16479: a = 0.590587616027; chi2_after = 629.829507227;
case 16482: a = 0.773003286149; chi2_after = 79.5323761349;
case 16483: a = 0.97177798996; chi2_after = 2.28002493027;
case 16488: a = 0.446477573505; chi2_after = 4177.30549053;
case 16489: a = 1.03539553254; chi2_after = 1.15970889519;
case 16491: a = 0.841734215014; chi2_after = 50.4535171943;
case 16492: a = 0.58119392284; chi2_after = 704.60500538;
case 16496: a = 0.380200685374; chi2_after = 7441.00799259;
case 16497: a = 0.857604044268; chi2_after = 103.292105248;
case 16498: a = 0.725875544929; chi2_after = 190.538262893;
case 16502: a = 0.846918279473; chi2_after = 23.6287194578;
case 16506: a = 0.668318370516; chi2_after = 296.640504179;
case 16508: a = 0.951225517006; chi2_after = 2.65405115629;
case 16509: a = 0.372757536119; chi2_after = 3504.80704921;
case 16634: a = 3.57870237024; chi2_after = 46231.2309186;
case 16784: a = 3.54079059228; chi2_after = 51367.1325463;
case 16785: a = 2.29505488501; chi2_after = 805.23814222;
case 16789: a = 2.10312571324; chi2_after = 390.642113739;
case 16818: a = 1.05928499487; chi2_after = 7.46538749539;
case 16819: a = 0.652694394079; chi2_after = 277.802397854;
case 16821: a = 0.49605876291; chi2_after = 851.689881889;
case 16825: a = 0.989059503561; chi2_after = 1.86377152278;
case 16826: a = 0.364780255824; chi2_after = 11622.2979826;
case 16827: a = 0.927111202314; chi2_after = 14.0287811975;
case 16865: a = 0.6988360987; chi2_after = 184.404804203;
case 16867: a = 0.640996302171; chi2_after = 333.459945611;
case 16886: a = 3.94435369815; chi2_after = 55359.9355792;
case 16907: a = 0.866805580284; chi2_after = 10.8431149405;
case 16910: a = 0.868158998531; chi2_after = 3.52537086683;
case 16911: a = 0.986830801742; chi2_after = 3.44798555593;
case 16927: a = 0.909704288465; chi2_after = 2.41967215729;
case 17091: a = 0.396896660456; chi2_after = 4136.05601625;
case 17097: a = 0.642649452599; chi2_after = 369.432034028;
case 17107: a = 0.620847208257; chi2_after = 57.7980336448;
case 17117: a = 1.884119558; chi2_after = 123.706689008;
case 17202: a = 3.01463554885; chi2_after = 4780.54960875;
case 17204: a = 4.61943382884; chi2_after = 6710.190102;
case 17206: a = 1.41539877065; chi2_after = 15.2864674827;
case 17237: a = 2.04853974655; chi2_after = 40.5744695891;
case 17295: a = 2.42212815907; chi2_after = 768.16792482;
case 17396: a = 0.719776018533; chi2_after = 14.234113405;
case 17452: a = 0.766535713502; chi2_after = 23.2655011592;
case 17463: a = 0.782806020627; chi2_after = 21.7596832248;
case 17537: a = 0.455579491892; chi2_after = 307.542646254;

	}
	return chi2_after < 100 && a < 1.2 ? a : 1.0;
#endif
}

namespace {

	bool bad_cell(unsigned int cell_id)
	{
		switch (cell_id) {

case 3: return true;
case 13: return true;
case 15: return true;
case 57: return true;
case 59: return true;
case 61: return true;
case 63: return true;
case 74: return true;
case 103: return true;
case 143: return true;
case 145: return true;
case 176: return true;
case 191: return true;
case 198: return true;
case 287: return true;
case 288: return true;
case 289: return true;
case 290: return true;
case 291: return true;
case 292: return true;
case 293: return true;
case 294: return true;
case 295: return true;
case 296: return true;
case 297: return true;
case 298: return true;
case 299: return true;
case 300: return true;
case 301: return true;
case 302: return true;
case 303: return true;
case 336: return true;
case 337: return true;
case 338: return true;
case 339: return true;
case 340: return true;
case 341: return true;
case 342: return true;
case 343: return true;
case 344: return true;
case 345: return true;
case 346: return true;
case 347: return true;
case 348: return true;
case 349: return true;
case 350: return true;
case 351: return true;
case 353: return true;
case 399: return true;
case 447: return true;
case 554: return true;
case 593: return true;
case 594: return true;
case 595: return true;
case 597: return true;
case 599: return true;
case 601: return true;
case 603: return true;
case 605: return true;
case 607: return true;
case 640: return true;
case 641: return true;
case 643: return true;
case 644: return true;
case 645: return true;
case 646: return true;
case 647: return true;
case 648: return true;
case 649: return true;
case 650: return true;
case 651: return true;
case 652: return true;
case 653: return true;
case 655: return true;
case 720: return true;
case 909: return true;
case 917: return true;
case 1002: return true;
case 1038: return true;
case 1061: return true;
case 1067: return true;
case 1204: return true;
case 1212: return true;
case 1222: return true;
case 1276: return true;
case 1288: return true;
case 1366: return true;
case 1376: return true;
case 1380: return true;
case 1384: return true;
case 1386: return true;
case 1414: return true;
case 1519: return true;
case 1534: return true;
case 1535: return true;
case 1704: return true;
case 1711: return true;
case 1738: return true;
case 1836: return true;
case 1837: return true;
case 1838: return true;
case 1844: return true;
case 1860: return true;
case 1892: return true;
case 1963: return true;
case 1967: return true;
case 2014: return true;
case 2022: return true;
case 2026: return true;
case 2047: return true;
case 2112: return true;
case 2113: return true;
case 2114: return true;
case 2115: return true;
case 2117: return true;
case 2118: return true;
case 2120: return true;
case 2123: return true;
case 2124: return true;
case 2125: return true;
case 2161: return true;
case 2196: return true;
case 2245: return true;
case 2298: return true;
case 2325: return true;
case 2341: return true;
case 2371: return true;
case 2395: return true;
case 2397: return true;
case 2399: return true;
case 2424: return true;
case 2487: return true;
case 2505: return true;
case 2506: return true;
case 2533: return true;
case 2534: return true;
case 2540: return true;
case 2575: return true;
case 2581: return true;
case 2586: return true;
case 2665: return true;
case 2688: return true;
case 2776: return true;
case 2778: return true;
case 2793: return true;
case 2797: return true;
case 2805: return true;
case 2823: return true;
case 2824: return true;
case 2825: return true;
case 2857: return true;
case 2884: return true;
case 2888: return true;
case 2891: return true;
case 2915: return true;
case 2921: return true;
case 2985: return true;
case 3039: return true;
case 3135: return true;
case 3160: return true;
case 3161: return true;
case 3196: return true;
case 3223: return true;
case 3236: return true;
case 3244: return true;
case 3259: return true;
case 3297: return true;
case 3339: return true;
case 3488: return true;
case 3503: return true;
case 3732: return true;
case 3740: return true;
case 3748: return true;
case 3754: return true;
case 3772: return true;
case 3776: return true;
case 3777: return true;
case 3778: return true;
case 3779: return true;
case 3780: return true;
case 3781: return true;
case 3782: return true;
case 3783: return true;
case 3784: return true;
case 3785: return true;
case 3786: return true;
case 3787: return true;
case 3788: return true;
case 3789: return true;
case 3790: return true;
case 3791: return true;
case 3796: return true;
case 3824: return true;
case 3825: return true;
case 3826: return true;
case 3827: return true;
case 3828: return true;
case 3829: return true;
case 3830: return true;
case 3831: return true;
case 3832: return true;
case 3833: return true;
case 3834: return true;
case 3835: return true;
case 3836: return true;
case 3837: return true;
case 3838: return true;
case 3839: return true;
case 3854: return true;
case 3906: return true;
case 3908: return true;
case 3914: return true;
case 3940: return true;
case 3955: return true;
case 3974: return true;
case 4011: return true;
case 4026: return true;
case 4027: return true;
case 4058: return true;
case 4129: return true;
case 4212: return true;
case 4230: return true;
case 4282: return true;
case 4320: return true;
case 4421: return true;
case 4516: return true;
case 4530: return true;
case 4531: return true;
case 4532: return true;
case 4538: return true;
case 4543: return true;
case 4596: return true;
case 4605: return true;
case 4613: return true;
case 4617: return true;
case 4621: return true;
case 4627: return true;
case 4637: return true;
case 4639: return true;
case 4817: return true;
case 4835: return true;
case 4837: return true;
case 4839: return true;
case 4847: return true;
case 5039: return true;
case 5135: return true;
case 5183: return true;
case 5217: return true;
case 5219: return true;
case 5221: return true;
case 5223: return true;
case 5225: return true;
case 5227: return true;
case 5229: return true;
case 5231: return true;
case 5264: return true;
case 5265: return true;
case 5266: return true;
case 5267: return true;
case 5268: return true;
case 5269: return true;
case 5270: return true;
case 5271: return true;
case 5272: return true;
case 5273: return true;
case 5274: return true;
case 5275: return true;
case 5276: return true;
case 5277: return true;
case 5278: return true;
case 5279: return true;
case 5329: return true;
case 5331: return true;
case 5333: return true;
case 5335: return true;
case 5337: return true;
case 5339: return true;
case 5341: return true;
case 5343: return true;
case 5430: return true;
case 5448: return true;
case 5608: return true;
case 5648: return true;
case 5650: return true;
case 5652: return true;
case 5654: return true;
case 5656: return true;
case 5658: return true;
case 5660: return true;
case 5662: return true;
case 5663: return true;
case 5714: return true;
case 5767: return true;
case 5824: return true;
case 5826: return true;
case 5832: return true;
case 5833: return true;
case 5850: return true;
case 6084: return true;
case 6092: return true;
case 6095: return true;
case 6104: return true;
case 6132: return true;
case 6136: return true;
case 6137: return true;
case 6140: return true;
case 6143: return true;
case 6184: return true;
case 6275: return true;
case 6331: return true;
case 6340: return true;
case 6481: return true;
case 6624: return true;
case 6625: return true;
case 6626: return true;
case 6627: return true;
case 6628: return true;
case 6629: return true;
case 6630: return true;
case 6631: return true;
case 6632: return true;
case 6633: return true;
case 6634: return true;
case 6635: return true;
case 6636: return true;
case 6637: return true;
case 6638: return true;
case 6639: return true;
case 6641: return true;
case 6644: return true;
case 6645: return true;
case 6647: return true;
case 6649: return true;
case 6650: return true;
case 6654: return true;
case 6672: return true;
case 6673: return true;
case 6674: return true;
case 6675: return true;
case 6676: return true;
case 6677: return true;
case 6678: return true;
case 6679: return true;
case 6680: return true;
case 6681: return true;
case 6682: return true;
case 6683: return true;
case 6684: return true;
case 6685: return true;
case 6686: return true;
case 6687: return true;
case 6701: return true;
case 6735: return true;
case 6811: return true;
case 6886: return true;
case 6943: return true;
case 6945: return true;
case 6991: return true;
case 7150: return true;
case 7371: return true;
case 7417: return true;
case 7425: return true;
case 7430: return true;
case 7491: return true;
case 7572: return true;
case 7751: return true;
case 7793: return true;
case 8047: return true;
case 8244: return true;
case 8260: return true;
case 8275: return true;
case 8340: return true;
case 8353: return true;
case 8356: return true;
case 8357: return true;
case 8399: return true;
case 8420: return true;
case 8453: return true;
case 8591: return true;
case 8610: return true;
case 8639: return true;
case 8807: return true;
case 8809: return true;
case 8813: return true;
case 8855: return true;
case 8904: return true;
case 8907: return true;
case 8911: return true;
case 8916: return true;
case 8938: return true;
case 8944: return true;
case 8947: return true;
case 9078: return true;
case 9200: return true;
case 9243: return true;
case 9269: return true;
case 9282: return true;
case 9286: return true;
case 9291: return true;
case 9338: return true;
case 9354: return true;
case 9361: return true;
case 9459: return true;
case 9533: return true;
case 9598: return true;
case 9703: return true;
case 9706: return true;
case 9769: return true;
case 9798: return true;
case 9801: return true;
case 9802: return true;
case 9807: return true;
case 9849: return true;
case 9874: return true;
case 9878: return true;
case 9927: return true;
case 9940: return true;
case 9942: return true;
case 9945: return true;
case 9951: return true;
case 9960: return true;
case 10091: return true;
case 10134: return true;
case 10135: return true;
case 10138: return true;
case 10143: return true;
case 10164: return true;
case 10171: return true;
case 10203: return true;
case 10326: return true;
case 10363: return true;
case 10451: return true;
case 10505: return true;
case 10560: return true;
case 10564: return true;
case 10565: return true;
case 10568: return true;
case 10570: return true;
case 10573: return true;
case 10590: return true;
case 10591: return true;
case 10610: return true;
case 10613: return true;
case 10617: return true;
case 10619: return true;
case 10638: return true;
case 10639: return true;
case 10678: return true;
case 10829: return true;
case 10831: return true;
case 10921: return true;
case 10981: return true;
case 10982: return true;
case 10988: return true;
case 11024: return true;
case 11034: return true;
case 11042: return true;
case 11044: return true;
case 11048: return true;
case 11050: return true;
case 11052: return true;
case 11054: return true;
case 11091: return true;
case 11093: return true;
case 11095: return true;
case 11097: return true;
case 11099: return true;
case 11102: return true;
case 11141: return true;
case 11148: return true;
case 11150: return true;
case 11197: return true;
case 11198: return true;
case 11246: return true;
case 11247: return true;
case 11294: return true;
case 11295: return true;
case 11411: return true;
case 11462: return true;
case 11525: return true;
case 11534: return true;
case 11535: return true;
case 11551: return true;
case 11630: return true;
case 11647: return true;
case 11738: return true;
case 11904: return true;
case 11918: return true;
case 12032: return true;
case 12033: return true;
case 12034: return true;
case 12035: return true;
case 12042: return true;
case 12043: return true;
case 12044: return true;
case 12045: return true;
case 12046: return true;
case 12047: return true;
case 12117: return true;
case 12126: return true;
case 12127: return true;
case 12147: return true;
case 12160: return true;
case 12162: return true;
case 12164: return true;
case 12166: return true;
case 12168: return true;
case 12174: return true;
case 12271: return true;
case 12276: return true;
case 12278: return true;
case 12280: return true;
case 12282: return true;
case 12286: return true;
case 12325: return true;
case 12329: return true;
case 12334: return true;
case 12335: return true;
case 12371: return true;
case 12373: return true;
case 12375: return true;
case 12381: return true;
case 12382: return true;
case 12383: return true;
case 12604: return true;
case 12607: return true;
case 12608: return true;
case 12610: return true;
case 12613: return true;
case 12614: return true;
case 12615: return true;
case 12618: return true;
case 12621: return true;
case 12622: return true;
case 12623: return true;
case 12783: return true;
case 12800: return true;
case 12801: return true;
case 12802: return true;
case 12803: return true;
case 12805: return true;
case 12806: return true;
case 12807: return true;
case 12809: return true;
case 12810: return true;
case 12811: return true;
case 12813: return true;
case 12814: return true;
case 12815: return true;
case 12818: return true;
case 12821: return true;
case 12822: return true;
case 12823: return true;
case 12827: return true;
case 12828: return true;
case 12831: return true;
case 12848: return true;
case 12850: return true;
case 12853: return true;
case 12856: return true;
case 12859: return true;
case 12860: return true;
case 12862: return true;
case 12870: return true;
case 12871: return true;
case 12874: return true;
case 12879: return true;
case 12914: return true;
case 12915: return true;
case 12918: return true;
case 12921: return true;
case 12923: return true;
case 12927: return true;
case 12988: return true;
case 13049: return true;
case 13055: return true;
case 13056: return true;
case 13058: return true;
case 13059: return true;
case 13064: return true;
case 13068: return true;
case 13120: return true;
case 13122: return true;
case 13123: return true;
case 13124: return true;
case 13126: return true;
case 13128: return true;
case 13133: return true;
case 13156: return true;
case 13248: return true;
case 13249: return true;
case 13280: return true;
case 13282: return true;
case 13284: return true;
case 13286: return true;
case 13287: return true;
case 13288: return true;
case 13290: return true;
case 13292: return true;
case 13293: return true;
case 13294: return true;
case 13392: return true;
case 13394: return true;
case 13395: return true;
case 13398: return true;
case 13399: return true;
case 13400: return true;
case 13402: return true;
case 13403: return true;
case 13406: return true;
case 13407: return true;
case 13433: return true;
case 13458: return true;
case 13459: return true;
case 13460: return true;
case 13462: return true;
case 13463: return true;
case 13465: return true;
case 13466: return true;
case 13467: return true;
case 13469: return true;
case 13470: return true;
case 13471: return true;
case 13477: return true;
case 13478: return true;
case 13479: return true;
case 13481: return true;
case 13485: return true;
case 13486: return true;
case 13487: return true;
case 13508: return true;
case 13509: return true;
case 13510: return true;
case 13512: return true;
case 13513: return true;
case 13517: return true;
case 13521: return true;
case 13524: return true;
case 13525: return true;
case 13553: return true;
case 13556: return true;
case 13558: return true;
case 13559: return true;
case 13574: return true;
case 13578: return true;
case 13589: return true;
case 13786: return true;
case 13798: return true;
case 13828: return true;
case 13835: return true;
case 13848: return true;
case 13863: return true;
case 13903: return true;
case 13921: return true;
case 13926: return true;
case 13933: return true;
case 13953: return true;
case 13968: return true;
case 13969: return true;
case 13972: return true;
case 13973: return true;
case 13975: return true;
case 13977: return true;
case 13980: return true;
case 13981: return true;
case 13983: return true;
case 14058: return true;
case 14060: return true;
case 14081: return true;
case 14116: return true;
case 14117: return true;
case 14119: return true;
case 14125: return true;
case 14153: return true;
case 14157: return true;
case 14193: return true;
case 14207: return true;
case 14228: return true;
case 14233: return true;
case 14234: return true;
case 14236: return true;
case 14237: return true;
case 14239: return true;
case 14249: return true;
case 14303: return true;
case 14306: return true;
case 14307: return true;
case 14308: return true;
case 14309: return true;
case 14311: return true;
case 14312: return true;
case 14314: return true;
case 14315: return true;
case 14317: return true;
case 14318: return true;
case 14319: return true;
case 14325: return true;
case 14362: return true;
case 14363: return true;
case 14373: return true;
case 14374: return true;
case 14400: return true;
case 14403: return true;
case 14406: return true;
case 14409: return true;
case 14442: return true;
case 14447: return true;
case 14498: return true;
case 14501: return true;
case 14502: return true;
case 14503: return true;
case 14505: return true;
case 14544: return true;
case 14545: return true;
case 14549: return true;
case 14550: return true;
case 14553: return true;
case 14554: return true;
case 14557: return true;
case 14558: return true;
case 14592: return true;
case 14594: return true;
case 14596: return true;
case 14597: return true;
case 14598: return true;
case 14602: return true;
case 14607: return true;
case 14608: return true;
case 14610: return true;
case 14611: return true;
case 14616: return true;
case 14617: return true;
case 14618: return true;
case 14620: return true;
case 14621: return true;
case 14623: return true;
case 14636: return true;
case 14638: return true;
case 14639: return true;
case 14640: return true;
case 14647: return true;
case 14660: return true;
case 14661: return true;
case 14662: return true;
case 14663: return true;
case 14668: return true;
case 14669: return true;
case 14672: return true;
case 14676: return true;
case 14678: return true;
case 14679: return true;
case 14688: return true;
case 14689: return true;
case 14692: return true;
case 14694: return true;
case 14695: return true;
case 14697: return true;
case 14700: return true;
case 14701: return true;
case 14702: return true;
case 14703: return true;
case 14704: return true;
case 14709: return true;
case 14711: return true;
case 14712: return true;
case 14713: return true;
case 14721: return true;
case 14724: return true;
case 14725: return true;
case 14726: return true;
case 14729: return true;
case 14731: return true;
case 14732: return true;
case 14733: return true;
case 14734: return true;
case 14735: return true;
case 14736: return true;
case 14737: return true;
case 14738: return true;
case 14739: return true;
case 14741: return true;
case 14742: return true;
case 14744: return true;
case 14746: return true;
case 14747: return true;
case 14750: return true;
case 14752: return true;
case 14754: return true;
case 14755: return true;
case 14756: return true;
case 14757: return true;
case 14758: return true;
case 14760: return true;
case 14761: return true;
case 14762: return true;
case 14763: return true;
case 14764: return true;
case 14765: return true;
case 14766: return true;
case 14767: return true;
case 14772: return true;
case 14773: return true;
case 14774: return true;
case 14775: return true;
case 14777: return true;
case 14779: return true;
case 14783: return true;
case 14839: return true;
case 14862: return true;
case 14867: return true;
case 14874: return true;
case 14895: return true;
case 14977: return true;
case 14980: return true;
case 14986: return true;
case 14991: return true;
case 15009: return true;
case 15012: return true;
case 15015: return true;
case 15020: return true;
case 15021: return true;
case 15025: return true;
case 15027: return true;
case 15033: return true;
case 15096: return true;
case 15137: return true;
case 15140: return true;
case 15143: return true;
case 15149: return true;
case 15152: return true;
case 15154: return true;
case 15155: return true;
case 15157: return true;
case 15161: return true;
case 15164: return true;
case 15214: return true;
case 15217: return true;
case 15266: return true;
case 15272: return true;
case 15276: return true;
case 15280: return true;
case 15281: return true;
case 15282: return true;
case 15285: return true;
case 15286: return true;
case 15287: return true;
case 15288: return true;
case 15289: return true;
case 15290: return true;
case 15291: return true;
case 15292: return true;
case 15293: return true;
case 15295: return true;
case 15296: return true;
case 15298: return true;
case 15299: return true;
case 15318: return true;
case 15323: return true;
case 15326: return true;
case 15327: return true;
case 15331: return true;
case 15332: return true;
case 15335: return true;
case 15336: return true;
case 15338: return true;
case 15339: return true;
case 15340: return true;
case 15342: return true;
case 15343: return true;
case 15345: return true;
case 15348: return true;
case 15350: return true;
case 15354: return true;
case 15359: return true;
case 15367: return true;
case 15371: return true;
case 15398: return true;
case 15445: return true;
case 15460: return true;
case 15462: return true;
case 15463: return true;
case 15469: return true;
case 15476: return true;
case 15493: return true;
case 15526: return true;
case 15624: return true;
case 15670: return true;
case 15690: return true;
case 15718: return true;
case 15719: return true;
case 15754: return true;
case 15777: return true;
case 15783: return true;
case 15785: return true;
case 15821: return true;
case 15861: return true;
case 15872: return true;
case 15877: return true;
case 15908: return true;
case 15936: return true;
case 15945: return true;
case 15947: return true;
case 15949: return true;
case 15955: return true;
case 15958: return true;
case 15963: return true;
case 15965: return true;
case 15970: return true;
case 15971: return true;
case 15972: return true;
case 15982: return true;
case 15983: return true;
case 16018: return true;
case 16019: return true;
case 16030: return true;
case 16031: return true;
case 16054: return true;
case 16088: return true;
case 16089: return true;
case 16090: return true;
case 16091: return true;
case 16094: return true;
case 16215: return true;
case 16258: return true;
case 16268: return true;
case 16304: return true;
case 16306: return true;
case 16307: return true;
case 16308: return true;
case 16311: return true;
case 16317: return true;
case 16349: return true;
case 16388: return true;
case 16395: return true;
case 16422: return true;
case 16428: return true;
case 16429: return true;
case 16431: return true;
case 16440: return true;
case 16442: return true;
case 16443: return true;
case 16446: return true;
case 16467: return true;
case 16468: return true;
case 16469: return true;
case 16473: return true;
case 16474: return true;
case 16475: return true;
case 16477: return true;
case 16478: return true;
case 16479: return true;
case 16482: return true;
case 16483: return true;
case 16488: return true;
case 16489: return true;
case 16491: return true;
case 16492: return true;
case 16493: return true;
case 16496: return true;
case 16497: return true;
case 16500: return true;
case 16501: return true;
case 16502: return true;
case 16506: return true;
case 16508: return true;
case 16510: return true;
case 16511: return true;
case 16542: return true;
case 16543: return true;
case 16590: return true;
case 16591: return true;
case 16622: return true;
case 16623: return true;
case 16634: return true;
case 16670: return true;
case 16671: return true;
case 16784: return true;
case 16785: return true;
case 16789: return true;
case 16818: return true;
case 16819: return true;
case 16821: return true;
case 16824: return true;
case 16825: return true;
case 16826: return true;
case 16827: return true;
case 16865: return true;
case 16867: return true;
case 16886: return true;
case 16907: return true;
case 16910: return true;
case 16911: return true;
case 16927: return true;
case 17046: return true;
case 17088: return true;
case 17090: return true;
case 17091: return true;
case 17097: return true;
case 17107: return true;
case 17117: return true;
case 17202: return true;
case 17204: return true;
case 17206: return true;
case 17233: return true;
case 17236: return true;
case 17237: return true;
case 17291: return true;
case 17295: return true;
case 17396: return true;
case 17452: return true;
case 17463: return true;
case 17537: return true;

#if 0
case 1060: return true;
case 1066: return true;
case 1275: return true;
case 1681: return true;
case 1684: return true;
case 1759: return true;
case 1960: return true;
case 2122: return true;
case 2777: return true;
case 3134: return true;
case 3484: return true;
case 3543: return true;
case 3848: return true;
case 3857: return true;
case 3877: return true;
case 3883: return true;
case 3954: return true;
case 4025: return true;
case 4046: return true;
case 4604: return true;
case 5327: return true;
case 5607: return true;
case 5825: return true;
case 5835: return true;
case 6110: return true;
case 6339: return true;
case 6648: return true;
case 6700: return true;
case 6799: return true;
case 6800: return true;
case 6801: return true;
case 6803: return true;
case 6804: return true;
case 6805: return true;
case 6806: return true;
case 6808: return true;
case 7088: return true;
case 7873: return true;
case 8810: return true;
case 8811: return true;
case 8848: return true;
case 8853: return true;
case 8902: return true;
case 8905: return true;
case 8909: return true;
case 8942: return true;
case 8945: return true;
case 8953: return true;
case 9336: return true;
case 9343: return true;
case 9796: return true;
case 9799: return true;
case 9805: return true;
case 9870: return true;
case 9876: return true;
case 10089: return true;
case 10132: return true;
case 10133: return true;
case 10136: return true;
case 10141: return true;
case 10821: return true;
case 10827: return true;
case 10979: return true;
case 10980: return true;
case 10986: return true;
case 11022: return true;
case 11032: return true;
case 11046: return true;
case 11086: return true;
case 11100: return true;
case 11139: return true;
case 11140: return true;
case 11142: return true;
case 11144: return true;
case 11146: return true;
case 11195: return true;
case 11409: return true;
case 11628: return true;
case 11645: return true;
case 11948: return true;
case 12061: return true;
case 12115: return true;
case 12124: return true;
case 12125: return true;
case 12605: return true;
case 12612: return true;
case 12616: return true;
case 12781: return true;
case 12799: return true;
case 12804: return true;
case 12808: return true;
case 12812: return true;
case 12816: return true;
case 12824: return true;
case 12846: return true;
case 12854: return true;
case 12868: return true;
case 12869: return true;
case 12872: return true;
case 12877: return true;
case 12919: return true;
case 12925: return true;
case 12986: return true;
case 13054: return true;
case 13062: return true;
case 13118: return true;
case 13154: return true;
case 13246: return true;
case 13247: return true;
case 13278: return true;
case 13291: return true;
case 13401: return true;
case 13404: return true;
case 13431: return true;
case 13474: return true;
case 13475: return true;
case 13483: return true;
case 13506: return true;
case 13507: return true;
case 13515: return true;
case 13519: return true;
case 13523: return true;
case 13557: return true;
case 13576: return true;
case 13784: return true;
case 13826: return true;
case 13833: return true;
case 13901: return true;
case 13919: return true;
case 13931: return true;
case 14056: return true;
case 14115: return true;
case 14123: return true;
case 14155: return true;
case 14304: return true;
case 14323: return true;
case 14360: return true;
case 14361: return true;
case 14371: return true;
case 14394: return true;
case 14401: return true;
case 14407: return true;
case 14499: return true;
case 14500: return true;
case 14547: return true;
case 14552: return true;
case 14590: return true;
case 14595: return true;
case 14603: return true;
case 14615: return true;
case 14658: return true;
case 14659: return true;
case 14666: return true;
case 14667: return true;
case 14670: return true;
case 14674: return true;
case 14686: return true;
case 14687: return true;
case 14707: return true;
case 14740: return true;
case 14748: return true;
case 14751: return true;
case 14759: return true;
case 14771: return true;
case 14837: return true;
case 14872: return true;
case 14975: return true;
case 14984: return true;
case 14989: return true;
case 15007: return true;
case 15013: return true;
case 15018: return true;
case 15023: return true;
case 15031: return true;
case 15092: return true;
case 15135: return true;
case 15138: return true;
case 15141: return true;
case 15147: return true;
case 15159: return true;
case 15162: return true;
case 15264: return true;
case 15270: return true;
case 15278: return true;
case 15297: return true;
case 15314: return true;
case 15316: return true;
case 15324: return true;
case 15329: return true;
case 15346: return true;
case 15352: return true;
case 15357: return true;
case 15396: return true;
case 15458: return true;
case 15467: return true;
case 15474: return true;
case 15491: return true;
case 15524: return true;
case 15668: return true;
case 15716: return true;
case 15775: return true;
case 15859: return true;
case 15868: return true;
case 15875: return true;
case 15906: return true;
case 15953: return true;
case 15976: return true;
case 15977: return true;
case 16024: return true;
case 16052: return true;
case 16086: return true;
case 16087: return true;
case 16092: return true;
case 16302: return true;
case 16338: return true;
case 16344: return true;
case 16386: return true;
case 16420: return true;
case 16426: return true;
case 16427: return true;
case 16438: return true;
case 16441: return true;
case 16444: return true;
case 16463: return true;
case 16466: return true;
case 16471: return true;
case 16480: return true;
case 16486: return true;
case 16490: return true;
case 16494: return true;
case 16504: return true;
case 16507: return true;
case 16817: return true;
case 16863: return true;
case 16905: return true;
case 16908: return true;
case 16914: return true;
case 17089: return true;
case 17095: return true;
case 17105: return true;
case 17394: return true;
case 17440: return true;
case 17450: return true;
case 17455: return true;
case 17459: return true;
case 17461: return true;
case 17535: return true;
#endif

		default: return false;
		}
	}

	bool bad_cluster(AliVCluster *c)
	{
		for (Int_t j = 0; j < c->GetNCells(); j++) {
			const Int_t cell_id = c->GetCellsAbsId()[j];
			if (bad_cell(cell_id)) {
				return true;
			}
		}
		return false;
	}

}

void AliAnalysisTaskCalibEmcal::init_temperature_calib(const Int_t run_number)
{
	AliOADBContainer *cont_rftd = new AliOADBContainer("");
  
	cont_rftd->InitFromFile(
		"$ALICE_PHYSICS/OADB/EMCAL/"
		"EMCALTemperatureCorrCalib.root",
		"AliEMCALRunDepTempCalibCorrections");
  
	TH1S *htd = (TH1S*)cont_rftd->GetObject(run_number); 
  
	// If it did not exist for this run, get closes one
	if (!htd) {
		AliWarning(Form("No TemperatureCorrCalib Objects for run: %d",run_number));

		// let's get the closest runnumber instead then..
		Int_t lower = 0;
		Int_t ic = 0;
		Int_t maxEntry = cont_rftd->GetNumberOfEntries();

		while ((ic < maxEntry) && (cont_rftd->UpperLimit(ic) < run_number)) {
			lower = ic;
			ic++;
		}
    
		Int_t closest = lower;
		if ((ic<maxEntry) &&
			(cont_rftd->LowerLimit(ic) - run_number) <
			(run_number - cont_rftd->UpperLimit(lower)) ) {
			closest = ic;
		}
    
		AliWarning(Form("TemperatureCorrCalib Objects found closest id %d from run: %d",
						closest, cont_rftd->LowerLimit(closest)));
      
		htd = (TH1S*) cont_rftd->GetObjectByIndex(closest);
	}

	// Fill parameters
	if (htd) {
		AliInfo("Recalibrate (Temperature) EMCAL");

		Int_t nSM = _emcal_geometry->GetNumberOfSuperModules();

		for (Int_t ism = 0; ism < nSM; ++ism) {
			for (Int_t icol = 0; icol < 48; ++icol) {
				for (Int_t irow = 0; irow < 24; ++irow) {
					Float_t factor =
						_reco_util->GetEMCALChannelRecalibrationFactor(ism,icol,irow);

					Int_t absID =
						_emcal_geometry->GetAbsCellIdFromCellIndexes(ism, irow, icol); // original calibration factor
          
					AliDebug(3, Form(" ism %d, icol %d, irow %d,absID %d - Calib factor %1.5f - ",ism, icol, irow, absID, factor));

					factor *= htd->GetBinContent(absID) / 10000.0; // correction dependent on T
          
					_reco_util->SetEMCALChannelRecalibrationFactor(ism, icol, irow, factor);
            
					AliDebug(3,Form("T  factor %1.5f - final factor %1.5f",
									htd->GetBinContent(absID) / 10000.0,
									_reco_util->GetEMCALChannelRecalibrationFactor(ism, icol, irow)));
				} // columns
			} // rows 
		} // SM loop
	}
	else AliInfo("Do NOT recalibrate EMCAL with T variations, no params TH1");
  
	delete cont_rftd; 
}


void AliAnalysisTaskCalibEmcal::UserExec(Option_t *option)
{
	AliVEvent *event = InputEvent();
    AliESDEvent *esd_event = dynamic_cast<AliESDEvent *>(event);
    AliAODEvent *aod_event = dynamic_cast<AliAODEvent *>(event);

	const Int_t run_number = event->GetRunNumber();
	const Int_t bunch_crossing = event->GetBunchCrossNumber();

#if 0
	static const char *method[] = {
		"V0M", "V0A", "V0C", "V0MEq", "V0AEq", "V0CEq", "FMD", "TRK",
		"TKL", "CL0", "CL1", "CND", "ZNA", "ZNC", "ZPA", "ZPC", "NPA",
		"V0MvsFMD", "TKLvsV0M", "ZEMvsZDC", "V0Mtrue", "V0Atrue",
		"V0Ctrue", "V0MEqtrue", "V0AEqtrue", "V0CEqtrue", "FMDtrue",
		"TRKtrue", "TKLtrue", "CL0true", "CL1true", "CNDtrue",
		"ZNAtrue", "ZNCtrue", "ZPAtrue", "ZPCtrue", NULL
	};

	for (const char **p = method; *p != NULL; p++) {
		fprintf(stdout, "%s:%d: %s %f\n", __FILE__, __LINE__, *p,
				event->GetCentrality()->GetCentralityPercentile(*p));
	}
#endif

	if (_emcal_geometry == NULL) {
		TGeoManager::Import(
			"$ALICE_PHYSICS/OADB/EMCAL/geometry_2015.root");
		_emcal_geometry =
			AliEMCALGeometry::GetInstance(_emcal_geometry_name);
		init_temperature_calib(run_number);
	}

#if 0
	if (_emcal_geometry != NULL) {
		for (size_t i = 0; i < _ncell; i++) {
			_emcal_geometry->PrintCellIndexes(i, 1, "");
		}
	}
#endif

	if (_emcal_pi0_calibration_factor.empty()) {
		TFile *file = TFile::Open(
			"alien:///alice/cern.ch/user/y/ylai/"
			"multiplyPi0CalibrationFactors_TextToHisto_"
			"EMCALcoeffs2012pizOnlyNoEspectra_"
			"DCALandThirdsAllOne.root");

		if (file != NULL) {
			for (size_t i = 0; i < 20; i++) {
				_emcal_pi0_calibration_factor.push_back(
					static_cast<TH2D *>(file->Get(
						Form("EMCALRecalFactors_SM%d", i))));
			}
		}
	}

    TString trigger_class =
		esd_event ? esd_event->GetFiredTriggerClasses() :
		aod_event ? aod_event->GetFiredTriggerClasses() :
		"";

#if 0
	fprintf(stdout, "%s:%d: trigger_class = %s\n", __FILE__, __LINE__, trigger_class.Data());
#endif
#if 0
	if (!trigger_class.Contains(_trigger_name)) {
		return;
	}
#endif

	_histogram_run_number_event_triggered->Fill(run_number);

	TRefArray *calo_cluster = new TRefArray();

	event->GetEMCALClusters(calo_cluster);

	AliVCaloCells *emcal_cell = event->GetEMCALCells();
	Double_t primary_vertex[3];

	event->GetPrimaryVertex()->GetXYZ(primary_vertex);

	for (Int_t i = 0; i < calo_cluster->GetEntriesFast(); i++) {
		AliVCluster *c = (AliVCluster *)calo_cluster->At(i);

#define RECAL_CLUSTER(c)										\
		_reco_util->RecalculateClusterPosition(					\
			_emcal_geometry, emcal_cell, (c));					\
		(c)->SetE(												\
			_reco_util->CorrectClusterEnergyLinearity(c));		\
		_reco_util->RecalibrateClusterEnergy(					\
			_emcal_geometry, (c), emcal_cell);					\
		_reco_util->RecalculateClusterShowerShapeParameters(	\
			_emcal_geometry, emcal_cell, (c));					\
		_reco_util->RecalculateClusterPID(c);

		RECAL_CLUSTER(c);
	}

	unsigned int multiplicity_0_5 = 0;

	for (Int_t i = 0; i < calo_cluster->GetEntriesFast(); i++) {
		AliVCluster *c = (AliVCluster *)calo_cluster->At(i);
		if ((c)->E() >= 0.5) {
			multiplicity_0_5++;
		}
	}
	_histogram_multiplicity_0_5->Fill(multiplicity_0_5);

#define CLUSTER_CUT_PASS(c)								\
		((c)->GetTOF() * 1e+9 >= _timing_min &&			\
			(c)->GetTOF() * 1e+9 < _timing_max &&		\
			(c)->E() >= _pi0_e_min &&					\
			(c)->E() < _pi0_e_max &&					\
			(c)->GetNCells() >= _pi0_ncell_min &&		\
			(c)->GetM02() >= _pi0_m02_min &&			\
			(c)->GetM02() < _pi0_m02_max)
 // &&				\
 // 			_reco_util->IsGoodCluster(					\
 // 				(c), _emcal_geometry, emcal_cell,		\
 // 				bunch_crossing) &&						\
 // 			!_reco_util->ClusterContainsBadChannel(	\
 // 				_emcal_geometry, (c)->GetCellsAbsId(),	\
 // 				(c)->GetNCells()) &&					\
 // 			_reco_util->CheckCellFiducialRegion(		\
 // 				_emcal_geometry, (c), emcal_cell))

	const unsigned int multiplicity_0_5_bin =
		multiplicity_0_5 >= 25 ? 21 :
		multiplicity_0_5 >= 20 ? 20 :
		multiplicity_0_5;

	_mixed_event[multiplicity_0_5_bin].push_back(
		std::vector<std::pair<TLorentzVector, double> >());
	for (Int_t i = 0; i < calo_cluster->GetEntriesFast(); i++) {
		AliVCluster *c = (AliVCluster *)calo_cluster->At(i);
		TLorentzVector p;

		c->GetMomentum(p, primary_vertex);
		if (CLUSTER_CUT_PASS(c) && (!bad_cluster(c))) {
			_mixed_event[multiplicity_0_5_bin].back().push_back(
				std::pair<TLorentzVector, double>(p, c->GetTOF()));
		}
	}

	std::vector<std::pair<TLorentzVector, double> > mixed_event =
		_mixed_event[multiplicity_0_5_bin].front();

	if (_mixed_event[multiplicity_0_5_bin].size() >= 20) {
		_mixed_event[multiplicity_0_5_bin].pop_front();
	}

	for (Int_t i = 0; i < calo_cluster->GetEntriesFast(); i++) {
		AliVCluster *c = (AliVCluster *)calo_cluster->At(i);
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
			const Double_t scale = emcal_scale(cell_id);

			_histogram_cell_id_amplitude->Fill(cell_id, amplitude);
			_histogram_cell_id_amplitude_scale->Fill(
				cell_id, amplitude * scale);
			if (multiplicity_0_5 < 200) {
				_histogram_cell_id_amplitude_low_mult->Fill(cell_id, amplitude);
			}
			if (multiplicity_0_5 >= 300) {
				_histogram_cell_id_amplitude_high_mult->Fill(cell_id, amplitude);
			}
		}

#define RECAL_ENERGY(_e_recal, _cell_id_max, _e_cell_max, _c)	  \
		Double_t _e_recal = 0;									  \
		Int_t _cell_id_max = -1;								  \
		Double_t _e_cell_max = -DBL_MAX;						  \
		for (Int_t cell_index = 0;								  \
			 cell_index < (c)->GetNCells(); cell_index++) {		  \
			const Int_t cell_id =								  \
				(c)->GetCellsAbsId()[cell_index];				  \
			const Double_t amplitude =							  \
				emcal_cell->GetCellAmplitude(cell_id);			  \
			const Double_t scale = emcal_scale(cell_id);		  \
			const Double_t e_cell = amplitude * scale;			  \
			if (e_cell > e_cell_max) {							  \
				_cell_id_max = cell_id;							  \
				_e_cell_max = e_cell;							  \
			}													  \
			_e_recal += e_cell;									  \
		}

		RECAL_ENERGY(e_recal, cell_id_max, e_cell_max, c);

		TLorentzVector p;

		c->GetMomentum(p, primary_vertex);

		const TLorentzVector p_recal = p * (e_recal / p.E());

		unsigned int sm;
		unsigned int ieta;
		unsigned int iphi;

		to_sm_ieta_iphi(sm, ieta, iphi, cell_id_max);

		if (CLUSTER_CUT_PASS(c)) {
			for (Int_t j = 0; j < i; j++) {
				AliVCluster *c2 = (AliVCluster *)calo_cluster->At(j);

				if (CLUSTER_CUT_PASS(c2) &&
					(!bad_cluster(c2)) &&
					fabs(c2->GetTOF() - c->GetTOF()) * 1e+9 <
					_pi0_dtiming_max) {
					RECAL_ENERGY(e2_recal, cell_id_max_2,
								 e_cell_max_2, c2);

					TLorentzVector p2;

					c2->GetMomentum(p2, primary_vertex);

					const TLorentzVector p2_recal =
						p2 * (e2_recal / p2.E());

					unsigned int sm2;
					unsigned int ieta2;
					unsigned int iphi2;

					to_sm_ieta_iphi(sm2, ieta2, iphi2, cell_id_max_2);

					// if (sm2 == sm) {

					const Double_t mass = (p + p2).M();

					_histogram_cell_id_pi0_mass->Fill(
						cell_id_max, mass);

					_histogram_sm_ieta_pi0_mass->Fill(
						flat_sm_ieta(sm, ieta), mass);
					// _histogram_sm_ieta_pi0_mass->Fill(
					// 	flat_sm_ieta(sm2, ieta2), mass);
					_histogram_sm_iphi_pi0_mass->Fill(
						flat_sm_iphi(sm, iphi), mass);
					// _histogram_sm_iphi_pi0_mass->Fill(
					// 	flat_sm_iphi(sm2, iphi2), mass);

					if (mass >= 0.1 && mass < 0.2) {
						_histogram_pi0_deta_dphi->Fill(
							p2.Eta() - p.Eta(),
							angular_range_reduce(p2.Phi() - p.Phi()));
						if (ieta % 2 == 0) {
							_histogram_pi0_deta_dphi_ieta_even->Fill(
								p2.Eta() - p.Eta(),
								angular_range_reduce(p2.Phi() - p.Phi()));
						}
						else {
							_histogram_pi0_deta_dphi_ieta_odd->Fill(
								p2.Eta() - p.Eta(),
								angular_range_reduce(p2.Phi() - p.Phi()));
						}
					}

					const double dr =
						sqrt(std::pow(p2.Eta() - p.Eta(), 2) +
							 std::pow(angular_range_reduce(
								p2.Phi() - p.Phi()), 2));
					const double e_geom_mean = sqrt(p.E() * p2.E());

					// if (e_geom_mean < 1.5) {
					_histogram_pi0_e_geom_mean_dr->Fill(
						e_geom_mean, dr);
					if (ieta % 2 == 0) {
						_histogram_pi0_e_geom_mean_dr_ieta_even->Fill(
							e_geom_mean, dr);
					}
					else {
						_histogram_pi0_e_geom_mean_dr_ieta_odd->Fill(
							e_geom_mean, dr);
					}
					// }
					if (e_geom_mean < 1.4 && dr >= 0.05 && dr < 0.07) {
						_histogram_pi0_correlated_cell_id_e->Fill(
							cell_id_max, p.E());
						_histogram_pi0_correlated_cell_id_e->Fill(
							cell_id_max_2, p2.E());
					}

					const Double_t mass_recal = (p_recal + p2_recal).M();

					_histogram_cell_id_pi0_mass_recal->Fill(
						cell_id_max, mass_recal);

					_histogram_sm_ieta_pi0_mass_recal->Fill(
						flat_sm_ieta(sm, ieta), mass_recal);
					// _histogram_sm_ieta_pi0_mass_recal->Fill(
					// 	flat_sm_ieta(sm2, ieta2), mass_recal);
					_histogram_sm_iphi_pi0_mass_recal->Fill(
						flat_sm_iphi(sm, iphi), mass_recal);
					// _histogram_sm_iphi_pi0_mass_recal->Fill(
					// 	flat_sm_iphi(sm2, iphi2), mass_recal);

					if (!sm_ieta_warm(sm, ieta) &&
						!sm_ieta_warm(sm2, ieta2) ||
						sm_ieta_warm(sm, ieta) &&
						sm_ieta_warm(sm2, ieta2)) {
						_histogram_sm_ieta_pi0_mass_same->Fill(
							flat_sm_ieta(sm, ieta), mass);
						_histogram_sm_ieta_pi0_mass_same->Fill(
							flat_sm_ieta(sm2, ieta2), mass);
						_histogram_sm_ieta_pi0_mass_same_recal->Fill(
							flat_sm_ieta(sm, ieta), mass_recal);
						_histogram_sm_ieta_pi0_mass_same_recal->Fill(
							flat_sm_ieta(sm2, ieta2), mass_recal);
					}
					else {
						_histogram_sm_ieta_pi0_mass_cross->Fill(
							flat_sm_ieta(sm, ieta), mass);
						_histogram_sm_ieta_pi0_mass_cross->Fill(
							flat_sm_ieta(sm2, ieta2), mass);
						_histogram_sm_ieta_pi0_mass_cross_recal->Fill(
							flat_sm_ieta(sm, ieta), mass_recal);
						_histogram_sm_ieta_pi0_mass_cross_recal->Fill(
							flat_sm_ieta(sm2, ieta2), mass_recal);
					}

					// fprintf(stdout, "%s:%d: %f %f\n", __FILE__,
					// 		__LINE__, mass, mass_recal);

					// }
				}
			}

			for (size_t j = 0; j < mixed_event.size(); j++) {
				std::pair<TLorentzVector, double> c2 = mixed_event[j];

				if (fabs(c2.second - c->GetTOF()) * 1e+9 <
					_pi0_dtiming_max) {
					const TLorentzVector p2 = c2.first;
					const Double_t mass = (p + p2).M();

					_histogram_cell_id_pi0_mass->Fill(
						cell_id_max, mass);

					_histogram_sm_ieta_pi0_mass_mixed->Fill(
						flat_sm_ieta(sm, ieta), mass);
					_histogram_sm_iphi_pi0_mass_mixed->Fill(
						flat_sm_iphi(sm, iphi), mass);
				}
			}

		}
	}
}

AliEMCALRecoUtils *AliAnalysisTaskCalibEmcal::GetEMCALRecoUtils(void)
{
	return _reco_util;
}
