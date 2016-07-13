#include <cstdio>

void to_sm_ieta_iphi_eta_phi(unsigned int &sm, unsigned int &ieta,
							 unsigned int &iphi, double &eta,
							 double &phi, unsigned int n)
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

	const double coeff_eta[20][2] = {
		{0.65846, -0.0138514},
		{-0.0074445, -0.0138514},
		{0.65846, -0.0138514},
		{-0.0074445, -0.0138514},
		{0.65846, -0.0138514},
		{-0.0074445, -0.0138514},
		{0.65846, -0.0138514},
		{-0.0074445, -0.0138514},
		{0.65846, -0.0138514},
		{-0.0074445, -0.0138514},
		{0.657297, -0.0138276},
		{-0.00739987, -0.0138276},
		{0.658282, -0.0138382},
		{-0.229297, -0.0138382},
		{0.658282, -0.0138382},
		{-0.229297, -0.0138382},
		{0.658282, -0.0138382},
		{-0.229297, -0.0138382},
		{0.657297, -0.0138276},
		{-0.00741356, -0.0138267}
	};
	const double coeff_phi[20][2] = {
		{1.41567, 0.0134896},
		{1.41567, 0.0134896},
		{1.76473, 0.0134893},
		{1.76473, 0.0134893},
		{2.1138, 0.0134895},
		{2.1138, 0.0134895},
		{2.46286, 0.0134896},
		{2.46286, 0.0134896},
		{2.81193, 0.0134892},
		{2.81193, 0.0134892},
		{-3.1219, 0.0133945},
		{-3.1219, 0.0133945},
		{-1.72585, 0.0134832},
		{-1.72585, 0.0134832},
		{-1.37679, 0.0134831},
		{-1.37679, 0.0134831},
		{-1.02772, 0.0134828},
		{-1.02772, 0.0134828},
		{-0.678436, 0.0133937},
		{-0.678437, 0.0133935}
	};

	eta = coeff_eta[sm][0] + coeff_eta[sm][1] * ieta;
	phi = coeff_phi[sm][0] + coeff_phi[sm][1] * iphi;
}

int main(void)
{
	FILE *fp = fopen("mapping.txt", "r");
	unsigned int flat_index;
	unsigned int index_supermodule;
	unsigned int index_azimuth;
	unsigned int index_pseudorapidity;
	double pseudorapidity;
	double azimuth;

	while (fscanf(fp, "%u%u%u%u%lf%lf", &flat_index,
				  &index_supermodule, &index_azimuth,
				  &index_pseudorapidity, &pseudorapidity,
				  &azimuth) == 6) {
		unsigned int sm;
		unsigned int ieta;
		unsigned int iphi;
		double eta;
		double phi;

		to_sm_ieta_iphi_eta_phi(sm, ieta, iphi, eta, phi,
								flat_index);
		fprintf(stderr, "%u %u %u %u %u %u %u %lf %lf %lf %lf\n",
				flat_index, index_supermodule, sm,
				index_azimuth, iphi,
				index_pseudorapidity, ieta,
				pseudorapidity, eta, azimuth, phi);
	}

	fclose(fp);
}
