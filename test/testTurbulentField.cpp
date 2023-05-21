#include <stdexcept>

#include "crpropa/Grid.h"
#include "crpropa/Random.h"
#include "crpropa/Units.h"
#include "crpropa/Common.h"
#include "crpropa/GridTools.h"
#include "crpropa/magneticField/turbulentField/TurbulentField.h"
#include "crpropa/magneticField/turbulentField/GridTurbulence.h"
#include "crpropa/magneticField/turbulentField/PlaneWaveTurbulence.h"
#include "crpropa/magneticField/turbulentField/SimpleGridTurbulence.h"

#include "gtest/gtest.h"
#include <numeric>

using namespace crpropa;

//check problems brought up in https://github.com/CRPropa/CRPropa3/issues/322
TEST(testTurbulenceSpectrum, constructor) {
	double sIndex = 5./3.;
	double qIndex = 4.;
	double bendOver = 1.;
	double lMin = 1.; 
	double lMax = 10.; 
	double brms = 1*muG;
	auto spectrum = TurbulenceSpectrum(brms, lMin, lMax);
	EXPECT_DOUBLE_EQ(spectrum.getBrms(), brms);
	EXPECT_DOUBLE_EQ(spectrum.getLmin(), lMin);
	EXPECT_DOUBLE_EQ(spectrum.getLmax(), lMax);
	EXPECT_DOUBLE_EQ(spectrum.getLbendover(), bendOver); //default
	EXPECT_DOUBLE_EQ(spectrum.getSindex(), sIndex); //default
	EXPECT_DOUBLE_EQ(spectrum.getQindex(), qIndex); //default
}

TEST(testTurbulenceSpectrum, PW) {
    const double r = 1e-5 * pc;
    const int n = static_cast<int>(std::pow(10, 4));
    std::vector<double> sq_diffs(n);

    const double l_min = 1e-7 * pc;
    const double l_max = 1e5 * pc;
    const double brms = 1 * muG;
    const int n_wavemodes = 1024;

    auto spectrum = SimpleTurbulenceSpectrum(brms, l_min, l_max);
    auto turbulence = PlaneWaveTurbulence(spectrum, n_wavemodes, 1);

    Random &R = Random::instance();
    for (int i = 0; i < n; i++) {
        Vector3d randomDirection = R.randVector();
        Vector3d offset = randomDirection * l_max;
        Vector3d displacement = randomDirection * r;

        Vector3d b0 = turbulence.getField(offset);
        Vector3d b = turbulence.getField(offset + displacement);
        
        double diff = (b0 - b).getR();
        double sq_diff = diff * diff;
        sq_diffs[i] = sq_diff;
    }

    std::transform(sq_diffs.begin(), sq_diffs.end(), sq_diffs.begin(), [brms](double& d) { return d / (brms * brms); });
    double mean_s2_diff = std::accumulate(sq_diffs.begin(), sq_diffs.end(), 0.0) / sq_diffs.size();
    
    double theory_s2_diff = std::pow((r / l_max * 5), (2.0 / 3.0)) * 2;
    std::cout << "Mean of squared differences: " << mean_s2_diff << std::endl;
    std::cout << "Theory of squared differences: " << theory_s2_diff << std::endl;
}

TEST(testTurbulenceSpectrum, correlationLength) {
	double lMin = 0.00001; // not used for Lc
	double lMax = 9999999; // not used for Lc
	double lBo = 100;
	auto spectrum = TurbulenceSpectrum(1*muG, lMin, lMax, lBo);
	auto Lc = spectrum.getCorrelationLength();
    EXPECT_NEAR(Lc, 0.498*lBo, 0.001*lBo);
}

#ifdef CRPROPA_HAVE_FFTW3F

TEST(testSimpleGridTurbulence, oldFunctionForCrrelationLength) { //TODO: remove in future
	double lMin = 1*kpc;
	double lMax = 1*Gpc;
	double alpha = -11/3.;
	auto Lc = turbulentCorrelationLength(lMin, lMax, alpha);
    EXPECT_NEAR(Lc, lMax/5, 1*Mpc);
}

TEST(testVectorFieldGrid, Turbulence_bmean_brms) {
	// Test for zero mean: <B> = 0
	size_t n = 64;
	double spacing = 10 * Mpc / n;
	double Brms = 1;
	double lMin = 2 * spacing;
	double lMax = 8 * spacing;

	auto spectrum = SimpleTurbulenceSpectrum(Brms, lMin, lMax);
	auto gp = GridProperties(Vector3d(0, 0, 0), n, spacing);
    auto tf = SimpleGridTurbulence(spectrum, gp);
	auto grid = tf.getGrid();

	double precision = 1e-7;
	Vector3f bMean = meanFieldVector(grid);
	EXPECT_NEAR(0, bMean.x, precision);
	EXPECT_NEAR(0, bMean.y, precision);
	EXPECT_NEAR(0, bMean.z, precision);
	EXPECT_NEAR(1, rmsFieldStrength(grid), precision);
}

TEST(testVectorFieldGrid, Turbulence_seed) {
	// Test if seeding produces 2 identical fields
	size_t n = 64;
	double spacing = 1 * Mpc;
	double Brms = 1;
	double lMin = 2 * spacing;
	double lMax = 8 * spacing;
	int seed = 753;
	
	auto spectrum = SimpleTurbulenceSpectrum(Brms, lMin, lMax);

	auto gp1 = GridProperties(Vector3d(0, 0, 0), n, spacing);
    auto tf1 = SimpleGridTurbulence(spectrum, gp1,  seed);

	auto gp2 = GridProperties(Vector3d(0, 0, 0), n, spacing);
    auto tf2 = SimpleGridTurbulence(spectrum, gp2, seed);

	Vector3d pos(22 * Mpc);
	EXPECT_FLOAT_EQ(tf1.getField(pos).x, tf2.getField(pos).x);
}

TEST(testVectorFieldGrid, turbulence_Exceptions) {
	// Test exceptions
	size_t n = 64;
	double spacing = 10 * Mpc / n;
	double brms = 1;
	ref_ptr<Grid3f> grid = new Grid3f(Vector3d(0, 0, 0), n, spacing);

	// should be fine
	EXPECT_NO_THROW(initTurbulence(grid, brms, 2 * spacing, 8 * spacing));
	// lMin too small
	EXPECT_THROW(initTurbulence(grid, brms, 1.5 * spacing, 8 * spacing),
			std::runtime_error);
	// lMin > lMax
	EXPECT_THROW(initTurbulence(grid, brms, 8.1 * spacing, 8 * spacing),
			std::runtime_error);
	// lMax too large
	EXPECT_THROW(initTurbulence(grid, brms, 2 * spacing, 65 * spacing),
			std::runtime_error);
}

TEST(testGridTurbulence, Turbulence_seed) {
	// Test if seeding produces 2 identical fields
	size_t n = 64;
	double spacing = 1 * Mpc;
	double Brms = 1;
	double lMin = 2 * spacing;
	double lMax = 8 * spacing;
	double lBo = lMax/6;
	int seed = 137;
	
	auto spectrum = TurbulenceSpectrum(Brms, lMin, lMax, lBo);

	auto gp1 = GridProperties(Vector3d(0, 0, 0), n, spacing);
    auto tf1 = GridTurbulence(spectrum, gp1, seed);

	auto gp2 = GridProperties(Vector3d(0, 0, 0), n, spacing);
    auto tf2 = GridTurbulence(spectrum, gp2, seed);

	Vector3d pos(22 * Mpc);
	EXPECT_FLOAT_EQ(tf1.getField(pos).x, tf2.getField(pos).x);
}
#endif // CRPROPA_HAVE_FFTW3F

int main(int argc, char **argv) {
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
