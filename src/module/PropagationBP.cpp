#include "crpropa/module/PropagationBP.h"
#include "crpropa/Random.h"

#include <sstream>
#include <stdexcept>
#include <vector>
#include <random>

namespace crpropa {
	void PropagationBP::tryStep(const Y &y, Y &out, Y &error, double h,
			ParticleState &particle, double z, double q, double m) const {
	}


	PropagationBP::Y PropagationBP::dY(Vector3d pos, Vector3d dir, double step,
			double z, double q, double m, Vector3d &B) const {
		// half leap frog step in the position
		pos += dir * step / 2.;

		// get B field at particle position
		B = getFieldAtPosition(pos, z);

		// Boris help vectors
		Vector3d t = B * q / 2 / m * step / c_light;
		Vector3d s = t * 2 / (1 + t.dot(t));
		Vector3d v_help;

		// Boris push
		v_help = dir + dir.cross(t);
		dir = dir + v_help.cross(s);

		// the other half leap frog step in the position
		pos += dir * step / 2.;
		return Y(pos, dir);
	}


	// with a fixed step size
	PropagationBP::PropagationBP(ref_ptr<MagneticField> field, double fixedStep) :
			minStep(0) {
		setField(field);
		setTolerance(0.42);
		setMaximumStep(fixedStep);
		setMinimumStep(fixedStep);
	}


	// with a fixed step size and additional MC scattering
	PropagationBP::PropagationBP(int test, ref_ptr<MagneticField> field, double fixedStep, double scatterRate, double minB) :
			minStep(0) {
		setField(field);
		setTolerance(0.42);
		setMaximumStep(fixedStep);
		setMinimumStep(fixedStep);
		setScatterRate(scatterRate);
		setMinB(minB);
	}


	// with adaptive step size
	PropagationBP::PropagationBP(ref_ptr<MagneticField> field, double tolerance, double minStep, double maxStep) :
			minStep(0) {
		setField(field);
		setTolerance(tolerance);
		setMaximumStep(maxStep);
		setMinimumStep(minStep);
	}


	void PropagationBP::process(Candidate *candidate) const {
		// save the new previous particle state
		ParticleState &current = candidate->current;
		candidate->previous = current;

		Y yIn(current.getPosition(), current.getDirection());

		// calculate charge of particle
		double q = current.getCharge();
		double step = maxStep;
		Vector3d B;

		// rectilinear propagation for neutral particles
		if (q == 0) {
			step = clip(candidate->getNextStep(), minStep, maxStep);
			current.setPosition(yIn.x + yIn.u * step);
			candidate->setCurrentStep(step);
			candidate->setNextStep(maxStep);
			return;
		}

		Y yOut, yErr;
		double newStep = step;
		double z = candidate->getRedshift();
		double m = current.getEnergy()/(c_light * c_light);

		// if minStep is the same as maxStep the adaptive algorithm with its error
		// estimation is not needed and the computation time can be saved:
		if (minStep == maxStep){
			yOut = dY(yIn.x, yIn.u, step, z, q, m, B);
		} 

		Vector3d dir = yOut.u.getUnitVector();

		// MC scatter only when field strength is above threshold
		if (B.getR() >= minB) {
			
			// Initialize random number generator
			int seed = 1;
			std::mt19937 gen(seed != 0 ? seed : std::time(nullptr));
			std::normal_distribution<double> gaussianDist(0.0, 1.0);

			double deltaPhi = sqrt(step * scatterRate / c_light) * gaussianDist(gen);
			Vector3d rv = Random::instance().randVector();
			
			Vector3d rotationAxis = dir.cross(rv);
			dir = dir.getRotated(rotationAxis, deltaPhi);
			
			current.setNrScatter(current.getNrScatter()+1);
		}

		current.setPosition(yOut.x);
		current.setDirection(dir);
		candidate->setCurrentStep(step);
		candidate->setNextStep(newStep);
	}


	void PropagationBP::setField(ref_ptr<MagneticField> f) {
		field = f;
	}


	ref_ptr<MagneticField> PropagationBP::getField() const {
		return field;
	}


	Vector3d PropagationBP::getFieldAtPosition(Vector3d pos, double z) const {
		Vector3d B(0, 0, 0);
		try {
			// check if field is valid and use the field vector at the
			// position pos with the redshift z
			if (field.valid())
				B = field->getField(pos, z);
		} catch (std::exception &e) {
			KISS_LOG_ERROR 	<< "PropagationBP: Exception in PropagationBP::getFieldAtPosition.\n"
					<< e.what();
		}	
		return B;
	}


	double PropagationBP::errorEstimation(const Vector3d x1, const Vector3d x2, double step) const {
		// compare the position after one step with the position after two steps with step/2.
		Vector3d diff = (x1 - x2);

		double S = diff.getR() / (step * (1 - 1/4.) );	// 1/4 = (1/2)²  number of steps for x1 divided by number of steps for x2 to the power of p (order)

		return S;
	}


	void PropagationBP::setTolerance(double tol) {
		if ((tol > 1) or (tol < 0))
			throw std::runtime_error(
					"PropagationBP: target error not in range 0-1");
		tolerance = tol;
	}


	void PropagationBP::setMinimumStep(double min) {
		if (min < 0)
			throw std::runtime_error("PropagationBP: minStep < 0 ");
		if (min > maxStep)
			throw std::runtime_error("PropagationBP: minStep > maxStep");
		minStep = min;
	}


	void PropagationBP::setMaximumStep(double max) {
		if (max < minStep)
			throw std::runtime_error("PropagationBP: maxStep < minStep");
		maxStep = max;
	}


	void PropagationBP::setScatterRate(double sRate) {
		scatterRate = sRate;
	}


	void PropagationBP::setMinB(double mB) {
		minB = mB;
	}


	double PropagationBP::getTolerance() const {
		return tolerance;
	}


	double PropagationBP::getMinimumStep() const {
		return minStep;
	}


	double PropagationBP::getMaximumStep() const {
		return maxStep;
	}


	std::string PropagationBP::getDescription() const {
		std::stringstream s;
		s << "Propagation in magnetic fields using the adaptive Boris push method.";
		s << " Target error: " << tolerance;
		s << ", Minimum Step: " << minStep / kpc << " kpc";
		s << ", Maximum Step: " << maxStep / kpc << " kpc";
		return s.str();
	}
} // namespace crpropa
