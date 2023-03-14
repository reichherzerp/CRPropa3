#include "crpropa/module/PropagationRS.h"

#include <sstream>
#include <stdexcept>
#include <vector>

namespace crpropa {
	void PropagationRS::tryStep(const Y &y, Y &out, Y &error, double h,
			ParticleState &particle, double z, double m, double q) const {
		out = dY(y.x, y.u, h, z, q, m);  // 1 step with h

		Y outHelp = dY(y.x, y.u, h/2, z, q, m);  // 2 steps with h/2
		Y outCompare = dY(outHelp.x, outHelp.u, h/2, z, q, m);

		error = errorEstimation(out.x , outCompare.x , h);
	}


	PropagationRS::Y PropagationRS::dY(Vector3d pos, Vector3d dir, double step,
			double z, double q, double m) const {
		// half leap frog step in the position
		pos += dir * step / 2.;

		// get B field at particle position
		Vector3d B = getFieldAtPosition(pos, z);

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
	PropagationRS::PropagationRS(ref_ptr<MagneticField> field, double fixedStep) :
			minStep(0) {
		setField(field);
		setTolerance(0.42);
		setMaximumStep(fixedStep);
		setMinimumStep(fixedStep);
	}


	// with adaptive step size
	PropagationRS::PropagationRS(ref_ptr<MagneticField> field, double tolerance, double minStep, double maxStep) :
			minStep(0) {
		setField(field);
		setTolerance(tolerance);
		setMaximumStep(maxStep);
		setMinimumStep(minStep);
	}


	void PropagationRS::process(Candidate *candidate) const {
		// save the new previous particle state
		ParticleState &current = candidate->current;
		candidate->previous = current;

		Y yIn(current.getPosition(), current.getDirection());

		// calculate charge of particle
		double q = current.getCharge();
		double step = maxStep;

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
			tryStep(yIn, yOut, yErr, step, current, z, m, q);
		} else {
			step = clip(candidate->getNextStep(), minStep, maxStep);
			newStep = step;
			double r = 42;  // arbitrary value

			// try performing step until the target error (tolerance) or the minimum/maximum step size has been reached
			while (true) {
				tryStep(yIn, yOut, yErr, step, current, z, m, q);
				r = yErr.u.getR() / tolerance;  // ratio of absolute direction error and tolerance
				if (r > 1) {  // large direction error relative to tolerance, try to decrease step size
					if (step == minStep)  // already minimum step size
						break;
					else {
						newStep = step * 0.95 * pow(r, -0.2);
						newStep = std::max(newStep, 0.1 * step); // limit step size decrease
						newStep = std::max(newStep, minStep); // limit step size to minStep
						step = newStep;
					}
				} else {  // small direction error relative to tolerance, try to increase step size
					if (step != maxStep) {  // only update once if maximum step size yet not reached
						newStep = step * 0.95 * pow(r, -0.2);
						newStep = std::min(newStep, 5 * step); // limit step size increase
						newStep = std::min(newStep, maxStep); // limit step size to maxStep
					}
					break;
				}
			}
		}

		current.setPosition(yOut.x);
		current.setDirection(yOut.u.getUnitVector());
		candidate->setCurrentStep(step);
		candidate->setNextStep(newStep);
	}


	void PropagationRS::setField(ref_ptr<MagneticField> f) {
		field = f;
	}


	ref_ptr<MagneticField> PropagationRS::getField() const {
		return field;
	}


	Vector3d PropagationRS::getFieldAtPosition(Vector3d pos, double z) const {
		Vector3d B(0, 0, 0);
		try {
			// check if field is valid and use the field vector at the
			// position pos with the redshift z
			if (field.valid())
				B = field->getField(pos, z);
		} catch (std::exception &e) {
			KISS_LOG_ERROR 	<< "PropagationRS: Exception in PropagationRS::getFieldAtPosition.\n"
					<< e.what();
		}	
		return B;
	}


	double PropagationRS::errorEstimation(const Vector3d x1, const Vector3d x2, double step) const {
		// compare the position after one step with the position after two steps with step/2.
		Vector3d diff = (x1 - x2);

		double S = diff.getR() / (step * (1 - 1/4.) );	// 1/4 = (1/2)²  number of steps for x1 divided by number of steps for x2 to the power of p (order)

		return S;
	}


	void PropagationRS::setTolerance(double tol) {
		if ((tol > 1) or (tol < 0))
			throw std::runtime_error(
					"PropagationRS: target error not in range 0-1");
		tolerance = tol;
	}


	void PropagationRS::setMinimumStep(double min) {
		if (min < 0)
			throw std::runtime_error("PropagationRS: minStep < 0 ");
		if (min > maxStep)
			throw std::runtime_error("PropagationRS: minStep > maxStep");
		minStep = min;
	}


	void PropagationRS::setMaximumStep(double max) {
		if (max < minStep)
			throw std::runtime_error("PropagationRS: maxStep < minStep");
		maxStep = max;
	}


	double PropagationRS::getTolerance() const {
		return tolerance;
	}


	double PropagationRS::getMinimumStep() const {
		return minStep;
	}


	double PropagationRS::getMaximumStep() const {
		return maxStep;
	}


	std::string PropagationRS::getDescription() const {
		std::stringstream s;
		s << "Propagation in magnetic fields using the adaptive Boris push method.";
		s << " Target error: " << tolerance;
		s << ", Minimum Step: " << minStep / kpc << " kpc";
		s << ", Maximum Step: " << maxStep / kpc << " kpc";
		return s.str();
	}
} // namespace crpropa
