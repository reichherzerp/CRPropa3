#include "crpropa/module/Scatter.h"
#include "crpropa/Random.h"

#include <sstream>
#include <stdexcept>
#include <vector>



namespace crpropa {
	
	/*
	Scatter::Y Scatter::dY(Vector3d pos, Vector3d dir, double step,
			double z, double q, double m) const {
		
		Random random;
		int seed = 1;
		if (seed != 0)
			random.seed(seed); // use given seed
		double phi = random.randUniform(0, 2 * M_PI);
		
		if (pitchAngleScattering) {
			// get B field at particle position
			Vector3d B = getFieldAtPosition(pos, z);

			// compute pitch angle
			double theta = dir.getAngleTo(B); // between 0 and pi
			double mu = cos(theta);
			
			double deltaT = step / c_light;
			// scatterRate < c/(4s)
			double deltaMu = 2 * sqrt(scatterRate * (1 - pow(mu, 2)) * deltaT) * sin(phi);

			if (deltaMu > 1 || deltaMu < -1) {
				throw std::runtime_error("Error: mu is larger than 1 or smaller than -1. Reduce scatter rate or step size.");
			}

			double deltaPhi = acos(deltaMu);

			Vector3d perp = dir.cross(B);
			dir = dir.getRotated(perp, deltaPhi);
			pos += dir * step;
		} else {
			double kappa = scatterRate;
			double deltaPhi = sqrt(2 * step * c_light / kappa); //* sin(phi);
			
			Vector3d rv = crpropa::Random::instance().randVector();
			Vector3d rotationAxis = dir.cross(rv);

			dir = dir.getRotated(rotationAxis, deltaPhi);
	
			pos += dir * step;
		}
		
		return Y(pos, dir);
	}*/

	Scatter::Scatter(double scatterRate) {
		setScatterRate(scatterRate);
	}

	void Scatter::process(Candidate *c) const {

		double step = c->getCurrentStep();
		Vector3d dir = c->current.getDirection();
		
		double deltaPhi = sqrt(2 * step * scatterRate / c_light);	
		Vector3d rv = crpropa::Random::instance().randVector();
		Vector3d rotationAxis = dir.cross(rv);
		dir = dir.getRotated(rotationAxis, deltaPhi);
		
		c->current.setDirection(dir);
		c->setNextStep(step);
	} 

	void Scatter::setScatterRate(double sRate) {
		scatterRate = sRate;
	}

	std::string Scatter::getDescription() const {
		std::stringstream s;
		s << "Scattering module to change the direction due to a general scattering process.";
		s << " Scatter rate: " << scatterRate << " 1/s";
		return s.str();
	}
} // namespace crpropa
