#include "crpropa/module/Scatter.h"
#include "crpropa/Random.h"

#include <sstream>
#include <stdexcept>
#include <vector>
#include <iostream>
#include <random>



namespace crpropa {
	
	Scatter::Scatter(double scatterRate, double scatterRateExternal, double expansionFactor, double expansionStart) {
		setScatterRate(scatterRate, scatterRateExternal);
		setExpansion(expansionFactor, expansionStart);
	}

	void Scatter::process(Candidate *c) const {
		
		double step = c->getCurrentStep();
		double trajectoryLength = c->getTrajectoryLength();
		Vector3d dir = c->current.getDirection();
		Vector3d pos = c->current.getPosition();

		double currentScatterRate = scatterRate;

		if (pos.getDistanceTo(Vector3d()) > trajectoryLength * expansionFactor + expansionStart) {
			// the particle is outside of the bubble
			// 3D expansion with ceter (0,0,0) and constant expansion speed
			currentScatterRate = scatterRateExternal;
		}

		// Initialize random number generator
		int seed = 1;
    	std::mt19937 gen(seed != 0 ? seed : std::time(nullptr));
    	std::normal_distribution<double> gaussianDist(0.0, 1.0);

		double deltaPhi = sqrt(step * currentScatterRate / c_light) * gaussianDist(gen);
		Vector3d rv = Random::instance().randVector();
		Vector3d rotationAxis = dir.cross(rv);
		dir = dir.getRotated(rotationAxis, deltaPhi);
		
		c->current.setDirection(dir);
		c->setNextStep(step);
	} 

	void Scatter::setScatterRate(double sRate, double sRateExternal) {
		scatterRate = sRate;
		scatterRateExternal = sRateExternal;
	}

	void Scatter::setExpansion(double expFactor, double expStart) {
		expansionFactor = expFactor;
		expansionStart = expStart;
	}

	std::string Scatter::getDescription() const {
		std::stringstream s;
		s << "Scattering module to change the direction due to a general scattering process.";
		s << "Scatter rate: " << scatterRate << " 1/s";
		return s.str();
	}
} // namespace crpropa
