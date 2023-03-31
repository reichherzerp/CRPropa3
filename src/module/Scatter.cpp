#include "crpropa/module/Scatter.h"
#include "crpropa/Random.h"

#include <sstream>
#include <stdexcept>
#include <vector>
#include <iostream>
#include <random>



namespace crpropa {
	
	Scatter::Scatter(double scatterRate) {
		setScatterRate(scatterRate);
	}

	void Scatter::process(Candidate *c) const {
		
		double step = c->getCurrentStep();
		Vector3d dir = c->current.getDirection();

		// Initialize random number generator
		int seed = 1;
    	std::mt19937 gen(seed != 0 ? seed : std::time(nullptr));
    	std::normal_distribution<double> gaussianDist(0.0, 1.0);

		double deltaPhi = sqrt(step * scatterRate / c_light) * gaussianDist(gen);
		Vector3d rv = Random::instance().randVector();
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
		s << "Scatter rate: " << scatterRate << " 1/s";
		return s.str();
	}
} // namespace crpropa
