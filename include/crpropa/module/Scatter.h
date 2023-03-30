#ifndef CRPROPA_Scatter_H
#define CRPROPA_Scatter_H

#include "crpropa/Module.h"
#include "crpropa/Units.h"
#include "crpropa/magneticField/MagneticField.h"
#include "kiss/logger.h"

namespace crpropa {
/**
 * \addtogroup Propagation
 * @{
 */

/**
 @class Scatter
 @brief Propagation through magnetic fields using the Boris method.

 This module solves the equations of motion of a relativistic charged particle when propagating through a magnetic field.\n
 It uses the Boris push integration method.\n
 It can be used with a fixed step size or an adaptive version which supports the step size control.
 The step size control tries to keep the relative error close to, but smaller than the designated tolerance.
 Additionally a minimum and maximum size for the steps can be set.
 For neutral particles a rectilinear propagation is applied and a next step of the maximum step size proposed.
 */
class Scatter: public Module {

private:
	double scatterRate;

public:
	/** Default constructor for the Boris push. It is constructed with a fixed step size.
	 * @param field
	 * @param scatterRate
	 */
	Scatter(ref_ptr<MagneticField> field = NULL, double scatterRate = 1);

	/** Propagates the particle. Is called once per iteration.
	 * @param candidate	 The Candidate is a passive object, that holds the information about the state of the cosmic ray and the simulation itself. */
	void process(Candidate *candidate) const;

	/** Set functions for the parameters of the class Scatter */
	void setScatterRate(double sRate);

	/** Get functions for the parameters of the class Scatter, similar to the set functions */
	std::string getDescription() const;
};
/** @}*/

} // namespace crpropa

#endif CRPROPA_Scatter_H
