// This file contains an implementation of a vectorized cosine, which
// is based in part on the implementations in the library "SLEEF" by
// Naoki Shibata. SLEEF was used under the Boost Software License,
// Version 1.0. The original source file contained the following
// copyright notice:
//
//   //          Copyright Naoki Shibata 2010 - 2018.
//   // Distributed under the Boost Software License, Version 1.0.
//   //    (See accompanying file LICENSE.txt or copy at
//   //          http://www.boost.org/LICENSE_1_0.txt)
//
// SLEEF was used under the following license, which is not necessarily the
// license that applies to this file:
//
//         Boost Software License - Version 1.0 - August 17th, 2003
//
//         Permission is hereby granted, free of charge, to any person or
//         organization obtaining a copy of the software and accompanying
//         documentation covered by this license (the "Software") to use,
//         reproduce, display, distribute, execute, and transmit the Software,
//         and to prepare derivative works of the Software, and to permit
//         third-parties to whom the Software is furnished to do so, all subject
//         to the following:
//
//         The copyright notices in the Software and this entire statement,
//         including the above license grant, this restriction and the following
//         disclaimer, must be included in all copies of the Software, in whole
//         or in part, and all derivative works of the Software, unless such
//         copies or derivative works are solely in the form of
//         machine-executable object code generated by a source language
//         processor.
//
//         THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
//         EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
//         MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE, TITLE AND
//         NON-INFRINGEMENT. IN NO EVENT SHALL THE COPYRIGHT HOLDERS OR ANYONE
//         DISTRIBUTING THE SOFTWARE BE LIABLE FOR ANY DAMAGES OR OTHER
//         LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT
//         OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
//         THE SOFTWARE.

#include "crpropa/magneticField/turbulentField/PlaneWaveTurbulenceCylinder.h"
#include "crpropa/GridTools.h"
#include "crpropa/Random.h"
#include "crpropa/Units.h"

#include "kiss/logger.h"

#include <iostream>

#if defined(FAST_WAVES)
#if defined(__SSE__) && defined(__SSE2__) && defined(__SSE3__) && defined(__SSE4_1__) && defined(__SSE4_2__) && defined(__AVX__)
#define ENABLE_FAST_WAVES
#else
#error "FAST_WAVES is enabled, but it appears that not all required SIMD extensions are enabled in the compiler. Without these extensions, the FAST_WAVES optimization cannot be used. Please make sure that the SIMD_EXTENSIONS option in cmake matches the capabilities of your target CPU, or (if your target CPU does not support the required extensions), disable the FAST_WAVES flag in cmake."
#endif
#endif


#ifdef ENABLE_FAST_WAVES
#include <immintrin.h>
#endif

namespace crpropa {
#ifdef ENABLE_FAST_WAVES
// see
// https://stackoverflow.com/questions/49941645/get-sum-of-values-stored-in-m256d-with-sse-avx
double hsum_double_avx(__m256d v) {
	__m128d vlow = _mm256_castpd256_pd128(v);
	__m128d vhigh = _mm256_extractf128_pd(v, 1); // high 128
	vlow = _mm_add_pd(vlow, vhigh);              // reduce down to 128

	__m128d high64 = _mm_unpackhi_pd(vlow, vlow);
	return _mm_cvtsd_f64(_mm_add_sd(vlow, high64)); // reduce to scalar
}
#endif // defined(ENABLE_FAST_WAVES)

PlaneWaveTurbulenceCylinder::PlaneWaveTurbulenceCylinder(const TurbulenceSpectrum &spectrum,
                                         int Nm, int seed, std::string tType, const Vector3d c, double r, double d, double dFactor, double l, bool cons)
    : TurbulentField(spectrum), Nm(Nm) {

#ifdef ENABLE_FAST_WAVES
	KISS_LOG_INFO << "PlaneWaveTurbulenceCylinder: Using SIMD TD13 implementation"
	              << std::endl;

	// There used to be a cpuid check here, to see if the cpu running
	// this code would support SIMD (SSE + AVX). However, the library providing
	// the relevant function is no longer being used, and doing this manually
	// might be a bit too much work.
#endif

	if (Nm <= 1) {
		throw std::runtime_error(
		    "PlaneWaveTurbulenceCylinder: Nm <= 1. Specify at least two wavemodes in "
		    "order to generate the k distribution properly.");
	}

	Random random;
	if (seed != 0)
		random.seed(seed);

	double kmax = 2 * M_PI / spectrum.getLmin();
	double kmin = 2 * M_PI / spectrum.getLmax();

	xi = std::vector<Vector3d>(Nm, Vector3d(0.));
	kappa = std::vector<Vector3d>(Nm, Vector3d(0.));
	phi = std::vector<double>(Nm, 0.);
	costheta = std::vector<double>(Nm, 0.);
	beta = std::vector<double>(Nm, 0.);
	Ak = std::vector<double>(Nm, 0.);
	k = std::vector<double>(Nm, 0.);
    R = r;
    center = c;
    delta = d;
    decayFactor = dFactor;
    turbType = tType;
    length = l;
    constant = cons;


	double delta = log10(kmax / kmin);
	for (int i = 0; i < Nm; i++) {
		k[i] = pow(10, log10(kmin) + ((double)i) / ((double)(Nm - 1)) * delta);
	}

	// * compute Ak *

	double delta_k0 =
	    (k[1] - k[0]) / k[1]; // multiply this by k[i] to get delta_k[i]
	// Note: this is probably unnecessary since it's just a factor
	// and will get normalized out anyways. It's not like this is
	// performance-sensitive code though, and I don't want to change
	// anything numerical now.

	// For this loop, the Ak array actually contains Gk*delta_k (ie
	// non-normalized Ak^2). Normalization happens in a second loop,
	// once the total is known.
	double Ak2_sum = 0; // sum of Ak^2 over all k
    std::cout << tType << std::endl;
	for (int i = 0; i < Nm; i++) {
		double k = this->k[i];
		double kHat = k * spectrum.getLbendover();
		double Gk = spectrum.energySpectrum(k) * (1 + kHat * kHat);	// correct different implementation in TD 13 (eq. 5, missing + 1 in the denuminators exponent)
		Ak[i] = Gk * delta_k0 * k;
		Ak2_sum += Ak[i];

		// phi, costheta, and sintheta are for drawing vectors with
		// uniform distribution on the unit sphere.
		// This is similar to Random::randVector(): their t is our phi,
		// z is costheta, and r is sintheta. Our kappa is equivalent to
		// the return value of randVector(); however, TD13 then reuse
		// these values to generate a random vector perpendicular to kappa.
		double phi = random.randUniform(-M_PI, M_PI);
		double costheta = 0.0;
		if (tType == "3D") {
			costheta = random.randUniform(-1., 1.);
		} else if (tType == "slab") {
			costheta = 1.0;
		} else if (tType == "cylindrical") {
			costheta = 1.0      ;
		} 
		double sintheta = sqrt(1 - costheta * costheta);

		double alpha = 0.0;
		if (tType != "slab" && tType != "cylindrical") {
			alpha = random.randUniform(0, 2 * M_PI);
		}
		double beta = random.randUniform(0, 2 * M_PI);

		Vector3d kappa =
		    Vector3d(sintheta * cos(phi), sintheta * sin(phi), costheta);

		// NOTE: all other variable names match the ones from the TD13 paper.
		// However, our xi is actually their psi, and their xi isn't used at
		// all. (Though both can be used for the polarization vector, according
		// to the paper.) The reason for this discrepancy is that this code
		// used to be based on the original GJ99 paper, which provided only a
		// xi vector, and this xi happens to be almost the same as TD13's psi.
        Vector3d xi = Vector3d(costheta * cos(phi) * cos(alpha) + sin(phi) * sin(alpha),
		             costheta * sin(phi) * cos(alpha) - cos(phi) * sin(alpha),
		             -sintheta * cos(alpha));
        this->xi[i] = xi;
		this->kappa[i] = kappa;
		this->phi[i] = phi;
		this->costheta[i] = costheta;
		this->beta[i] = beta;
	}

	// Only in this loop are the actual Ak computed and stored.
	// This two-step process is necessary in order to normalize the values
	// properly.
	for (int i = 0; i < Nm; i++) {
		Ak[i] = sqrt(2 * Ak[i] / Ak2_sum) * spectrum.getBrms();
	}

#ifdef ENABLE_FAST_WAVES
	// * copy data into AVX-compatible arrays *
	//
	// AVX requires all data to be aligned to 256 bit, or 32 bytes, which is the
	// same as 4 double precision floating point numbers. Since support for
	// alignments this big seems to be somewhat tentative in C++ allocators,
	// we're aligning them manually by allocating a normal double array, and
	// then computing the offset to the first value with the correct alignment.
	// This is a little bit of work, so instead of doing it separately for each
	// of the individual data arrays, we're doing it once for one big array that
	// all of the component arrays get packed into.
	//
	// The other thing to keep in mind is that AVX always reads in units of 256
	// bits, or 4 doubles. This means that our number of wavemodes must be
	// divisible by 4. If it isn't, we simply pad it out with zeros. Since the
	// final step of the computation of each wavemode is multiplication by the
	// amplitude, which will be set to 0, these padding wavemodes won't affect
	// the result.

	avx_Nm = ((Nm + 4 - 1) / 4) * 4; // round up to next larger multiple of 4:
	                                 // align is 256 = 4 * sizeof(double) bit
	avx_data = std::vector<double>(itotal * avx_Nm + 3, 0.);

	// get the first 256-bit aligned element
	size_t size = avx_data.size() * sizeof(double);
	void *pointer = avx_data.data();
	align_offset =
	    (double *)std::align(32, 32, pointer, size) - avx_data.data();

	// copy into the AVX arrays
	for (int i = 0; i < Nm; i++) {
		avx_data[i + align_offset + avx_Nm * iAxi0] = Ak[i] * xi[i].x;
		avx_data[i + align_offset + avx_Nm * iAxi1] = Ak[i] * xi[i].y;
		avx_data[i + align_offset + avx_Nm * iAxi2] = Ak[i] * xi[i].z;

		// The cosine implementation computes cos(pi*x), so we'll divide out the
		// pi here.
		avx_data[i + align_offset + avx_Nm * ikkappa0] =
		    k[i] / M_PI * kappa[i].x;
		avx_data[i + align_offset + avx_Nm * ikkappa1] =
		    k[i] / M_PI * kappa[i].y;
		avx_data[i + align_offset + avx_Nm * ikkappa2] =
		    k[i] / M_PI * kappa[i].z;

		// We also need to divide beta by pi, since that goes into the argument
		// of the cosine as well.
		avx_data[i + align_offset + avx_Nm * ibeta] = beta[i] / M_PI;
	}
#endif // ENABLE_FAST_WAVES
}

Vector3d PlaneWaveTurbulenceCylinder::getField(const Vector3d &pos) const {
    Vector3d posPlane = Vector3d(pos.x, pos.y, 0);
    double scaling_factor = 1.0;
    if (length > 0.0 && pos.z > length) {
        return Vector3d(0.0);
    } else {
        if (!constant && length > 0.0) {
            scaling_factor = pos.z / length;
        }
    }

    
    double dist = posPlane.getDistanceTo(center);
    double smoothing_factor = 1.0;

    if (dist > R && turbType != "cylindrical") {
        double transition = 1 / (1 + std::exp(-(dist - R) / delta));
        smoothing_factor = (1 - transition) * std::exp(-(dist - R) / decayFactor);
    }

    

#ifndef ENABLE_FAST_WAVES
    Vector3d B(0.);
    for (int i = 0; i < Nm; i++) {
        double z_ = pos.dot(kappa[i]);
        B += xi[i] * Ak[i] * cos(k[i] * z_ + beta[i]);
    }

    if (turbType == "cylindrical") {
        Vector3d center_to_pos = posPlane - center;
        double x = center_to_pos.x;
        double y = center_to_pos.y;
        
        // solenoidal vector components
        double r = center_to_pos.getR(); // same as dist
        double eps = pow(10,-8);
        double a = 1.0;
        // v_x = (-y * k / (sqrt(x^2 + y^2) + eps)) * (1 - tanh((sqrt(x^2 + y^2) - R)/a))
        // v_y = (x * k / (sqrt(x^2 + y^2) + eps)) * (1 - tanh((sqrt(x^2 + y^2) - R)/a))
        double v_x = (-y / (r + eps)) * (1 - tanh((r - R)/a));
        double v_y = (x / (r + eps)) * (1 - tanh((r - R)/a));
        Vector3d b_solenoidal(v_x, v_y, 0);

        // Scale the rotated vector by the current magnitude of B
        double B_magnitude = B.getR();
        B.x = b_solenoidal.x * B_magnitude;
        B.y = b_solenoidal.y * B_magnitude;
    }

    return smoothing_factor * B * scaling_factor;

#else  // ENABLE_FAST_WAVES

	// Initialize accumulators
	//
	// There is one accumulator per component of the result vector.
	// Note that each accumulator contains four numbers. At the end of
	// the loop, each of these numbers will contain the sum of every
	// fourth wavemode, starting at a different offset. In the end, each
	// of the accumulator's numbers are added together (using
	// hsum_double_avx), resulting in the total sum for that component.

	__m256d acc0 = _mm256_setzero_pd();
	__m256d acc1 = _mm256_setzero_pd();
	__m256d acc2 = _mm256_setzero_pd();

	// broadcast position into AVX registers
	__m256d pos0 = _mm256_set1_pd(pos.x);
	__m256d pos1 = _mm256_set1_pd(pos.y);
	__m256d pos2 = _mm256_set1_pd(pos.z);

	for (int i = 0; i < avx_Nm; i += 4) {

		// Load data from memory into AVX registers:
		//  - the three components of the vector A * xi
		__m256d Axi0 =
		    _mm256_load_pd(avx_data.data() + i + align_offset + avx_Nm * iAxi0);
		__m256d Axi1 =
		    _mm256_load_pd(avx_data.data() + i + align_offset + avx_Nm * iAxi1);
		__m256d Axi2 =
		    _mm256_load_pd(avx_data.data() + i + align_offset + avx_Nm * iAxi2);

		//  - the three components of the vector k * kappa
		__m256d kkappa0 = _mm256_load_pd(avx_data.data() + i + align_offset +
		                                 avx_Nm * ikkappa0);
		__m256d kkappa1 = _mm256_load_pd(avx_data.data() + i + align_offset +
		                                 avx_Nm * ikkappa1);
		__m256d kkappa2 = _mm256_load_pd(avx_data.data() + i + align_offset +
		                                 avx_Nm * ikkappa2);

		//  - the phase beta.
		__m256d beta =
		    _mm256_load_pd(avx_data.data() + i + align_offset + avx_Nm * ibeta);

		// Then, do the computation.

		// This is the scalar product between k*kappa and pos:
		__m256d z = _mm256_add_pd(_mm256_mul_pd(pos0, kkappa0),
		                          _mm256_add_pd(_mm256_mul_pd(pos1, kkappa1),
		                                        _mm256_mul_pd(pos2, kkappa2)));

		// Here, the phase is added on. This is the argument of the cosine.
		__m256d cos_arg = _mm256_add_pd(z, beta);

		// ********
		// * Computing the cosine
		// * Part 1: Argument reduction
		//
		//  To understand the computation of the cosine, first note that the
		//  cosine is periodic and we thus only need to model its behavior
		//  between 0 and 2*pi to be able compute the function anywhere. In
		//  fact, by mirroring the function along the x and y axes, even the
		//  range between 0 and pi/2 is sufficient for this purpose. In this
		//  range, the cosine can be efficiently evaluated with high precision
		//  by using a polynomial approximation. Thus, to compute the cosine,
		//  the input value is first reduced so that it lies within this range.
		//  Then, the polynomial approximation is evaluated. Finally, if
		//  necessary, the sign of the result is flipped (mirroring the function
		//  along the x axis).
		//
		//  The actual computation is slightly more involved. First, argument
		//  reduction can be simplified drastically by computing cos(pi*x),
		//  such that the values are reduced to the range [0, 0.5) instead of
		//  [0, pi/2). Since the cosine is even (independent of the sign), we
		//  can first reduce values to [-0.5, 0.5) – that is, a simple rounding
		//  operation – and then neutralize the sign. In fact, precisely because
		//  the cosine is even, all terms of the polynomial are powers of x^2,
		//  so the value of x^2 (computed as x*x) forms the basis for the
		//  polynomial approximation. If I understand things correctly, then (in
		//  IEEE-754 floating point) x*x and (-x)*(-x) will always result in the
		//  exact same value, which means that any error bound over [0, 0.5)
		//  automatically applies to (-0.5, 0] as well.

		// First, compute round(x), and store it in q. If this value is odd,
		// we're looking at the negative half-wave of the cosine, and thus
		// will have to invert the sign of the result.
		__m256d q = _mm256_round_pd(
		    cos_arg, (_MM_FROUND_TO_NEAREST_INT | _MM_FROUND_NO_EXC));

		// Since we're computing cos(pi*x), round(x) always yields the center of
		// a half-wave (where cos(pi*x) achieves an extremum). This point
		// logically corresponds to x=0. Therefore, we subtract this center from
		// the actual input argument to find the corresponding point on the
		// half-wave that is centered around zero.
		__m256d s = _mm256_sub_pd(cos_arg, q);

		// We now want to check whether q (the index of our half-wave) is even
		// or odd, since all of the odd-numbered half-waves are negative, so
		// we'll have to flip the final result. On an int, this is as simple as
		// checking the 0th bit. Idea: manipulate the double in such a way that
		// we can do this. So, we add 2^52, such that the last digit of the
		// mantissa is actually in the ones' position. Since q may be negative,
		// we'll also add 2^51 to make sure it's positive. Note that 2^51 is
		// even and thus leaves evenness invariant, which is the only thing we
		// care about here.
		//
		// This is based on the int extraction process described here:
		// https://stackoverflow.com/questions/41144668/how-to-efficiently-perform-double-int64-conversions-with-sse-avx/41223013
		//
		// We assume -2^51 <= q < 2^51 for this, which is unproblematic, as
		// double precision has decayed far enough at that point that the
		// usefulness of the cosine becomes limited.
		//
		// Explanation: The mantissa of a double-precision float has 52 bits
		// (excluding the implicit first bit, which is always one). If |q| >
		// 2^51, this implicit first bit has a place value of at least 2^51,
		// while the first stored bit of the mantissa has a place value of at
		// least 2^50. This means that the LSB of the mantissa has a place value
		// of at least 2^(-1), or 0.5. For a cos(pi*x), this corresponds to a
		// quarter of a cycle (pi/2), so at this point the precision of the
		// input argument is so low that going from one representable number to
		// the next causes the result to jump by +/-1.

		q = _mm256_add_pd(q, _mm256_set1_pd(0x0018000000000000));

		// Unfortunately, integer comparisons were only introduced in AVX2, so
		// we'll have to make do with a floating point comparison to check
		// whether the last bit is set. However, masking out all but the last
		// bit will result in a denormal float, which may either result in
		// performance problems or just be rounded down to zero, neither of
		// which is what we want here. To fix this, we'll mask in not only bit
		// 0, but also the exponent (and sign, but that doesn't matter) of q.
		// Luckily, the exponent of q is guaranteed to have the fixed value of
		// 1075 (corresponding to 2^52) after our addition.

		__m256d invert = _mm256_and_pd(
		    q, _mm256_castsi256_pd(_mm256_set1_epi64x(0xfff0000000000001)));

		// If we did have a one in bit 0, our result will be equal to 2^52 + 1.
		invert = _mm256_cmp_pd(
		    invert, _mm256_castsi256_pd(_mm256_set1_epi64x(0x4330000000000001)),
		    _CMP_EQ_OQ);

		// Now we know whether to flip the sign of the result. However, remember
		// that we're working on multiple values at a time, so an if statement
		// won't be of much use here (plus it might perform badly). Instead,
		// we'll make use of the fact that the result of the comparison is all
		// ones if the comparison was true (i.e. q is odd and we need to flip
		// the result), and all zeroes otherwise. If we now mask out all bits
		// except the sign bit, we get something that, when xor'ed into our
		// final result, will flip the sign exactly when q is odd.
		invert = _mm256_and_pd(invert, _mm256_set1_pd(-0.0));
		// (Note that the binary representation of -0.0 is all 0 bits, except
		// for the sign bit, which is set to 1.)

		// TODO: clamp floats between 0 and 1? This would ensure that we never
		// see inf's, but maybe we want that, so that things dont just fail
		// silently...

		// * end of argument reduction
		// *******

		// ******
		// * Evaluate the cosine using a polynomial approximation for the zeroth
		// half-wave.
		// * The coefficients for this were generated using sleefs gencoef.c.
		// * These coefficients are probably far from optimal; however, they
		// should be sufficient for this case.
		s = _mm256_mul_pd(s, s);

		__m256d u = _mm256_set1_pd(+0.2211852080653743946e+0);

		u = _mm256_add_pd(_mm256_mul_pd(u, s),
		                  _mm256_set1_pd(-0.1332560668688523853e+1));
		u = _mm256_add_pd(_mm256_mul_pd(u, s),
		                  _mm256_set1_pd(+0.4058509506474178075e+1));
		u = _mm256_add_pd(_mm256_mul_pd(u, s),
		                  _mm256_set1_pd(-0.4934797516664651162e+1));
		u = _mm256_add_pd(_mm256_mul_pd(u, s), _mm256_set1_pd(1.));

		// Then, flip the sign of each double for which invert is not zero.
		// Since invert has only zero bits except for a possible one in bit 63,
		// we can xor it onto our result to selectively invert the 63rd (sign)
		// bit in each double where invert is set.
		u = _mm256_xor_pd(u, invert);

		// * end computation of cosine
		// **********

		// Finally, Ak*xi is multiplied on. Since this is a vector, the
		// multiplication needs to be done for each of the three
		// components, so it happens separately.
		acc0 = _mm256_add_pd(_mm256_mul_pd(u, Axi0), acc0);
		acc1 = _mm256_add_pd(_mm256_mul_pd(u, Axi1), acc1);
		acc2 = _mm256_add_pd(_mm256_mul_pd(u, Axi2), acc2);
	}

	return Vector3d(hsum_double_avx(acc0) * smoothing_factor * scaling_factor, 
                    hsum_double_avx(acc1) * smoothing_factor * scaling_factor,
                    hsum_double_avx(acc2) * smoothing_factor * scaling_factor);
#endif // ENABLE_FAST_WAVES
}

} // namespace crpropa
