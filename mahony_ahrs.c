//=====================================================================================================
// MahonyAHRS.c
//=====================================================================================================
//
// Madgwick's implementation of Mayhony's AHRS algorithm.
// See: http://www.x-io.co.uk/node/8#open_source_ahrs_and_imu_algorithms
//
// Date			Author			Notes
// 29/09/2011	SOH Madgwick    Initial release
// 02/10/2011	SOH Madgwick	Optimised for reduced CPU load
// 08/23/2020	Disi A	        Object-oriented fashion
// 08/31/2020	Disi A	        Variable sample frequence implementation
//
//=====================================================================================================

//---------------------------------------------------------------------------------------------------
// Header files

#include "mahony_ahrs.h"
#include <stdlib.h>
#include <math.h>

//---------------------------------------------------------------------------------------------------
// Variable definitions
MahonyAHRS* create_mahony_ahrs(MA_PRECISION sample_rate){
	if(sample_rate <= 0) return NULL;
	MahonyAHRS* workspace = (MahonyAHRS *) malloc(sizeof(MahonyAHRS));
	workspace->sample_rate = sample_rate;
	// quaternion of sensor frame relative to auxiliary frame
	workspace->q0 = 1.0f;
	workspace->q1 = 0.0f;
	workspace->q2 = 0.0f;
	workspace->q3 = 0.0f;
	// integral error terms scaled by Ki
	workspace->integralFBx = 0.0f;
	workspace->integralFBy = 0.0f;
	workspace->integralFBz = 0.0f;
	return workspace;
}

void mahony_ahrs_update_sample_rate(MahonyAHRS* workspace, MA_PRECISION sample_rate) {
	if(workspace == NULL) return;
	if(sample_rate <= 0) return;
	if(sample_rate != workspace -> sample_rate) {// Reset the parameters when sample rates are different;
		workspace->sample_rate = sample_rate;
		workspace->q0 = 1.0f;
		workspace->q1 = 0.0f;
		workspace->q2 = 0.0f;
		workspace->q3 = 0.0f;
		workspace->integralFBx = 0.0f;
		workspace->integralFBy = 0.0f;
		workspace->integralFBz = 0.0f;
	}
}

void free_mahony_ahrs(MahonyAHRS* workspace){
	free(workspace);
}

//====================================================================================================
// Functions

//---------------------------------------------------------------------------------------------------
// IMU algorithm update
void mahony_ahrs_update_imu(MahonyAHRS* workspace, MA_PRECISION gx, MA_PRECISION gy, MA_PRECISION gz, MA_PRECISION ax, MA_PRECISION ay, MA_PRECISION az)
{
	// temp vars
	MA_PRECISION recipNorm;
	MA_PRECISION halfvx, halfvy, halfvz;
	MA_PRECISION halfex, halfey, halfez;
	MA_PRECISION qa, qb, qc;

	// Compute feedback only if accelerometer measurement valid (avoids NaN in accelerometer normalisation)
	if (!((ax == 0.0f) && (ay == 0.0f) && (az == 0.0f)))
	{

		// Normalise accelerometer measurement
		recipNorm = inv_sqrt(ax * ax + ay * ay + az * az);
		ax *= recipNorm;
		ay *= recipNorm;
		az *= recipNorm;

		// Estimated direction of gravity and vector perpendicular to magnetic flux
		halfvx = workspace->q1 * workspace->q3 - workspace->q0 * workspace->q2;
		halfvy = workspace->q0 * workspace->q1 + workspace->q2 * workspace->q3;
		halfvz = workspace->q0 * workspace->q0 - 0.5f + workspace->q3 * workspace->q3;

		// Error is sum of cross product between estimated and measured direction of gravity
		halfex = (ay * halfvz - az * halfvy);
		halfey = (az * halfvx - ax * halfvz);
		halfez = (ax * halfvy - ay * halfvx);

		// Compute and apply integral feedback if enabled
		if (TWO_KI > 0.0f)
		{
			workspace->integralFBx += TWO_KI * halfex * (1.0f / workspace->sample_rate); // integral error scaled by Ki
			workspace->integralFBy += TWO_KI * halfey * (1.0f / workspace->sample_rate);
			workspace->integralFBz += TWO_KI * halfez * (1.0f / workspace->sample_rate);
			gx += workspace->integralFBx; // apply integral feedback
			gy += workspace->integralFBy;
			gz += workspace->integralFBz;
		} else {
			workspace->integralFBx = 0.0f; // prevent integral windup
			workspace->integralFBy = 0.0f;
			workspace->integralFBz = 0.0f;
		}

		// Apply proportional feedback
		gx += TWO_KP * halfex;
		gy += TWO_KP * halfey;
		gz += TWO_KP * halfez;
	}

	// Integrate rate of change of quaternion
	gx *= (0.5f / workspace->sample_rate); // pre-multiply common factors
	gy *= (0.5f / workspace->sample_rate);
	gz *= (0.5f / workspace->sample_rate);
	qa = workspace->q0;
	qb = workspace->q1;
	qc = workspace->q2;
	workspace->q0 += (-qb * gx - qc * gy - workspace->q3 * gz);
	workspace->q1 += (qa * gx + qc * gz - workspace->q3 * gy);
	workspace->q2 += (qa * gy - qb * gz + workspace->q3 * gx);
	workspace->q3 += (qa * gz + qb * gy - qc * gx);

	// Normalise quaternion
	recipNorm = inv_sqrt(workspace->q0 * workspace->q0 + workspace->q1 * workspace->q1 + workspace->q2 * workspace->q2 + workspace->q3 * workspace->q3);
	workspace->q0 *= recipNorm;
	workspace->q1 *= recipNorm;
	workspace->q2 *= recipNorm;
	workspace->q3 *= recipNorm;
	COMPUTE_EULER_ANGLE(workspace);
}

//---------------------------------------------------------------------------------------------------
// AHRS algorithm update
void mahony_ahrs_update(MahonyAHRS* workspace, MA_PRECISION gx, MA_PRECISION gy, MA_PRECISION gz, MA_PRECISION ax, MA_PRECISION ay, MA_PRECISION az, MA_PRECISION mx, MA_PRECISION my, MA_PRECISION mz)
{
	MA_PRECISION recipNorm;
	MA_PRECISION q0q0, q0q1, q0q2, q0q3, q1q1, q1q2, q1q3, q2q2, q2q3, q3q3;
	MA_PRECISION hx, hy, bx, bz;
	MA_PRECISION halfvx, halfvy, halfvz, halfwx, halfwy, halfwz;
	MA_PRECISION halfex, halfey, halfez;
	MA_PRECISION qa, qb, qc;

	// Use IMU algorithm if magnetometer measurement invalid (avoids NaN in magnetometer normalisation)
	if ((mx == 0.0f) && (my == 0.0f) && (mz == 0.0f))
	{
		mahony_ahrs_update_imu(workspace, gx, gy, gz, ax, ay, az);
		return;
	}

	// Compute feedback only if accelerometer measurement valid (avoids NaN in accelerometer normalisation)
	if (!((ax == 0.0f) && (ay == 0.0f) && (az == 0.0f)))
	{

		// Normalise accelerometer measurement
		recipNorm = inv_sqrt(ax * ax + ay * ay + az * az);
		ax *= recipNorm;
		ay *= recipNorm;
		az *= recipNorm;

		// Normalise magnetometer measurement
		recipNorm = inv_sqrt(mx * mx + my * my + mz * mz);
		mx *= recipNorm;
		my *= recipNorm;
		mz *= recipNorm;

		// Auxiliary variables to avoid repeated arithmetic
		q0q0 = workspace->q0 * workspace->q0;
		q0q1 = workspace->q0 * workspace->q1;
		q0q2 = workspace->q0 * workspace->q2;
		q0q3 = workspace->q0 * workspace->q3;
		q1q1 = workspace->q1 * workspace->q1;
		q1q2 = workspace->q1 * workspace->q2;
		q1q3 = workspace->q1 * workspace->q3;
		q2q2 = workspace->q2 * workspace->q2;
		q2q3 = workspace->q2 * workspace->q3;
		q3q3 = workspace->q3 * workspace->q3;

		// Reference direction of Earth's magnetic field
		hx = 2.0f * (mx * (0.5f - q2q2 - q3q3) + my * (q1q2 - q0q3) + mz * (q1q3 + q0q2));
		hy = 2.0f * (mx * (q1q2 + q0q3) + my * (0.5f - q1q1 - q3q3) + mz * (q2q3 - q0q1));
		bx = sqrt(hx * hx + hy * hy);
		bz = 2.0f * (mx * (q1q3 - q0q2) + my * (q2q3 + q0q1) + mz * (0.5f - q1q1 - q2q2));
		
		// Estimated direction of gravity and magnetic field
		halfvx = q1q3 - q0q2;
		halfvy = q0q1 + q2q3;
		halfvz = q0q0 - 0.5f + q3q3;
		halfwx = bx * (0.5f - q2q2 - q3q3) + bz * (q1q3 - q0q2);
		halfwy = bx * (q1q2 - q0q3) + bz * (q0q1 + q2q3);
		halfwz = bx * (q0q2 + q1q3) + bz * (0.5f - q1q1 - q2q2);

		// Error is sum of cross product between estimated direction and measured direction of field vectors
		halfex = (ay * halfvz - az * halfvy) + (my * halfwz - mz * halfwy);
		halfey = (az * halfvx - ax * halfvz) + (mz * halfwx - mx * halfwz);
		halfez = (ax * halfvy - ay * halfvx) + (mx * halfwy - my * halfwx);

		// Compute and apply integral feedback if enabled
		if (TWO_KI > 0.0f)
		{
			workspace->integralFBx += TWO_KI * halfex * (1.0f / workspace->sample_rate); // integral error scaled by Ki
			workspace->integralFBy += TWO_KI * halfey * (1.0f / workspace->sample_rate);
			workspace->integralFBz += TWO_KI * halfez * (1.0f / workspace->sample_rate);
			gx += workspace->integralFBx; // apply integral feedback
			gy += workspace->integralFBy;
			gz += workspace->integralFBz;
		} else {
			workspace->integralFBx = 0.0f; // prevent integral windup
			workspace->integralFBy = 0.0f;
			workspace->integralFBz = 0.0f;
		}

		// Apply proportional feedback
		gx += TWO_KP * halfex;
		gy += TWO_KP * halfey;
		gz += TWO_KP * halfez;
	}

	// Integrate rate of change of quaternion
	gx *= (0.5f / workspace->sample_rate); // pre-multiply common factors
	gy *= (0.5f / workspace->sample_rate);
	gz *= (0.5f / workspace->sample_rate);
	qa = workspace->q0;
	qb = workspace->q1;
	qc = workspace->q2;
	workspace->q0 += (-qb * gx - qc * gy - workspace->q3 * gz);
	workspace->q1 += (qa * gx + qc * gz - workspace->q3 * gy);
	workspace->q2 += (qa * gy - qb * gz + workspace->q3 * gx);
	workspace->q3 += (qa * gz + qb * gy - qc * gx);

	// Normalise quaternion
	recipNorm = inv_sqrt(workspace->q0 * workspace->q0 + workspace->q1 * workspace->q1 + workspace->q2 * workspace->q2 + workspace->q3 * workspace->q3);
	workspace->q0 *= recipNorm;
	workspace->q1 *= recipNorm;
	workspace->q2 *= recipNorm;
	workspace->q3 *= recipNorm;
	
	COMPUTE_EULER_ANGLE(workspace);
}

