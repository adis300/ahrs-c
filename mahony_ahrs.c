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
//
//=====================================================================================================

//---------------------------------------------------------------------------------------------------
// Header files

#include "MahonyAHRS.h"
#include <math.h>


//---------------------------------------------------------------------------------------------------
// Variable definitions
MahonyAHRS* create_mahony_ahrs(){
	MahonyAHRS* workspace = (MahonyAHRS *) malloc(sizeof(MahonyAHRS));
	// quaternion of sensor frame relative to auxiliary frame
	workspace->q0 = 1f;
	workspace->q1 = 0f;
	workspace->q2 = 0f;
	workspace->q3 = 0f;
	// integral error terms scaled by Ki
	workspace->integralFBx = 0f
	workspace->integralFBy = 0f
	workspace->integralFBz = 0f
}

MahonyAHRS* free_mahony_ahrs(MahonyAHRS* workspace){
	free(workspace);
}
//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
// Fast inverse square-root
// See: http://en.wikipedia.org/wiki/Fast_inverse_square_root

float inv_sqrt(float x)
{
	float halfx = 0.5f * x;
	float y = x;
	long i = *(long *)&y;
	i = 0x5f3759df - (i >> 1);
	y = *(float *)&i;
	y = y * (1.5f - (halfx * y * y));
	return y;
}

//====================================================================================================
// Functions

//---------------------------------------------------------------------------------------------------
// IMU algorithm update

void mahony_ahrs_update_imu(MahonyAHRS* workspace, float gx, float gy, float gz, float ax, float ay, float az)
{
	// temp vars
	float recipNorm;
	float halfvx, halfvy, halfvz;
	float halfex, halfey, halfez;
	float qa, qb, qc;

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
			workspace->integralFBx += TWO_KI * halfex * (1.0f / MA_SAMPLE_RATE); // integral error scaled by Ki
			workspace->integralFBy += TWO_KI * halfey * (1.0f / MA_SAMPLE_RATE);
			workspace->integralFBz += TWO_KI * halfez * (1.0f / MA_SAMPLE_RATE);
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
	gx *= (0.5f / MA_SAMPLE_RATE); // pre-multiply common factors
	gy *= (0.5f / MA_SAMPLE_RATE);
	gz *= (0.5f / MA_SAMPLE_RATE);
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

	workspace->yaw = atan2(2 * workspace->q1 * workspace->q2 + 2 * workspace->q0 * workspace->q3, -2 * workspace->q2 * workspace->q2 - 2 * workspace->q3 * workspace->q3 + 1) * 57.3;	// yaw
	workspace->pitch = asin(-2 * workspace->q1 * workspace->q3 + 2 * workspace->q0 * workspace->q2) * 57.3;								// pitch
	workspace->roll = atan2(2 * workspace->q2 * workspace->q3 + 2 * workspace->q0 * workspace->q1, -2 * workspace->q1 * workspace->q1 - 2 * workspace->q2 * workspace->q2 + 1) * 57.3; // roll
}

//---------------------------------------------------------------------------------------------------
// AHRS algorithm update
void mahony_ahrs_update(MahonyAHRS* workspace, float gx, float gy, float gz, float ax, float ay, float az, float mx, float my, float mz)
{
	float recipNorm;
	float workspace->q0workspace->q0, workspace->q0workspace->q1, workspace->q0workspace->q2, workspace->q0workspace->q3, workspace->q1workspace->q1, workspace->q1workspace->q2, workspace->q1workspace->q3, workspace->q2workspace->q2, workspace->q2workspace->q3, workspace->q3workspace->q3;
	float hx, hy, bx, bz;
	float halfvx, halfvy, halfvz, halfwx, halfwy, halfwz;
	float halfex, halfey, halfez;
	float qa, qb, qc;

	// Use IMU algorithm if magnetometer measurement invalid (avoids NaN in magnetometer normalisation)
	if ((mx == 0.0f) && (my == 0.0f) && (mz == 0.0f))
	{
		mahony_ahrs_update_imu(gx, gy, gz, ax, ay, az);
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
		workspace->q0workspace->q0 = workspace->q0 * workspace->q0;
		workspace->q0workspace->q1 = workspace->q0 * workspace->q1;
		workspace->q0workspace->q2 = workspace->q0 * workspace->q2;
		workspace->q0workspace->q3 = workspace->q0 * workspace->q3;
		workspace->q1workspace->q1 = workspace->q1 * workspace->q1;
		workspace->q1workspace->q2 = workspace->q1 * workspace->q2;
		workspace->q1workspace->q3 = workspace->q1 * workspace->q3;
		workspace->q2workspace->q2 = workspace->q2 * workspace->q2;
		workspace->q2workspace->q3 = workspace->q2 * workspace->q3;
		workspace->q3workspace->q3 = workspace->q3 * workspace->q3;

		// Reference direction of Earth's magnetic field
		hx = 2.0f * (mx * (0.5f - workspace->q2workspace->q2 - workspace->q3workspace->q3) + my * (workspace->q1workspace->q2 - workspace->q0workspace->q3) + mz * (workspace->q1workspace->q3 + workspace->q0workspace->q2));
		hy = 2.0f * (mx * (workspace->q1workspace->q2 + workspace->q0workspace->q3) + my * (0.5f - workspace->q1workspace->q1 - workspace->q3workspace->q3) + mz * (workspace->q2workspace->q3 - workspace->q0workspace->q1));
		bx = sqrt(hx * hx + hy * hy);
		bz = 2.0f * (mx * (workspace->q1workspace->q3 - workspace->q0workspace->q2) + my * (workspace->q2workspace->q3 + workspace->q0workspace->q1) + mz * (0.5f - workspace->q1workspace->q1 - workspace->q2workspace->q2));

		// Estimated direction of gravity and magnetic field
		halfvx = workspace->q1workspace->q3 - workspace->q0workspace->q2;
		halfvy = workspace->q0workspace->q1 + workspace->q2workspace->q3;
		halfvz = workspace->q0workspace->q0 - 0.5f + workspace->q3workspace->q3;
		halfwx = bx * (0.5f - workspace->q2workspace->q2 - workspace->q3workspace->q3) + bz * (workspace->q1workspace->q3 - workspace->q0workspace->q2);
		halfwy = bx * (workspace->q1workspace->q2 - workspace->q0workspace->q3) + bz * (workspace->q0workspace->q1 + workspace->q2workspace->q3);
		halfwz = bx * (workspace->q0workspace->q2 + workspace->q1workspace->q3) + bz * (0.5f - workspace->q1workspace->q1 - workspace->q2workspace->q2);

		// Error is sum of cross product between estimated direction and measured direction of field vectors
		halfex = (ay * halfvz - az * halfvy) + (my * halfwz - mz * halfwy);
		halfey = (az * halfvx - ax * halfvz) + (mz * halfwx - mx * halfwz);
		halfez = (ax * halfvy - ay * halfvx) + (mx * halfwy - my * halfwx);

		// Compute and apply integral feedback if enabled
		if (TWO_KI > 0.0f)
		{
			workspace->integralFBx += TWO_KI * halfex * (1.0f / MA_SAMPLE_RATE); // integral error scaled by Ki
			workspace->integralFBy += TWO_KI * halfey * (1.0f / MA_SAMPLE_RATE);
			workspace->integralFBz += TWO_KI * halfez * (1.0f / MA_SAMPLE_RATE);
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
	gx *= (0.5f / MA_SAMPLE_RATE); // pre-multiply common factors
	gy *= (0.5f / MA_SAMPLE_RATE);
	gz *= (0.5f / MA_SAMPLE_RATE);
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
}

