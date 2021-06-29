//=====================================================================================================
// madgwick_ahrs.h
//=====================================================================================================
//
// Madgwick's implementation of Madgwick's AHRS algorithm.
// See: http://www.x-io.co.uk/node/8#open_source_ahrs_and_imu_algorithms
//
// Date			Author			Notes
// 29/09/2011	SOH Madgwick    Initial release
// 02/10/2011	SOH Madgwick	Optimised for reduced CPU load
// 19/02/2012	SOH Madgwick	Magnetometer measurement is normalised
// 08/31/2020	Disi A	        Object-oriented fashion
//
//=====================================================================================================// Header files

#include "madgwick_ahrs.h"
#include <stdlib.h>
#include <math.h>


#if MA_DOUBLE_PRECISION
#define ATAN2 atan2
#define ASIN asin
#define SQRT sqrt
#define MAGIC_R 0x5fe6eb50c7b537a9
#else
#define ATAN2 atan2f
#define ASIN asinf
#define SQRT sqrtf
#define MAGIC_R 0x5f3759df
#endif

//---------------------------------------------------------------------------------------------------
// Variable definitions
MadgwickAHRS* create_madgwick_ahrs(MA_PRECISION sample_rate){
	if(sample_rate <= 0) return NULL;
	MadgwickAHRS* workspace = (MadgwickAHRS *) malloc(sizeof(MadgwickAHRS));
	workspace->sample_rate = sample_rate;
	// quaternion of sensor frame relative to auxiliary frame
	workspace->q0 = 1.0f;
	workspace->q1 = 0.0f;
	workspace->q2 = 0.0f;
	workspace->q3 = 0.0f;
	return workspace;
}

void madgwick_ahrs_update_sample_rate(MadgwickAHRS* workspace, MA_PRECISION sample_rate) {
	if(workspace == NULL) return;
	if(sample_rate <= 0) return;
	if(sample_rate != workspace -> sample_rate) {// Reset the parameters when sample rates are different;
		workspace->sample_rate = sample_rate;
		workspace->q0 = 1.0f;
		workspace->q1 = 0.0f;
		workspace->q2 = 0.0f;
		workspace->q3 = 0.0f;
	}
}

void free_madgwick_ahrs(MadgwickAHRS* workspace){
	free(workspace);
}
//---------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------------------------------
// Fast inverse square-root
// See: http://en.wikipedia.org/wiki/Fast_inverse_square_root

MA_PRECISION inv_sqrt(MA_PRECISION x)
{
	MA_PRECISION halfx = 0.5f * x;
	MA_PRECISION y = x;
	long i = *(long *)&y;
	i = MAGIC_R - (i >> 1);
	y = *(MA_PRECISION *)&i;
	y = y * (1.5f - (halfx * y * y));
	return y;
}

void compute_euler_angle(MadgwickAHRS* workspace){
	//*yaw = atan2f((2*q.q2*q.q3 - 2*q.q1*q.q4), (2*q.q1*q.q1 + 2*q.q2*q.q2 -1));  // equation (7)
    //*pitch = -asinf(2*q.q2*q.q4 + 2*q.q1*q.q3);                                  // equatino (8)
    //*roll  = atan2f((2*q.q3*q.q4 - 2*q.q1*q.q2), (2*q.q1*q.q1 + 2*q.q4*q.q4 -1));
	workspace->yaw = atan2(2 * workspace->q1 * workspace->q2 + 2 * workspace->q0 * workspace->q3, -2 * workspace->q2 * workspace->q2 - 2 * workspace->q3 * workspace->q3 + 1) * 57.3;	// yaw
	workspace->pitch = asin(-2 * workspace->q1 * workspace->q3 + 2 * workspace->q0 * workspace->q2) * 57.3;								// pitch
	workspace->roll = atan2(2 * workspace->q2 * workspace->q3 + 2 * workspace->q0 * workspace->q1, -2 * workspace->q1 * workspace->q1 - 2 * workspace->q2 * workspace->q2 + 1) * 57.3; // roll
}

//====================================================================================================
// Functions

//---------------------------------------------------------------------------------------------------
// IMU algorithm update

void madgwick_ahrs_update_imu(MadgwickAHRS* workspace, MA_PRECISION gx, MA_PRECISION gy, MA_PRECISION gz, MA_PRECISION ax, MA_PRECISION ay, MA_PRECISION az)
{
	// temp vars
	MA_PRECISION recipNorm;
	MA_PRECISION s0, s1, s2, s3;
	MA_PRECISION qDot1, qDot2, qDot3, qDot4;
	MA_PRECISION _2q0, _2q1, _2q2, _2q3, _4q0, _4q1, _4q2 ,_8q1, _8q2, q0q0, q1q1, q2q2, q3q3;

	// Rate of change of quaternion from gyroscope
	qDot1 = 0.5f * (-workspace->q1 * gx - workspace->q2 * gy - workspace->q3 * gz);
	qDot2 = 0.5f * (workspace->q0 * gx + workspace->q2 * gz - workspace->q3 * gy);
	qDot3 = 0.5f * (workspace->q0 * gy - workspace->q1 * gz + workspace->q3 * gx);
	qDot4 = 0.5f * (workspace->q0 * gz + workspace->q1 * gy - workspace->q2 * gx);

	// Compute feedback only if accelerometer measurement valid (avoids NaN in accelerometer normalisation)
	if(!((ax == 0.0f) && (ay == 0.0f) && (az == 0.0f))) {

		// Normalise accelerometer measurement
		recipNorm = inv_sqrt(ax * ax + ay * ay + az * az);
		ax *= recipNorm;
		ay *= recipNorm;
		az *= recipNorm;   

		// Auxiliary variables to avoid repeated arithmetic
		_2q0 = 2.0f * workspace->q0;
		_2q1 = 2.0f * workspace->q1;
		_2q2 = 2.0f * workspace->q2;
		_2q3 = 2.0f * workspace->q3;
		_4q0 = 4.0f * workspace->q0;
		_4q1 = 4.0f * workspace->q1;
		_4q2 = 4.0f * workspace->q2;
		_8q1 = 8.0f * workspace->q1;
		_8q2 = 8.0f * workspace->q2;
		q0q0 = workspace->q0 * workspace->q0;
		q1q1 = workspace->q1 * workspace->q1;
		q2q2 = workspace->q2 * workspace->q2;
		q3q3 = workspace->q3 * workspace->q3;

		// Gradient decent algorithm corrective step
		s0 = _4q0 * q2q2 + _2q2 * ax + _4q0 * q1q1 - _2q1 * ay;
		s1 = _4q1 * q3q3 - _2q3 * ax + 4.0f * q0q0 * workspace->q1 - _2q0 * ay - _4q1 + _8q1 * q1q1 + _8q1 * q2q2 + _4q1 * az;
		s2 = 4.0f * q0q0 * workspace->q2 + _2q0 * ax + _4q2 * q3q3 - _2q3 * ay - _4q2 + _8q2 * q1q1 + _8q2 * q2q2 + _4q2 * az;
		s3 = 4.0f * q1q1 * workspace->q3 - _2q1 * ax + 4.0f * q2q2 * workspace->q3 - _2q2 * ay;
		recipNorm = inv_sqrt(s0 * s0 + s1 * s1 + s2 * s2 + s3 * s3); // normalise step magnitude
		s0 *= recipNorm;
		s1 *= recipNorm;
		s2 *= recipNorm;
		s3 *= recipNorm;

		// Apply feedback step
		qDot1 -= BETA * s0;
		qDot2 -= BETA * s1;
		qDot3 -= BETA * s2;
		qDot4 -= BETA * s3;
	}

	// Integrate rate of change of quaternion to yield quaternion
	workspace->q0 += qDot1 * (1.0f / workspace->sample_rate);
	workspace->q1 += qDot2 * (1.0f / workspace->sample_rate);
	workspace->q2 += qDot3 * (1.0f / workspace->sample_rate);
	workspace->q3 += qDot4 * (1.0f / workspace->sample_rate);

	// Normalise quaternion
	recipNorm = inv_sqrt(workspace->q0 * workspace->q0 + workspace->q1 * workspace->q1 + workspace->q2 * workspace->q2 + workspace->q3 * workspace->q3);
	workspace->q0 *= recipNorm;
	workspace->q1 *= recipNorm;
	workspace->q2 *= recipNorm;
	workspace->q3 *= recipNorm;
	compute_euler_angle(workspace);
}

//---------------------------------------------------------------------------------------------------
// AHRS algorithm update
void madgwick_ahrs_update(MadgwickAHRS* workspace, MA_PRECISION gx, MA_PRECISION gy, MA_PRECISION gz, MA_PRECISION ax, MA_PRECISION ay, MA_PRECISION az, MA_PRECISION mx, MA_PRECISION my, MA_PRECISION mz)
{
	MA_PRECISION recipNorm;
	MA_PRECISION s0, s1, s2, s3;
	MA_PRECISION qDot1, qDot2, qDot3, qDot4;
	MA_PRECISION hx, hy;
	MA_PRECISION _2q0mx, _2q0my, _2q0mz, _2q1mx, _2bx, _2bz, _4bx, _4bz, _2q0, _2q1, _2q2, _2q3, _2q0q2, _2q2q3, q0q0, q0q1, q0q2, q0q3, q1q1, q1q2, q1q3, q2q2, q2q3, q3q3;

	// Use IMU algorithm if magnetometer measurement invalid (avoids NaN in magnetometer normalisation)
	if ((mx == 0.0f) && (my == 0.0f) && (mz == 0.0f))
	{
		madgwick_ahrs_update_imu(workspace, gx, gy, gz, ax, ay, az);
		return;
	}

    // Rate of change of quaternion from gyroscope
	qDot1 = 0.5f * (-workspace->q1 * gx - workspace->q2 * gy - workspace->q3 * gz);
	qDot2 = 0.5f * (workspace->q0 * gx + workspace->q2 * gz - workspace->q3 * gy);
	qDot3 = 0.5f * (workspace->q0 * gy - workspace->q1 * gz + workspace->q3 * gx);
	qDot4 = 0.5f * (workspace->q0 * gz + workspace->q1 * gy - workspace->q2 * gx);

	// Compute feedback only if accelerometer measurement valid (avoids NaN in accelerometer normalisation)
	if(!((ax == 0.0f) && (ay == 0.0f) && (az == 0.0f))) {
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
		_2q0mx = 2.0f * workspace->q0 * mx;
		_2q0my = 2.0f * workspace->q0 * my;
		_2q0mz = 2.0f * workspace->q0 * mz;
		_2q1mx = 2.0f * workspace->q1 * mx;
		_2q0 = 2.0f * workspace->q0;
		_2q1 = 2.0f * workspace->q1;
		_2q2 = 2.0f * workspace->q2;
		_2q3 = 2.0f * workspace->q3;
		_2q0q2 = 2.0f * workspace->q0 * workspace->q2;
		_2q2q3 = 2.0f * workspace->q2 * workspace->q3;
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
		hx = mx * q0q0 - _2q0my * workspace->q3 + _2q0mz * workspace->q2 + mx * q1q1 + _2q1 * my * workspace->q2 + _2q1 * mz * workspace->q3 - mx * q2q2 - mx * q3q3;
		hy = _2q0mx * workspace->q3 + my * q0q0 - _2q0mz * workspace->q1 + _2q1mx * workspace->q2 - my * q1q1 + my * q2q2 + _2q2 * mz * workspace->q3 - my * q3q3;
		_2bx = SQRT(hx * hx + hy * hy);
		_2bz = -_2q0mx * workspace->q2 + _2q0my * workspace->q1 + mz * q0q0 + _2q1mx * workspace->q3 - mz * q1q1 + _2q2 * my * workspace->q3 - mz * q2q2 + mz * q3q3;
		_4bx = 2.0f * _2bx;
		_4bz = 2.0f * _2bz;

        // Gradient decent algorithm corrective step
		s0 = -_2q2 * (2.0f * q1q3 - _2q0q2 - ax) + _2q1 * (2.0f * q0q1 + _2q2q3 - ay) - _2bz * workspace->q2 * (_2bx * (0.5f - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (-_2bx * workspace->q3 + _2bz * workspace->q1) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + _2bx * workspace->q2 * (_2bx * (q0q2 + q1q3) + _2bz * (0.5f - q1q1 - q2q2) - mz);
		s1 = _2q3 * (2.0f * q1q3 - _2q0q2 - ax) + _2q0 * (2.0f * q0q1 + _2q2q3 - ay) - 4.0f * workspace->q1 * (1 - 2.0f * q1q1 - 2.0f * q2q2 - az) + _2bz * workspace->q3 * (_2bx * (0.5f - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (_2bx * workspace->q2 + _2bz * workspace->q0) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + (_2bx * workspace->q3 - _4bz * workspace->q1) * (_2bx * (q0q2 + q1q3) + _2bz * (0.5f - q1q1 - q2q2) - mz);
		s2 = -_2q0 * (2.0f * q1q3 - _2q0q2 - ax) + _2q3 * (2.0f * q0q1 + _2q2q3 - ay) - 4.0f * workspace->q2 * (1 - 2.0f * q1q1 - 2.0f * q2q2 - az) + (-_4bx * workspace->q2 - _2bz * workspace->q0) * (_2bx * (0.5f - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (_2bx * workspace->q1 + _2bz * workspace->q3) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + (_2bx * workspace->q0 - _4bz * workspace->q2) * (_2bx * (q0q2 + q1q3) + _2bz * (0.5f - q1q1 - q2q2) - mz);
		s3 = _2q1 * (2.0f * q1q3 - _2q0q2 - ax) + _2q2 * (2.0f * q0q1 + _2q2q3 - ay) + (-_4bx * workspace->q3 + _2bz * workspace->q1) * (_2bx * (0.5f - q2q2 - q3q3) + _2bz * (q1q3 - q0q2) - mx) + (-_2bx * workspace->q0 + _2bz * workspace->q2) * (_2bx * (q1q2 - q0q3) + _2bz * (q0q1 + q2q3) - my) + _2bx * workspace->q1 * (_2bx * (q0q2 + q1q3) + _2bz * (0.5f - q1q1 - q2q2) - mz);
		recipNorm = inv_sqrt(s0 * s0 + s1 * s1 + s2 * s2 + s3 * s3); // normalise step magnitude
		s0 *= recipNorm;
		s1 *= recipNorm;
		s2 *= recipNorm;
		s3 *= recipNorm;

		// Apply feedback step
		qDot1 -= BETA * s0;
		qDot2 -= BETA * s1;
		qDot3 -= BETA * s2;
		qDot4 -= BETA * s3;
	}

	// Integrate rate of change of quaternion to yield quaternion
	workspace->q0 += qDot1 * (1.0f / workspace->sample_rate);
	workspace->q1 += qDot2 * (1.0f / workspace->sample_rate);
	workspace->q2 += qDot3 * (1.0f / workspace->sample_rate);
	workspace->q3 += qDot4 * (1.0f / workspace->sample_rate);

	// Normalise quaternion
	recipNorm = inv_sqrt(workspace->q0 * workspace->q0 + workspace->q1 * workspace->q1 + workspace->q2 * workspace->q2 + workspace->q3 * workspace->q3);
	workspace->q0 *= recipNorm;
	workspace->q1 *= recipNorm;
	workspace->q2 *= recipNorm;
	workspace->q3 *= recipNorm;

	compute_euler_angle(workspace);
}

