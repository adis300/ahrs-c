//=====================================================================================================
// madgwick.h
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
//=====================================================================================================
#ifndef __MADGWICK_AHRS_H__
#define __MADGWICK_AHRS_H__

#include "arhs.h"

typedef struct {
    MA_PRECISION sample_rate;
    MA_PRECISION q0;
    MA_PRECISION q1;
    MA_PRECISION q2;
    MA_PRECISION q3;
    
    // Result variables
    MA_PRECISION yaw;
    MA_PRECISION pitch;
    MA_PRECISION roll;
} MadgwickAHRS;
//----------------------------------------------------------------------------------------------------
// Variable declaration
#define BETA 0.033f   // 2 * proportional gain

MadgwickAHRS* create_madgwick_ahrs(MA_PRECISION sample_rate);
void free_madgwick_ahrs(MadgwickAHRS* workspace);

//---------------------------------------------------------------------------------------------------
// Function declarations
void madgwick_ahrs_update_sample_rate(MadgwickAHRS* workspace, MA_PRECISION sample_rate);

void madgwick_ahrs_update_imu(MadgwickAHRS* workspace, MA_PRECISION gx, MA_PRECISION gy, MA_PRECISION gz, MA_PRECISION ax, MA_PRECISION ay, MA_PRECISION az);

void madgwick_ahrs_update(MadgwickAHRS* workspace, MA_PRECISION gx, MA_PRECISION gy, MA_PRECISION gz, MA_PRECISION ax, MA_PRECISION ay, MA_PRECISION az, MA_PRECISION mx, MA_PRECISION my, MA_PRECISION mz);

#endif /* __MADGWICK_AHRS_H__ */