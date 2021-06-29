//=====================================================================================================
// mahony_ahrs.h
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
#ifndef mahony_ahrs_h
#define mahony_ahrs_h

#include "arhs.h"

typedef struct {
    MA_PRECISION sample_rate;
    MA_PRECISION q0;
    MA_PRECISION q1;
    MA_PRECISION q2;
    MA_PRECISION q3;
    // Intermediate variables
    MA_PRECISION integralFBx;
    MA_PRECISION integralFBy;
    MA_PRECISION integralFBz;
    
    // Result variables
    MA_PRECISION yaw;
    MA_PRECISION pitch;
    MA_PRECISION roll;
} MahonyAHRS;
//----------------------------------------------------------------------------------------------------
// Variable declaration
#define TWO_KP (2.0f * 0.5f)   // 2 * proportional gain (Kp)
#define TWO_KI (2.0f * 0.0f)     // 2 * integral gain (Ki)

MahonyAHRS* create_mahony_ahrs(MA_PRECISION sample_rate);
void free_mahony_ahrs(MahonyAHRS* workspace);

//---------------------------------------------------------------------------------------------------
// Function declarations
void mahony_ahrs_update_sample_rate(MahonyAHRS* workspace, MA_PRECISION sample_rate);

void mahony_ahrs_update_imu(MahonyAHRS* workspace, MA_PRECISION gx, MA_PRECISION gy, MA_PRECISION gz, MA_PRECISION ax, MA_PRECISION ay, MA_PRECISION az);

void mahony_ahrs_update(MahonyAHRS* workspace, MA_PRECISION gx, MA_PRECISION gy, MA_PRECISION gz, MA_PRECISION ax, MA_PRECISION ay, MA_PRECISION az, MA_PRECISION mx, MA_PRECISION my, MA_PRECISION mz);

#endif /* mahony_ahrs_h */