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
//
//=====================================================================================================
#ifndef mahony_ahrs_h
#define mahony_ahrs_h

#define MA_DOUBLE_PRECISION 0

#if MA_DOUBLE_PRECISION
#define MA_PRECISION double
#else
#define MA_PRECISION float
#endif

typedef struct {
    MA_PRECISION q0;
    MA_PRECISION q1;
    MA_PRECISION q2;
    MA_PRECISION q3;
    // Intermediate variables
    MA_PRECISION integralFBx;
    MA_PRECISION integralFBy;
    MA_PRECISION integralFBz;
    MA_PRECISION yaw;
    MA_PRECISION pitch;
    MA_PRECISION roll;
} MahonyAHRS;
//----------------------------------------------------------------------------------------------------
// Variable declaration
#define MA_SAMPLE_RATE 512.0
#define TWO_KP (2.0f * 0.5f)   // 2 * proportional gain (Kp)
#define TWO_KI (2.0f * 0f)     // 2 * integral gain (Ki)

MahonyAHRS* create_mahony_ahrs();
MahonyAHRS* free_mahony_ahrs();

//---------------------------------------------------------------------------------------------------
// Function declarations

void mahony_ahrs_update_imu(MahonyAHRS* workspace, float gx, float gy, float gz, float ax, float ay, float az);

void mahony_ahrs_update(MahonyAHRS* workspace, float gx, float gy, float gz, float ax, float ay, float az, float mx, float my, float mz);

#endif /* mahony_ahrs_h */