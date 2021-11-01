#include "Common.h"
#include "QuadControl.h"

#include "Utility/SimpleConfig.h"

#include "Utility/StringUtils.h"
#include "Trajectory.h"
#include "BaseController.h"
#include "Math/Mat3x3F.h"

#ifdef __PX4_NUTTX
#include <systemlib/param/param.h>
#endif

void QuadControl::Init()
{
  BaseController::Init();

  // variables needed for integral control
  integratedAltitudeError = 0;
    
#ifndef __PX4_NUTTX
  // Load params from simulator parameter system
  ParamsHandle config = SimpleConfig::GetInstance();
   
  // Load parameters (default to 0)
  kpPosXY = config->Get(_config+".kpPosXY", 0);
  kpPosZ = config->Get(_config + ".kpPosZ", 0);
  KiPosZ = config->Get(_config + ".KiPosZ", 0);
     
  kpVelXY = config->Get(_config + ".kpVelXY", 0);
  kpVelZ = config->Get(_config + ".kpVelZ", 0);

  kpBank = config->Get(_config + ".kpBank", 0);
  kpYaw = config->Get(_config + ".kpYaw", 0);

  kpPQR = config->Get(_config + ".kpPQR", V3F());

  maxDescentRate = config->Get(_config + ".maxDescentRate", 100);
  maxAscentRate = config->Get(_config + ".maxAscentRate", 100);
  maxSpeedXY = config->Get(_config + ".maxSpeedXY", 100);
  maxAccelXY = config->Get(_config + ".maxHorizAccel", 100);

  maxTiltAngle = config->Get(_config + ".maxTiltAngle", 100);

  minMotorThrust = config->Get(_config + ".minMotorThrust", 0);
  maxMotorThrust = config->Get(_config + ".maxMotorThrust", 100);
#else
  // load params from PX4 parameter system
  //TODO
  param_get(param_find("MC_PITCH_P"), &Kp_bank);
  param_get(param_find("MC_YAW_P"), &Kp_yaw);
#endif
}

VehicleCommand QuadControl::GenerateMotorCommands(float collThrustCmd, V3F momentCmd)
{
  // Convert a desired 3-axis moment and collective thrust command to 
  //   individual motor thrust commands
  // INPUTS: 
  //   collThrustCmd: desired collective thrust [N]
  //   momentCmd: desired rotation moment about each axis [N m]
  // OUTPUT:
  //   set class member variable cmd (class variable for graphing) where
  //   cmd.desiredThrustsN[0..3]: motor commands, in [N]

  // HINTS: 
  // - you can access parts of momentCmd via e.g. momentCmd.x
  // You'll need the arm length parameter L, and the drag/thrust ratio kappa

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

  //  cmd.desiredThrustsN[0] = mass * 9.81f / 4.f; // front left
  //  cmd.desiredThrustsN[1] = mass * 9.81f / 4.f; // front right
  //  cmd.desiredThrustsN[2] = mass * 9.81f / 4.f; // rear left
  //  cmd.desiredThrustsN[3] = mass * 9.81f / 4.f; // rear right

    float l = L / sqrt(2);

    float fsum_tau_x = momentCmd.x / l ;      // fsum_tau_x =  f1 - f2 -f4 + f3
    float fsum_tau_y = momentCmd.y / l ;      // fsum_tau_y =  f1 + f2 -f4 - f3
    float fsum_tau_z = momentCmd.z / kappa ;  // fsum_tau_z = -f1 + f2 -f4 + f3
                                                // collThrustCmd = f1 + f2 +f4 + f3

    cmd.desiredThrustsN[0] = (collThrustCmd + fsum_tau_x + fsum_tau_y - fsum_tau_z) / 4; // f1: front left
    cmd.desiredThrustsN[1] = (collThrustCmd - fsum_tau_x + fsum_tau_y + fsum_tau_z) / 4; // f2: front right
    cmd.desiredThrustsN[2] = (collThrustCmd + fsum_tau_x - fsum_tau_y + fsum_tau_z) / 4; // f3: rear left
    cmd.desiredThrustsN[3] = (collThrustCmd - fsum_tau_x - fsum_tau_y - fsum_tau_z) / 4; // f4: rear right

  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return cmd;
}

V3F QuadControl::BodyRateControl(V3F pqrCmd, V3F pqr)
{
  // Calculate a desired 3-axis moment given a desired and current body rate
  // INPUTS: 
  //   pqrCmd: desired body rates [rad/s]
  //   pqr: current or estimated body rates [rad/s]
  // OUTPUT:
  //   return a V3F containing the desired moments for each of the 3 axes

  // HINTS: 
  //  - you can use V3Fs just like scalars: V3F a(1,1,1), b(2,3,4), c; c=a-b;
  //  - you'll need parameters for moments of inertia Ixx, Iyy, Izz
  //  - you'll also need the gain parameter kpPQR (it's a V3F)


  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////
    V3F rate_error = pqrCmd - pqr;
    V3F momentCmd = V3F(Ixx, Iyy, Izz) * kpPQR * rate_error;
  
  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return momentCmd;
}

// returns a desired roll and pitch rate 
V3F QuadControl::RollPitchControl(V3F accelCmd, Quaternion<float> attitude, float collThrustCmd)
{
  // Calculate a desired pitch and roll angle rates based on a desired global
  //   lateral acceleration, the current attitude of the quad, and desired
  //   collective thrust command
  // INPUTS: 
  //   accelCmd: desired acceleration in global XY coordinates [m/s2]
  //   attitude: current or estimated attitude of the vehicle
  //   collThrustCmd: desired collective thrust of the quad [N]
  // OUTPUT:
  //   return a V3F containing the desired pitch and roll rates. The Z
  //     element of the V3F should be left at its default value (0)

  // HINTS: 
  //  - we already provide rotation matrix R: to get element R[1,2] (python) use R(1,2) (C++)
  //  - you'll need the roll/pitch gain kpBank
  //  - collThrustCmd is a force in Newtons! You'll likely want to convert it to acceleration first

  V3F pqrCmd;
  Mat3x3F R = attitude.RotationMatrix_IwrtB();

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////
  
  if (collThrustCmd > 0)
  {
      float b_x = R(0, 2);
      float c = -collThrustCmd / mass;
      float b_x_com = accelCmd.x / c;
      float b_x_com_constrained = CONSTRAIN(b_x_com, -maxTiltAngle, maxTiltAngle);
      float b_x_error = b_x_com_constrained - b_x;
      float b_x_p = kpBank * b_x_error;

      float b_y = R(1, 2);
      float b_y_com = accelCmd.y / c;
      float b_y_com_constrained = CONSTRAIN(b_y_com, -maxTiltAngle, maxTiltAngle);
      float b_y_error = b_y_com_constrained - b_y;
      float b_y_p = kpBank * b_y_error;


      pqrCmd.x = (R(1, 0) * b_x_p - R(0, 0) * b_y_p) / R(2, 2); //p_commanded
      pqrCmd.y = (R(1, 1) * b_x_p - R(0, 1) * b_y_p) / R(2, 2); //q_commanded

  }
  else {
      pqrCmd.x = 0;
      pqrCmd.y = 0;
  }
    /////////////////////////////// END STUDENT CODE ////////////////////////////

  return pqrCmd;
}

float QuadControl::AltitudeControl(float posZCmd, float velZCmd, float posZ, float velZ, Quaternion<float> attitude, float accelZCmd, float dt)
{
  // Calculate desired quad thrust based on altitude setpoint, actual altitude,
  //   vertical velocity setpoint, actual vertical velocity, and a vertical 
  //   acceleration feed-forward command
  // INPUTS: 
  //   posZCmd, velZCmd: desired vertical position and velocity in NED [m]
  //   posZ, velZ: current vertical position and velocity in NED [m]
  //   accelZCmd: feed-forward vertical acceleration in NED [m/s2]
  //   dt: the time step of the measurements [seconds]
  // OUTPUT:
  //   return a collective thrust command in [N]

  // HINTS: 
  //  - we already provide rotation matrix R: to get element R[1,2] (python) use R(1,2) (C++)
  //  - you'll need the gain parameters kpPosZ and kpVelZ
  //  - maxAscentRate and maxDescentRate are maximum vertical speeds. Note they're both >=0!
  //  - make sure to return a force, not an acceleration
  //  - remember that for an upright quad in NED, thrust should be HIGHER if the desired Z acceleration is LOWER

  Mat3x3F R = attitude.RotationMatrix_IwrtB();
  float thrust = 0;

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////
  ////float velZ_com = kpPosZ * (posZCmd - posZ) + velZCmd;  sil     !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! ! ! 1  
  ////////velZ_com = CONSTRAIN(velZ_com, -maxDescentRate, maxAscentRate);
  velZCmd = CONSTRAIN(velZCmd, -maxDescentRate, maxAscentRate);
  integratedAltitudeError += (posZCmd - posZ)*dt;
  float u_1_bar = accelZCmd + kpPosZ * (posZCmd - posZ) + kpVelZ * (velZCmd - velZ) + KiPosZ * integratedAltitudeError;
  float acceleration_cmd = (u_1_bar - 9.81) / R(2, 2);
  acceleration_cmd = CONSTRAIN(acceleration_cmd, -maxAscentRate / dt, maxAscentRate / dt);
  thrust = -mass * acceleration_cmd  ;
  /////////////////////////////// END STUDENT CODE ////////////////////////////
  
  return thrust;
}

// returns a desired acceleration in global frame
V3F QuadControl::LateralPositionControl(V3F posCmd, V3F velCmd, V3F pos, V3F vel, V3F accelCmdFF)
{
  // Calculate a desired horizontal acceleration based on 
  //  desired lateral position/velocity/acceleration and current pose
  // INPUTS: 
  //   posCmd: desired position, in NED [m]
  //   velCmd: desired velocity, in NED [m/s]
  //   pos: current position, NED [m]
  //   vel: current velocity, NED [m/s]
  //   accelCmdFF: feed-forward acceleration, NED [m/s2]
  // OUTPUT:
  //   return a V3F with desired horizontal accelerations. 
  //     the Z component should be 0
  // HINTS: 
  //  - use the gain parameters kpPosXY and kpVelXY
  //  - make sure you limit the maximum horizontal velocity and acceleration
  //    to maxSpeedXY and maxAccelXY

  // make sure we don't have any incoming z-component
  accelCmdFF.z = 0;
  velCmd.z = 0;
  posCmd.z = pos.z;

  // we initialize the returned desired acceleration to the feed-forward value.
  // Make sure to _add_, not simply replace, the result of your controller
  // to this variable
  V3F accelCmd = accelCmdFF;

  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////

  if (velCmd.mag() > maxSpeedXY) {
      velCmd = velCmd.norm() * maxSpeedXY;
  }

  V3F kpPosXYZ = V3F(kpPosXY, kpPosXY, 0);
  V3F kpVelXYZ = V3F(kpVelXY, kpVelXY, 0);
  V3F pos_error = posCmd - pos;
  V3F vel_error = velCmd - vel;
  
  accelCmd = kpPosXYZ * pos_error + kpVelXYZ * vel_error + accelCmdFF ;

  if (accelCmd.mag() > maxAccelXY) {
      accelCmd = accelCmd.norm() * maxAccelXY;
  }
  accelCmd.z = 0;


  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return accelCmd;
}

// returns desired yaw rate
float QuadControl::YawControl(float yawCmd, float yaw)
{
  // Calculate a desired yaw rate to control yaw to yawCmd
  // INPUTS: 
  //   yawCmd: commanded yaw [rad]
  //   yaw: current yaw [rad]
  // OUTPUT:
  //   return a desired yaw rate [rad/s]
  // HINTS: 
  //  - use fmodf(foo,b) to unwrap a radian angle measure float foo to range [0,b]. 
  //  - use the yaw control gain parameter kpYaw

  float yawRateCmd=0;
  ////////////////////////////// BEGIN STUDENT CODE ///////////////////////////
  yawCmd = fmodf(yawCmd, 2.0 * M_PI);  // Ensure the target is within range of 0 to 2*pi

  float psi_error = yawCmd - yaw;
  psi_error = fmodf(psi_error, 2.0 * M_PI);
  yawRateCmd = kpYaw * psi_error;

  /////////////////////////////// END STUDENT CODE ////////////////////////////

  return yawRateCmd;

}

VehicleCommand QuadControl::RunControl(float dt, float simTime)
{
  curTrajPoint = GetNextTrajectoryPoint(simTime);

  float collThrustCmd = AltitudeControl(curTrajPoint.position.z, curTrajPoint.velocity.z, estPos.z, estVel.z, estAtt, curTrajPoint.accel.z, dt);

  // reserve some thrust margin for angle control
  float thrustMargin = .1f*(maxMotorThrust - minMotorThrust);
  collThrustCmd = CONSTRAIN(collThrustCmd, (minMotorThrust+ thrustMargin)*4.f, (maxMotorThrust-thrustMargin)*4.f);
  
  V3F desAcc = LateralPositionControl(curTrajPoint.position, curTrajPoint.velocity, estPos, estVel, curTrajPoint.accel);
  
  V3F desOmega = RollPitchControl(desAcc, estAtt, collThrustCmd);
  desOmega.z = YawControl(curTrajPoint.attitude.Yaw(), estAtt.Yaw());

  V3F desMoment = BodyRateControl(desOmega, estOmega);

  return GenerateMotorCommands(collThrustCmd, desMoment);
}
